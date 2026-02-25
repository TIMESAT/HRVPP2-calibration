#!/usr/bin/env python3
"""
Site-first summarization for *_GPP_VPP-style CSVs.

This script computes metrics in the correct order:
  1) Compute per-site metrics (each site is one statistical unit)
  2) Summarize across sites (equal weight per site), overall and optionally by land cover

Input format (wide settings columns):
  site, lc_id, id, GPP, settings 1, settings 2, ..., settings 6000 ...
Where:
  - id looks like "2017_s1_SOSD" and the VPP name is the last token ("SOSD", "EOSD", "AMPL", ...)

Metrics:
  - For DIRECT_VPPS = {SOSD, EOSD, MAXD, LENGTH}: MAE, RMSE, Bias (per site, over years/rows)
  - For all other VPPs: Pearson R^2 (per site, over years/rows)

Outputs:
  - Site-level table: one row per (site, Setting) [and lc_id if by-lc]
  - Summary table: mean/median across sites per Setting [and lc_id if by-lc], plus n_sites

Usage:
  python summarize_vpp_site_first.py \
    --input PPI_GPP_VPP.csv \
    --out-site-all output/PPI_SITE__ALL.csv \
    --out-sum-all  output/PPI_SUMMARY__ALL.csv \
    --out-site-by-lc output/PPI_SITE__BY_LC.csv \
    --out-sum-by-lc  output/PPI_SUMMARY__BY_LC.csv \
    --workers 24 --batch-size 1024 --engine pyarrow

Notes:
  - This is "site-first": it never pools rows across sites when computing RMSE/R2.
  - Batch-size trades memory vs overhead. For 6000 settings, 512–2048 is a good start.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable, List, Dict, Tuple, Optional

import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

DIRECT_VPPS = {"SOSD", "EOSD", "MAXD", "LENGTH"}


# -------------------------------
# Helpers
# -------------------------------

def detect_settings_columns(df: pd.DataFrame) -> List[str]:
    cols = [c for c in df.columns if c.startswith("settings ")]
    if not cols:
        raise ValueError("No columns starting with 'settings ' found.")
    return cols

def ensure_vpp_column(df: pd.DataFrame) -> pd.DataFrame:
    if "vpp" in df.columns:
        return df
    out = df.copy()
    out["vpp"] = out["id"].astype(str).str.split("_").str[-1]
    return out

def chunk_list(xs: List[str], n: int) -> List[List[str]]:
    return [xs[i:i+n] for i in range(0, len(xs), n)]


# -------------------------------
# Core metric kernels (site-level)
# -------------------------------

def _site_direct_stats(ref: np.ndarray, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Per-column MAE/RMSE/Bias for DIRECT_VPPS inside one site, masking NaNs.
    ref: (n_rows,), X: (n_rows, n_cols)
    Returns: mae, rmse, bias, n (all (n_cols,))
    """
    y = pd.to_numeric(pd.Series(ref), errors="coerce").to_numpy(dtype=float)
    X = np.asarray(X, dtype=float)

    m = np.isfinite(y)[:, None] & np.isfinite(X)
    n = m.sum(axis=0).astype(float)

    n_safe = np.where(n == 0, np.nan, n)
    diffs = np.where(m, X - y[:, None], 0.0)

    mae = np.sum(np.abs(diffs), axis=0) / n_safe
    rmse = np.sqrt(np.sum(diffs * diffs, axis=0) / n_safe)
    bias = np.sum(diffs, axis=0) / n_safe

    return mae, rmse, bias, n

def _site_r2(ref: np.ndarray, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Per-column Pearson r^2 inside one site, masking NaNs independently per column.
    Returns: r2, n (both (n_cols,))
    """
    y = pd.to_numeric(pd.Series(ref), errors="coerce").to_numpy(dtype=float)
    X = np.asarray(X, dtype=float)

    m = np.isfinite(y)[:, None] & np.isfinite(X)
    n = m.sum(axis=0).astype(float)

    r2 = np.full(X.shape[1], np.nan, dtype=float)
    for j in range(X.shape[1]):
        mj = m[:, j]
        if mj.sum() < 3:
            continue
        xj = X[mj, j]
        yj = y[mj]
        if np.nanstd(xj) == 0 or np.nanstd(yj) == 0:
            continue
        r = np.corrcoef(xj, yj)[0, 1]
        r2[j] = r * r

    return r2, n


# -------------------------------
# Worker: compute site-level block for a batch of settings columns
# -------------------------------

def compute_site_block(
    df_slice: pd.DataFrame,
    setting_cols: List[str],
    by_lc: bool
) -> pd.DataFrame:
    """
    Returns site-level metrics for this block:
      keys: site, Setting (+ lc_id if by_lc)
      columns: <VPP>_{MAE,RMSE,Bias,R2,N} for VPPs present in this slice
    """
    df_slice = ensure_vpp_column(df_slice)

    key_site = (["lc_id"] if by_lc else []) + ["site"]
    rows: List[Dict[str, float]] = []

    for site_key, sub_site in df_slice.groupby(key_site, dropna=False):
        # Unpack keys
        if by_lc:
            lc_val, site_val = site_key
        else:
            lc_val, site_val = None, site_key

        # For each VPP variable at this site
        for vpp_name, part in sub_site.groupby("vpp", dropna=False):
            ref = part["GPP"].to_numpy()
            X = part[setting_cols].to_numpy(dtype=float)

            if str(vpp_name) in DIRECT_VPPS:
                mae, rmse, bias, n = _site_direct_stats(ref, X)
                for j, s in enumerate(setting_cols):
                    rec = {
                        "site": site_val,
                        "Setting": s,
                        f"{vpp_name}_MAE": float(mae[j]),
                        f"{vpp_name}_RMSE": float(rmse[j]),
                        f"{vpp_name}_Bias": float(bias[j]),
                        f"{vpp_name}_N": float(n[j]),
                    }
                    if by_lc:
                        rec["lc_id"] = lc_val
                    rows.append(rec)
            else:
                r2, n = _site_r2(ref, X)
                for j, s in enumerate(setting_cols):
                    rec = {
                        "site": site_val,
                        "Setting": s,
                        f"{vpp_name}_R2": float(r2[j]),
                        f"{vpp_name}_N": float(n[j]),
                    }
                    if by_lc:
                        rec["lc_id"] = lc_val
                    rows.append(rec)

    if not rows:
        key_cols = (["lc_id"] if by_lc else []) + ["site", "Setting"]
        return pd.DataFrame(columns=key_cols)

    out = pd.DataFrame(rows)

    # Collapse: one row per (site, Setting) by taking first non-null across partial metric rows
    key_cols = (["lc_id"] if by_lc else []) + ["site", "Setting"]
    out = out.groupby(key_cols, dropna=False, as_index=False).first()
    return out


# -------------------------------
# Orchestrator: site-level then summaries
# -------------------------------

def build_site_table(df: pd.DataFrame, by_lc: bool, workers: int, batch_size: int) -> pd.DataFrame:
    df = ensure_vpp_column(df)
    settings = detect_settings_columns(df)

    keep = ["site", "GPP", "vpp"] + (["lc_id"] if by_lc else []) + settings
    df = df[keep]

    batches = chunk_list(settings, batch_size)
    parts: List[pd.DataFrame] = []

    if workers <= 1:
        for cols in batches:
            df_slice = df[["site", "GPP", "vpp"] + (["lc_id"] if by_lc else []) + cols]
            parts.append(compute_site_block(df_slice, cols, by_lc))
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = []
            for cols in batches:
                df_slice = df[["site", "GPP", "vpp"] + (["lc_id"] if by_lc else []) + cols]
                futs.append(ex.submit(compute_site_block, df_slice, cols, by_lc))
            for f in as_completed(futs):
                parts.append(f.result())

    if not parts:
        key_cols = (["lc_id"] if by_lc else []) + ["site", "Setting"]
        return pd.DataFrame(columns=key_cols)

    # Stitch blocks: concat then groupby-first (fills metrics from different blocks)
    key_cols = (["lc_id"] if by_lc else []) + ["site", "Setting"]
    out = pd.concat(parts, ignore_index=True)
    out = out.groupby(key_cols, dropna=False, as_index=False).first()

    # Stable ordering
    sort_cols = (["lc_id"] if by_lc else []) + ["site", "Setting"]
    out.sort_values(sort_cols, inplace=True, kind="mergesort")

    return out

def summarize_across_sites(site_tbl: pd.DataFrame, by_lc: bool) -> pd.DataFrame:
    key_cols = (["lc_id"] if by_lc else []) + ["Setting"]

    metric_cols = [c for c in site_tbl.columns if c not in key_cols + ["site"]]
    if not metric_cols:
        out = site_tbl[key_cols].drop_duplicates().copy()
        out["n_sites"] = site_tbl.groupby(key_cols, dropna=False)["site"].nunique().values
        return out

    agg = {c: ["mean", "median"] for c in metric_cols}
    out = site_tbl.groupby(key_cols, dropna=False).agg(agg)
    out.columns = ["__".join(x) for x in out.columns.to_flat_index()]
    out = out.reset_index()

    out["n_sites"] = site_tbl.groupby(key_cols, dropna=False)["site"].nunique().values

    # Stable ordering
    sort_cols = (["lc_id"] if by_lc else []) + ["Setting"]
    out.sort_values(sort_cols, inplace=True, kind="mergesort")

    return out


# -------------------------------
# IO
# -------------------------------

def read_csv_fast(path: Path, engine: Optional[str]) -> pd.DataFrame:
    read_kwargs = {}
    if engine:
        read_kwargs["engine"] = engine
    else:
        # prefer pyarrow if available
        try:
            read_kwargs["engine"] = "pyarrow"
        except Exception:
            pass
    return pd.read_csv(path, **read_kwargs)

def write_parent(p: Path) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)

def run(
    inp: Path,
    out_site_all: Path,
    out_sum_all: Path,
    out_site_by_lc: Optional[Path],
    out_sum_by_lc: Optional[Path],
    workers: int,
    batch_size: int,
    engine: Optional[str],
) -> None:
    df = read_csv_fast(inp, engine=engine)

    # Overall site table and summary
    site_all = build_site_table(df, by_lc=False, workers=workers, batch_size=batch_size)
    sum_all = summarize_across_sites(site_all, by_lc=False)

    write_parent(out_site_all)
    write_parent(out_sum_all)
    site_all.to_csv(out_site_all, index=False)
    sum_all.to_csv(out_sum_all, index=False)

    # By land cover (optional)
    if out_site_by_lc and out_sum_by_lc:
        site_lc = build_site_table(df, by_lc=True, workers=workers, batch_size=batch_size)
        sum_lc = summarize_across_sites(site_lc, by_lc=True)

        write_parent(out_site_by_lc)
        write_parent(out_sum_by_lc)
        site_lc.to_csv(out_site_by_lc, index=False)
        sum_lc.to_csv(out_sum_by_lc, index=False)

    print(f"Wrote site-level (ALL): {out_site_all.resolve()}")
    print(f"Wrote summary    (ALL): {out_sum_all.resolve()}")
    if out_site_by_lc and out_sum_by_lc:
        print(f"Wrote site-level (BY_LC): {out_site_by_lc.resolve()}")
        print(f"Wrote summary    (BY_LC): {out_sum_by_lc.resolve()}")


# -------------------------------
# CLI
# -------------------------------

def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Site-first summarization for *_GPP_VPP.csv files.")
    ap.add_argument("--input", required=False, default="PPI_GPP_VPP.csv", help="Input *_GPP_VPP.csv")
    ap.add_argument("--out-site-all", required=False, default="PPI_GPP_VPP_SITE_ALL.csv", help="Output site-level CSV (overall)")
    ap.add_argument("--out-sum-all", required=False, default="PPI_GPP_VPP_SUM_ALL.csv", help="Output summary CSV (overall, aggregated across sites)")

    ap.add_argument("--out-site-by-lc", default=None, help="Optional output site-level CSV grouped by lc_id")
    ap.add_argument("--out-sum-by-lc", default=None, help="Optional output summary CSV grouped by lc_id")

    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 1) - 1),
                    help="Process workers (default: cpu_count-1). Use 1 to disable multiprocessing.")
    ap.add_argument("--batch-size", type=int, default=1024,
                    help="Number of settings columns per batch (e.g., 512–2048).")
    ap.add_argument("--engine", choices=["pyarrow", "c", "python"], default=None,
                    help="CSV engine (default: try pyarrow).")

    args = ap.parse_args(argv)

    inp = Path(args.input)
    out_site_all = Path(args.out_site_all)
    out_sum_all = Path(args.out_sum_all)

    out_site_by_lc = Path(args.out_site_by_lc) if args.out_site_by_lc else None
    out_sum_by_lc = Path(args.out_sum_by_lc) if args.out_sum_by_lc else None

    # If one lc output is set, require both to avoid confusion
    if (out_site_by_lc is None) ^ (out_sum_by_lc is None):
        raise SystemExit("Provide both --out-site-by-lc and --out-sum-by-lc, or neither.")

    run(
        inp=inp,
        out_site_all=out_site_all,
        out_sum_all=out_sum_all,
        out_site_by_lc=out_site_by_lc,
        out_sum_by_lc=out_sum_by_lc,
        workers=args.workers,
        batch_size=args.batch_size,
        engine=args.engine,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
