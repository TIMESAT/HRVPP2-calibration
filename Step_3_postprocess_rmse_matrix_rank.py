#!/usr/bin/env python3
"""
Compute per-setting win rates from binary matrices and create joint SOS+EOS binary.

Inputs:
  - SOS_RMSE_binary.csv  (rows=site, cols=settings, values 0/1 or NA)
  - EOS_RMSE_binary.csv  (rows=site, cols=settings, values 0/1 or NA)

Outputs:
  1) SOS_setting_winrate.csv        (near-optimal frequency by setting)
  2) EOS_setting_winrate.csv
  3) JOINT_SOS_EOS_binary.csv       (AND by default)
  4) JOINT_setting_winrate.csv

Notes:
  - "win rate" here means fraction of sites where binary==1 (near-optimal frequency),
    not necessarily strict rank-1 unless you use the optional strict mode.
  - Missing values are ignored per setting when computing rates (configurable).

Usage:
  python postprocess_binary_wins.py \
    --sos SOS_RMSE_binary.csv \
    --eos EOS_RMSE_binary.csv \
    --out-dir wins \
    --joint-op and \
    --min-sites 10
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def _read_binary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    # ensure numeric 0/1 with possible NA
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    df.index = df.index.astype(str)
    return df


def _align_matrices(a: pd.DataFrame, b: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # align on common sites and common settings columns
    common_sites = a.index.intersection(b.index)
    common_cols = a.columns.intersection(b.columns)
    if len(common_sites) == 0 or len(common_cols) == 0:
        raise ValueError("No overlapping sites or settings between SOS and EOS matrices.")
    a2 = a.loc[common_sites, common_cols].copy()
    b2 = b.loc[common_sites, common_cols].copy()
    return a2, b2


def winrate_from_binary(
    bin_df: pd.DataFrame,
    min_sites: int = 1,
    na_policy: str = "ignore",
) -> pd.DataFrame:
    """
    Compute near-optimal frequency per setting from a binary matrix.

    na_policy:
      - "ignore": rate = sum(1) / count(non-NA)
      - "zero":   treat NA as 0, rate = sum(1) / total_sites
    """
    if na_policy not in {"ignore", "zero"}:
        raise ValueError("na_policy must be 'ignore' or 'zero'")

    n_sites_total = bin_df.shape[0]

    if na_policy == "zero":
        x = bin_df.fillna(0).astype(int)
        wins = x.sum(axis=0)
        denom = pd.Series(n_sites_total, index=bin_df.columns)
    else:
        wins = bin_df.sum(axis=0, skipna=True)
        denom = bin_df.notna().sum(axis=0)

    rate = (wins / denom).astype(float)
    out = pd.DataFrame({
        "Setting": bin_df.columns.astype(str),
        "win_sites": wins.astype(float).values,
        "n_sites_used": denom.astype(float).values,
        "win_rate": rate.values
    })

    # filter weakly-supported settings
    out = out[out["n_sites_used"] >= float(min_sites)].copy()
    out.sort_values(["win_rate", "win_sites"], ascending=[False, False], inplace=True)
    out.reset_index(drop=True, inplace=True)
    return out


def joint_binary(
    sos: pd.DataFrame,
    eos: pd.DataFrame,
    op: str = "and"
) -> pd.DataFrame:
    """
    Create joint binary matrix from SOS and EOS binaries.

    op:
      - "and": 1 if SOS==1 AND EOS==1
      - "or":  1 if SOS==1 OR  EOS==1
    NA handling:
      - if either is NA, result is NA (conservative)
    """
    if op not in {"and", "or"}:
        raise ValueError("op must be 'and' or 'or'")

    # Keep NA if either is NA
    na_mask = sos.isna() | eos.isna()

    if op == "and":
        jb = (sos.fillna(0).astype(int) & eos.fillna(0).astype(int)).astype("Int64")
    else:
        jb = (sos.fillna(0).astype(int) | eos.fillna(0).astype(int)).astype("Int64")

    jb = jb.mask(na_mask)  # restore NA where information missing
    return jb


def main() -> None:
    ap = argparse.ArgumentParser(description="Compute per-setting win rates and joint SOS+EOS binary.")
    ap.add_argument("--sos", required=False, default="VPP_RMSE_matrix_SP10000/SOS_RMSE_binary.csv", help="SOS_RMSE_binary.csv")
    ap.add_argument("--eos", required=False, default="VPP_RMSE_matrix_SP10000/EOS_RMSE_binary.csv",help="EOS_RMSE_binary.csv")
    ap.add_argument("--out-dir", required=False, default="VPP_RMSE_matrix_SP10000/", help="Output directory")
    ap.add_argument("--joint-op", choices=["and", "or"], default="and", help="Joint rule (default: and)")
    ap.add_argument("--na-policy", choices=["ignore", "zero"], default="ignore",
                    help="How to treat NA when computing win rates (default: ignore)")
    ap.add_argument("--min-sites", type=int, default=1,
                    help="Require at least this many sites contributing to a setting's rate (default: 1)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sos = _read_binary(Path(args.sos))
    eos = _read_binary(Path(args.eos))
    sos, eos = _align_matrices(sos, eos)

    # 1) win rates separately
    wr_sos = winrate_from_binary(sos, min_sites=args.min_sites, na_policy=args.na_policy)
    wr_eos = winrate_from_binary(eos, min_sites=args.min_sites, na_policy=args.na_policy)

    wr_sos.to_csv(out_dir / "SOS_setting_winrate.csv", index=False)
    wr_eos.to_csv(out_dir / "EOS_setting_winrate.csv", index=False)

    # 2) joint binary
    jb = joint_binary(sos, eos, op=args.joint_op)
    jb.to_csv(out_dir / "JOINT_SOS_EOS_binary.csv")

    # 3) joint win rates
    wr_joint = winrate_from_binary(jb, min_sites=args.min_sites, na_policy=args.na_policy)
    wr_joint.to_csv(out_dir / "JOINT_setting_winrate.csv", index=False)

    print("Wrote:")
    print(f"  {out_dir / 'SOS_setting_winrate.csv'}")
    print(f"  {out_dir / 'EOS_setting_winrate.csv'}")
    print(f"  {out_dir / 'JOINT_SOS_EOS_binary.csv'}")
    print(f"  {out_dir / 'JOINT_setting_winrate.csv'}")


if __name__ == "__main__":
    main()
