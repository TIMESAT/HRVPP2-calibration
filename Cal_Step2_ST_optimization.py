#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import Step_2_state_function
import Step_2_ts_function


def parse_site_lc_from_ppi(path):
    name = Path(path).name
    m = re.match(
        r"^PPI_\d{4}_\d{2}_\d{2}_\d{4}_\d{2}_\d{2}_(?P<site>[^_]+)_(?P<lc>LC\d+)_mean\.csv$",
        name,
    )
    if not m:
        raise ValueError(f"Cannot parse PPI filename: {name}")
    return m.group("site"), m.group("lc")


def parse_site_from_gpp(path):
    name = Path(path).name
    m = re.search(r"_([^_]+)_DAILY_GPP", name)
    if not m:
        raise ValueError(f"Cannot parse site from GPP file: {name}")
    return m.group(1)


def collect_file_paths(ppi_dir, gpp_dir):
    files = defaultdict(lambda: {"ppi": None, "gpp": None})

    for f in Path(ppi_dir).glob("*.csv"):
        try:
            site, lc = parse_site_lc_from_ppi(f)
            files[(site, lc)]["ppi"] = str(f)
        except ValueError:
            continue

    gpp_map = {}
    for f in Path(gpp_dir).glob("*.csv"):
        try:
            site = parse_site_from_gpp(f)
            gpp_map[site] = str(f)
        except ValueError:
            continue

    for site, lc in list(files.keys()):
        if site in gpp_map:
            files[(site, lc)]["gpp"] = gpp_map[site]

    return files


def read_ppi_daily(path):
    df = pd.read_csv(path)
    df["t"] = pd.to_datetime(df["t"]).dt.date
    df = df.rename(columns={"PPI_mean": "PPI"})[["t", "PPI"]]
    df["PPI"] = df["PPI"].clip(-1, 5)
    df["PPI"] = df["PPI"].clip(lower=0)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df


def read_gpp_daily(path):
    df = pd.read_csv(path)
    df["t"] = pd.to_datetime(df["t"]).dt.date
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df[["t", "GPP"]]


def diagnose_vi(merged, vi, gpp_col="GPP"):
    df = merged[[vi, gpp_col]].dropna()
    n = len(df)
    if n < 5:
        return None

    rho_all, _ = spearmanr(df[vi], df[gpp_col])

    gpp = df[gpp_col].values
    vi_values = df[vi].values
    threshold = np.nanpercentile(gpp, 5)
    gs_mask = gpp > threshold

    if gs_mask.sum() >= 5:
        rho_gs, _ = spearmanr(gpp[gs_mask], vi_values[gs_mask])
    else:
        rho_gs = np.nan

    return {
        "VI": vi,
        "n": n,
        "rho_all": rho_all,
        "rho_gs": rho_gs,
    }


def merge_ppi_gpp_for_site_lc(site, lc, info):
    if not info["ppi"]:
        raise ValueError(f"Missing PPI for {site}, {lc}")
    if not info["gpp"]:
        raise ValueError(f"Missing GPP for {site}, {lc}")

    ppi = read_ppi_daily(info["ppi"])
    gpp = read_gpp_daily(info["gpp"])

    merged = pd.merge(ppi, gpp, on="t", how="outer").sort_values("t")
    merged.replace([np.inf, -np.inf], np.nan, inplace=True)

    merged_cut = merged.dropna(subset=["PPI", "GPP"], how="any")

    diag = diagnose_vi(merged_cut, "PPI", gpp_col="GPP")
    if diag is None:
        diag_df = pd.DataFrame(columns=["VI", "n", "rho_all", "rho_gs", "site", "lc"])
    else:
        diag["site"] = site
        diag["lc"] = lc
        diag_df = pd.DataFrame([diag])

    return merged[["t", "PPI", "GPP"]], merged_cut[["t", "PPI", "GPP"]], diag_df


def summarize_raw_diags(raw_diags, out_vi_dir):
    if not raw_diags:
        return

    df_diags_all = pd.concat(raw_diags, ignore_index=True)
    out_csv = Path(out_vi_dir) / "Raw_PPI_GPP_diag_ALL_sites.csv"
    df_diags_all.to_csv(out_csv, index=False)
    print(f"[INFO] wrote {out_csv}")

    summary = (
        df_diags_all.groupby("VI")
        .agg(
            mean_rho_gs=("rho_gs", "mean"),
            median_rho_gs=("rho_gs", "median"),
            num_sites=("site", "nunique"),
        )
        .reset_index()
    )
    out_csv = Path(out_vi_dir) / "Raw_PPI_GPP_diag_summary.csv"
    summary.to_csv(out_csv, index=False)
    print(f"[INFO] wrote {out_csv}")


def summarize_st_diags(st_diags, out_vi_dir):
    if not st_diags:
        return

    df_diags_st_all = pd.concat(st_diags, ignore_index=True)

    out_csv_long = Path(out_vi_dir) / "ST_PPI_GPP_diag_ALL_sites_long.csv"
    df_diags_st_all.to_csv(out_csv_long, index=False)
    print(f"[INFO] wrote {out_csv_long}")

    df_unique = (
        df_diags_st_all.sort_values("n", ascending=False)
        .drop_duplicates(subset=["site", "lc", "VI"])
    )

    wide = df_unique.pivot_table(
        index=["site", "lc"],
        columns="VI",
        values="rho_gs",
    )

    n_per_site_lc = (
        df_unique.groupby(["site", "lc"])["n"]
        .min()
        .rename("n")
    )

    wide = wide.join(n_per_site_lc)
    cols = ["n"] + [c for c in wide.columns if c != "n"]
    wide = wide[cols].reset_index()
    wide.columns.name = None

    out_csv_wide = Path(out_vi_dir) / "ST_PPI_GPP_diag_ALL_sites.csv"
    wide.to_csv(out_csv_wide, index=False)
    print(f"[INFO] wrote {out_csv_wide}")

    out_csv_mask = Path(out_vi_dir) / "ST_PPI_GPP_diag_ALL_sites_bestmask.csv"
    Step_2_state_function.make_best_vi_mask(out_csv_wide, out_csv_mask, alpha=0.05)


def main():
    p = argparse.ArgumentParser(
        description="TIMESAT optimization for PPI only, evaluated against measured GPP."
    )
    p.add_argument(
        "--ppi_dir",
        type=str,
        default="/Users/zzcai/Library/CloudStorage/OneDrive-LundUniversity/HR-VPP2/calval/VI/PPI/flux_cal/V1/csv_lc",
    )
    p.add_argument(
        "--gpp_dir",
        type=str,
        default="/Users/zzcai/Library/CloudStorage/OneDrive-LundUniversity/HR-VPP2/calval/VI/GPP_NT_VUT_MEAN",
    )
    p.add_argument("--raw_vi", type=int, default=0)
    p.add_argument("--st_vi", type=int, default=1)
    p.add_argument("--vpp_vi", type=int, default=0)
    p.add_argument("--out_vi_dir", type=str, default="output/Cal/PPI_GPP")
    p.add_argument("--out_plot_dir", type=str, default="output/Cal/PPI_GPP")
    args = p.parse_args()

    Path(args.out_vi_dir).mkdir(parents=True, exist_ok=True)
    Path(args.out_plot_dir).mkdir(parents=True, exist_ok=True)

    files = collect_file_paths(args.ppi_dir, args.gpp_dir)
    print(f"Found {len(files)} site/LC groups")

    rawvi_diags = []
    st_diags = []

    for (site, lc), info in sorted(files.items()):
        try:
            merged, merged_cut, diag_raw = merge_ppi_gpp_for_site_lc(site, lc, info)
        except Exception as e:
            print(f"[WARN] {site},{lc}: {e}")
            continue

        if args.raw_vi:
            out_csv = Path(args.out_vi_dir) / f"Raw_PPI_GPP_data_{site}_{lc}.csv"
            merged.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")

            out_csv = Path(args.out_vi_dir) / f"Raw_PPI_GPP_diag_{site}_{lc}.csv"
            diag_raw.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")

            if not diag_raw.empty:
                rawvi_diags.append(diag_raw)

        if args.st_vi or args.vpp_vi:
            print(f"[INFO] site={site}, lc={lc}, vi=PPI")
            ppi_slice = merged[["t", "PPI"]].copy()
            yfit_data, vpp_data, settings_df = Step_2_ts_function._ts_run_(ppi_slice, args.vpp_vi)

            if isinstance(yfit_data, int) or yfit_data is None or len(yfit_data) == 0:
                print(f"[WARN] No ST output for {site},{lc}")
                continue

            yfit_data = yfit_data.copy()
            yfit_data["t"] = pd.to_datetime(yfit_data["t"])

            gpp_data = merged[["t", "GPP"]].copy()
            gpp_data["t"] = pd.to_datetime(gpp_data["t"])

            merged_st = pd.merge(yfit_data, gpp_data, on="t", how="outer").sort_values("t")

            out_csv = Path(args.out_vi_dir) / f"Fitted_PPI_GPP_data_{site}_{lc}.csv"
            merged_st.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")

            out_csv = Path(args.out_vi_dir) / f"Fitted_PPI_GPP_settings_{site}_{lc}.csv"
            settings_df.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")

            if args.vpp_vi and not isinstance(vpp_data, int):
                out_csv = Path(args.out_vi_dir) / f"Fitted_PPI_GPP_vpp_{site}_{lc}.csv"
                vpp_data.to_csv(out_csv, index=False)
                print(f"[INFO] wrote {out_csv}")

            st_diags_list = []
            st_cols = [
                c for c in merged_st.columns
                if c not in {"t", "GPP"} and not str(c).endswith("_qa")
            ]

            for vi in st_cols:
                d = diagnose_vi(merged_st, vi, gpp_col="GPP")
                if d is None:
                    continue
                d["site"] = site
                d["lc"] = lc
                st_diags_list.append(d)

            diag_st = pd.DataFrame(st_diags_list)
            out_csv = Path(args.out_vi_dir) / f"Fitted_PPI_GPP_diag_{site}_{lc}.csv"
            diag_st.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")

            if not diag_st.empty:
                st_diags.append(diag_st)

    summarize_raw_diags(rawvi_diags, args.out_vi_dir)
    summarize_st_diags(st_diags, args.out_vi_dir)

    if not rawvi_diags and not st_diags:
        print("[WARN] no diagnostics found")


if __name__ == "__main__":
    main()
