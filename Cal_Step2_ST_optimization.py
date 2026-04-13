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


def parse_site_lc_from_ppi(path):
    name = Path(path).name
    m = re.match(
        r"^PPI_\d{4}_\d{2}_\d{2}_\d{4}_\d{2}_\d{2}_(?P<site>[^_]+)_(?P<lc>LC\d+)_mean\.csv$",
        name,
    )
    if not m:
        raise ValueError(f"Cannot parse PPI filename: {name}")
    return m.group("site"), m.group("lc")


def parse_site_lc_from_nbar_ppi(path):
    name = Path(path).name
    m = re.match(
        r"^S2_NBAR_SZA11_(?P<site>[^_]+)_(?P<lc>LC\d+)_NDVI_NIRv_EVI2_mean\.csv$",
        name,
    )
    if not m:
        raise ValueError(f"Cannot parse NBAR PPI filename: {name}")
    return m.group("site"), m.group("lc")


def parse_site_from_gpp(path):
    name = Path(path).name
    m = re.search(r"_([^_]+)_DAILY_GPP", name)
    if not m:
        raise ValueError(f"Cannot parse site from GPP file: {name}")
    return m.group(1)


def collect_file_paths(ppi_dir, nbar_ppi_dir, gpp_dir):
    files = defaultdict(lambda: {"ppi": None, "nbar_ppi": None, "gpp": None})

    for f in Path(ppi_dir).glob("*.csv"):
        try:
            site, lc = parse_site_lc_from_ppi(f)
            files[(site, lc)]["ppi"] = str(f)
        except ValueError:
            continue

    for f in Path(nbar_ppi_dir).glob("*.csv"):
        try:
            site, lc = parse_site_lc_from_nbar_ppi(f)
            files[(site, lc)]["nbar_ppi"] = str(f)
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


def read_nbar_ppi_daily(path):
    df = pd.read_csv(path)
    df["t"] = pd.to_datetime(df["t"]).dt.date
    df = df.rename(columns={"PPI": "NBAR_PPI"})[["t", "NBAR_PPI"]]
    df["NBAR_PPI"] = df["NBAR_PPI"].clip(-1, 5)
    df["NBAR_PPI"] = df["NBAR_PPI"].clip(lower=0)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df


def read_gpp_daily(path):
    df = pd.read_csv(path)
    df["t"] = pd.to_datetime(df["t"]).dt.date
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df[["t", "GPP"]]


def diagnose_pair(merged, x_col, y_col, label):
    df = merged[[x_col, y_col]].dropna()
    n = len(df)
    if n < 5:
        return None

    rho_all, _ = spearmanr(df[x_col], df[y_col])

    y_values = df[y_col].values
    x_values = df[x_col].values
    threshold = np.nanpercentile(y_values, 5)
    gs_mask = y_values > threshold

    if gs_mask.sum() >= 5:
        rho_gs, _ = spearmanr(x_values[gs_mask], y_values[gs_mask])
    else:
        rho_gs = np.nan

    return {
        "VI": label,
        "x_col": x_col,
        "y_col": y_col,
        "n": n,
        "rho_all": rho_all,
        "rho_gs": rho_gs,
    }


def diagnose_vi(merged, vi, gpp_col="GPP"):
    return diagnose_pair(merged, vi, gpp_col, vi)


def merge_ppi_gpp_for_site_lc(site, lc, info):
    if not info["ppi"]:
        raise ValueError(f"Missing PPI for {site}, {lc}")
    if not info["nbar_ppi"]:
        raise ValueError(f"Missing NBAR PPI for {site}, {lc}")
    if not info["gpp"]:
        raise ValueError(f"Missing GPP for {site}, {lc}")

    ppi = read_ppi_daily(info["ppi"])
    nbar_ppi = read_nbar_ppi_daily(info["nbar_ppi"])
    gpp = read_gpp_daily(info["gpp"])

    merged = pd.merge(ppi, nbar_ppi, on="t", how="outer")
    merged = pd.merge(merged, gpp, on="t", how="outer").sort_values("t")
    merged.replace([np.inf, -np.inf], np.nan, inplace=True)

    diag_specs = [
        ("PPI", "NBAR_PPI", "PPI_vs_NBAR_PPI"),
        ("PPI", "GPP", "GPP_vs_PPI"),
        ("NBAR_PPI", "GPP", "GPP_vs_NBAR_PPI"),
    ]
    diag_rows = []
    for x_col, y_col, label in diag_specs:
        diag = diagnose_pair(merged, x_col, y_col, label)
        if diag is None:
            continue
        diag["site"] = site
        diag["lc"] = lc
        diag_rows.append(diag)

    diag_df = pd.DataFrame(
        diag_rows,
        columns=["VI", "x_col", "y_col", "n", "rho_all", "rho_gs", "site", "lc"],
    )

    return merged[["t", "PPI", "NBAR_PPI", "GPP"]], diag_df


def summarize_raw_diags(raw_diags, out_vi_dir):
    if not raw_diags:
        return

    df_diags_all = pd.concat(raw_diags, ignore_index=True)
    out_csv = Path(out_vi_dir) / "Raw_PPI_NBARPPI_GPP_diag_ALL_sites.csv"
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
    out_csv = Path(out_vi_dir) / "Raw_PPI_NBARPPI_GPP_diag_summary.csv"
    summary.to_csv(out_csv, index=False)
    print(f"[INFO] wrote {out_csv}")


def summarize_st_diags(st_diags, out_vi_dir):
    if not st_diags:
        return

    df_diags_st_all = pd.concat(st_diags, ignore_index=True)

    out_csv_long = Path(out_vi_dir) / "ST_PPI_NBARPPI_GPP_diag_ALL_sites_long.csv"
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

    out_csv_wide = Path(out_vi_dir) / "ST_PPI_NBARPPI_GPP_diag_ALL_sites.csv"
    wide.to_csv(out_csv_wide, index=False)
    print(f"[INFO] wrote {out_csv_wide}")

    out_csv_mask = Path(out_vi_dir) / "ST_PPI_NBARPPI_GPP_diag_ALL_sites_bestmask.csv"
    Step_2_state_function.make_best_vi_mask(out_csv_wide, out_csv_mask, alpha=0.05)


def rename_timesat_columns(df, old_prefix, new_prefix):
    rename_map = {}
    for col in df.columns:
        col_str = str(col)
        if col_str == old_prefix or col_str.startswith(f"{old_prefix}-") or col_str.startswith(f"{old_prefix}_"):
            rename_map[col] = col_str.replace(old_prefix, new_prefix, 1)
    return df.rename(columns=rename_map)


def run_timesat_for_series(merged, input_col, output_prefix, run_vpp):
    import Step_2_ts_function

    indata = merged[["t", input_col]].copy().rename(columns={input_col: "PPI"})
    yfit_data, vpp_data, settings_df = Step_2_ts_function._ts_run_(indata, run_vpp)

    if isinstance(yfit_data, int) or yfit_data is None or len(yfit_data) == 0:
        return None, None, None

    yfit_data = yfit_data.copy()
    yfit_data["t"] = pd.to_datetime(yfit_data["t"])
    yfit_data = rename_timesat_columns(yfit_data, "PPI", output_prefix)

    settings_df = settings_df.copy()
    if "settings_id" in settings_df.columns:
        settings_df["settings_id"] = settings_df["settings_id"].astype(str).str.replace(
            r"^PPI",
            output_prefix,
            n=1,
            regex=True,
        )

    if not isinstance(vpp_data, int):
        vpp_data = vpp_data.copy()
        if "id" in vpp_data.columns:
            vpp_data = rename_timesat_columns(vpp_data, "PPI", output_prefix)

    return yfit_data, vpp_data, settings_df


def summarize_best_settings_by_class(out_vi_dir, example_settings_csv, top_n=3):
    from Cal_Step2_best_settings_by_class import build_best_settings_by_class

    diag_long_csv = Path(out_vi_dir) / "ST_PPI_NBARPPI_GPP_diag_ALL_sites_long.csv"
    if not diag_long_csv.exists():
        print(f"[WARN] missing {diag_long_csv}, skip class summary")
        return

    if not example_settings_csv or not Path(example_settings_csv).exists():
        print("[WARN] missing example settings CSV, skip class summary")
        return

    build_best_settings_by_class(
        diag_long_csv=diag_long_csv,
        settings_csv=example_settings_csv,
        out_dir=out_vi_dir,
        top_n=top_n,
    )


def main():
    p = argparse.ArgumentParser(
        description="TIMESAT optimization for PPI, with raw comparisons against NBAR PPI and measured GPP."
    )
    p.add_argument(
        "--ppi_dir",
        type=str,
        default="Data/VI/PPI/flux_cal/V1/csv_lc",
    )
    p.add_argument(
        "--nbar_ppi_dir",
        type=str,
        default="Data/VI/DVI_NIRv_EVI2_nbar_sza11/flux_cal/V1/csv_lc",
    )
    p.add_argument(
        "--gpp_dir",
        type=str,
        default="Data/VI/GPP_NT_VUT_MEAN",
    )
    p.add_argument("--raw_vi", type=int, default=0)
    p.add_argument("--st_vi", type=int, default=1)
    p.add_argument("--vpp_vi", type=int, default=0)
    p.add_argument("--out_vi_dir", type=str, default="output/Cal/PPI_GPP")
    p.add_argument("--out_plot_dir", type=str, default="output/Cal/PPI_GPP")
    args = p.parse_args()

    Path(args.out_vi_dir).mkdir(parents=True, exist_ok=True)
    Path(args.out_plot_dir).mkdir(parents=True, exist_ok=True)

    files = collect_file_paths(args.ppi_dir, args.nbar_ppi_dir, args.gpp_dir)
    print(f"Found {len(files)} site/LC groups")

    rawvi_diags = []
    st_diags = []
    first_settings_csv = None

    for (site, lc), info in sorted(files.items()):
        try:
            merged, diag_raw = merge_ppi_gpp_for_site_lc(site, lc, info)
        except Exception as e:
            print(f"[WARN] {site},{lc}: {e}")
            continue

        if args.raw_vi:
            out_csv = Path(args.out_vi_dir) / f"Raw_PPI_NBARPPI_GPP_data_{site}_{lc}.csv"
            merged.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")

            out_csv = Path(args.out_vi_dir) / f"Raw_PPI_NBARPPI_GPP_diag_{site}_{lc}.csv"
            diag_raw.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")

            if not diag_raw.empty:
                rawvi_diags.append(diag_raw)

        if args.st_vi or args.vpp_vi:
            print(f"[INFO] site={site}, lc={lc}, vi=PPI")
            ppi_fit, ppi_vpp, ppi_settings = run_timesat_for_series(
                merged, "PPI", "ST_PPI", args.vpp_vi
            )
            if ppi_fit is None:
                print(f"[WARN] No ST output for PPI at {site},{lc}")
                continue

            print(f"[INFO] site={site}, lc={lc}, vi=NBAR_PPI")
            nbar_fit, nbar_vpp, nbar_settings = run_timesat_for_series(
                merged, "NBAR_PPI", "ST_NBAR_PPI", args.vpp_vi
            )
            if nbar_fit is None:
                print(f"[WARN] No ST output for NBAR_PPI at {site},{lc}")
                continue

            gpp_data = merged[["t", "GPP"]].copy()
            gpp_data["t"] = pd.to_datetime(gpp_data["t"])

            merged_st = pd.merge(ppi_fit, nbar_fit, on="t", how="outer")
            merged_st = pd.merge(merged_st, gpp_data, on="t", how="outer").sort_values("t")

            out_csv = Path(args.out_vi_dir) / f"Fitted_PPI_NBARPPI_GPP_data_{site}_{lc}.csv"
            merged_st.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")

            settings_all = pd.concat([ppi_settings, nbar_settings], ignore_index=True)
            out_csv = Path(args.out_vi_dir) / f"Fitted_PPI_NBARPPI_GPP_settings_{site}_{lc}.csv"
            settings_all.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")
            if first_settings_csv is None:
                first_settings_csv = str(out_csv)

            if args.vpp_vi and not isinstance(ppi_vpp, int) and not isinstance(nbar_vpp, int):
                vpp_all = pd.merge(ppi_vpp, nbar_vpp, on="id", how="outer")
                out_csv = Path(args.out_vi_dir) / f"Fitted_PPI_NBARPPI_GPP_vpp_{site}_{lc}.csv"
                vpp_all.to_csv(out_csv, index=False)
                print(f"[INFO] wrote {out_csv}")

            st_diags_list = []
            diag_specs = []
            st_ppi_cols = [
                c for c in merged_st.columns
                if str(c).startswith("ST_PPI") and not str(c).endswith("_qa")
            ]
            st_nbar_cols = [
                c for c in merged_st.columns
                if str(c).startswith("ST_NBAR_PPI") and not str(c).endswith("_qa")
            ]

            for vi in st_ppi_cols:
                diag_specs.append((vi, "GPP", "GPP_vs_ST_PPI"))
            for vi in st_nbar_cols:
                diag_specs.append((vi, "GPP", "GPP_vs_ST_NBAR_PPI"))

            for ppi_col in st_ppi_cols:
                suffix = str(ppi_col).replace("ST_PPI", "", 1)
                nbar_col = f"ST_NBAR_PPI{suffix}"
                if nbar_col in merged_st.columns:
                    diag_specs.append((ppi_col, nbar_col, "ST_PPI_vs_ST_NBAR_PPI"))

            for x_col, y_col, label in diag_specs:
                d = diagnose_pair(merged_st, x_col, y_col, label)
                if d is None:
                    continue
                d["site"] = site
                d["lc"] = lc
                st_diags_list.append(d)

            diag_st = pd.DataFrame(st_diags_list)
            out_csv = Path(args.out_vi_dir) / f"Fitted_PPI_NBARPPI_GPP_diag_{site}_{lc}.csv"
            diag_st.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")

            if not diag_st.empty:
                st_diags.append(diag_st)

    summarize_raw_diags(rawvi_diags, args.out_vi_dir)
    summarize_st_diags(st_diags, args.out_vi_dir)

    if st_diags:
        try:
            summarize_best_settings_by_class(args.out_vi_dir, first_settings_csv, top_n=3)
        except Exception as e:
            print(f"[WARN] failed to summarize best settings by class: {e}")

    if not rawvi_diags and not st_diags:
        print("[WARN] no diagnostics found")


if __name__ == "__main__":
    main()
