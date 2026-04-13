#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import timesat

import Step_2_ts_function


CLASS_GROUPS = {
    "Forest": {
        "landcover_codes": [7, 8, 16],
        "selected_setting": "ST_NBAR_PPI-15",
        "smooth": 10000.0,
        "p25_rho_gs": 0.793725,
    },
    "Grassland": {
        "landcover_codes": [2, 5, 9, 10, 14, 15],
        "selected_setting": "ST_NBAR_PPI-14",
        "smooth": 3000.0,
        "p25_rho_gs": 0.733175,
    },
    "Cropland": {
        "landcover_codes": [11, 12, 13],
        "selected_setting": "ST_NBAR_PPI-12",
        "smooth": 300.0,
        "p25_rho_gs": 0.771976,
    },
}


def build_lc_meta():
    lc_to_meta = {}
    for class_name, meta in CLASS_GROUPS.items():
        for code in meta["landcover_codes"]:
            lc_to_meta[f"LC{code}"] = {
                "class": class_name,
                "selected_setting": meta["selected_setting"],
                "smooth": meta["smooth"],
                "p25_rho_gs": meta["p25_rho_gs"],
            }
    return lc_to_meta


LC_TO_META = build_lc_meta()

DAYLIKE_VPPS = {"SOSD", "EOSD"}


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


def parse_site_lc_from_fitted(path):
    name = Path(path).name
    m = re.match(r"^Fitted_PPI_NBARPPI_GPP_data_(?P<site>[^_]+)_(?P<lc>LC\d+)\.csv$", name)
    if not m:
        raise ValueError(f"Cannot parse fitted filename: {name}")
    return m.group("site"), m.group("lc")


def collect_file_paths(nbar_ppi_dir, gpp_dir, fitted_dir):
    files = defaultdict(lambda: {"nbar_ppi": None, "gpp": None, "fitted": None})

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

    for f in Path(fitted_dir).glob("Fitted_PPI_NBARPPI_GPP_data_*.csv"):
        try:
            site, lc = parse_site_lc_from_fitted(f)
            files[(site, lc)]["fitted"] = str(f)
        except ValueError:
            continue

    for site, lc in list(files.keys()):
        if site in gpp_map:
            files[(site, lc)]["gpp"] = gpp_map[site]

    return files


def read_nbar_ppi_daily(path):
    df = pd.read_csv(path)
    df["t"] = pd.to_datetime(df["t"], errors="coerce").dt.date
    df = df.dropna(subset=["t"]).copy()
    df = df.rename(columns={"PPI": "NBAR_PPI"})[["t", "NBAR_PPI"]]
    df["NBAR_PPI"] = pd.to_numeric(df["NBAR_PPI"], errors="coerce").clip(-1, 5).clip(lower=0)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df


def read_gpp_daily(path):
    df = pd.read_csv(path)
    df["t"] = pd.to_datetime(df["t"], errors="coerce").dt.date
    df = df.dropna(subset=["t"]).copy()
    df["GPP"] = pd.to_numeric(df["GPP"], errors="coerce")
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df[["t", "GPP"]]


def make_selected_daily_context(raw_nbar, gpp, fitted_csv, selected_setting):
    context = pd.merge(raw_nbar, gpp, on="t", how="outer")

    if fitted_csv and Path(fitted_csv).exists():
        fitted_df = pd.read_csv(fitted_csv)
        if selected_setting in fitted_df.columns:
            fitted_df["t"] = pd.to_datetime(fitted_df["t"]).dt.date
            fitted_df = fitted_df[["t", selected_setting]].copy()
            context = pd.merge(context, fitted_df, on="t", how="outer")

    return context.sort_values("t")


def _fixed_st_vpp_run(indata, smooth_value):
    viname = [c for c in indata.columns if c != "t"][0]
    if viname != "NBAR_PPI":
        raise ValueError(f"Expected NBAR_PPI input, got {viname}")

    p_a = [[-10000.0, 10000.0, 1.0], [-10000.0, 10000.0, 1.0], [-10000.0, 10000.0, 1.0]]
    p_ignoreday = 366
    p_ylu = [0.0, 5.0]
    p_printflag = 0
    p_nodata = -9999
    p_davailwin = 45
    p_outlier = 0
    p_hrvppformat = 1
    p_nclasses = 1

    landuse = np.zeros(255, dtype="uint8")
    p_fitmethod = np.zeros(255, dtype="uint8")
    p_smooth_vec = np.zeros(255, dtype="double")
    p_nenvi_vec = np.ones(255, dtype="uint8")
    p_wfactnum_vec = np.ones(255, dtype="double")
    p_startmethod_vec = np.ones(255, dtype="uint8")
    p_startcutoff_vec = np.full((255, 2), 0.5, order="F", dtype="double")
    p_low_percentile_vec = np.full(255, 0.05, dtype="double")
    p_fillbase_vec = np.ones(255, dtype="uint8")
    p_seasonmethod_vec = np.zeros(255, dtype="uint8")
    p_seapar_vec = np.zeros(255, dtype="double")
    landuse[0] = 1

    timevector, tv_yyyymmdd, vi, qa, lc, yr, yrstart, yrend, npt, col = Step_2_ts_function.read_table_data(
        indata, p_a
    )

    p_outindex = np.arange(1, yr * 365 + 1)
    seasons = ["s1", "s2"]
    vppname = [
        "SOSD",
        "SOSV",
        "LSLOPE",
        "EOSD",
        "EOSV",
        "RSLOPE",
        "LENGTH",
        "MINV",
        "MAXD",
        "MAXV",
        "AMPL",
        "TPROD",
        "SPROD",
    ]
    vpplist = [
        f"{year}_{season}_{vpp}"
        for year in range(yrstart, yrend + 1)
        for season in seasons
        for vpp in vppname
    ]

    p_fitmethod.fill(2)
    p_smooth_vec.fill(float(smooth_value))

    seasonmethod_options = [1, 2]
    seapar_options = [0.0, 0.25, 0.5, 0.75, 1.0]
    cutoff_options = np.arange(0.05, 0.55, 0.05)

    settings_rows = []
    vpp_columns = []
    counter = 0

    for seasonmethod_value in seasonmethod_options:
        p_seasonmethod_vec.fill(seasonmethod_value)

        for seapar_value in seapar_options:
            p_seapar_vec.fill(seapar_value)

            for soscutoff_value in cutoff_options:
                p_startcutoff_vec[0, 0] = soscutoff_value

                for eoscutoff_value in cutoff_options:
                    p_startcutoff_vec[1, 0] = eoscutoff_value

                    counter += 1
                    settings_id = f"{viname}_VPP-{counter}"
                    settings_rows.append(
                        {
                            "settings_id": settings_id,
                            "fitmethod": 2,
                            "method_name": "SP",
                            "smooth": float(smooth_value),
                            "seasonmethod": int(seasonmethod_value),
                            "seapar": float(seapar_value),
                            "sos_cutoff": float(soscutoff_value),
                            "eos_cutoff": float(eoscutoff_value),
                        }
                    )

                    vpp, vppqa, nseason, yfit, yfitqa, seasonfit, tseq = timesat.tsfprocess(
                        yr,
                        vi,
                        qa,
                        timevector,
                        lc,
                        p_nclasses,
                        landuse,
                        p_outindex,
                        p_ignoreday,
                        p_ylu,
                        p_printflag,
                        p_fitmethod,
                        p_smooth_vec,
                        p_nodata,
                        p_davailwin,
                        p_outlier,
                        p_nenvi_vec,
                        p_wfactnum_vec,
                        p_startmethod_vec,
                        p_startcutoff_vec,
                        p_low_percentile_vec,
                        p_fillbase_vec,
                        p_hrvppformat,
                        p_seasonmethod_vec,
                        p_seapar_vec,
                        1,
                    )

                    vpp = np.squeeze(np.moveaxis(vpp, -1, 0), axis=-1)
                    vpp_columns.append((settings_id, np.asarray(vpp, dtype=np.float32)))

    settings_df = pd.DataFrame(settings_rows)
    vpp_data = pd.DataFrame({"id": vpplist})
    if vpp_columns:
        vpp_array = np.column_stack([col_vals for _, col_vals in vpp_columns])
        vpp_cols_names = [sid for sid, _ in vpp_columns]
        vpp_df = pd.DataFrame(vpp_array, columns=vpp_cols_names)
        vpp_data = pd.concat([vpp_data, vpp_df], axis=1)
    vpp_data.replace(-9999, np.nan, inplace=True)
    return vpp_data, settings_df


def add_metadata(df, site, lc, class_name, selected_setting, p25_rho_gs):
    df = df.copy()
    df["site"] = site
    df["lc"] = lc
    df["class"] = class_name
    df["selected_setting"] = selected_setting
    df["selected_setting_p25_rho_gs"] = p25_rho_gs
    return df


def parse_vpp_id(vpp_id):
    year, season, vpp = str(vpp_id).split("_", 2)
    return int(year), season, vpp


def to_seq_day(value):
    if pd.isna(value):
        return np.nan
    try:
        value = float(value)
    except Exception:
        return np.nan
    int_part = int(np.floor(value))
    frac_part = value - int_part
    s = str(int_part).zfill(7)
    year = int(s[:4])
    doy = int(s[-3:])
    return (year - 2017) * 365 + doy + frac_part


def convert_day_value(value, seasonmethod, vpp_name):
    if vpp_name not in DAYLIKE_VPPS:
        return np.nan
    if pd.isna(value):
        return np.nan
    value = float(value)
    if int(seasonmethod) == 1:
        if value > 10000:
            return to_seq_day(value)
        return value
    return value


def summarize_best_vpp_settings(
    gpp_vpp_long,
    gpp_settings_all,
    selected_vpp_long,
    selected_settings_all,
    out_vi_dir,
    max_sosd_abs_diff=60.0,
    max_eosd_abs_diff=60.0,
):
    gpp_meta = gpp_settings_all[["site", "lc", "settings_id", "seasonmethod"]].drop_duplicates()
    gpp_df = gpp_vpp_long.merge(gpp_meta, on=["site", "lc", "settings_id"], how="left")
    gpp_df[["year", "season", "vpp"]] = gpp_df["id"].apply(lambda s: pd.Series(parse_vpp_id(s)))
    gpp_df = gpp_df[gpp_df["vpp"].isin(DAYLIKE_VPPS)].copy()
    gpp_df["gpp_day"] = gpp_df.apply(
        lambda r: convert_day_value(r["value"], r["seasonmethod"], r["vpp"]), axis=1
    )
    gpp_ref = gpp_df[["site", "lc", "id", "vpp", "gpp_day"]].drop_duplicates()

    sel_meta_cols = [
        "site",
        "lc",
        "class",
        "selected_setting",
        "selected_setting_p25_rho_gs",
        "settings_id",
        "seasonmethod",
        "seapar",
        "sos_cutoff",
        "eos_cutoff",
        "smooth",
        "fitmethod",
        "method_name",
    ]
    sel_meta = selected_settings_all[sel_meta_cols].drop_duplicates()
    sel_df = selected_vpp_long.merge(sel_meta, on=["site", "lc", "settings_id"], how="left")
    sel_df[["year", "season", "vpp"]] = sel_df["id"].apply(lambda s: pd.Series(parse_vpp_id(s)))
    sel_df = sel_df[sel_df["vpp"].isin(DAYLIKE_VPPS)].copy()
    sel_df["selected_day"] = sel_df.apply(
        lambda r: convert_day_value(r["value"], r["seasonmethod"], r["vpp"]), axis=1
    )

    comp = sel_df.merge(gpp_ref, on=["site", "lc", "id", "vpp"], how="inner")
    if "class_x" in comp.columns and "class" not in comp.columns:
        comp = comp.rename(columns={"class_x": "class"})
    if "selected_setting_x" in comp.columns and "selected_setting" not in comp.columns:
        comp = comp.rename(columns={"selected_setting_x": "selected_setting"})
    if (
        "selected_setting_p25_rho_gs_x" in comp.columns
        and "selected_setting_p25_rho_gs" not in comp.columns
    ):
        comp = comp.rename(
            columns={"selected_setting_p25_rho_gs_x": "selected_setting_p25_rho_gs"}
        )
    comp["site_lc"] = comp["site"].astype(str) + "|" + comp["lc"].astype(str)
    comp["abs_diff_days"] = np.abs(comp["selected_day"] - comp["gpp_day"])
    comp["signed_diff_days"] = comp["selected_day"] - comp["gpp_day"]
    comp = comp[np.isfinite(comp["abs_diff_days"])].copy()

    out_csv = Path(out_vi_dir) / "Selected_ST_NBARPPI_VPP_vs_GPP_SOSD_EOSD_ALL_sites_long.csv"
    comp.to_csv(out_csv, index=False)
    print(f"[INFO] wrote {out_csv}")

    # Keep only seasons where both SOSD and EOSD can be reasonably matched to GPP.
    season_wide = (
        comp.pivot_table(
            index=["site", "lc", "class", "settings_id", "year", "season"],
            columns="vpp",
            values="signed_diff_days",
            aggfunc="first",
        )
        .reset_index()
    )
    for col in ["SOSD", "EOSD"]:
        if col not in season_wide.columns:
            season_wide[col] = np.nan

    season_wide["keep_match"] = (
        np.isfinite(season_wide["SOSD"])
        & np.isfinite(season_wide["EOSD"])
        & (np.abs(season_wide["SOSD"]) <= max_sosd_abs_diff)
        & (np.abs(season_wide["EOSD"]) <= max_eosd_abs_diff)
    )

    out_csv = Path(out_vi_dir) / "Selected_ST_NBARPPI_VPP_vs_GPP_SOSD_EOSD_match_filter.csv"
    season_wide.to_csv(out_csv, index=False)
    print(f"[INFO] wrote {out_csv}")

    keep_keys = season_wide.loc[
        season_wide["keep_match"],
        ["site", "lc", "class", "settings_id", "year", "season"],
    ].copy()
    keep_keys["__keep__"] = 1

    comp = comp.merge(
        keep_keys,
        on=["site", "lc", "class", "settings_id", "year", "season"],
        how="inner",
    )
    comp = comp.drop(columns="__keep__")

    out_csv = Path(out_vi_dir) / "Selected_ST_NBARPPI_VPP_vs_GPP_SOSD_EOSD_ALL_sites_long_filtered.csv"
    comp.to_csv(out_csv, index=False)
    print(f"[INFO] wrote {out_csv}")

    required_cols = [
        "class",
        "settings_id",
        "fitmethod",
        "method_name",
        "smooth",
        "seasonmethod",
        "seapar",
        "sos_cutoff",
        "eos_cutoff",
    ]
    missing_cols = [c for c in required_cols if c not in comp.columns]
    if missing_cols:
        raise KeyError(
            f"Missing required columns in Step3 comparison table: {missing_cols}. "
            f"Available columns: {list(comp.columns)}"
        )

    by_setting = (
        comp.groupby(
            required_cols
        )
        .agg(
            mean_abs_diff_days=("abs_diff_days", "mean"),
            median_abs_diff_days=("abs_diff_days", "median"),
            rmsd_days=("signed_diff_days", lambda s: float(np.sqrt(np.mean(np.square(s))))),
            bias_days=("signed_diff_days", "mean"),
            n_obs=("abs_diff_days", "count"),
            n_site_lc=("site_lc", "nunique"),
        )
        .reset_index()
    )

    best_rows = []
    for class_name, sub in by_setting.groupby("class"):
        sub = sub.sort_values(
            ["mean_abs_diff_days", "median_abs_diff_days", "rmsd_days", "n_obs"],
            ascending=[True, True, True, False],
        )
        best_rows.append(sub.iloc[0].to_dict())

    best_df = pd.DataFrame(best_rows).sort_values("class")

    out_csv = Path(out_vi_dir) / "Selected_ST_NBARPPI_VPP_SOSD_EOSD_summary_by_setting.csv"
    by_setting.to_csv(out_csv, index=False)
    print(f"[INFO] wrote {out_csv}")

    out_csv = Path(out_vi_dir) / "Selected_ST_NBARPPI_VPP_SOSD_EOSD_best_settings_by_class.csv"
    best_df.to_csv(out_csv, index=False)
    print(f"[INFO] wrote {out_csv}")


def main():
    p = argparse.ArgumentParser(
        description=(
            "Step3 VPP optimization using class-based selected ST_NBAR_PPI smoothing "
            "and full phenology parameter sweep."
        )
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
    p.add_argument(
        "--fitted_dir",
        type=str,
        default="output/Cal/PPI_GPP",
    )
    p.add_argument(
        "--out_vi_dir",
        type=str,
        default="output/Cal/Step3_ST_NBARPPI_VPP",
    )
    p.add_argument(
        "--max_sosd_abs_diff",
        type=float,
        default=60.0,
        help="Keep only matched seasons with |SOSD_selected - SOSD_GPP| <= this many days.",
    )
    p.add_argument(
        "--max_eosd_abs_diff",
        type=float,
        default=60.0,
        help="Keep only matched seasons with |EOSD_selected - EOSD_GPP| <= this many days.",
    )
    args = p.parse_args()

    out_vi_dir = Path(args.out_vi_dir)
    out_vi_dir.mkdir(parents=True, exist_ok=True)

    files = collect_file_paths(args.nbar_ppi_dir, args.gpp_dir, args.fitted_dir)
    print(f"Found {len(files)} site/LC groups")

    lookup_rows = []
    all_gpp_settings = []
    all_selected_settings = []
    all_gpp_vpp_long = []
    all_selected_vpp_long = []

    for (site, lc), info in sorted(files.items()):
        if lc not in LC_TO_META:
            print(f"[WARN] {site},{lc}: unsupported LC class, skip")
            continue
        if not info["nbar_ppi"]:
            print(f"[WARN] {site},{lc}: missing NBAR_PPI")
            continue
        if not info["gpp"]:
            print(f"[WARN] {site},{lc}: missing GPP")
            continue

        meta = LC_TO_META[lc]
        class_name = meta["class"]
        selected_setting = meta["selected_setting"]
        smooth_value = meta["smooth"]
        p25_rho_gs = meta["p25_rho_gs"]

        print(
            f"[INFO] site={site}, lc={lc}, class={class_name}, "
            f"selected_setting={selected_setting}, smooth={smooth_value}"
        )

        raw_nbar = read_nbar_ppi_daily(info["nbar_ppi"])
        gpp = read_gpp_daily(info["gpp"])

        context_df = make_selected_daily_context(raw_nbar, gpp, info["fitted"], selected_setting)
        out_csv = out_vi_dir / f"Selected_ST_NBARPPI_GPP_data_{site}_{lc}.csv"
        context_df.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")

        gpp_slice = gpp[["t", "GPP"]].copy()
        _, gpp_vpp_data, gpp_settings_df = Step_2_ts_function._ts_run_(gpp_slice, 1)
        if isinstance(gpp_vpp_data, int) or gpp_vpp_data is None or len(gpp_vpp_data) == 0:
            print(f"[WARN] {site},{lc}: no GPP VPP output")
            continue

        gpp_settings_df = add_metadata(
            gpp_settings_df, site, lc, class_name, selected_setting, p25_rho_gs
        )
        out_csv = out_vi_dir / f"GPP_reference_VPP_{site}_{lc}.csv"
        gpp_vpp_data.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")

        out_csv = out_vi_dir / f"GPP_reference_VPP_settings_{site}_{lc}.csv"
        gpp_settings_df.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")

        selected_vpp_data, selected_settings_df = _fixed_st_vpp_run(
            raw_nbar.rename(columns={"NBAR_PPI": "NBAR_PPI"})[["t", "NBAR_PPI"]].copy(),
            smooth_value=smooth_value,
        )

        selected_settings_df["source_setting"] = selected_setting
        selected_settings_df = add_metadata(
            selected_settings_df, site, lc, class_name, selected_setting, p25_rho_gs
        )

        out_csv = out_vi_dir / f"Selected_ST_NBARPPI_VPP_{site}_{lc}.csv"
        selected_vpp_data.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")

        out_csv = out_vi_dir / f"Selected_ST_NBARPPI_VPP_settings_{site}_{lc}.csv"
        selected_settings_df.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")

        lookup_rows.append(
            {
                "site": site,
                "lc": lc,
                "class": class_name,
                "selected_setting": selected_setting,
                "fitmethod": 2,
                "method_name": "SP",
                "smooth": smooth_value,
                "selected_setting_p25_rho_gs": p25_rho_gs,
            }
        )

        all_gpp_settings.append(gpp_settings_df)
        all_selected_settings.append(selected_settings_df)

        gpp_long = gpp_vpp_data.melt(id_vars="id", var_name="settings_id", value_name="value")
        gpp_long = add_metadata(gpp_long, site, lc, class_name, selected_setting, p25_rho_gs)
        all_gpp_vpp_long.append(gpp_long)

        selected_long = selected_vpp_data.melt(id_vars="id", var_name="settings_id", value_name="value")
        selected_long = add_metadata(selected_long, site, lc, class_name, selected_setting, p25_rho_gs)
        all_selected_vpp_long.append(selected_long)

    if lookup_rows:
        out_csv = out_vi_dir / "Selected_ST_NBARPPI_setting_lookup.csv"
        pd.DataFrame(lookup_rows).to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")

    if all_gpp_settings:
        out_csv = out_vi_dir / "GPP_reference_VPP_settings_ALL_sites.csv"
        gpp_settings_all = pd.concat(all_gpp_settings, ignore_index=True)
        gpp_settings_all.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")
    else:
        gpp_settings_all = None

    if all_selected_settings:
        out_csv = out_vi_dir / "Selected_ST_NBARPPI_VPP_settings_ALL_sites.csv"
        selected_settings_all = pd.concat(all_selected_settings, ignore_index=True)
        selected_settings_all.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")
    else:
        selected_settings_all = None

    if all_gpp_vpp_long:
        out_csv = out_vi_dir / "GPP_reference_VPP_ALL_sites_long.csv"
        gpp_vpp_long = pd.concat(all_gpp_vpp_long, ignore_index=True)
        gpp_vpp_long.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")
    else:
        gpp_vpp_long = None

    if all_selected_vpp_long:
        out_csv = out_vi_dir / "Selected_ST_NBARPPI_VPP_ALL_sites_long.csv"
        selected_vpp_long = pd.concat(all_selected_vpp_long, ignore_index=True)
        selected_vpp_long.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")
    else:
        selected_vpp_long = None

    if (
        gpp_vpp_long is not None
        and gpp_settings_all is not None
        and selected_vpp_long is not None
        and selected_settings_all is not None
    ):
        summarize_best_vpp_settings(
            gpp_vpp_long,
            gpp_settings_all,
            selected_vpp_long,
            selected_settings_all,
            out_vi_dir,
            max_sosd_abs_diff=args.max_sosd_abs_diff,
            max_eosd_abs_diff=args.max_eosd_abs_diff,
        )


if __name__ == "__main__":
    main()
