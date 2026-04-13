#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import pandas as pd


CLASS_GROUPS = {
    "Forest": {
        "landcover_codes": [7, 8, 16],
    },
    "Grassland": {
        "landcover_codes": [2, 5, 9, 10, 14, 15],
    },
    "Cropland": {
        "landcover_codes": [11, 12, 13],
    },
}


def build_lc_to_class():
    mapping = {}
    for class_name, meta in CLASS_GROUPS.items():
        for code in meta["landcover_codes"]:
            mapping[f"LC{code}"] = class_name
    return mapping


LC_TO_CLASS = build_lc_to_class()

SELECTED_NBAR_SETTINGS = {
    "Cropland": "ST_NBAR_PPI-12",
    "Grassland": "ST_NBAR_PPI-14",
    "Forest": "ST_NBAR_PPI-15",
}


def _prepare_diag_long(diag_long_csv):
    diag = pd.read_csv(diag_long_csv)
    diag["class"] = diag["lc"].map(LC_TO_CLASS)
    diag = diag.dropna(subset=["class"]).copy()

    families = {
        "PPI": "GPP_vs_ST_PPI",
        "NBAR_PPI": "GPP_vs_ST_NBAR_PPI",
    }

    diag["vi_family"] = pd.NA
    for family, label in families.items():
        diag.loc[diag["VI"] == label, "vi_family"] = family

    diag = diag.dropna(subset=["vi_family", "x_col", "rho_gs"]).copy()
    return diag


def _summarize_family(diag, settings, vi_family, out_dir, top_n):
    sub = diag.loc[diag["vi_family"] == vi_family].copy()
    if sub.empty:
        raise ValueError(f"No rows found for {vi_family} in diagnostics table.")

    best_rows = []
    top_rows = []

    for class_name, class_df in sub.groupby("class"):
        summary = (
            class_df.groupby("x_col")
            .agg(
                mean_rho_gs=("rho_gs", "mean"),
                median_rho_gs=("rho_gs", "median"),
                p25_rho_gs=("rho_gs", lambda s: s.quantile(0.25)),
                valid_sites=("site", "nunique"),
                valid_site_lc=("site", "count"),
            )
            .sort_values(["mean_rho_gs", "median_rho_gs", "valid_sites"], ascending=[False, False, False])
            .reset_index()
            .rename(columns={"x_col": "settings_id"})
        )

        if summary.empty:
            continue

        class_settings = settings.loc[settings["settings_id"].isin(summary["settings_id"])].copy()
        merged = summary.merge(class_settings, on="settings_id", how="left")

        best = merged.iloc[0].to_dict()
        best_rows.append(
            {
                "class": class_name,
                "vi_family": vi_family,
                "selection_metric": "highest mean rho_gs across site-lc pairs in class",
                **best,
            }
        )

        for rank, (_, row) in enumerate(merged.head(top_n).iterrows(), start=1):
            top_rows.append(
                {
                    "class": class_name,
                    "vi_family": vi_family,
                    "rank": rank,
                    "selection_metric": "highest mean rho_gs across site-lc pairs in class",
                    **row.to_dict(),
                }
            )

    best_df = pd.DataFrame(best_rows).sort_values(["class", "vi_family"])
    top_df = pd.DataFrame(top_rows).sort_values(["class", "vi_family", "rank"])

    best_out = out_dir / f"ST_{vi_family}_best_settings_by_class.csv"
    top_out = out_dir / f"ST_{vi_family}_top{top_n}_settings_by_class.csv"

    best_df.to_csv(best_out, index=False)
    top_df.to_csv(top_out, index=False)

    print(f"[INFO] wrote {best_out}")
    print(f"[INFO] wrote {top_out}")

    return best_df, top_df


def summarize_selected_nbar_settings(diag, settings, out_dir):
    sub = diag.loc[diag["vi_family"] == "NBAR_PPI"].copy()
    rows = []

    for class_name, settings_id in SELECTED_NBAR_SETTINGS.items():
        class_df = sub[(sub["class"] == class_name) & (sub["x_col"] == settings_id)].copy()
        if class_df.empty:
            continue

        settings_row = settings.loc[settings["settings_id"] == settings_id]
        if settings_row.empty:
            raise ValueError(f"Missing settings row for {settings_id}")

        row = {
            "class": class_name,
            "vi_family": "NBAR_PPI",
            "selection_type": "preselected operational setting",
            "settings_id": settings_id,
            "mean_rho_gs": class_df["rho_gs"].mean(),
            "median_rho_gs": class_df["rho_gs"].median(),
            "p25_rho_gs": class_df["rho_gs"].quantile(0.25),
            "valid_sites": class_df["site"].nunique(),
            "valid_site_lc": len(class_df),
        }
        row.update(settings_row.iloc[0].to_dict())
        rows.append(row)

    selected_df = pd.DataFrame(rows).sort_values("class")
    selected_out = out_dir / "ST_NBAR_PPI_selected_settings_by_class.csv"
    selected_df.to_csv(selected_out, index=False)
    print(f"[INFO] wrote {selected_out}")
    return selected_df


def build_best_settings_by_class(diag_long_csv, settings_csv, out_dir, top_n=3):
    diag_long_csv = Path(diag_long_csv)
    settings_csv = Path(settings_csv)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    diag = _prepare_diag_long(diag_long_csv)
    settings = pd.read_csv(settings_csv)

    wanted_prefixes = ("ST_PPI-", "ST_NBAR_PPI-")
    settings = settings.loc[
        settings["settings_id"].astype(str).str.startswith(wanted_prefixes)
    ].copy()

    if settings.empty:
        raise ValueError(f"No ST_PPI / ST_NBAR_PPI settings found in {settings_csv}")

    best_frames = []
    top_frames = []
    for vi_family in ["PPI", "NBAR_PPI"]:
        best_df, top_df = _summarize_family(diag, settings, vi_family, out_dir, top_n)
        best_frames.append(best_df)
        top_frames.append(top_df)

    summarize_selected_nbar_settings(diag, settings, out_dir)

    combined_best = pd.concat(best_frames, ignore_index=True)
    combined_top = pd.concat(top_frames, ignore_index=True)

    combined_best_out = out_dir / "ST_PPI_NBARPPI_best_settings_by_class.csv"
    combined_top_out = out_dir / f"ST_PPI_NBARPPI_top{top_n}_settings_by_class.csv"

    combined_best.to_csv(combined_best_out, index=False)
    combined_top.to_csv(combined_top_out, index=False)

    print(f"[INFO] wrote {combined_best_out}")
    print(f"[INFO] wrote {combined_top_out}")
    print("\n[Best settings by class]")
    print(combined_best.to_string(index=False))


def main():
    p = argparse.ArgumentParser(
        description=(
            "Summarize best ST_PPI and ST_NBAR_PPI smoothing settings for Forest, "
            "Grassland, and Cropland using mean rho_gs against GPP."
        )
    )
    p.add_argument(
        "--diag_long_csv",
        default="output/Cal/PPI_GPP/ST_PPI_NBARPPI_GPP_diag_ALL_sites_long.csv",
        help="Long diagnostics table produced by Cal_Step2_ST_optimization.py",
    )
    p.add_argument(
        "--settings_csv",
        default="output/Cal/PPI_GPP/Fitted_PPI_NBARPPI_GPP_settings_AT-Neu_LC10.csv",
        help="Any per-site settings CSV; settings_id to parameter mapping is shared across sites.",
    )
    p.add_argument(
        "--out_dir",
        default="output/Cal/PPI_GPP",
        help="Directory for summary outputs.",
    )
    p.add_argument(
        "--top_n",
        type=int,
        default=3,
        help="How many top settings to keep for each class.",
    )
    args = p.parse_args()

    build_best_settings_by_class(
        diag_long_csv=args.diag_long_csv,
        settings_csv=args.settings_csv,
        out_dir=args.out_dir,
        top_n=args.top_n,
    )


if __name__ == "__main__":
    main()
