#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import pandas as pd


LC_TO_CLASS = {
    "LC7": "Forest",
    "LC8": "Forest",
    "LC16": "Forest",
    "LC9": "Grassland",
    "LC10": "Grassland",
    "LC14": "Grassland",
    "LC15": "Grassland",
    "LC11": "Cropland",
    "LC12": "Cropland",
    "LC13": "Cropland",
}


def build_best_settings_by_class(diag_csv, settings_csv, out_dir, top_n=3):
    diag_csv = Path(diag_csv)
    settings_csv = Path(settings_csv)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    diag = pd.read_csv(diag_csv)
    settings = pd.read_csv(settings_csv)

    diag["class"] = diag["lc"].map(LC_TO_CLASS)
    diag = diag.dropna(subset=["class"]).copy()

    setting_cols = [c for c in diag.columns if c.startswith("PPI-")]
    if not setting_cols:
        raise ValueError(f"No PPI setting columns found in {diag_csv}")

    best_rows = []
    top_rows = []

    for cls, sub in diag.groupby("class"):
        means = sub[setting_cols].mean(axis=0, skipna=True).sort_values(ascending=False)
        medians = sub[setting_cols].median(axis=0, skipna=True)
        counts = sub[setting_cols].notna().sum(axis=0)

        best_id = means.index[0]
        best_setting = settings.loc[settings["settings_id"] == best_id]
        if best_setting.empty:
            raise ValueError(f"Missing settings row for {best_id} in {settings_csv}")

        best_row = {
            "class": cls,
            "n_sites": sub["site"].nunique(),
            "best_setting_id": best_id,
            "selection_metric": "highest mean rho_gs across sites in class",
            "mean_rho_gs": means[best_id],
            "median_rho_gs": medians[best_id],
            "valid_sites_for_best": int(counts[best_id]),
        }
        best_row.update(best_setting.iloc[0].to_dict())
        best_rows.append(best_row)

        for rank, (settings_id, mean_val) in enumerate(means.head(top_n).items(), start=1):
            settings_row = settings.loc[settings["settings_id"] == settings_id]
            if settings_row.empty:
                raise ValueError(f"Missing settings row for {settings_id} in {settings_csv}")

            top_row = {
                "class": cls,
                "rank": rank,
                "selection_metric": "highest mean rho_gs across sites in class",
                "settings_id": settings_id,
                "mean_rho_gs": mean_val,
                "median_rho_gs": medians[settings_id],
                "valid_sites": int(counts[settings_id]),
            }
            top_row.update(settings_row.iloc[0].to_dict())
            top_rows.append(top_row)

    best_df = pd.DataFrame(best_rows).sort_values("class")
    top_df = pd.DataFrame(top_rows).sort_values(["class", "rank"])

    best_out = out_dir / "ST_PPI_GPP_best_settings_by_class.csv"
    top_out = out_dir / f"ST_PPI_GPP_top{top_n}_settings_by_class.csv"

    best_df.to_csv(best_out, index=False)
    top_df.to_csv(top_out, index=False)

    print(f"[INFO] wrote {best_out}")
    print(f"[INFO] wrote {top_out}")
    print("\n[Best settings by class]")
    print(best_df.to_string(index=False))


def main():
    p = argparse.ArgumentParser(
        description=(
            "Summarize best PPI smoothing settings for Forest, Grassland, and Cropland "
            "using the mean rho_gs across sites in each class."
        )
    )
    p.add_argument(
        "--diag_csv",
        default="output/Cal/PPI_GPP/ST_PPI_GPP_diag_ALL_sites.csv",
        help="Wide diagnostics table produced by Cal_Step2_ST_optimization.py",
    )
    p.add_argument(
        "--settings_csv",
        default="output/Cal/PPI_GPP/Fitted_PPI_GPP_settings_AT-Neu_LC10.csv",
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
        diag_csv=args.diag_csv,
        settings_csv=args.settings_csv,
        out_dir=args.out_dir,
        top_n=args.top_n,
    )


if __name__ == "__main__":
    main()
