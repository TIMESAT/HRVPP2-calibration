#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


CLASS_GROUPS = {
    "Forest": {
        "landcover_codes": [7, 8, 16],
        "setting": "ST_NBAR_PPI-15",
    },
    "Grassland": {
        "landcover_codes": [2, 5, 9, 10, 14, 15],
        "setting": "ST_NBAR_PPI-14",
    },
    "Cropland": {
        "landcover_codes": [11, 12, 13],
        "setting": "ST_NBAR_PPI-12",
    },
}


def build_lc_to_class_and_setting():
    lc_to_class = {}
    lc_to_setting = {}
    for class_name, meta in CLASS_GROUPS.items():
        for code in meta["landcover_codes"]:
            lc = f"LC{code}"
            lc_to_class[lc] = class_name
            lc_to_setting[lc] = meta["setting"]
    return lc_to_class, lc_to_setting


LC_TO_CLASS, LC_TO_SETTING = build_lc_to_class_and_setting()


def parse_site_lc_from_fitted(path):
    name = Path(path).name
    m = re.match(r"^Fitted_PPI_NBARPPI_GPP_data_(?P<site>[^_]+)_(?P<lc>LC\d+)\.csv$", name)
    if not m:
        raise ValueError(f"Cannot parse fitted filename: {name}")
    return m.group("site"), m.group("lc")


def raw_nbar_path(raw_nbar_dir, site, lc):
    return Path(raw_nbar_dir) / f"S2_NBAR_SZA11_{site}_{lc}_NDVI_NIRv_EVI2_mean.csv"


def plot_one_site(site, lc, raw_nbar_csv, fitted_csv, out_plot_dir):
    class_name = LC_TO_CLASS.get(lc)
    setting = LC_TO_SETTING.get(lc)
    if class_name is None or setting is None:
        raise ValueError(f"Unsupported land-cover class for {site}, {lc}")

    raw_df = pd.read_csv(raw_nbar_csv)
    fit_df = pd.read_csv(fitted_csv)

    if "PPI" not in raw_df.columns:
        raise ValueError(f"Missing raw NBAR PPI column in {raw_nbar_csv}")
    if setting not in fit_df.columns:
        raise ValueError(f"Missing selected ST column {setting} in {fitted_csv}")
    if "GPP" not in fit_df.columns:
        raise ValueError(f"Missing GPP column in {fitted_csv}")

    raw_df = raw_df.rename(columns={"PPI": "NBAR_PPI"})[["t", "NBAR_PPI"]].copy()
    raw_df["t"] = pd.to_datetime(raw_df["t"])
    raw_df["NBAR_PPI"] = pd.to_numeric(raw_df["NBAR_PPI"], errors="coerce")

    fit_df = fit_df[["t", setting, "GPP"]].copy()
    fit_df["t"] = pd.to_datetime(fit_df["t"])
    fit_df[setting] = pd.to_numeric(fit_df[setting], errors="coerce")
    fit_df["GPP"] = pd.to_numeric(fit_df["GPP"], errors="coerce")

    merged = pd.merge(raw_df, fit_df, on="t", how="outer").sort_values("t")

    dates = merged["t"]
    raw_mask = np.isfinite(merged["NBAR_PPI"])
    st_mask = np.isfinite(merged[setting])
    gpp_mask = np.isfinite(merged["GPP"])

    rho_df = merged[[setting, "GPP"]].dropna()
    if len(rho_df) >= 5:
        rho_gs = spearmanr(rho_df[setting], rho_df["GPP"])[0]
    else:
        rho_gs = np.nan

    fig, ax_vi = plt.subplots(figsize=(14, 4.8))
    ax_gpp = ax_vi.twinx()

    ax_gpp.plot(
        dates[gpp_mask],
        merged.loc[gpp_mask, "GPP"],
        color="black",
        linewidth=1.0,
        alpha=0.85,
        label="GPP",
        zorder=1,
    )

    ax_vi.scatter(
        dates[raw_mask],
        merged.loc[raw_mask, "NBAR_PPI"],
        s=10,
        color="0.8",
        alpha=0.9,
        label="Raw NBAR_PPI",
        zorder=3,
    )

    ax_vi.plot(
        dates[st_mask],
        merged.loc[st_mask, setting],
        color="#d62728",
        linewidth=1.5,
        label=setting,
        zorder=4,
    )

    ax_vi.set_ylabel("NBAR_PPI / ST_NBAR_PPI", color="#444444")
    ax_gpp.set_ylabel("GPP", color="black")
    ax_vi.set_xlabel("Date")
    ax_gpp.tick_params(axis="y", labelcolor="black")
    ax_vi.grid(True, axis="y", color="0.92", linewidth=0.8)
    title = f"{site} ({lc}, {class_name}) - {setting}"
    if np.isfinite(rho_gs):
        title += f" | rho_gs={rho_gs:.3f}"
    ax_vi.set_title(title)

    ax_vi.set_zorder(ax_gpp.get_zorder() + 1)
    ax_vi.patch.set_visible(False)

    handles_vi, labels_vi = ax_vi.get_legend_handles_labels()
    handles_gpp, labels_gpp = ax_gpp.get_legend_handles_labels()
    ax_vi.legend(
        handles_gpp + handles_vi,
        labels_gpp + labels_vi,
        loc="upper left",
        frameon=False,
        ncol=3,
    )

    if np.isfinite(rho_gs):
        ax_vi.text(
            0.99,
            0.96,
            f"rho_gs = {rho_gs:.3f}",
            transform=ax_vi.transAxes,
            ha="right",
            va="top",
            fontsize=10,
            color="#d62728",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.85", alpha=0.9),
        )

    fig.autofmt_xdate()
    fig.tight_layout()

    out_plot_dir = Path(out_plot_dir)
    out_plot_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_plot_dir / f"ST_NBARPPI_GPP_timeseries_{site}_{lc}.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] wrote {out_png}")


def main():
    p = argparse.ArgumentParser(
        description=(
            "Plot GPP (black), raw NBAR PPI (light-grey points), and selected "
            "ST NBAR PPI (red) for each site."
        )
    )
    p.add_argument(
        "--raw_nbar_dir",
        type=str,
        default="Data/VI/DVI_NIRv_EVI2_nbar_sza11/flux_cal/V1/csv_lc",
    )
    p.add_argument(
        "--fitted_dir",
        type=str,
        default="output/Cal/PPI_GPP",
    )
    p.add_argument(
        "--out_plot_dir",
        type=str,
        default="output/Cal/PPI_GPP/ST_NBARPPI_selected_timeseries",
    )
    args = p.parse_args()

    fitted_files = sorted(Path(args.fitted_dir).glob("Fitted_PPI_NBARPPI_GPP_data_*.csv"))
    print(f"Found {len(fitted_files)} fitted site files")

    for fitted_csv in fitted_files:
        try:
            site, lc = parse_site_lc_from_fitted(fitted_csv)
            raw_csv = raw_nbar_path(args.raw_nbar_dir, site, lc)
            if not raw_csv.exists():
                print(f"[WARN] missing raw NBAR PPI file for {site}, {lc}: {raw_csv}")
                continue
            plot_one_site(site, lc, raw_csv, fitted_csv, args.out_plot_dir)
        except Exception as e:
            print(f"[WARN] failed for {fitted_csv.name}: {e}")


if __name__ == "__main__":
    main()
