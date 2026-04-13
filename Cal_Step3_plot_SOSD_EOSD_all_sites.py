#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DAYLIKE_VPPS = {"SOSD", "EOSD"}


def parse_site_lc_from_data(path):
    name = Path(path).name
    m = re.match(r"^Selected_ST_NBARPPI_GPP_data_(?P<site>[^_]+)_(?P<lc>LC\d+)\.csv$", name)
    if not m:
        raise ValueError(f"Cannot parse site/lc from {name}")
    return m.group("site"), m.group("lc")


def parse_vpp_id(vpp_id):
    year, season, vpp = str(vpp_id).split("_", 2)
    return int(year), season, vpp


def to_seq_day_from_yyyydoy(value):
    if pd.isna(value):
        return np.nan
    value = float(value)
    int_part = int(np.floor(value))
    frac_part = value - int_part
    s = str(int_part).zfill(7)
    year = int(s[:4])
    doy = int(s[-3:])
    return (year - 2017) * 365 + doy + frac_part


def seq_day_to_datetime(value):
    if pd.isna(value):
        return pd.NaT
    value = float(value)
    base = datetime(2017, 1, 1)
    return base + timedelta(days=value - 1)


def convert_to_datetime(value, seasonmethod, vpp_name):
    if vpp_name not in DAYLIKE_VPPS or pd.isna(value):
        return pd.NaT
    value = float(value)
    if int(seasonmethod) == 1 and value > 10000:
        value = to_seq_day_from_yyyydoy(value)
    return seq_day_to_datetime(value)


def extract_events(vpp_csv, settings_csv, chosen_settings_id):
    vpp_df = pd.read_csv(vpp_csv)
    settings_df = pd.read_csv(settings_csv)

    settings_row = settings_df.loc[settings_df["settings_id"] == chosen_settings_id]
    if settings_row.empty:
        raise ValueError(f"Missing settings row {chosen_settings_id} in {settings_csv}")
    seasonmethod = int(settings_row.iloc[0]["seasonmethod"])

    sub = vpp_df[["id", chosen_settings_id]].copy()
    sub[["year", "season", "vpp"]] = sub["id"].apply(lambda s: pd.Series(parse_vpp_id(s)))
    sub = sub[sub["vpp"].isin(DAYLIKE_VPPS)].copy()
    sub["event_time"] = sub.apply(
        lambda r: convert_to_datetime(r[chosen_settings_id], seasonmethod, r["vpp"]),
        axis=1,
    )
    return sub[["id", "year", "season", "vpp", "event_time"]]


def pick_selected_settings_id(best_df, lc):
    lc_code = int(str(lc).replace("LC", ""))
    if lc_code in {11, 12, 13}:
        cls = "Cropland"
    elif lc_code in {7, 8, 16}:
        cls = "Forest"
    elif lc_code in {2, 5, 9, 10, 14, 15}:
        cls = "Grassland"
    else:
        raise ValueError(f"Unsupported LC {lc}")

    row = best_df.loc[best_df["class"] == cls]
    if row.empty:
        raise ValueError(f"Missing best setting for class {cls}")
    return cls, row.iloc[0]["settings_id"]


def plot_one_site(data_csv, gpp_vpp_csv, gpp_settings_csv, selected_vpp_csv, selected_settings_csv, best_df, out_dir):
    site, lc = parse_site_lc_from_data(data_csv)
    cls, selected_settings_id = pick_selected_settings_id(best_df, lc)

    df = pd.read_csv(data_csv)
    df["t"] = pd.to_datetime(df["t"])

    st_col = [c for c in df.columns if c.startswith("ST_NBAR_PPI-")]
    if len(st_col) != 1:
        raise ValueError(f"Expected exactly one ST_NBAR_PPI column in {data_csv}, got {st_col}")
    st_col = st_col[0]

    gpp_settings_df = pd.read_csv(gpp_settings_csv)
    gpp_settings_id = gpp_settings_df.iloc[0]["settings_id"]

    gpp_events = extract_events(gpp_vpp_csv, gpp_settings_csv, gpp_settings_id)
    selected_events = extract_events(selected_vpp_csv, selected_settings_csv, selected_settings_id)

    fig, ax_vi = plt.subplots(figsize=(14, 5.2))
    ax_gpp = ax_vi.twinx()

    gpp_mask = np.isfinite(df["GPP"])
    raw_mask = np.isfinite(df["NBAR_PPI"])
    st_mask = np.isfinite(df[st_col])

    ax_gpp.plot(
        df.loc[gpp_mask, "t"],
        df.loc[gpp_mask, "GPP"],
        color="black",
        linewidth=1.0,
        alpha=0.85,
        label="GPP",
        zorder=1,
    )
    ax_vi.scatter(
        df.loc[raw_mask, "t"],
        df.loc[raw_mask, "NBAR_PPI"],
        s=10,
        color="0.8",
        alpha=0.85,
        label="Raw NBAR_PPI",
        zorder=3,
    )
    ax_vi.plot(
        df.loc[st_mask, "t"],
        df.loc[st_mask, st_col],
        color="#d62728",
        linewidth=1.5,
        label=st_col,
        zorder=4,
    )

    for _, row in gpp_events.iterrows():
        if pd.isna(row["event_time"]):
            continue
        color = "black"
        linestyle = "--" if row["vpp"] == "SOSD" else ":"
        ax_vi.axvline(row["event_time"], color=color, linestyle=linestyle, linewidth=0.9, alpha=0.75, zorder=2)

    for _, row in selected_events.iterrows():
        if pd.isna(row["event_time"]):
            continue
        color = "#d62728"
        linestyle = "--" if row["vpp"] == "SOSD" else ":"
        ax_vi.axvline(row["event_time"], color=color, linestyle=linestyle, linewidth=0.9, alpha=0.85, zorder=5)

    ax_vi.set_ylabel("NBAR_PPI / ST_NBAR_PPI")
    ax_gpp.set_ylabel("GPP", color="black")
    ax_gpp.tick_params(axis="y", labelcolor="black")
    ax_vi.set_xlabel("Date")
    ax_vi.grid(True, axis="y", color="0.92", linewidth=0.8)
    ax_vi.set_title(f"{site} ({lc}, {cls}) | selected phenology setting: {selected_settings_id}")

    ax_vi.text(
        0.99,
        0.98,
        "black dashed/dotted: GPP SOSD/EOSD\nred dashed/dotted: NBAR_PPI SOSD/EOSD",
        transform=ax_vi.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color="0.25",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.85", alpha=0.9),
    )

    ax_vi.set_zorder(ax_gpp.get_zorder() + 1)
    ax_vi.patch.set_visible(False)

    handles_vi, labels_vi = ax_vi.get_legend_handles_labels()
    handles_gpp, labels_gpp = ax_gpp.get_legend_handles_labels()
    ax_vi.legend(handles_gpp + handles_vi, labels_gpp + labels_vi, loc="upper left", frameon=False, ncol=3)

    fig.autofmt_xdate()
    fig.tight_layout()

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"Step3_SOSD_EOSD_check_{site}_{lc}.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] wrote {out_png}")


def main():
    p = argparse.ArgumentParser(
        description="Plot all sites to inspect whether Step3 SOSD/EOSD seasons are defined correctly."
    )
    p.add_argument("--step3_dir", type=str, default="output/Cal/Step3_ST_NBARPPI_VPP")
    p.add_argument(
        "--best_settings_csv",
        type=str,
        default="output/Cal/Step3_ST_NBARPPI_VPP/Selected_ST_NBARPPI_VPP_SOSD_EOSD_best_settings_by_class.csv",
    )
    p.add_argument(
        "--out_plot_dir",
        type=str,
        default="output/Cal/Step3_ST_NBARPPI_VPP/plots_SOSD_EOSD_check",
    )
    args = p.parse_args()

    step3_dir = Path(args.step3_dir)
    best_df = pd.read_csv(args.best_settings_csv)

    data_files = sorted(step3_dir.glob("Selected_ST_NBARPPI_GPP_data_*.csv"))
    print(f"Found {len(data_files)} site files")

    for data_csv in data_files:
        try:
            site, lc = parse_site_lc_from_data(data_csv)
            plot_one_site(
                data_csv,
                step3_dir / f"GPP_reference_VPP_{site}_{lc}.csv",
                step3_dir / f"GPP_reference_VPP_settings_{site}_{lc}.csv",
                step3_dir / f"Selected_ST_NBARPPI_VPP_{site}_{lc}.csv",
                step3_dir / f"Selected_ST_NBARPPI_VPP_settings_{site}_{lc}.csv",
                best_df,
                args.out_plot_dir,
            )
        except Exception as e:
            print(f"[WARN] failed for {data_csv.name}: {e}")


if __name__ == "__main__":
    main()
