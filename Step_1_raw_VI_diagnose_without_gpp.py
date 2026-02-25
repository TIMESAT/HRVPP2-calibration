#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge VI time series (NDVI, EVI2, PPI, LAI, FAPAR) for all sites and
generate land-cover specific scatterplot matrices with Spearman correlations.

Inputs (defaults from your description):
1) PPI CSVs:
   /home/hongxiao/fs7/proj/HRVPP2/calval/VI/PPI/flux_cal/V1/csv_lc/
   e.g. PPI_2017_01_01_2024_12_31_AT-Neu_LC10_mean.csv
      t,PPI_mean

2) NDVI + EVI2 CSVs:
   /home/hongxiao/fs7/proj/HRVPP2/calval/VI/NDVI_EVI2/flux_cal/V1/csv_lc/
   e.g. S2_2017_01_01_2024_12_31_AT-Neu_LC10_NDVI_EVI2_mean.csv
      t,NDVI,EVI2

3) FAPAR CSVs:
   /home/hongxiao/fs7/proj/HRVPP2/calval/VI/FAPAR/V5/timeseries/csv_lc_200m/
   e.g. FAPAR_2017_01_01_2024_12_31_AT-Neu_LC10_mean.csv
      t,FAPAR_mean

4) LAI CSVs:
   /home/hongxiao/fs7/proj/HRVPP2/calval/VI/LAI/flux_cal/V1/csv_lc/
   e.g. LAI_2017_01_01_2024_12_31_AT-Neu_LC10_mean.csv
      t,LAI_mean

Outputs:
- Per-site+LC merged VI CSVs:
    <out_vi_dir>/VIs_2017_01_01_2024_12_31_<SITE>_<LC>_mean.csv
    columns: t, NDVI, EVI2, PPI, LAI, FAPAR

- Per-LC scatterplot matrix PNGs:
    <out_plot_dir>/scatter_matrix_<LC>.png

- Per-LC Spearman correlation matrices:
    <out_plot_dir>/spearman_corr_<LC>.csv

Usage example:
    python vi_scatter_matrix_by_lc.py \
        --out_vi_dir /home/hongxiao/fs7/proj/HRVPP2/calval/VI/PNG \
        --out_plot_dir /home/hongxiao/fs7/proj/HRVPP2/calval/PNG

"""


import re
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

LC_PALETTE = {
    0:  "Other (not defined)",
    1:  "Water",
    2:  "Wetland",
    3:  "Non-vegetated",
    4:  "Snow & Ice",
    5:  "Lichens & Mosses",
    6:  "Sealed",
    7:  "Broadleaf Forest",
    8:  "Coniferous Forest",
    9:  "Periodically Grassland",
    10: "Permanent Grassland",
    11: "Annual Cropland",
    12: "Permanent Cropland",
    13: "Rice",
    14: "Permanent herbaceous",
    15: "Low Woody",
    16: "Evergreen forest",
    254: "Outside Area",
}

# =============================================================================
# 1. FILENAME PARSING
# =============================================================================

def parse_site_lc_from_ndvi_nirv_evi2(path):
    """
    Example:
    S2_2017_01_01_2024_12_31_AT-Neu_LC10_NDVI_NIRv_EVI2_mean.csv
    """
    name = Path(path).name
    m = re.match(
        r'^S2_\d{4}_\d{2}_\d{2}_\d{4}_\d{2}_\d{2}_(?P<site>[^_]+)_(?P<lc>LC\d+)_NDVI_NIRv_EVI2_mean\.csv$',
        name
    )
    if not m:
        raise ValueError(f"Cannot parse NDVI/NIRv/EVI2 filename: {name}")
    return m.group("site"), m.group("lc")


def parse_site_lc_generic(path, product):
    """
    Handles:
    PPI_YYYY_MM_DD_YYYY_MM_DD_SITE_LCxx_mean.csv
    FAPAR_...
    LAI_...
    """
    name = Path(path).name
    pattern = rf'^{product}_\d{{4}}_\d{{2}}_\d{{2}}_\d{{4}}_\d{{2}}_\d{{2}}_(?P<site>[^_]+)_(?P<lc>LC\d+)_mean\.csv$'
    m = re.match(pattern, name)
    if not m:
        raise ValueError(f"Cannot parse {product} filename: {name}")
    return m.group("site"), m.group("lc")


def parse_site_from_gpp(path):
    """
    New measured GPP format:
    FLUX_EO_cal_data_AT-Neu_DAILY_GPP_NT_VUT_MEAN.csv
    Extract SITE = AT-Neu
    """
    name = Path(path).name
    m = re.search(r'_([^_]+)_DAILY_GPP', name)
    if not m:
        raise ValueError(f"Cannot parse site from GPP file: {name}")
    return m.group(1)


def parse_site_from_gpp_lue(path):
    """
    LUE GPP format:
    FLUX_EO_cal_data_AT-Neu_DAILY_GPP_LUE.csv
    """
    name = Path(path).name
    m = re.search(r'_([^_]+)_DAILY_GPP_LUE', name)
    if not m:
        raise ValueError(f"Cannot parse site from LUE GPP filename: {name}")
    return m.group(1)


# =============================================================================
# 2. COLLECT PATHS
# =============================================================================

def collect_file_paths(ppi_dir, ndvi_dir, fapar_dir, lai_dir, gpp_dir, gpp_lue_dir):

    files = defaultdict(lambda: {
        "ndvi_nirv_evi2": None,
        "ppi": None,
        "fapar": None,
        "lai": None,
        "gpp": None,
        "gpp_m": None,
    })

    # NDVI + NIRv + EVI2
    for f in Path(ndvi_dir).glob("*.csv"):
        try:
            site, lc = parse_site_lc_from_ndvi_nirv_evi2(f)
            files[(site, lc)]["ndvi_nirv_evi2"] = str(f)
        except ValueError:
            continue

    # PPI
    for f in Path(ppi_dir).glob("*.csv"):
        try:
            site, lc = parse_site_lc_generic(f, "PPI")
            files[(site, lc)]["ppi"] = str(f)
        except ValueError:
            continue

    # FAPAR
    for f in Path(fapar_dir).glob("*.csv"):
        try:
            site, lc = parse_site_lc_generic(f, "FAPAR")
            files[(site, lc)]["fapar"] = str(f)
        except ValueError:
            continue

    # LAI
    for f in Path(lai_dir).glob("*.csv"):
        try:
            site, lc = parse_site_lc_generic(f, "LAI")
            files[(site, lc)]["lai"] = str(f)
        except ValueError:
            continue

    # Measured GPP (new NT_VUT_MEAN)
    gpp_map = {}
    for f in Path(gpp_dir).glob("*.csv"):
        site = parse_site_from_gpp(f)
        gpp_map[site] = str(f)

    # Modelled GPP LUE
    gpp_lue_map = {}
    for f in Path(gpp_lue_dir).glob("*.csv"):
        site = parse_site_from_gpp_lue(f)
        gpp_lue_map[site] = str(f)

    # Attach GPP & GPP_m
    for (site, lc) in files:
        if site in gpp_map:
            files[(site, lc)]["gpp"] = gpp_map[site]
        if site in gpp_lue_map:
            files[(site, lc)]["gpp_m"] = gpp_lue_map[site]

    return files


# =============================================================================
# 3. READ HELPERS
# =============================================================================

def read_gpp_daily(path):
    """
    New measured GPP format:
    t,GPP
    2017-01-01,
    2017-01-05,0.43
    """
    df = pd.read_csv(path)
    df["t"] = pd.to_datetime(df["t"]).dt.date
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df[["t", "GPP"]]


def read_gpp_lue_daily(path):
    df = pd.read_csv(path)
    df["t"] = pd.to_datetime(df["t"]).dt.date
    df = df.rename(columns={"GPP": "GPP_m"})[["t", "GPP_m"]]
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df


# =============================================================================
# 4. MERGE PER SITE + LC
# =============================================================================

def diagnose_vi(merged, vi, gpp_col="GPP"):
    """
    merged: merged DataFrame for a single site
    vi: VI column name (e.g., "NDVI", "EVI2", "PPI")
    gpp_col: ground GPP col (default: "GPP")
    """

    df = merged[[vi, gpp_col]].dropna()
    n = len(df)
    if n < 5:
        return None

    # 1) Spearman correlation (all data)
    rho_all, p_all = spearmanr(df[vi], df[gpp_col])

    # 2) Growing-season (only if enough data)
    gpp = df[gpp_col].values
    vi_values = df[vi].values
    gpp_max = np.nanmax(gpp)
    threshold = np.nanpercentile(gpp, 5)  # 或 5th percentile
    gs_mask = gpp > threshold
    
    if gs_mask.sum() >= 5:
        rho_gs, p_gs = spearmanr(gpp[gs_mask], vi_values[gs_mask])
    else:
        rho_gs, p_gs = np.nan, np.nan

    return {
        "VI": vi,
        "n": n,
        "rho_all": rho_all,
        "rho_gs": rho_gs
    }


def merge_vis_for_site_lc(site, lc, info):

    dfs = []

    # NDVI + NIRv + EVI2
    if info["ndvi_nirv_evi2"]:
        df = pd.read_csv(info["ndvi_nirv_evi2"])
        df["t"] = pd.to_datetime(df["t"]).dt.date
        dfs.append(df)

    # PPI
    if info["ppi"]:
        df = pd.read_csv(info["ppi"])
        df["t"] = pd.to_datetime(df["t"]).dt.date
        df = df.rename(columns={"PPI_mean": "PPI"})[["t", "PPI"]]
        df["PPI"] = df["PPI"].clip(-1, 5)
        # NEW: ReLU – force negative PPI to 0
        df["PPI"] = df["PPI"].clip(lower=0)
        dfs.append(df)

    # FAPAR
    if info["fapar"]:
        df = pd.read_csv(info["fapar"])
        df["t"] = pd.to_datetime(df["t"]).dt.date
        df = df.rename(columns={"FAPAR_mean": "FAPAR"})[["t", "FAPAR"]]
        dfs.append(df)

    # LAI
    if info["lai"]:
        df = pd.read_csv(info["lai"])
        df["t"] = pd.to_datetime(df["t"]).dt.date
        df = df.rename(columns={"LAI_mean": "LAI"})[["t", "LAI"]]
        dfs.append(df)

    # Measured GPP (new NT_VUT_MEAN)
    if info["gpp"]:
        df = read_gpp_daily(info["gpp"])
        dfs.append(df)

    # Modelled GPP from LUE
    if info["gpp_m"]:
        df = read_gpp_lue_daily(info["gpp_m"])
        dfs.append(df)

    if not dfs:
        raise ValueError(f"No data for {site}, {lc}")

    merged = dfs[0]
    for df in dfs[1:]:
        merged = pd.merge(merged, df, on="t", how="outer")

    merged = merged.sort_values("t")

    vis_cols = ["NDVI","NIRv","EVI2","PPI","LAI","FAPAR","GPP_m","GPP"]

    # ensure all columns present
    for col in vis_cols:
        if col not in merged:
            merged[col] = np.nan

    merged.replace([np.inf, -np.inf], np.nan, inplace=True)

    # --- NEW: drop any row that has a missing value in *any* VI column ---
    merged_cut = merged.dropna(subset=vis_cols, how="any")

    # --- Perform VI diagnostics ---
    diags = []
    for vi in vis_cols:
        if vi == "GPP":  # skip ground GPP
            continue
        d = diagnose_vi(merged_cut, vi, gpp_col="GPP")
        if d is not None:
            d["site"] = site
            d["lc"] = lc
            diags.append(d)


    return merged[["t"] + vis_cols], merged_cut[["t"] + vis_cols], pd.DataFrame(diags)


# =============================================================================
# 5. ROBUST LIMITS
# =============================================================================

def robust_limits(a, q1=0.01, q2=0.99, expand=1.1):
    a = a[np.isfinite(a)]
    if len(a) < 20:
        return None
    lo = np.nanquantile(a, q1)
    hi = np.nanquantile(a, q2)
    return expand*lo, expand*hi


# =============================================================================
# 6. SCATTER MATRIX
# =============================================================================
VI_ORDER_FIXED = ["GPP_m","PPI","LAI","EVI2","NIRv","FAPAR","NDVI"]
VI_ORDER_FIXED = ["NDVI","EVI2","NIRv","PPI","LAI","FAPAR"]

VI_COLORS = {
    "NDVI": "#1f77b4",
    "EVI2": "#ff7f0e",
    "NIRv": "#2ca02c",
    "PPI":  "#d62728",
    "LAI":  "#9467bd",
    "FAPAR": "#8c564b",
    "GPP_m": "#e377c2"
}


def make_scatter_matrix_for_site(site, lc, df_site, out_plot_dir):

    vis = ["NDVI","EVI2","NIRv","PPI","LAI","FAPAR","GPP_m","GPP"]
    n = len(vis)

    fig, axes = plt.subplots(n, n, figsize=(1.45*n, 1.45*n), squeeze=False)

    plt.subplots_adjust(
        left=0.07, right=0.93,
        top=0.93, bottom=0.07,
        wspace=0.02, hspace=0.02
    )

    corr_mat = pd.DataFrame(np.nan, index=vis, columns=vis)

    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            y = vis[i]
            x = vis[j]

            # 对角：直方图
            if i == j:
                arr = df_site[y].to_numpy()
                arr = arr[np.isfinite(arr)]
                if len(arr) > 0:
                    lim = robust_limits(arr)
                    if lim:
                        ax.hist(arr, bins=30, range=lim, edgecolor='black', alpha=0.7)
                        ax.set_xlim(lim)
                    else:
                        ax.hist(arr, bins=30, edgecolor='black', alpha=0.7)

                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
                ax.set_ylabel(y, fontsize=10)
                ax.set_xlabel("")

            # 下三角：散点 + Spearman
            elif i > j:
                xx = df_site[x].to_numpy()
                yy = df_site[y].to_numpy()
                mask = np.isfinite(xx) & np.isfinite(yy)
                if np.sum(mask) > 0:
                    ax.scatter(xx[mask], yy[mask], s=6, alpha=0.4)

                if np.sum(mask) > 2:
                    rho = spearmanr(xx[mask], yy[mask])[0]
                else:
                    rho = np.nan
                N = int(np.sum(mask))

                corr_mat.loc[y, x] = rho
                corr_mat.loc[x, y] = rho

                ax.text(
                    0.05, 0.95,
                    f"ρ={rho:.2f}\nN={N}",
                    transform=ax.transAxes,
                    ha="left", va="top",
                    fontsize=9,
                    bbox=dict(boxstyle="round",
                              facecolor="white",
                              alpha=0.7,
                              edgecolor='none')
                )

                if j != 0:
                    ax.set_yticklabels([])
                else:
                    ax.set_ylabel(y, fontsize=9)

                if i != n - 1:
                    ax.set_xticklabels([])
                else:
                    ax.set_xlabel(x, fontsize=9)

                limx = robust_limits(xx)
                if limx:
                    ax.set_xlim(limx)
                limy = robust_limits(yy)
                if limy:
                    ax.set_ylim(limy)

            # 上三角空
            else:
                ax.axis("off")

    fig.suptitle(f"VI Scatter Matrix (site={site}, LC={lc})", y=0.995, fontsize=10)

    out_plot_dir = Path(out_plot_dir)
    out_plot_dir.mkdir(parents=True, exist_ok=True)

    fig.savefig(out_plot_dir / f"scatter_matrix_{site}_{lc}.png", dpi=350)
    plt.close(fig)

    # corr_mat.to_csv(out_plot_dir / f"spearman_corr_{site}_{lc}.csv", float_format="%.4f")


def plot_box_overall(df_diags_all, out_plot_dir):
    """
    图 1：所有 site 混合，一个 VI 一个 box，
    按固定 VI 顺序排序。
    """
    out_plot_dir = Path(out_plot_dir)
    out_plot_dir.mkdir(parents=True, exist_ok=True)

    # 只保留数据中存在的 VI
    vi_order = [vi for vi in VI_ORDER_FIXED if vi in df_diags_all["VI"].unique()]

    data = [
        df_diags_all.loc[df_diags_all["VI"] == vi, "rho_gs"].dropna().values
        for vi in vi_order
    ]

    fig, ax = plt.subplots(figsize=(max(7, 1.2*len(vi_order)), 4))

    bp = ax.boxplot(
        data,
        labels=vi_order,
        showfliers=True,
        patch_artist=True
    )

    # 设置颜色
    for patch, vi in zip(bp["boxes"], vi_order):
        color = VI_COLORS.get(vi, "#999999")
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xlabel("VI")
    ax.set_ylabel("Spearman ρ (Growing Season)")
    ax.set_title("Distribution of VI–GPP Spearman ρ (All Sites)")

    fig.tight_layout()
    fig.savefig(out_plot_dir / "Raw_VIs_rho_gs_box_overall.png", dpi=300)
    plt.close(fig)



def plot_box_by_lc(df_diags_all, out_plot_dir):
    """
    图 2：按 LC 分组，组内为不同 VI 的 box，
    按固定 VI 顺序排列，颜色固定。
    """
    out_plot_dir = Path(out_plot_dir)
    out_plot_dir.mkdir(parents=True, exist_ok=True)

    vi_order = [vi for vi in VI_ORDER_FIXED if vi in df_diags_all["VI"].unique()]
    lcs = sorted(df_diags_all["lc"].unique())

    n_lc = len(lcs)
    n_vi = len(vi_order)

    fig, ax = plt.subplots(figsize=(max(9, 1.6*n_lc), 6))

    group_width = 0.8
    box_width = group_width / n_vi

    vi_handles = {}

    for i_lc, lc in enumerate(lcs):
        base = i_lc

        for i_vi, vi in enumerate(vi_order):
            sub = df_diags_all[
                (df_diags_all["lc"] == lc) &
                (df_diags_all["VI"] == vi)
            ]["rho_gs"].dropna()

            if len(sub) == 0:
                continue

            pos = base + (i_vi - (n_vi - 1) / 2) * box_width

            bp = ax.boxplot(
                [sub.values],
                positions=[pos],
                widths=box_width * 0.8,
                patch_artist=True
            )

            color = VI_COLORS.get(vi, "#999999")
            bp["boxes"][0].set_facecolor(color)
            bp["boxes"][0].set_alpha(0.7)

            if vi not in vi_handles:
                vi_handles[vi] = bp["boxes"][0]

    ax.set_xticks(range(n_lc))
    ax.set_xticklabels(lcs, rotation=45, ha="right")

    ax.set_xlabel("Land-cover (LC)")
    ax.set_ylabel("Spearman ρ (Growing Season)")
    ax.set_title("VI–GPP Spearman ρ by LC")

    handles = [vi_handles[vi] for vi in vi_order if vi in vi_handles]
    labels = [vi for vi in vi_order if vi in vi_handles]
    ax.legend(handles, labels, title="VI", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_plot_dir / "Raw_VIs_rho_gs_box_by_LC.png", dpi=300)
    plt.close(fig)

def plot_timeseries_for_site(site, lc, df_site, out_plot_dir):
    """
    对单个 site+LC 画时间序列：
    - 每个 VI 一个子图（按照 VI_ORDER_FIXED 顺序）
    - 每个子图左轴是该 VI，右轴是参考 GPP 作为灰色背景
    """

    out_plot_dir = Path(out_plot_dir)
    out_plot_dir.mkdir(parents=True, exist_ok=True)

    df = df_site.sort_values("t").copy()
    dates = pd.to_datetime(df["t"])

    # 如果没有 GPP，就不画
    if "GPP" not in df.columns:
        return

    gpp = df["GPP"].to_numpy()
    gpp_mask = np.isfinite(gpp)
    if not gpp_mask.any():
        return

    # 确定要画的 VI 列表（按固定顺序，只保留存在的列）
    vis_list = [vi for vi in VI_ORDER_FIXED if vi in df.columns]

    if len(vis_list) == 0:
        return

    n_rows = len(vis_list)
    fig, axes = plt.subplots(
        n_rows, 1,
        figsize=(10, 2.0 * n_rows),
        sharex=True
    )

    # 如果只有一个 VI，axes 不是 list，需要统一成 list
    if n_rows == 1:
        axes = [axes]

    for i, vi in enumerate(vis_list):
        ax = axes[i]

        # 左轴：VI
        y_vi = df[vi].to_numpy()
        mask_vi = np.isfinite(y_vi)

        color_vi = VI_COLORS.get(vi, None)
        if color_vi is None:
            ax.plot(dates[mask_vi], y_vi[mask_vi],
                    "-o", markersize=2, linewidth=0.8)
        else:
            ax.plot(dates[mask_vi], y_vi[mask_vi],
                    "-o", markersize=2, linewidth=0.8, color=color_vi)

        ax.set_ylabel(vi)

        # 右轴：GPP 灰色背景
        ax2 = ax.twinx()
        ax2.plot(dates[gpp_mask], gpp[gpp_mask],
                 "-", linewidth=0.8, color="0.7", alpha=0.7)
        ax2.set_ylabel("GPP", color="0.5")
        ax2.tick_params(axis="y", labelcolor="0.5")

        # 让 VI 线更“前景”，GPP 在“背景”
        ax.set_zorder(ax2.get_zorder() + 1)
        ax.patch.set_visible(False)

        if i == 0:
            ax.set_title(f"Time series – {site} ({lc})")

    axes[-1].set_xlabel("Date")
    fig.autofmt_xdate()
    fig.tight_layout()

    out_png = out_plot_dir / f"time_series_{site}_{lc}.png"
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_timeseries_and_scatter_for_site(site, lc, df_site, out_plot_dir):
    """
    Combined diagnostic plot for one site + LC:
    - Left: VI time series (left axis = VI, right axis = GPP in grey)
    - Right: VI vs GPP scatter
    - Rows follow VI_ORDER_FIXED
    """

    out_plot_dir = Path(out_plot_dir)
    out_plot_dir.mkdir(parents=True, exist_ok=True)

    # Required columns
    if "GPP" not in df_site.columns:
        return

    df = df_site.sort_values("t").copy()
    dates = pd.to_datetime(df["t"])

    gpp = df["GPP"].to_numpy()
    gpp_mask = np.isfinite(gpp)
    if not gpp_mask.any():
        return

    # Keep only VIs that exist in data
    vis_list = [vi for vi in VI_ORDER_FIXED if vi in df.columns]
    if len(vis_list) == 0:
        return

    n_rows = len(vis_list)

    fig, axes = plt.subplots(
        n_rows, 2,
        figsize=(12, 2.2 * n_rows),
        sharex="col",
        gridspec_kw={"width_ratios": [2.2, 1.3]}
    )

    if n_rows == 1:
        axes = np.array([axes])

    for i, vi in enumerate(vis_list):

        # =====================
        # LEFT: time series
        # =====================
        ax_ts = axes[i, 0]

        y_vi = df[vi].to_numpy()
        mask_vi = np.isfinite(y_vi)

        color_vi = VI_COLORS.get(vi, None)
        ax_ts.plot(
            dates[mask_vi], y_vi[mask_vi],
            "-o", markersize=2, linewidth=0.8,
            color=color_vi
        )

        ax_ts.set_ylabel(vi)

        # GPP background (right axis)
        ax_ts2 = ax_ts.twinx()
        ax_ts2.plot(
            dates[gpp_mask], gpp[gpp_mask],
            "-", linewidth=0.8, color="0.7", alpha=0.7
        )
        ax_ts2.set_ylabel("GPP", color="0.5")
        ax_ts2.tick_params(axis="y", labelcolor="0.5")

        ax_ts.set_zorder(ax_ts2.get_zorder() + 1)
        ax_ts.patch.set_visible(False)

        if i == 0:
            ax_ts.set_title("Time series")

        # =====================
        # RIGHT: scatter
        # =====================
        ax_sc = axes[i, 1]

        xx = gpp          # x-axis: GPP
        yy = y_vi         # y-axis: VI
        mask = np.isfinite(xx) & np.isfinite(yy)

        if mask.sum() > 0:
            ax_sc.scatter(
                xx[mask], yy[mask],
                s=8, alpha=0.4, color=color_vi
            )

        if mask.sum() > 2:
            rho = spearmanr(xx[mask], yy[mask])[0]
        else:
            rho = np.nan

        N = int(mask.sum())


        ax_sc.text(
            0.05, 0.95,
            f"ρ={rho:.2f}\nN={N}",
            transform=ax_sc.transAxes,
            ha="left", va="top",
            fontsize=9,
            bbox=dict(
                boxstyle="round",
                facecolor="white",
                alpha=0.7,
                edgecolor="none"
            )
        )

        limx = robust_limits(xx)
        if limx:
            ax_sc.set_xlim(limx)

        limy = robust_limits(yy)
        if limy:
            ax_sc.set_ylim(limy)

        ax_sc.set_xlabel("GPP")
        ax_sc.set_ylabel(vi)

        if i == 0:
            ax_sc.set_title("VI vs GPP")

    axes[-1, 0].set_xlabel("Date")
    axes[-1, 1].set_xlabel("VI value")

    fig.suptitle(f"{site} – {lc}", y=0.995)
    fig.autofmt_xdate()
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    out_png = out_plot_dir / f"time_series_scatter_{site}_{lc}.png"
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

# =============================================================================
# 7. MAIN
# =============================================================================

def main():

    p = argparse.ArgumentParser()
    p.add_argument("--ppi_dir", type=str, default="/projects/eko/fs7/proj/HRVPP2/calval/VI/PPI/flux_cal/V1/csv_lc/")
    p.add_argument("--ndvi_dir", type=str, default="/projects/eko/fs7/proj/HRVPP2/calval/VI/NDVI_EVI2/flux_cal/V1/csv_lc/")
    p.add_argument("--fapar_dir", type=str, default="/projects/eko/fs7/proj/HRVPP2/calval/VI/FAPAR/V5/timeseries/csv_lc_200m/")
    p.add_argument("--lai_dir", type=str, default="/projects/eko/fs7/proj/HRVPP2/calval/VI/LAI/flux_cal/V1/csv_lc/")
    p.add_argument("--gpp_dir", type=str, default="/projects/eko/fs7/proj/HRVPP2/calval/VI/GPP_NT_VUT_MEAN/")
    p.add_argument("--gpp_lue_dir", type=str, default="/projects/eko/fs7/proj/HRVPP2/calval/VI/GPP_LUE/")
    p.add_argument("--out_vi_dir", type=str, default="output/Cal/rawvi")
    p.add_argument("--out_plot_dir", type=str, default="output/Cal/rawvi")

    args = p.parse_args()

    Path(args.out_vi_dir).mkdir(parents=True, exist_ok=True)

    files = collect_file_paths(
        args.ppi_dir, args.ndvi_dir, args.fapar_dir, args.lai_dir,
        args.gpp_dir, args.gpp_lue_dir
    )

    print(f"Found {len(files)} site/LC groups")

    lc_to_dfs = defaultdict(list)
    all_diags = []   # 新增：收集所有站点的诊断结果


    for (site,lc),info in sorted(files.items()):
        try:
            merged, merged_cut, diag_vi = merge_vis_for_site_lc(site,lc,info)
        except Exception as e:
            print(f"[WARN] {site},{lc}: {e}")
            continue

        if True:
            plot_timeseries_and_scatter_for_site(
                site, lc, merged, args.out_plot_dir
            )


                

        # out_csv = Path(args.out_vi_dir) / f"Raw_VIs_diag_{site}_{lc}.csv"
        # diag_vi.to_csv(out_csv, index=False)
        # print(f"[INFO] wrote {out_csv}")
        
        if 0: # if want to save png
            # 每个 site 单独画 scatter matrix
            make_scatter_matrix_for_site(site, lc, merged_cut, args.out_plot_dir)

            # --- 每个 site+LC 画时间序列图 ---
            plot_timeseries_for_site(site, lc, merged, args.out_plot_dir)
            # save merged CSV
            out_csv = Path(args.out_vi_dir) / f"Raw_VIs_data_{site}_{lc}.csv"
            merged.to_csv(out_csv, index=False)
            print(f"[INFO] wrote {out_csv}")
            merged["site"] = site
            merged["lc"] = lc

        # 如果你暂时不需要按 LC 聚合的 matrix，可以注释掉下一行
        # lc_to_dfs[lc].append(merged)

        # 新增：把诊断结果收集起来（做跨-site诊断）
        all_diags.append(diag_vi)

    # ==============================
    # 汇总所有站点的诊断结果
    # ==============================
    if all_diags:
        df_diags_all = pd.concat(all_diags, ignore_index=True)
        out_csv = Path(args.out_vi_dir) / "Raw_VIs_diag_ALL_sites.csv"
        df_diags_all.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")
    else:
        print("[WARN] no diagnostics found")
        df_diags_all = None


    # for lc,dfs in sorted(lc_to_dfs.items()):
    #     df_all = pd.concat(dfs, ignore_index=True)

    #     # skip empty
    #     if df_all[["NDVI","NIRv","EVI2","PPI","LAI","FAPAR","GPP_m","GPP"]].isna().all().all():
    #         print(f"[WARN] LC={lc}: all NaN")
    #         continue

    #     print(f"[INFO] Scatter matrix LC={lc}, N={len(df_all)}")
    #     make_scatter_matrix_for_lc(lc, df_all, args.out_plot_dir)

    if df_diags_all is not None:

        # --- Convert LC tag "LC10" → "Permanent Grassland" ---
        def lc_to_name(lc_str):
            match = re.match(r"LC(\d+)", lc_str)
            if match:
                code = int(match.group(1))
                return LC_PALETTE.get(code, lc_str)
            else:
                return lc_str

        df_diags_all["lc"] = df_diags_all["lc"].apply(lc_to_name)


        # 按VI做一个总的 summary
        summary_vi = (
            df_diags_all
            .groupby("VI")
            .agg(
                mean_rho_gs=("rho_gs", "mean"),
                median_rho_gs=("rho_gs", "median"),
                num_sites=("site", "nunique")
            )
            .reset_index()
        )

        out_csv = Path(args.out_vi_dir) / "Raw_VIs_diag_summary_by_VI.csv"
        summary_vi.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")

        print("\n[SUMMARY by VI]")
        print(summary_vi)


        summary_vi_lc = (
            df_diags_all
            .groupby(["lc", "VI"])
            .agg(
                mean_rho_gs=("rho_gs", "mean"),
                median_rho_gs=("rho_gs", "median"),
                num_sites=("site", "nunique")
            )
            .reset_index()
        )

        out_csv = Path(args.out_vi_dir) / "Raw_VIs_diag_summary_by_LC_VI.csv"
        summary_vi_lc.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")


        # ============================
        # 每个 VI 的最差 & 最好站点（按 rho_gs）
        # ============================
        valid = df_diags_all.dropna(subset=["rho_gs"])

        if valid.empty:
            print("[WARN] No valid rho_gs values for best/worst site analysis")
        else:
            # ---- Worst (min rho_gs) ----
            idx_min = valid.groupby("VI")["rho_gs"].idxmin()
            worst_sites = (
                valid.loc[idx_min, ["VI", "site", "lc", "rho_gs", "n"]]
                .sort_values("VI")
                .reset_index(drop=True)
            )

            out_csv_min = Path(args.out_vi_dir) / "Raw_VIs_worst_sites_by_VI.csv"
            worst_sites.to_csv(out_csv_min, index=False)
            print(f"[INFO] wrote {out_csv_min}")

            print("\n[Worst sites per VI (lowest rho_gs)]")
            print(worst_sites)

            # ---- Best (max rho_gs) ----
            idx_max = valid.groupby("VI")["rho_gs"].idxmax()
            best_sites = (
                valid.loc[idx_max, ["VI", "site", "lc", "rho_gs", "n"]]
                .sort_values("VI")
                .reset_index(drop=True)
            )

            out_csv_max = Path(args.out_vi_dir) / "Raw_VIs_best_sites_by_VI.csv"
            best_sites.to_csv(out_csv_max, index=False)
            print(f"[INFO] wrote {out_csv_max}")

            print("\n[Best sites per VI (highest rho_gs)]")
            print(best_sites)



        # 画两个 boxplot 图
        plot_box_overall(df_diags_all, args.out_plot_dir)
        plot_box_by_lc(df_diags_all, args.out_plot_dir)





if __name__ == "__main__":
    main()
