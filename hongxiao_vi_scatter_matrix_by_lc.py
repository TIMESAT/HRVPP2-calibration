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

    # ensure all columns present
    for col in ["NDVI","NIRv","EVI2","PPI","LAI","FAPAR","GPP_m","GPP"]:
        if col not in merged:
            merged[col] = np.nan

    merged.replace([np.inf, -np.inf], np.nan, inplace=True)

    return merged[["t","NDVI","NIRv","EVI2","PPI","LAI","FAPAR","GPP_m","GPP"]]


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

def make_scatter_matrix_for_lc(lc, df_all, out_plot_dir):

    vis = ["NDVI","EVI2","NIRv","PPI","LAI","FAPAR","GPP_m","GPP"]
    n = len(vis)

    fig, axes = plt.subplots(n, n, figsize=(1.45*n, 1.45*n), squeeze=False)

    plt.subplots_adjust(
        left=0.07, right=0.93,
        top=0.93, bottom=0.07,
        wspace=0.02, hspace=0.02
    )

    sites = sorted(df_all["site"].unique())
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    markers = ['o','s','^','D','v','P','X','*','h','8','<','>']

    site_to_style = {
        s: (colors[i % len(colors)], markers[(i // len(colors)) % len(markers)])
        for i,s in enumerate(sites)
    }

    corr_mat = pd.DataFrame(np.nan, index=vis, columns=vis)

    for i in range(n):
        for j in range(n):

            ax = axes[i,j]
            y = vis[i]
            x = vis[j]

            # DIAGONAL HISTOGRAM
            if i == j:
                arr = df_all[y].to_numpy()
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

            # LOWER TRIANGLE SCATTERS
            elif i > j:

                for s,(c,m) in site_to_style.items():
                    sub = df_all[df_all["site"]==s]
                    xx = sub[x].to_numpy()
                    yy = sub[y].to_numpy()
                    mask = np.isfinite(xx) & np.isfinite(yy)
                    if np.sum(mask)>0:
                        ax.scatter(xx[mask], yy[mask], s=6, alpha=0.25, c=[c], marker=m)

                # Spearman
                xx = df_all[x].to_numpy()
                yy = df_all[y].to_numpy()
                mask = np.isfinite(xx) & np.isfinite(yy)
                N = np.sum(mask)
                rho = spearmanr(xx[mask], yy[mask])[0] if N>2 else np.nan

                corr_mat.loc[y,x] = rho
                corr_mat.loc[x,y] = rho

                ax.text(
                    0.05, 0.95,
                    f"ρ={rho:.2f}\nN={N}",
                    transform=ax.transAxes,
                    ha="left", va="top",
                    fontsize=10,
                    bbox=dict(boxstyle="round",
                              facecolor="white",
                              alpha=0.7,
                              edgecolor='none')
                )

                if j != 0:
                    ax.set_yticklabels([])
                else:
                    ax.set_ylabel(y, fontsize=10)

                if i != n - 1:
                    ax.set_xticklabels([])
                else:
                    ax.set_xlabel(x, fontsize=10)

                limx = robust_limits(xx)
                if limx: ax.set_xlim(limx)
                limy = robust_limits(yy)
                if limy: ax.set_ylim(limy)

            # UPPER TRIANGLE EMPTY
            else:
                ax.axis("off")

    # LEGEND
    legend_handles = [
        Line2D([0],[0], marker=m, linestyle='',
               markerfacecolor=c, markeredgecolor='none',
               markersize=7, label=s)
        for s,(c,m) in site_to_style.items()
    ]
    axes[0,n-1].legend(handles=legend_handles, title="Site", loc="upper right", fontsize=10)

    fig.suptitle(f"VI Scatter Matrix (LC={lc})", y=0.995, fontsize=10)

    out_plot_dir = Path(out_plot_dir)
    out_plot_dir.mkdir(parents=True, exist_ok=True)

    fig.savefig(out_plot_dir / f"scatter_matrix_{lc}.png", dpi=350)
    plt.close(fig)

    corr_mat.to_csv(out_plot_dir / f"spearman_corr_{lc}.csv", float_format="%.4f")


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
    p.add_argument("--out_vi_dir", type=str, default="temp/")
    p.add_argument("--out_plot_dir", type=str, default="temp/")

    args = p.parse_args()

    Path(args.out_vi_dir).mkdir(parents=True, exist_ok=True)

    files = collect_file_paths(
        args.ppi_dir, args.ndvi_dir, args.fapar_dir, args.lai_dir,
        args.gpp_dir, args.gpp_lue_dir
    )

    print(f"Found {len(files)} site/LC groups")

    lc_to_dfs = defaultdict(list)

    for (site,lc),info in sorted(files.items()):
        try:
            merged = merge_vis_for_site_lc(site,lc,info)
        except Exception as e:
            print(f"[WARN] {site},{lc}: {e}")
            continue

        # derive date-range tag
        if info["ndvi_nirv_evi2"]:
            name = Path(info["ndvi_nirv_evi2"]).name
            m = re.match(r'^S2_(\d{4}_\d{2}_\d{2}_\d{4}_\d{2}_\d{2})_', name)
            dr = m.group(1) if m else "merged"
        else:
            dr = "merged"

        # save merged CSV
        out_csv = Path(args.out_vi_dir) / f"VIs_{dr}_{site}_{lc}_with_GPP_GPPm.csv"
        merged.to_csv(out_csv, index=False)
        print(f"[INFO] wrote {out_csv}")

        merged["site"] = site
        merged["lc"] = lc
        lc_to_dfs[lc].append(merged)

    for lc,dfs in sorted(lc_to_dfs.items()):
        df_all = pd.concat(dfs, ignore_index=True)

        # skip empty
        if df_all[["NDVI","NIRv","EVI2","PPI","LAI","FAPAR","GPP_m","GPP"]].isna().all().all():
            print(f"[WARN] LC={lc}: all NaN")
            continue

        print(f"[INFO] Scatter matrix LC={lc}, N={len(df_all)}")
        make_scatter_matrix_for_lc(lc, df_all, args.out_plot_dir)


if __name__ == "__main__":
    main()
