import pandas as pd
import numpy as np
from scipy.stats import friedmanchisquare, wilcoxon
import matplotlib.pyplot as plt

# 可选：用于多重比较校正（Holm）
try:
    from statsmodels.stats.multitest import multipletests
    _HAS_STATSMODELS = True
except Exception:
    _HAS_STATSMODELS = False


CSV_PATH = "/projects/eko/fs7/proj/HRVPP2/calval/cal/HRVPP2-calibration/output/Cal/rawvi/Raw_VIs_diag_ALL_sites.csv"  # 若不在当前目录，请写完整路径
METRIC = "rho_gs"                        # 你说用 rho_gs
HIGHER_IS_BETTER = True                  # Spearman rho：越大越好


def holm_correction(pvals: np.ndarray) -> np.ndarray:
    """
    Holm-Bonferroni 校正的简易实现（当 statsmodels 不可用时备用）
    返回校正后的 p 值（step-down）
    """
    pvals = np.asarray(pvals, dtype=float)
    m = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]

    adj = np.empty(m, dtype=float)
    for i, p in enumerate(ranked):
        adj[i] = (m - i) * p

    # step-down：保证单调不减
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)

    out = np.empty(m, dtype=float)
    out[order] = adj
    return out


def prepare_wide(df_long: pd.DataFrame, metric: str) -> tuple[pd.DataFrame, pd.Series]:
    """
    long -> wide:
      wide: index=site, columns=VI, values=metric
      lc:   Series with site index
    """
    required = {"site", "VI", metric}
    missing = required - set(df_long.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found: {list(df_long.columns)}")

    wide = df_long.pivot_table(index="site", columns="VI", values=metric, aggfunc="first")

    lc = None
    if "lc" in df_long.columns:
        lc = df_long.drop_duplicates("site").set_index("site")["lc"]
        lc = lc.reindex(wide.index)

    return wide, lc


def summarize_performance(wide: pd.DataFrame, higher_is_better: bool = True) -> pd.DataFrame:
    """
    输出每个 VI 的：
      - mean_rho
      - median_rho
      - mean_rank（越小越好）
      - win_rate（第一名比例，平局算共享胜利）
    """
    # 完整配对：只保留所有 VI 都有值的 site
    wide_cc = wide.dropna(axis=0, how="any")
    if wide_cc.empty:
        raise ValueError("No complete-case sites after dropping NA.")

    # rank：最好=1
    ranks = wide_cc.rank(axis=1, ascending=not higher_is_better, method="average")
    mean_rank = ranks.mean(axis=0)

    mean_rho = wide_cc.mean(axis=0)
    median_rho = wide_cc.median(axis=0)

    # win rate：每个 site 的最大值对应的 VI（平局共享）
    best_per_site = wide_cc.max(axis=1) if higher_is_better else wide_cc.min(axis=1)
    winners = wide_cc.eq(best_per_site, axis=0)
    win_rate = winners.sum(axis=0) / len(wide_cc)

    out = pd.DataFrame({
        "mean_rho": mean_rho,
        "median_rho": median_rho,
        "mean_rank(1=best)": mean_rank,
        "win_rate": win_rate
    }).sort_values(["mean_rank(1=best)", "win_rate", "median_rho"], ascending=[True, False, False])

    out.index.name = "VI"
    return out


def friedman_test(wide: pd.DataFrame) -> tuple[float, float, pd.DataFrame]:
    """
    Friedman 检验（配对、多个方法/VI）
    返回：stat, p, summary_table
    """
    wide_cc = wide.dropna(axis=0, how="any")
    cols = list(wide_cc.columns)

    # scipy 需要每个“处理条件”传一个数组
    arrays = [wide_cc[c].to_numpy() for c in cols]
    stat, p = friedmanchisquare(*arrays)

    summary = summarize_performance(wide_cc, higher_is_better=HIGHER_IS_BETTER)
    return stat, p, summary


def posthoc_best_vs_others(wide: pd.DataFrame, best_vi: str) -> pd.DataFrame:
    """
    事后比较：最佳 VI vs 其它 VI（配对 Wilcoxon）
    默认做单侧检验：best_vi 的 rho_gs 是否更大
    并进行 Holm 校正
    """
    wide_cc = wide.dropna(axis=0, how="any")
    cols = [c for c in wide_cc.columns if c != best_vi]

    results = []
    pvals = []

    for vi in cols:
        x = wide_cc[best_vi].to_numpy()
        y = wide_cc[vi].to_numpy()

        # 单侧：best > other
        alt = "greater" if HIGHER_IS_BETTER else "less"
        stat, p = wilcoxon(x, y, alternative=alt, zero_method="wilcox")

        diff = x - y
        results.append({
            "comparison": f"{best_vi} vs {vi}",
            "n_sites": len(diff),
            "median_diff": float(np.median(diff)),
            "mean_diff": float(np.mean(diff)),
            "wilcoxon_stat": float(stat),
            "p_raw": float(p)
        })
        pvals.append(p)

    pvals = np.array(pvals, dtype=float)
    if _HAS_STATSMODELS:
        _, p_adj, _, _ = multipletests(pvals, method="holm")
    else:
        p_adj = holm_correction(pvals)

    for r, pa in zip(results, p_adj):
        r["p_holm"] = float(pa)

    out = pd.DataFrame(results).sort_values("p_holm")
    return out


def plot_box_and_paired(wide: pd.DataFrame, title_prefix: str = "") -> None:
    """
    两张图：
      1) rho_gs boxplot（描述性）
      2) paired lines（每个 site 一条线，体现配对结构）
    """
    wide_cc = wide.dropna(axis=0, how="any")
    cols = list(wide_cc.columns)

    # 1) boxplot
    plt.figure()
    plt.boxplot([wide_cc[c].to_numpy() for c in cols], labels=cols, showfliers=True)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel(METRIC)
    plt.title(f"{title_prefix}{METRIC} distribution (site-level, descriptive)")
    plt.tight_layout()
    plt.show()

    # 2) paired lines
    plt.figure()
    x = np.arange(len(cols))
    for _, row in wide_cc.iterrows():
        plt.plot(x, row[cols].to_numpy(), alpha=0.25)
    plt.xticks(x, cols, rotation=45, ha="right")
    plt.ylabel(METRIC)
    plt.title(f"{title_prefix}Paired trajectories across VIs (each line = one site)")
    plt.tight_layout()
    plt.show()


def run_overall_and_by_lc(df_long: pd.DataFrame) -> None:
    wide, lc = prepare_wide(df_long, METRIC)

    print("=== OVERALL (all sites, paired by site) ===")
    stat, p, summary = friedman_test(wide)
    print(f"Friedman chi2 = {stat:.4f}, p = {p:.4g}")
    print("\nPerformance summary (lower mean_rank is better):")
    print(summary.to_string(float_format=lambda v: f"{v:.4f}"))

    best_vi = summary.index[0]
    print(f"\nBest VI by mean_rank: {best_vi}")

    posthoc = posthoc_best_vs_others(wide, best_vi)
    print("\nPost-hoc: best vs others (paired Wilcoxon, one-sided), Holm-corrected:")
    print(posthoc.to_string(index=False, float_format=lambda v: f"{v:.4g}" if abs(v) < 0.001 else f"{v:.4f}"))

    plot_box_and_paired(wide, title_prefix="OVERALL: ")

    # 按 lc 分组（如果有）
    if lc is not None and lc.notna().any():
        print("\n=== BY LAND-COVER (lc) ===")
        for cls in lc.dropna().unique():
            sites = lc[lc == cls].index
            wide_sub = wide.loc[sites]
            wide_sub = wide_sub.dropna(axis=0, how="any")
            if len(wide_sub) < 5:
                print(f"\n[lc={cls}] skipped (too few complete-case sites: {len(wide_sub)})")
                continue

            stat, p, summary = friedman_test(wide_sub)
            best_vi = summary.index[0]

            print(f"\n[lc={cls}] n_sites={len(wide_sub)} | Friedman chi2={stat:.4f}, p={p:.4g} | best={best_vi}")
            print(summary.to_string(float_format=lambda v: f"{v:.4f}"))

            posthoc = posthoc_best_vs_others(wide_sub, best_vi)
            print("\nPost-hoc (Holm):")
            print(posthoc.to_string(index=False, float_format=lambda v: f"{v:.4g}" if abs(v) < 0.001 else f"{v:.4f}"))

            plot_box_and_paired(wide_sub, title_prefix=f"lc={cls}: ")


if __name__ == "__main__":
    df = pd.read_csv(CSV_PATH)
    # 可选：如果你只想用 growing-season 的点数 n 过滤（比如 n>=某阈值）
    # df = df[df["n"] >= 60].copy()
    # exclude GPP_m
    df = df[df["VI"] != "GPP_m"].copy()
    
    run_overall_and_by_lc(df)
