import os
import pandas as pd
import matplotlib.pyplot as plt

# ---------------- paths ----------------
PPI_PATH = "output/PPI_scores_by_Setting_x_lc.csv"
SETTINGS_PATH = "output/settings_vpp.csv"
OUT_DIR = "output/boxplots_by_lc"
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------- load data ----------------
ppi = pd.read_csv(PPI_PATH)
settings = pd.read_csv(SETTINGS_PATH)

# Normalize column names
ppi.columns = ppi.columns.astype(str).str.strip()
settings.columns = settings.columns.str.strip()

# ---------------- merge (IMPORTANT FIX) ----------------
# Setting  <->  settings_id
df = ppi.merge(
    settings,
    left_on="Setting",
    right_on="settings_id",
    how="left"
)

# Required grouping columns
group_cols = ["method_name", "smooth", "seasonmethod", "seapar"]
missing = [c for c in group_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# ---------------- build x-axis group label ----------------
df["group"] = (
    df["method_name"].astype(str) + " | "
    + "smooth=" + df["smooth"].astype(str) + " | "
    + "season=" + df["seasonmethod"].astype(str) + " | "
    + "seapar=" + df["seapar"].astype(str)
)

# ---------------- LC columns ----------------
lc_cols = [c for c in ppi.columns if c != "Setting"]

# ---------------- plotting ----------------
for lc in lc_cols:
    plot_df = df[["group", lc]].dropna()

    if plot_df.empty:
        continue

    groups = sorted(plot_df["group"].unique())
    data = [plot_df.loc[plot_df["group"] == g, lc] for g in groups]

    plt.figure(figsize=(max(12, 0.6 * len(groups)), 6))
    plt.boxplot(
        data,
        labels=groups,
        showfliers=False
    )

    plt.xticks(rotation=45, ha="right")
    plt.ylabel("PPI (VPP_SCORE)")
    plt.title(f"PPI Score Distribution by Method Group — LC {lc}")
    plt.tight_layout()

    out_path = os.path.join(OUT_DIR, f"PPI_boxplot_LC_{lc}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"Wrote: {out_path}")


# ---------------- NEW: export best 100 settings per LC ----------------
TOP_DIR = "output/top100_by_lc"
os.makedirs(TOP_DIR, exist_ok=True)

meta_cols = [
    "settings_id",
    "method_name",
    "smooth",
    "seasonmethod",
    "seapar",
    "sos_cutoff",
    "eos_cutoff",
]

# sanity checks
missing_meta = [c for c in meta_cols if c not in df.columns]
if missing_meta:
    raise ValueError(f"Missing required metadata columns in merged df: {missing_meta}")

for lc in lc_cols:
    if lc not in df.columns:
        continue

    tmp = df[meta_cols + [lc]].copy()
    tmp = tmp.rename(columns={lc: "scorevalue"})

    # ensure numeric scores; drop missing
    tmp["scorevalue"] = pd.to_numeric(tmp["scorevalue"], errors="coerce")
    tmp = tmp.dropna(subset=["scorevalue"])

    # sort best-to-worst and take top 100
    tmp = tmp.sort_values("scorevalue", ascending=False).head(100)

    out_path = os.path.join(TOP_DIR, f"best100_settings_LC_{lc}.csv")
    tmp.to_csv(out_path, index=False)
    print(f"Wrote: {out_path}")
