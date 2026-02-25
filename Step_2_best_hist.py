import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# Path to your CSV file
csv_path = "/projects/eko/fs7/proj/HRVPP2/calval/cal/HRVPP2-calibration/output/Cal/rawvi/ST_VIs_diag_ALL_sites.csv"   # <-- change this

# Read CSV
df = pd.read_csv(csv_path)

# Required columns
for col in ["PPI-5", "PPI-15"]:
    if col not in df.columns:
        raise ValueError(f"Column '{col}' not found in CSV")

# Extract data
ppi_5 = df["PPI-5"].dropna().values
ppi_15 = df["PPI-15"].dropna().values

# Define bins: 0–1 with width 0.05
bins = np.arange(0, 1.05, 0.05)

# Plot
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

axes[0].hist(ppi_5, bins=bins, edgecolor="black")
axes[0].set_title("PPI–GPP Spearman ρ (Setting 5)")
axes[0].set_xlabel("Spearman ρ (Growing Season)")
axes[0].set_ylabel("Frequency")
axes[0].set_xlim(0, 1)
axes[0].yaxis.set_major_locator(MaxNLocator(integer=True))

axes[1].hist(ppi_15, bins=bins, edgecolor="black")
axes[1].set_title("PPI–GPP Spearman ρ (Setting 15)")
axes[1].set_xlabel("Spearman ρ (Growing Season)")
axes[1].set_xlim(0, 1)
axes[1].yaxis.set_major_locator(MaxNLocator(integer=True))

plt.tight_layout()
# Save figure
out_png = "best_ST_Hist.png"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close(fig)

mean_ppi_5 = np.mean(ppi_5)
mean_ppi_15 = np.mean(ppi_15)




print(f"Average Spearman ρ (PPI-5):  {mean_ppi_5:.3f}")
print(f"Average Spearman ρ (PPI-15): {mean_ppi_15:.3f}")





# Boxplot
fig, ax = plt.subplots(figsize=(6, 4))

ax.boxplot(
    [ppi_5, ppi_15],
    labels=["HR-VPP1", "HR-VPP2"],
    showmeans=True,
    meanline=True
)

ax.set_ylabel("Spearman ρ (Growing Season)")
ax.set_title("Distribution of PPI–GPP Spearman ρ")
ax.set_ylim(0.3, 1)

plt.tight_layout()

# Save figure
out_box = "best_ST_Boxplot.png"
fig.savefig(out_box, dpi=300, bbox_inches="tight")
plt.close(fig)
