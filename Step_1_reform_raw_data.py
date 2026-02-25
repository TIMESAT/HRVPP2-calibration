import pandas as pd

# Input and output paths
in_csv = "/projects/eko/fs7/proj/HRVPP2/calval/cal/HRVPP2-calibration/output/Cal/rawvi/Raw_VIs_diag_ALL_sites.csv"     # <-- change if needed
out_csv = "/projects/eko/fs7/proj/HRVPP2/calval/cal/HRVPP2-calibration/output/Cal/rawvi/Raw_VIs_diag_PPI_only.csv"

# Read CSV
df = pd.read_csv(in_csv)

# Check required columns
required_cols = ["VI", "n", "rho_gs", "site", "lc"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in input CSV: {missing}")

# Filter PPI rows and select columns
df_ppi = df.loc[df["VI"] == "PPI", ["n", "rho_gs", "site", "lc"]]

# Save new CSV
df_ppi.to_csv(out_csv, index=False)

print(f"[INFO] wrote {out_csv}")
