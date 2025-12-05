import pandas as pd
import matplotlib.pyplot as plt

# 1. Read the CSV file
df = pd.read_csv("ST_VIs_diag_ALL_sites.csv")

# 2. Select all columns starting from the 4th column (index 3)
set_cols = df.columns[3:]   # columns 4,5,6,... to the end

# 3. Convert these columns to numeric (coercing non-numeric values)
df[set_cols] = df[set_cols].apply(pd.to_numeric, errors='coerce')

# 4. Reshape into long format
long_df = df.melt(
    value_vars=set_cols,
    var_name="Setting",
    value_name="Value"
)

# 5. Drop rows that could not convert to numeric
long_df = long_df.dropna(subset=["Value"])

# 6. Create box plot
plt.figure(figsize=(12, 6))
long_df.boxplot(column="Value", by="Setting", grid=False)

plt.title("Box Plot for Spearman's rank correlation coefficient")
plt.suptitle("")  # remove pandas auto-title
plt.xlabel("Setting")
plt.ylabel("Spearman's rank correlation coefficient")
plt.xticks(rotation=45)
plt.tight_layout()

plt.show()
