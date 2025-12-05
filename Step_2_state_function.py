def make_best_vi_mask(
    in_csv,
    out_csv,
    alpha=0.05,
    id_cols=("site", "lc", "n"),
):
    """
    Build mask tables from ST_VIs_diag_ALL_sites.csv:
      1) Full mask (all VI columns)
      2) Mask excluding any VI column containing 'GPP'
      3) Mask including ONLY VI columns containing 'FAPAR'
    """

    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy.stats import norm

    in_csv = Path(in_csv)
    out_csv = Path(out_csv)

    out_csv_noGPP = out_csv.with_name(out_csv.stem + "_noGPP.csv")
    out_csv_FAPAR = out_csv.with_name(out_csv.stem + "_FAPARonly.csv")

    df = pd.read_csv(in_csv)

    # Identify VI columns
    all_vi_cols = [c for c in df.columns if c not in id_cols]
    vi_cols_noGPP = [c for c in all_vi_cols if "GPP" not in c]
    vi_cols_FAPAR = [c for c in all_vi_cols if "FAPAR" in c]

    # Small notice: If no FAPAR column exists, we still produce a valid CSV with only id columns.

    # ---------------------------
    # COMMON MASK CREATOR
    # ---------------------------
    def build_mask(df_in, vi_cols, out_path):
        mask_df = df_in.copy()

        # Initialize mask columns with 0
        for col in vi_cols:
            mask_df[col] = 0  

        # z critical for Fisher test
        z_crit = norm.ppf(1 - alpha / 2.0)

        def fisher_z(r):
            r = np.clip(r, -0.999999, 0.999999)
            return 0.5 * np.log((1 + r) / (1 - r))

        # ------------------------------------------------
        # Fill mask
        # ------------------------------------------------
        for idx, row in df_in.iterrows():
            n = row["n"]
            if not np.isfinite(n) or n <= 3:
                continue

            rhos = row[vi_cols].astype(float).to_numpy()
            if np.all(~np.isfinite(rhos)):
                continue

            try:
                best_idx = int(np.nanargmax(rhos))
            except ValueError:
                continue

            best_rho = rhos[best_idx]
            if not np.isfinite(best_rho):
                continue

            z_best = fisher_z(best_rho)
            z_all = fisher_z(rhos)

            se = np.sqrt(2.0 / (n - 3.0))
            z_diff = np.abs(z_all - z_best) / se

            not_sig = (z_diff < z_crit) & np.isfinite(rhos)

            for j, flag in enumerate(not_sig):
                if flag:
                    mask_df.at[idx, vi_cols[j]] = 1

        # ===========================================
        # Add SUM row
        # ===========================================
        result = mask_df.copy()

        sum_row = {}
        sum_row["site"] = "SUM"
        sum_row["lc"] = ""
        sum_row["n"] = len(result)

        for col in vi_cols:
            sum_row[col] = int(result[col].sum())

        result = pd.concat([result, pd.DataFrame([sum_row])], ignore_index=True)

        # ===========================================
        # Add LC-SUM rows
        # ===========================================
        lcs = df_in["lc"].unique()
        for lc in lcs:
            sub = result[result["lc"] == lc]
            if sub.empty:
                continue

            lc_row = {}
            lc_row["site"] = "LC-SUM"
            lc_row["lc"] = lc
            lc_row["n"] = len(sub)

            for col in vi_cols:
                lc_row[col] = int(sub[col].sum())

            result = pd.concat([result, pd.DataFrame([lc_row])], ignore_index=True)

        # write CSV
        result.to_csv(out_path, index=False)
        print(f"[INFO] wrote {out_path}")

    # ------------------------------------------------
    # Produce masks
    # ------------------------------------------------
    build_mask(df, all_vi_cols, out_csv)
    build_mask(df, vi_cols_noGPP, out_csv_noGPP)
    build_mask(df, vi_cols_FAPAR, out_csv_FAPAR)
