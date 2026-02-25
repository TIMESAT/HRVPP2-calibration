#!/usr/bin/env python3
"""
Create SOS/EOS RMSE matrices and site-wise near-optimal binary matrices
from a site-level VPP CSV.

Input (site-level CSV):
  - Must contain columns: site, Setting, SOSD_RMSE, EOSD_RMSE
  - Produced by your site-first summarization step:
      summarize_vpp_site_first.py  ->  PPI_GPP_VPP_SITE_ALL.csv (example)

Outputs (CSV; rows=site, cols=Setting):
  1) SOS_RMSE_matrix.csv
  2) EOS_RMSE_matrix.csv
  3) SOS_RMSE_binary.csv  (1 if RMSE <= best(site)*(1+tol), else 0)
  4) EOS_RMSE_binary.csv

Usage:
  python postprocess_rmse_matrix.py \
    --input PPI_GPP_VPP_SITE_ALL.csv \
    --out-dir VPP_RMSE_matrix \
    --tol 0.05
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def export_rmse_and_binary_matrices(
    site_csv: Path,
    out_dir: Path,
    tol: float,
    setting_min: int | None = None,
    setting_max: int | None = None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(site_csv)

    # --------------------------------------------------
    # NEW: restrict settings range if requested
    # --------------------------------------------------
    if setting_min is not None or setting_max is not None:
        # extract numeric index from "settings XXXX"
        setting_id = (
            df["Setting"]
            .astype(str)
            .str.replace("settings", "", regex=False)
            .str.strip()
            .astype(int)
        )

        if setting_min is not None:
            df = df[setting_id >= setting_min]
        if setting_max is not None:
            df = df[setting_id <= setting_max]

    required = {"site", "Setting", "SOSD_RMSE", "EOSD_RMSE"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {site_csv}: {missing}")

    # Pivot to matrices (rows=site, cols=Setting)
    sos_rmse = df.pivot(index="site", columns="Setting", values="SOSD_RMSE")
    eos_rmse = df.pivot(index="site", columns="Setting", values="EOSD_RMSE")

    sos_rmse.to_csv(out_dir / "SOS_RMSE_matrix.csv")
    eos_rmse.to_csv(out_dir / "EOS_RMSE_matrix.csv")

    # Site-wise near-optimal binaries (conservative NA handling: NA stays NA)
    sos_min = sos_rmse.min(axis=1, skipna=True)
    eos_min = eos_rmse.min(axis=1, skipna=True)

    sos_binary = sos_rmse.le(sos_min * (1.0 + tol), axis=0).astype("Int64")
    eos_binary = eos_rmse.le(eos_min * (1.0 + tol), axis=0).astype("Int64")

    sos_binary.to_csv(out_dir / "SOS_RMSE_binary.csv")
    eos_binary.to_csv(out_dir / "EOS_RMSE_binary.csv")

    print("Successfully written:")
    print(f"  {out_dir / 'SOS_RMSE_matrix.csv'}")
    print(f"  {out_dir / 'EOS_RMSE_matrix.csv'}")
    print(f"  {out_dir / 'SOS_RMSE_binary.csv'}")
    print(f"  {out_dir / 'EOS_RMSE_binary.csv'}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Post-process site-level VPP CSV into SOS/EOS RMSE and binary matrices."
    )
    parser.add_argument(
        "--input",
        required=False,
        default="PPI_GPP_VPP_SITE_ALL.csv",
        help="Site-level CSV (default: PPI_GPP_VPP_SITE_ALL.csv)",
    )
    parser.add_argument(
        "--out-dir",
        required=False,
        default="VPP_RMSE_matrix_SP10000",
        help="Output directory for RMSE/binary matrices (default: VPP_RMSE_matrix)",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=0.05,
        help="Relative tolerance for near-optimal binary (default: 0.05 = 5%%)",
    )
    parser.add_argument(
        "--setting-min",
        type=int,
        default=5001,
        help="Minimum setting index to process (e.g. 5001)",
    )
    parser.add_argument(
        "--setting-max",
        type=int,
        default=6000,
        help="Maximum setting index to process (e.g. 6000)",
    )

    args = parser.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_path.resolve()}")

    export_rmse_and_binary_matrices(
        site_csv=in_path,
        out_dir=Path(args.out_dir),
        tol=float(args.tol),
        setting_min=args.setting_min,
        setting_max=args.setting_max,
    )


if __name__ == "__main__":
    main()
