# Step 3 Calibration Report

## Final fixed rule

This document records the final selected Step 3 calibration rule for `ST_NBAR_PPI` phenology extraction.

The final rule is **fixed**, not re-optimized by class after this step.

Important note:

- `output/Cal/Step3_ST_NBARPPI_VPP/Selected_ST_NBARPPI_VPP_SOSD_EOSD_best_settings_by_class.csv` is the automatic best-result table from the unrestricted search
- it is **not** the final adopted product rule
- the final adopted rule in this report corresponds to the fixed settings summarized in:
  `output/Cal/Step3_ST_NBARPPI_VPP/Selected_ST_NBARPPI_VPP_fixed_rule_summary.csv`

The final fixed constraints are:

- `Cropland`
  - `smooth = 300`
  - `seasonmethod = 1`
  - `seapar = 0.5`
  - `SOS cutoff = 0.35`
  - `EOS cutoff = 0.15`

- `Grassland`
  - `smooth = 3000`
  - `seasonmethod = 1`
  - `seapar = 1.0`
  - `SOS cutoff = 0.35`
  - `EOS cutoff = 0.15`

- `Forest`
  - `smooth = 10000`
  - `seasonmethod = 1`
  - `seapar = 1.0`
  - `SOS cutoff = 0.35`
  - `EOS cutoff = 0.15`


## Step 3 objective

The purpose of Step 3 was to extract phenology from the class-specific selected `ST_NBAR_PPI` signal and compare it to a `GPP` reference phenology.

Only two phenology dates were used for the final optimization and evaluation:

- `SOSD`
- `EOSD`


## Input background

### Smoothing settings inherited from Step 2

The Step 2 class-based smoothing settings were fixed before Step 3:

- `Cropland -> ST_NBAR_PPI-12 -> SP, smooth=300`
- `Grassland -> ST_NBAR_PPI-14 -> SP, smooth=3000`
- `Forest -> ST_NBAR_PPI-15 -> SP, smooth=10000`

These smoothing settings were not re-optimized in the final Step 3 rule.

### Site coverage

Step 3 started from `71` site-LC groups found in the inputs.

Processed successfully:

- `68` site-LC groups

Skipped due to missing matching information:

- `FI-Sod LC8`
- `IT-Niv LC10`
- `IT-Noe LC15`

Processed site-LC counts by class:

| Class | n_site_lc |
|---|---:|
| Cropland | 11 |
| Grassland | 29 |
| Forest | 28 |

After the season-mismatch filter was applied under the final fixed rule, the number of site-LC groups contributing at least one valid matched season became:

| Class | n_site_lc_with_valid_season |
|---|---:|
| Cropland | 11 |
| Grassland | 28 |
| Forest | 28 |

This means one processed grassland site-LC group did not contribute any valid retained season under the final fixed rule.


## Phenology extraction logic

### Reference phenology

`GPP` phenology was extracted as the reference phenology.

### Candidate phenology

`ST_NBAR_PPI` phenology was extracted from the smoothed time series using TIMESAT.

### Date conversion rule

The output date interpretation depended on `seasonmethod`.

When `seasonmethod = 1`:

- day-like variables such as `SOSD` and `EOSD` were treated as `YYYYDOY`
- these were converted to a continuous day count
- `2017-01-01 = 1`
- every year was treated as `365` days

When `seasonmethod = 2`:

- `SOSD`, `EOSD`, and `MAXD` were already treated as sequential day counts from the first day of the first year
- no extra conversion was applied

The final selected fixed rule uses only:

- `seasonmethod = 1`

So the final reported calibration is based on the `YYYYDOY -> sequential day` conversion.


## Season mismatch filtering

Not all site-year seasons were retained in the final calibration.

To remove obviously mismatched phenology seasons, a season was kept only if:

- `|SOSD_selected - SOSD_GPP| <= 60 days`
- `|EOSD_selected - EOSD_GPP| <= 60 days`

Only seasons passing both conditions were used in the final statistics.

This filtering was applied automatically and consistently across classes.


## Final selected fixed settings

Under the final fixed rule, the corresponding settings in the Step 3 search table are:

| Class | settings_id | smooth | seasonmethod | seapar | sos_cutoff | eos_cutoff |
|---|---|---:|---:|---:|---:|---:|
| Cropland | `NBAR_PPI_VPP-263` | 300 | 1 | 0.5 | 0.35 | 0.15 |
| Grassland | `NBAR_PPI_VPP-463` | 3000 | 1 | 1.0 | 0.35 | 0.15 |
| Forest | `NBAR_PPI_VPP-463` | 10000 | 1 | 1.0 | 0.35 | 0.15 |


## Final combined performance

The final combined `SOSD + EOSD` performance for the fixed rule is:

| Class | settings_id | mean_abs_diff_days | median_abs_diff_days | rmsd_days | bias_days | n_obs | n_site_lc |
|---|---|---:|---:|---:|---:|---:|---:|
| Cropland | `NBAR_PPI_VPP-263` | 10.480 | 7.500 | 14.163 | 0.941 | 102 | 11 |
| Grassland | `NBAR_PPI_VPP-463` | 13.921 | 10.625 | 18.218 | -1.021 | 202 | 28 |
| Forest | `NBAR_PPI_VPP-463` | 16.593 | 13.375 | 21.649 | -6.991 | 282 | 28 |

Interpretation:

- `bias_days > 0`: `ST_NBAR_PPI` phenology is later than `GPP`
- `bias_days < 0`: `ST_NBAR_PPI` phenology is earlier than `GPP`
- `rmsd_days`: overall day mismatch relative to `GPP`


## Final SOSD and EOSD statistics

The final fixed rule was also evaluated separately for `SOSD` and `EOSD`.

| Class | settings_id | Metric | bias_days | abs_bias_days | rmsd_days | mean_abs_diff_days | median_abs_diff_days | n_obs |
|---|---|---|---:|---:|---:|---:|---:|---:|
| Cropland | `NBAR_PPI_VPP-263` | SOSD | 5.708 | 5.708 | 17.643 | 14.179 | 10.625 | 51 |
| Cropland | `NBAR_PPI_VPP-263` | EOSD | -3.826 | 3.826 | 9.481 | 6.782 | 4.875 | 51 |
| Grassland | `NBAR_PPI_VPP-463` | SOSD | 1.391 | 1.391 | 14.841 | 11.010 | 9.250 | 101 |
| Grassland | `NBAR_PPI_VPP-463` | EOSD | -3.433 | 3.433 | 21.061 | 16.832 | 13.875 | 101 |
| Forest | `NBAR_PPI_VPP-463` | SOSD | 0.075 | 0.075 | 19.160 | 14.019 | 10.250 | 141 |
| Forest | `NBAR_PPI_VPP-463` | EOSD | -14.057 | 14.057 | 23.879 | 19.167 | 17.250 | 141 |

Main interpretation:

- `Cropland` remained the best-performing class overall
- `Grassland` was acceptable, although `EOSD` remained relatively noisy
- `Forest` remained the most difficult class
- the main remaining problem in `Forest` was still `EOSD`


## Valid and ignored seasons

After the mismatch filter, the final valid season counts were:

| Class | settings_id | total_seasons | valid_seasons | ignored_seasons | retention_pct |
|---|---|---:|---:|---:|---:|
| Cropland | `NBAR_PPI_VPP-263` | 60 | 51 | 9 | 85.000 |
| Grassland | `NBAR_PPI_VPP-463` | 119 | 101 | 18 | 84.874 |
| Forest | `NBAR_PPI_VPP-463` | 166 | 141 | 25 | 84.940 |

The final valid seasons actually used in the calibration were therefore:

- `Cropland`: `51`
- `Grassland`: `101`
- `Forest`: `141`


## Ignored seasons by site

### Cropland

| Site | LC | ignored_seasons |
|---|---|---:|
| DE-RuS | LC11 | 2 |
| FR-Lam | LC11 | 2 |
| IT-Lsn | LC12 | 2 |
| BE-Lon | LC11 | 1 |
| DE-Hdn | LC11 | 1 |
| FR-EM2 | LC11 | 1 |

### Grassland

| Site | LC | ignored_seasons |
|---|---|---:|
| FR-Tou | LC9 | 5 |
| BE-Dor | LC10 | 3 |
| CH-Cha | LC10 | 2 |
| CH-Oe2 | LC9 | 2 |
| DE-RuR | LC10 | 2 |
| ES-LMa | LC10 | 2 |
| DE-SfS | LC2 | 1 |
| FR-Lus | LC10 | 1 |

### Forest

| Site | LC | ignored_seasons |
|---|---|---:|
| FR-FBn | LC8 | 4 |
| FR-Pue | LC16 | 3 |
| IT-Ren | LC8 | 3 |
| CZ-Lnz | LC7 | 2 |
| FR-Bil | LC8 | 2 |
| IT-BFt | LC7 | 2 |
| IT-Cp2 | LC16 | 2 |
| CH-Dav | LC8 | 1 |
| CZ-RAJ | LC8 | 1 |
| CZ-Stn | LC7 | 1 |
| DE-RuW | LC8 | 1 |
| FI-Hyy | LC8 | 1 |
| IT-TrF | LC8 | 1 |
| SE-Nor | LC8 | 1 |


## Full ignored season list

The following seasons were excluded by the final mismatch filter.

| Class | settings_id | Site | LC | Year | Season | SOSD diff (days) | EOSD diff (days) |
|---|---|---|---|---:|---|---:|---:|
| Cropland | `NBAR_PPI_VPP-263` | BE-Lon | LC11 | 2021 | s1 | 25.500 | -116.250 |
| Cropland | `NBAR_PPI_VPP-263` | DE-Hdn | LC11 | 2020 | s1 | -147.875 | 15.250 |
| Cropland | `NBAR_PPI_VPP-263` | DE-RuS | LC11 | 2017 | s1 | -61.875 | -71.250 |
| Cropland | `NBAR_PPI_VPP-263` | DE-RuS | LC11 | 2021 | s1 | -93.625 | -162.875 |
| Cropland | `NBAR_PPI_VPP-263` | FR-EM2 | LC11 | 2017 | s1 | -12.875 | -62.375 |
| Cropland | `NBAR_PPI_VPP-263` | FR-Lam | LC11 | 2021 | s1 | -115.375 | -151.125 |
| Cropland | `NBAR_PPI_VPP-263` | FR-Lam | LC11 | 2023 | s1 | -120.875 | 5.500 |
| Cropland | `NBAR_PPI_VPP-263` | IT-Lsn | LC12 | 2018 | s1 | 16.500 | -65.250 |
| Cropland | `NBAR_PPI_VPP-263` | IT-Lsn | LC12 | 2022 | s1 | 5.250 | -75.125 |
| Grassland | `NBAR_PPI_VPP-463` | BE-Dor | LC10 | 2017 | s1 | -13.000 | -107.125 |
| Grassland | `NBAR_PPI_VPP-463` | BE-Dor | LC10 | 2018 | s1 | -1.625 | 115.375 |
| Grassland | `NBAR_PPI_VPP-463` | BE-Dor | LC10 | 2023 | s1 | 25.625 | -123.500 |
| Grassland | `NBAR_PPI_VPP-463` | CH-Cha | LC10 | 2018 | s1 | -16.000 | -79.000 |
| Grassland | `NBAR_PPI_VPP-463` | CH-Cha | LC10 | 2021 | s1 | -73.250 | 68.750 |
| Grassland | `NBAR_PPI_VPP-463` | CH-Oe2 | LC9 | 2017 | s1 | -166.375 | -5.125 |
| Grassland | `NBAR_PPI_VPP-463` | CH-Oe2 | LC9 | 2018 | s1 | -156.875 | -45.125 |
| Grassland | `NBAR_PPI_VPP-463` | DE-RuR | LC10 | 2020 | s1 | 6.250 | -89.500 |
| Grassland | `NBAR_PPI_VPP-463` | DE-RuR | LC10 | 2023 | s1 | -164.625 | -222.000 |
| Grassland | `NBAR_PPI_VPP-463` | DE-SfS | LC2 | 2020 | s1 | 35.500 | 136.500 |
| Grassland | `NBAR_PPI_VPP-463` | ES-LMa | LC10 | 2022 | s1 | -79.750 | -6.500 |
| Grassland | `NBAR_PPI_VPP-463` | ES-LMa | LC10 | 2023 | s1 | 54.625 | -79.250 |
| Grassland | `NBAR_PPI_VPP-463` | FR-Lus | LC10 | 2023 | s1 | -8.750 | -65.625 |
| Grassland | `NBAR_PPI_VPP-463` | FR-Tou | LC9 | 2019 | s1 | -68.875 | -14.625 |
| Grassland | `NBAR_PPI_VPP-463` | FR-Tou | LC9 | 2020 | s1 | -75.500 | -3.375 |
| Grassland | `NBAR_PPI_VPP-463` | FR-Tou | LC9 | 2021 | s1 | -86.875 | -81.750 |
| Grassland | `NBAR_PPI_VPP-463` | FR-Tou | LC9 | 2022 | s1 | -153.000 | -22.625 |
| Grassland | `NBAR_PPI_VPP-463` | FR-Tou | LC9 | 2023 | s1 | -68.750 | -36.375 |
| Forest | `NBAR_PPI_VPP-463` | CH-Dav | LC8 | 2018 | s1 | 35.125 | -63.000 |
| Forest | `NBAR_PPI_VPP-463` | CZ-Lnz | LC7 | 2020 | s1 | -2.500 | -60.375 |
| Forest | `NBAR_PPI_VPP-463` | CZ-Lnz | LC7 | 2023 | s1 | 4.250 | -62.375 |
| Forest | `NBAR_PPI_VPP-463` | CZ-RAJ | LC8 | 2018 | s1 | 34.000 | -75.375 |
| Forest | `NBAR_PPI_VPP-463` | CZ-Stn | LC7 | 2020 | s1 | -5.125 | -68.375 |
| Forest | `NBAR_PPI_VPP-463` | DE-RuW | LC8 | 2021 | s1 | 80.000 | 110.375 |
| Forest | `NBAR_PPI_VPP-463` | FI-Hyy | LC8 | 2021 | s1 | -73.500 | -19.000 |
| Forest | `NBAR_PPI_VPP-463` | FR-Bil | LC8 | 2018 | s1 | -13.250 | 64.875 |
| Forest | `NBAR_PPI_VPP-463` | FR-Bil | LC8 | 2020 | s1 | 63.750 | 16.750 |
| Forest | `NBAR_PPI_VPP-463` | FR-FBn | LC8 | 2017 | s1 | 39.250 | 69.750 |
| Forest | `NBAR_PPI_VPP-463` | FR-FBn | LC8 | 2019 | s1 | 14.625 | 86.125 |
| Forest | `NBAR_PPI_VPP-463` | FR-FBn | LC8 | 2020 | s1 | 37.750 | 66.375 |
| Forest | `NBAR_PPI_VPP-463` | FR-FBn | LC8 | 2023 | s1 | 83.500 | 71.875 |
| Forest | `NBAR_PPI_VPP-463` | FR-Pue | LC16 | 2017 | s1 | 45.750 | 61.375 |
| Forest | `NBAR_PPI_VPP-463` | FR-Pue | LC16 | 2022 | s1 | 76.875 | 75.250 |
| Forest | `NBAR_PPI_VPP-463` | IT-BFt | LC7 | 2020 | s1 | -1.125 | -69.875 |
| Forest | `NBAR_PPI_VPP-463` | IT-BFt | LC7 | 2022 | s1 | -5.125 | -92.125 |
| Forest | `NBAR_PPI_VPP-463` | IT-Cp2 | LC16 | 2020 | s1 | 29.500 | -91.875 |
| Forest | `NBAR_PPI_VPP-463` | IT-Cp2 | LC16 | 2022 | s1 | 23.500 | -73.625 |
| Forest | `NBAR_PPI_VPP-463` | IT-Ren | LC8 | 2021 | s1 | -183.375 | -141.000 |
| Forest | `NBAR_PPI_VPP-463` | IT-Ren | LC8 | 2022 | s1 | 21.250 | -69.875 |
| Forest | `NBAR_PPI_VPP-463` | IT-Ren | LC8 | 2023 | s1 | -80.875 | -51.000 |
| Forest | `NBAR_PPI_VPP-463` | IT-TrF | LC8 | 2021 | s1 | -171.500 | -125.500 |
| Forest | `NBAR_PPI_VPP-463` | SE-Nor | LC8 | 2023 | s1 | -82.125 | 49.375 |


## Interpretation

### Cropland

Cropland remained the strongest class under the fixed final rule.

- overall bias was close to zero
- `EOSD` was clearly better constrained than `SOSD`
- `SOSD` still showed some spread, but the class was overall acceptable

### Grassland

Grassland remained intermediate.

- `SOSD` behaved reasonably well
- `EOSD` was still noisier than desired
- the class remained usable, but with higher uncertainty than Cropland

### Forest

Forest remained the most difficult class.

- `SOSD` was not strongly biased, but still had relatively high spread
- the dominant problem was `EOSD`
- `EOSD` stayed too early relative to GPP and remained the weakest component of the final calibration


## Final conclusion

The final selected Step 3 rule is:

- `Cropland`: `smooth=300`, `seasonmethod=1`, `seapar=0.5`, `SOS=0.35`, `EOS=0.15`
- `Grassland`: `smooth=3000`, `seasonmethod=1`, `seapar=1.0`, `SOS=0.35`, `EOS=0.15`
- `Forest`: `smooth=10000`, `seasonmethod=1`, `seapar=1.0`, `SOS=0.35`, `EOS=0.15`

The final calibration was based on:

- `68` processed site-LC groups
- `51` valid Cropland seasons
- `101` valid Grassland seasons
- `141` valid Forest seasons

The ignored seasons were filtered automatically using explicit mismatch thresholds against the `GPP` reference phenology.

This final version should be regarded as:

- fixed
- reproducible
- transparent
- conservative with respect to obvious season mismatches

and suitable as the final selected calibration rule for the current workflow.
