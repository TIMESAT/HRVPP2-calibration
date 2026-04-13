# Step 3 Calibration Report

## Final selected calibration rule

This document records the final Step 3 calibration workflow and the final selected rule for extracting phenology from `ST_NBAR_PPI`.

The final version is:

- `Cropland`
  - `smooth = 300`
  - `seasonmethod = 1`
  - `seapar = 0.5`
- `Grassland`
  - `smooth = 3000`
  - `seasonmethod = 1`
  - `seapar = 1.0`
- `Forest`
  - `smooth = 10000`
  - `seasonmethod = 1`
  - `seapar = 1.0`

Within each class, `SOSD` and `EOSD` were optimized by searching all:

- `sos_cutoff = 0.05, 0.10, ..., 0.50`
- `eos_cutoff = 0.05, 0.10, ..., 0.50`

The final selected settings under the rule above are:

| Class | settings_id | smooth | seasonmethod | seapar | sos_cutoff | eos_cutoff |
|---|---|---:|---:|---:|---:|---:|
| Cropland | `NBAR_PPI_VPP-251` | 300 | 1 | 0.5 | 0.30 | 0.05 |
| Grassland | `NBAR_PPI_VPP-451` | 3000 | 1 | 1.0 | 0.30 | 0.05 |
| Forest | `NBAR_PPI_VPP-471` | 10000 | 1 | 1.0 | 0.40 | 0.05 |


## Workflow summary

### Inputs

Step 3 used:

- class-based selected `ST_NBAR_PPI` smoothing from Step 2
- `GPP` phenology as reference
- only `SOSD` and `EOSD` for the final optimization target

### Date handling

Special care was taken to make `SOSD` and `EOSD` comparable between different season methods.

- When `seasonmethod = 1`, TIMESAT day-like outputs were interpreted as `YYYYDOY`
- These values were converted into sequential day counts starting from `2017-01-01 = 1`
- Every year was treated as `365` days
- When `seasonmethod = 2`, day-like outputs were already treated as sequential day counts and were kept as they were

### Reference and candidate phenology

- `GPP` phenology was extracted as the reference phenology
- `ST_NBAR_PPI` phenology was extracted for all candidate parameter combinations
- Comparison was based on the date differences:
  - `selected_day - gpp_day`

### Filtering of mismatched seasons

Some seasons were clearly mismatched between `GPP` and `ST_NBAR_PPI`, with unrealistic `SOSD` or `EOSD` differences.

To avoid contaminating the final calibration, a season was kept only if:

- `|SOSD_selected - SOSD_GPP| <= 60 days`
- `|EOSD_selected - EOSD_GPP| <= 60 days`

Only seasons that passed both conditions were used in the final summary tables.


## Site coverage

Step 3 started from `71` site-LC groups found in the input directories.

- processed successfully: `68`
- missing / skipped: `3`

The skipped site-LC groups were:

- `FI-Sod LC8`
- `IT-Niv LC10`
- `IT-Noe LC15`

The final processed site-LC counts by class were:

| Class | n_site_lc |
|---|---:|
| Cropland | 11 |
| Grassland | 29 |
| Forest | 28 |


## Final class-level performance

The final best settings under the selected rule gave the following class-level combined performance:

| Class | settings_id | mean_abs_diff_days | median_abs_diff_days | rmsd_days | bias_days | n_obs | n_site_lc |
|---|---|---:|---:|---:|---:|---:|---:|
| Cropland | `NBAR_PPI_VPP-251` | 9.649 | 7.375 | 13.454 | -0.651 | 100 | 11 |
| Grassland | `NBAR_PPI_VPP-451` | 13.412 | 10.688 | 17.589 | -2.501 | 200 | 28 |
| Forest | `NBAR_PPI_VPP-471` | 16.232 | 12.812 | 21.144 | -3.773 | 280 | 28 |

Interpretation:

- `bias_days > 0`: `ST_NBAR_PPI` phenology is later than `GPP`
- `bias_days < 0`: `ST_NBAR_PPI` phenology is earlier than `GPP`
- `rmsd_days`: overall timing mismatch


## Final SOSD and EOSD performance

To better understand the behavior of each class, `SOSD` and `EOSD` were also assessed separately.

| Class | settings_id | Metric | bias_days | abs_bias_days | rmsd_days | mean_abs_diff_days | median_abs_diff_days | n_obs |
|---|---|---|---:|---:|---:|---:|---:|---:|
| Cropland | `NBAR_PPI_VPP-251` | SOSD | 2.737 | 2.737 | 16.394 | 12.248 | 10.375 | 50 |
| Cropland | `NBAR_PPI_VPP-251` | EOSD | -4.040 | 4.040 | 9.658 | 7.050 | 5.312 | 50 |
| Grassland | `NBAR_PPI_VPP-451` | SOSD | -2.016 | 2.016 | 13.929 | 10.306 | 8.188 | 100 |
| Grassland | `NBAR_PPI_VPP-451` | EOSD | -2.985 | 2.985 | 20.609 | 16.518 | 13.688 | 100 |
| Forest | `NBAR_PPI_VPP-471` | SOSD | 6.307 | 6.307 | 18.241 | 13.464 | 8.562 | 140 |
| Forest | `NBAR_PPI_VPP-471` | EOSD | -13.854 | 13.854 | 23.694 | 19.000 | 17.062 | 140 |

Main findings:

- `Cropland` performed best overall
- `Grassland` was acceptable but `EOSD` remained slightly unstable
- `Forest` remained the most difficult class, especially for `EOSD`
- the main challenge in `Forest` was not `SOSD`, but the early and unstable `EOSD`


## Valid and ignored seasons

After applying the mismatch filter, the final season retention was:

| Class | settings_id | total_seasons | valid_seasons | ignored_seasons | retention_pct |
|---|---|---:|---:|---:|---:|
| Cropland | `NBAR_PPI_VPP-251` | 61 | 50 | 11 | 81.967 |
| Grassland | `NBAR_PPI_VPP-451` | 119 | 100 | 19 | 84.034 |
| Forest | `NBAR_PPI_VPP-471` | 166 | 140 | 26 | 84.337 |

This means the final calibration was based on:

- `Cropland`: `50` valid seasons
- `Grassland`: `100` valid seasons
- `Forest`: `140` valid seasons


## Ignored season summary by site

Some sites contributed multiple ignored seasons, which indicates that their `GPP` and `ST_NBAR_PPI` season matching was unstable for the final selected rule.

### Cropland

| Site | LC | ignored_seasons |
|---|---|---:|
| DE-RuS | LC11 | 2 |
| FR-EM2 | LC11 | 2 |
| FR-Lam | LC11 | 2 |
| IT-Lsn | LC12 | 2 |
| BE-Lon | LC11 | 1 |
| DE-Hdn | LC11 | 1 |
| FR-Gri | LC11 | 1 |

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
| DK-Skj | LC2 | 1 |
| FR-Lus | LC10 | 1 |

### Forest

| Site | LC | ignored_seasons |
|---|---|---:|
| FR-FBn | LC8 | 4 |
| FR-Pue | LC16 | 4 |
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

The following site-season combinations were excluded from the final calibration because they failed the matching filter:

| Class | settings_id | Site | LC | Year | Season | SOSD diff (days) | EOSD diff (days) |
|---|---|---|---|---:|---|---:|---:|
| Cropland | `NBAR_PPI_VPP-251` | BE-Lon | LC11 | 2021 | s1 | 24.750 | -116.250 |
| Cropland | `NBAR_PPI_VPP-251` | DE-Hdn | LC11 | 2020 | s1 | -150.875 | 15.250 |
| Cropland | `NBAR_PPI_VPP-251` | DE-RuS | LC11 | 2017 | s1 | -64.875 | -71.250 |
| Cropland | `NBAR_PPI_VPP-251` | DE-RuS | LC11 | 2021 | s1 | -94.625 | -162.875 |
| Cropland | `NBAR_PPI_VPP-251` | FR-EM2 | LC11 | 2017 | s1 | -15.375 | -62.375 |
| Cropland | `NBAR_PPI_VPP-251` | FR-EM2 | LC11 | 2023 | s1 | -128.375 | -2.500 |
| Cropland | `NBAR_PPI_VPP-251` | FR-Gri | LC11 | 2023 | s1 | -63.375 | 0.125 |
| Cropland | `NBAR_PPI_VPP-251` | FR-Lam | LC11 | 2021 | s1 | -119.125 | -151.125 |
| Cropland | `NBAR_PPI_VPP-251` | FR-Lam | LC11 | 2023 | s1 | -126.625 | 5.500 |
| Cropland | `NBAR_PPI_VPP-251` | IT-Lsn | LC12 | 2018 | s1 | 11.625 | -65.250 |
| Cropland | `NBAR_PPI_VPP-251` | IT-Lsn | LC12 | 2022 | s1 | 0.875 | -75.125 |
| Grassland | `NBAR_PPI_VPP-451` | BE-Dor | LC10 | 2017 | s1 | -15.375 | -107.125 |
| Grassland | `NBAR_PPI_VPP-451` | BE-Dor | LC10 | 2018 | s1 | -17.750 | 115.375 |
| Grassland | `NBAR_PPI_VPP-451` | BE-Dor | LC10 | 2023 | s1 | 22.875 | -123.500 |
| Grassland | `NBAR_PPI_VPP-451` | CH-Cha | LC10 | 2018 | s1 | -18.500 | -79.000 |
| Grassland | `NBAR_PPI_VPP-451` | CH-Cha | LC10 | 2021 | s1 | -75.875 | 68.750 |
| Grassland | `NBAR_PPI_VPP-451` | CH-Oe2 | LC9 | 2017 | s1 | -168.625 | -5.125 |
| Grassland | `NBAR_PPI_VPP-451` | CH-Oe2 | LC9 | 2018 | s1 | -161.250 | -45.125 |
| Grassland | `NBAR_PPI_VPP-451` | DE-RuR | LC10 | 2020 | s1 | 4.000 | -89.500 |
| Grassland | `NBAR_PPI_VPP-451` | DE-RuR | LC10 | 2023 | s1 | -169.500 | -222.000 |
| Grassland | `NBAR_PPI_VPP-451` | DE-SfS | LC2 | 2020 | s1 | 20.250 | 136.500 |
| Grassland | `NBAR_PPI_VPP-451` | DK-Skj | LC2 | 2020 | s1 | -64.125 | -48.250 |
| Grassland | `NBAR_PPI_VPP-451` | ES-LMa | LC10 | 2022 | s1 | -88.500 | -6.500 |
| Grassland | `NBAR_PPI_VPP-451` | ES-LMa | LC10 | 2023 | s1 | 51.625 | -79.250 |
| Grassland | `NBAR_PPI_VPP-451` | FR-Lus | LC10 | 2023 | s1 | -11.250 | -65.625 |
| Grassland | `NBAR_PPI_VPP-451` | FR-Tou | LC9 | 2019 | s1 | -117.000 | -14.625 |
| Grassland | `NBAR_PPI_VPP-451` | FR-Tou | LC9 | 2020 | s1 | -85.750 | -3.375 |
| Grassland | `NBAR_PPI_VPP-451` | FR-Tou | LC9 | 2021 | s1 | -94.625 | -81.750 |
| Grassland | `NBAR_PPI_VPP-451` | FR-Tou | LC9 | 2022 | s1 | -156.750 | -22.625 |
| Grassland | `NBAR_PPI_VPP-451` | FR-Tou | LC9 | 2023 | s1 | -74.375 | -36.375 |
| Forest | `NBAR_PPI_VPP-471` | CH-Dav | LC8 | 2018 | s1 | 38.000 | -63.000 |
| Forest | `NBAR_PPI_VPP-471` | CZ-Lnz | LC7 | 2020 | s1 | -0.375 | -60.375 |
| Forest | `NBAR_PPI_VPP-471` | CZ-Lnz | LC7 | 2023 | s1 | 7.250 | -62.375 |
| Forest | `NBAR_PPI_VPP-471` | CZ-RAJ | LC8 | 2018 | s1 | 38.125 | -75.375 |
| Forest | `NBAR_PPI_VPP-471` | CZ-Stn | LC7 | 2020 | s1 | -2.875 | -68.375 |
| Forest | `NBAR_PPI_VPP-471` | DE-RuW | LC8 | 2021 | s1 | 89.250 | 110.375 |
| Forest | `NBAR_PPI_VPP-471` | FI-Hyy | LC8 | 2021 | s1 | -65.875 | -19.000 |
| Forest | `NBAR_PPI_VPP-471` | FR-Bil | LC8 | 2018 | s1 | -6.625 | 64.875 |
| Forest | `NBAR_PPI_VPP-471` | FR-Bil | LC8 | 2020 | s1 | 71.250 | 16.750 |
| Forest | `NBAR_PPI_VPP-471` | FR-FBn | LC8 | 2017 | s1 | 46.250 | 69.750 |
| Forest | `NBAR_PPI_VPP-471` | FR-FBn | LC8 | 2019 | s1 | 20.125 | 86.125 |
| Forest | `NBAR_PPI_VPP-471` | FR-FBn | LC8 | 2020 | s1 | 44.625 | 66.375 |
| Forest | `NBAR_PPI_VPP-471` | FR-FBn | LC8 | 2023 | s1 | 89.875 | 71.875 |
| Forest | `NBAR_PPI_VPP-471` | FR-Pue | LC16 | 2017 | s1 | 50.875 | 61.375 |
| Forest | `NBAR_PPI_VPP-471` | FR-Pue | LC16 | 2020 | s1 | 63.000 | -42.500 |
| Forest | `NBAR_PPI_VPP-471` | FR-Pue | LC16 | 2022 | s1 | 81.375 | 75.250 |
| Forest | `NBAR_PPI_VPP-471` | FR-Pue | LC16 | 2023 | s1 | 38.875 | 71.375 |
| Forest | `NBAR_PPI_VPP-471` | IT-BFt | LC7 | 2020 | s1 | 1.375 | -69.875 |
| Forest | `NBAR_PPI_VPP-471` | IT-BFt | LC7 | 2022 | s1 | -2.875 | -92.125 |
| Forest | `NBAR_PPI_VPP-471` | IT-Cp2 | LC16 | 2020 | s1 | 35.875 | -91.875 |
| Forest | `NBAR_PPI_VPP-471` | IT-Cp2 | LC16 | 2022 | s1 | 28.000 | -73.625 |
| Forest | `NBAR_PPI_VPP-471` | IT-Ren | LC8 | 2021 | s1 | -177.375 | -141.000 |
| Forest | `NBAR_PPI_VPP-471` | IT-Ren | LC8 | 2022 | s1 | 26.125 | -69.875 |
| Forest | `NBAR_PPI_VPP-471` | IT-Ren | LC8 | 2023 | s1 | -77.500 | -51.000 |
| Forest | `NBAR_PPI_VPP-471` | IT-TrF | LC8 | 2021 | s1 | -168.375 | -125.500 |
| Forest | `NBAR_PPI_VPP-471` | SE-Nor | LC8 | 2023 | s1 | -75.875 | 49.375 |


## Interpretation

### Cropland

Cropland gave the best overall performance among the three classes.

- combined `bias_days` was close to zero
- `EOSD` was more stable than `SOSD`
- `SOSD` still had moderate spread, but the overall result was acceptable

### Grassland

Grassland was intermediate.

- `SOSD` behaved reasonably well
- `EOSD` was still somewhat unstable
- the class remained usable, but with larger uncertainty than Cropland

### Forest

Forest remained the most difficult class.

- `SOSD` was not catastrophic, but still weaker than in Cropland and Grassland
- the main problem was `EOSD`
- `EOSD` was systematically early relative to GPP and had the largest spread

This is consistent with the visual and statistical impression that season end is intrinsically harder to define in forests, especially when using a common optical-season framework against GPP reference phenology.


## Final conclusion

The final selected Step 3 rule is:

- `Cropland`: `smooth=300`, `seasonmethod=1`, `seapar=0.5`, `SOS=0.30`, `EOS=0.05`
- `Grassland`: `smooth=3000`, `seasonmethod=1`, `seapar=1.0`, `SOS=0.30`, `EOS=0.05`
- `Forest`: `smooth=10000`, `seasonmethod=1`, `seapar=1.0`, `SOS=0.40`, `EOS=0.05`

The calibration was based on:

- `68` processed site-LC groups
- `50` valid Cropland seasons
- `100` valid Grassland seasons
- `140` valid Forest seasons

The ignored seasons were not removed manually; they were excluded using explicit mismatch thresholds against GPP reference phenology.

This final version should therefore be considered:

- reproducible
- traceable
- conservative with respect to obvious season mismatches

while still retaining the majority of seasons in all three classes.
