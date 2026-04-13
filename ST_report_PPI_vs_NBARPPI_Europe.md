# PPI vs NBAR PPI for Europe: calibration summary and recommended settings

## 1. Purpose

This note summarizes the current calibration results from the `Cal_Step2_ST_optimization.py` workflow after adding:

- raw comparisons among `PPI`, `NBAR_PPI`, and `GPP`
- TIMESAT-smoothed comparisons:
  - `GPP vs ST_PPI`
  - `GPP vs ST_NBAR_PPI`
  - `ST_PPI vs ST_NBAR_PPI`
- class-based best-setting selection for `PPI` and `NBAR_PPI`

The main practical question is:

> For Europe-wide production, should we prefer `PPI` or `NBAR_PPI`, and which TIMESAT setting should be used for each land-cover group?


## 2. Current grouping setup

The land-cover grouping used in the calibration is:

```json
{
  "Forest": {
    "landcover_codes": [7, 8, 16]
  },
  "Grassland": {
    "landcover_codes": [2, 5, 9, 10, 14, 15]
  },
  "Cropland": {
    "landcover_codes": [11, 12, 13]
  }
}
```

These codes come from WorldCover-based class labels attached to each site as `LCxx`.

Important interpretation note:

- the current workflow treats land cover as **static**
- the grouping is not time-varying
- a site assigned to `LC10`, for example, remains in the same class for the whole time series
- therefore, the grouping does **not move year by year**

This matters especially for `Grassland`, because a static 2021 class may not fully represent 2017-2024 management or surface-state variability.


## 3. What is being compared

The current ST comparison uses:

- `ST_PPI` or `ST_NBAR_PPI` after TIMESAT smoothing
- against **raw observed GPP**

So the comparison is:

- smoothed VI vs raw GPP

not:

- smoothed VI vs smoothed GPP

This is a reasonable design if the goal is to evaluate which VI better tracks observed GPP while keeping GPP as the reference target.


## 4. Main class-wise results

The current best settings by class are:

| Class | Best `PPI` setting | Mean rho_gs | Median rho_gs | n | Best `NBAR_PPI` setting | Mean rho_gs | Median rho_gs | n |
|---|---|---:|---:|---:|---|---:|---:|---:|
| Cropland | `ST_PPI-12` | 0.809603 | 0.810680 | 11 | `ST_NBAR_PPI-11` | 0.804655 | 0.808364 | 11 |
| Forest | `ST_PPI-15` | 0.794639 | 0.830255 | 28 | `ST_NBAR_PPI-15` | 0.810465 | 0.855738 | 28 |
| Grassland | `ST_PPI-15` | 0.806021 | 0.827918 | 29 | `ST_NBAR_PPI-15` | 0.802007 | 0.826100 | 29 |

### 4.1 Class-wise difference between best `NBAR_PPI` and best `PPI`

| Class | `NBAR_PPI - PPI` in mean rho_gs | `NBAR_PPI - PPI` in median rho_gs | Interpretation |
|---|---:|---:|---|
| Cropland | -0.004948 | -0.002316 | `PPI` slightly better |
| Forest | +0.015827 | +0.025483 | `NBAR_PPI` clearly better |
| Grassland | -0.004015 | -0.001818 | `PPI` slightly better |

### 4.2 First takeaway

- `Forest`: `NBAR_PPI` is meaningfully better than `PPI`
- `Cropland`: `PPI` is slightly better, but the difference is small
- `Grassland`: `PPI` is also slightly better, but again the difference is small

If one wants a single Europe-wide product with stronger physical consistency across time and illumination, `NBAR_PPI` remains a defensible choice even though it is not numerically best in every class.


## 5. Why `NBAR_PPI` is still a strong Europe-wide candidate

Even though `PPI` slightly outperforms `NBAR_PPI` in Cropland and Grassland in the current calibration, there are still good reasons to prefer `NBAR_PPI` for Europe-wide production:

1. `NBAR_PPI` is physically more standardized for broad-area application.
   It should be less sensitive to angular and illumination-related variation than raw `PPI`.

2. The strongest class-level advantage in the current calibration is for Forest, where `NBAR_PPI` is clearly better.

3. The disadvantages of `NBAR_PPI` in Cropland and Grassland are small in absolute magnitude.

4. For continental mapping, temporal consistency and transferability can matter as much as small calibration gains at tower sites.

In short:

- if the goal is **best tower-calibration score only**, a mixed strategy may favor `PPI` in some classes
- if the goal is **stable Europe-wide production**, `NBAR_PPI` is a strong default choice


## 6. Candidate class-based `NBAR_PPI` settings

A practical Europe-wide `NBAR_PPI` class-based scheme under discussion is:

- `Cropland -> ST_NBAR_PPI-12`
- `Grassland -> ST_NBAR_PPI-14`
- `Forest -> ST_NBAR_PPI-15`

This is not exactly the pure class-wise top-mean solution, because:

- Cropland best by mean is `ST_NBAR_PPI-11`
- Grassland best by mean is `ST_NBAR_PPI-15`
- Forest best by mean is `ST_NBAR_PPI-15`

However, the proposed scheme can still be defended as a **robust and conservative operational choice**.


## 7. Cropland: why `ST_NBAR_PPI-12` can be defended over `ST_NBAR_PPI-11`

Top 3 Cropland `NBAR_PPI` settings:

| Rank | Setting | Mean rho_gs | Median rho_gs | n |
|---|---|---:|---:|---:|
| 1 | `ST_NBAR_PPI-11` | 0.804655 | 0.808364 | 11 |
| 2 | `ST_NBAR_PPI-12` | 0.804241 | 0.814905 | 11 |
| 3 | `ST_NBAR_PPI-13` | 0.802081 | 0.816986 | 11 |

Direct comparison between `ST_NBAR_PPI-11` and `ST_NBAR_PPI-12`:

- number of Cropland site-lc pairs: `11`
- mean difference (`12 - 11`): `-0.000413`
- median difference (`12 - 11`): `+0.001027`
- station wins:
  - `12` wins at `6` sites
  - `11` wins at `5` sites

Interpretation:

- the difference in mean performance is extremely small
- `ST_NBAR_PPI-12` has a higher median rho_gs
- the site-by-site win count is almost balanced, with a slight edge to `12`
- `12` is therefore a plausible **robust** choice rather than a purely mean-maximizing choice

Parameter interpretation:

- `ST_NBAR_PPI-11`: `SP`, `smooth = 100`
- `ST_NBAR_PPI-12`: `SP`, `smooth = 300`

Why `12` is defendable:

1. It is nearly indistinguishable from `11` in mean performance.
2. It gives a slightly better median performance.
3. It is a bit smoother, which can be useful for managed agricultural systems with noisier temporal dynamics.
4. Cropland sample size is small, so a slightly more conservative choice is reasonable.

Recommended wording:

> `ST_NBAR_PPI-12` was selected for Cropland because its mean performance was nearly identical to the top-ranked `ST_NBAR_PPI-11`, while its median performance was slightly higher. Given the limited Cropland sample size and the higher temporal noisiness expected in managed agricultural systems, the slightly smoother setting was preferred as a more robust operational choice.


## 8. Grassland: why `ST_NBAR_PPI-14` can be defended over `ST_NBAR_PPI-15`

Top 3 Grassland `NBAR_PPI` settings:

| Rank | Setting | Mean rho_gs | Median rho_gs | n |
|---|---|---:|---:|---:|
| 1 | `ST_NBAR_PPI-15` | 0.802007 | 0.826100 | 29 |
| 2 | `ST_NBAR_PPI-3`  | 0.800653 | 0.835101 | 29 |
| 3 | `ST_NBAR_PPI-14` | 0.800266 | 0.818770 | 29 |

Direct comparison between `ST_NBAR_PPI-14` and `ST_NBAR_PPI-15`:

- number of Grassland site-lc pairs: `29`
- mean difference (`15 - 14`): `+0.001740`
- median difference (`15 - 14`): `+0.002203`
- station wins:
  - `15` wins at `20` sites
  - `14` wins at `9` sites

Interpretation:

- `15` is slightly better by the calibration metric
- but the difference is still small in absolute magnitude

Parameter interpretation:

- `ST_NBAR_PPI-14`: `SP`, `smooth = 3000`
- `ST_NBAR_PPI-15`: `SP`, `smooth = 10000`

Why `14` can still be defended:

1. Grassland is the most heterogeneous class in the current grouping.
   It combines `LC2`, `LC5`, `LC9`, `LC10`, `LC14`, and `LC15`.

2. The class definition is static, based on WorldCover labels rather than time-varying land cover.

3. Grassland systems can be strongly affected by mowing, grazing, wetness, and seasonal disturbance.

4. Under these conditions, a slightly less aggressive smoothing level can be argued to be more conservative and less likely to suppress real temporal variability.

So the defense for `14` is **not** that it outperforms `15` in the current calibration. It does not. The defense is:

- the difference is small
- Grassland class uncertainty and heterogeneity are relatively high
- a slightly less smoothed setting may be more robust for operational use

Recommended wording:

> Although `ST_NBAR_PPI-15` achieved the highest mean rho_gs for Grassland, the margin over `ST_NBAR_PPI-14` was small. Because the Grassland group is heterogeneous and based on a static WorldCover-derived class assignment, `ST_NBAR_PPI-14` can be justified as a slightly more conservative smoothing setting for large-scale application.


## 9. Forest: why `ST_NBAR_PPI-15` is the easiest choice

For Forest, the evidence is the cleanest:

- best `PPI`: `ST_PPI-15`, mean rho_gs `0.794639`
- best `NBAR_PPI`: `ST_NBAR_PPI-15`, mean rho_gs `0.810465`
- difference: `+0.015827` in favor of `NBAR_PPI`

This is the strongest class-wise advantage observed for `NBAR_PPI`.

Forest recommendation:

- `Forest -> ST_NBAR_PPI-15`

This choice is both numerically best and conceptually consistent with using NBAR-corrected input for broad-area monitoring.


## 10. Should the proposed class-based NBAR scheme be used for Europe?

### 10.1 Short answer

Yes, it is reasonable to use:

- `Cropland -> ST_NBAR_PPI-12`
- `Grassland -> ST_NBAR_PPI-14`
- `Forest -> ST_NBAR_PPI-15`

for Europe-wide production, as long as it is presented as a **pragmatic, robust, class-based operational choice** rather than a statistically final optimum.

### 10.2 Why this is reasonable

1. All three settings are among the best-performing options for their groups.
2. The differences from the nominal top settings are small.
3. The selected settings form a coherent smoothness gradient:
   - Cropland: moderate smoothing
   - Grassland: stronger smoothing
   - Forest: strongest smoothing
4. This gradient matches the intuition that canopy structure and temporal noise differ by class.

### 10.3 Why some caution is still needed

1. Cropland sample size is limited (`n = 11`).
2. The class assignment is static and based on WorldCover.
3. No formal paired significance test has yet been used to compare candidate settings within each class.


## 11. What would count as “statistically significant” here

There is no universal threshold such as “difference > 0.01 means significant”.

Statistical significance depends on:

- sample size
- site-to-site variability
- paired structure
- the chosen statistical test

For this application, the right way to test significance is:

1. within each class, compare two settings at the **same site-lc pairs**
2. compute paired differences in `rho_gs`
3. run a paired test such as:
   - Wilcoxon signed-rank test
   - or paired t-test

Until that is done, statements should remain cautious:

- “slightly better”
- “practically similar”
- “a more conservative choice”

rather than:

- “significantly better”


## 12. Recommended final interpretation

### Preferred scientific interpretation

- `NBAR_PPI` is the more suitable basis for Europe-wide application because it offers stronger physical consistency and clearly improves Forest performance.
- `PPI` remains slightly stronger in Cropland and Grassland in the current tower calibration, but these advantages are small.
- A class-based `NBAR_PPI` configuration is therefore a reasonable compromise between calibration performance and large-scale product consistency.

### Recommended operational class-based settings

- `Cropland -> ST_NBAR_PPI-12`
- `Grassland -> ST_NBAR_PPI-14`
- `Forest -> ST_NBAR_PPI-15`

### How to phrase this conservatively

> The selected class-based NBAR_PPI settings are not always the top-ranked options by mean rho_gs, but they remain among the best-performing candidates and provide a pragmatic, robust configuration for Europe-wide application. In Cropland and Grassland, the selected settings differ only marginally from the numerical optimum, while in Forest the selected setting is also the clear top performer.


## 13. Suggested next step

To strengthen the argument further, the next analysis step should be:

- paired within-class comparison between candidate settings
- especially:
  - `Cropland: ST_NBAR_PPI-11 vs ST_NBAR_PPI-12`
  - `Grassland: ST_NBAR_PPI-14 vs ST_NBAR_PPI-15`

The output should include:

- mean and median paired delta
- number of winning sites
- Wilcoxon p-value

That would allow the current practical recommendation to be upgraded into a statistically better-supported conclusion.
