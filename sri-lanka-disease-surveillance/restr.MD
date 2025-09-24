Here’s a structured, thorough data dictionary for your `model_df`. I grouped fields by function and gave unit/type + how they’re derived. For repeated patterns (lags/rolls/anoms/z), see the “Suffix conventions” section—so we don’t repeat the same explanation 80 times.

---
  
  # Data dictionary for `model_df`
  
  ## 1) Identifiers & time axes
  
  * **district** · *factor/chr*
  Administrative district name (normalized). Key for all group-by ops.

* **date\_start** · *Date*
  Start date of the epidemiological week parsed from WER.

* **date\_end** · *Date*
  End date of the epidemiological week parsed from WER.

* **date\_mid** · *Date*
  Midpoint of the epidemiological week (used to align climate & other features).

* **wx\_date** · *Date*
  Climate-alignment date used during joins (should match `date_mid` after merge).

* **date** · *Date*
  Canonical date column for modeling. If present, this equals `date_mid`.

* **week\_of\_year** · *int* (1–53)
Week index as parsed/derived from WER or from the date. (Week 53 is coerced to 52 in cyclic splines.)

* **week** · *int* (1–52)
Week index standardized for cyclic-seasonality models (53 → 52).

* **week\_id** · *int*
  Serial week index (project-specific utility column; monotonically increasing).

* **yday** · *int* (1–366)
Day-of-year (for optional seasonal work outside of ISO week).

* **year, year2merge** · *int*
  Calendar year; `year2merge` is a helper for mid-year population merges.

* **n\_days\_week** · *int* (≤7)
Count of days represented in the weekly aggregation (accounts for partial weeks at edges).

## 2) Population & land cover

* **poptot** · *numeric* (persons)
District total population used as exposure. Must be >0 for offsets.

* **BuiltUp, Cropland, Grass, Paddy, Shrub, Water, Wetland, NoData** · *numeric* (proportion 0–1)
Land-cover shares (areal proportions) per district from the LC raster aggregation.

## 3) Disease counts (WER) and derived rates

> In WER tables, each disease typically has two columns. We parse them as:
  >
  > * `*_A` = **weekly new cases** in the current epidemiological week
> * `*_B` = **cumulative cases** (YTD) as reported in that issue
>   (If any issue deviates from this convention, refer to that WER’s header legend.)

Weekly A/B for the following (subset shown):
  
  * **dengue\_A, dengue\_B**
  * **dysentery\_A, dysentery\_B**
  * **encephalitis\_A, encephalitis\_B**
  * **enteric\_fever\_A, enteric\_fever\_B**
  * **food\_poisoning\_A, food\_poisoning\_B**
  * **leptospirosis\_A, leptospirosis\_B**
  * **typhus\_f\_A, typhus\_f\_B**
  * **viral\_hep\_A, viral\_hep\_B**
  * **rabies\_A, rabies\_B**
  * **chickenpox\_A, chickenpox\_B**
  * **meningitis\_A, meningitis\_B**
  * **leishmania\_A, leishmania\_B**
  * **tuberculosis\_A, tuberculosis\_B**
  * **wrcd\_A, wrcd\_B** *(WER label as-is; keep code mapping in your parsing notes)*
  
  Convenience fields:
  
  * **lepto** · *int*
  Alias for weekly leptospirosis cases (= `leptospirosis_A`).

* **dengue** · *int*
  Alias for weekly dengue cases (= `dengue_A`).

Rates per 100k (computed when `poptot` available):
  
  * **lepto\_100k, dengue\_100k, …, wrcd\_rate\_100k** · *numeric* (cases per 100,000)
`1e5 * weekly_A / poptot` created for multiple diseases:
  `dengue_100k, dysentery_100k, encephalitis_100k, enteric_fever_100k, food_poisoning_100k, leptospirosis_100k, typhus_f_100k, viral_hep_100k, rabies_100k, chickenpox_100k, meningitis_100k, leishmania_100k, tuberculosis_100k, wrcd_100k`.

## 4) Station/Open-Meteo style daily/weekly weather (joined earlier)

*(These are your pre-ERA5 fields used in earlier versions; units follow Open-Meteo defaults.)*
  
  * **temperature\_2m\_max / min / mean** · *°C*
  * **apparent\_temperature\_max / min / mean** · *°C*
  * **shortwave\_radiation\_sum** · *MJ/m²* (daily/weekly sum)
* **precipitation\_sum, rain\_sum** · *mm*
  * **precipitation\_hours** · *hours* with precipitation > 0
* **windspeed\_10m\_max** · *m/s*
  * **windgusts\_10m\_max** · *m/s*
  * **et0\_fao\_evapotranspiration** · *mm* (reference ET0)

> Use these only if your current analysis includes station-derived features; your ERA5-only scripts typically exclude them.

## 5) ERA5 weekly features (district-area-weighted)

Weekly aggregates aligned to WER weeks (units align to the underlying ERA5 variable definitions used in your feature builder).

### Temperature & heat/moisture

* **tmax\_mean, tmin\_mean, tmean\_mean** · *°C* (weekly means of daily maxima/minima/means)
* **tmax\_p90, tmax\_p95** · *°C* (weekly 90th/95th percentile of daily max temp)
* **tmax\_range** · *°C* (weekly mean daily diurnal range)
* **wbgt\_mean\_week** · *°C WBGT* (weekly mean wet-bulb globe temperature, modelled index)
* **rh\_mean\_week** · *%* (weekly mean relative humidity)
* **vpd\_mean\_week** · *kPa* (weekly mean vapor pressure deficit)
* **ssrd\_MJ\_mean\_week** · *MJ/m²* (weekly mean of daily shortwave radiation; total daily ssrd converted to MJ/m² first)

### Precipitation (two sources where provided)

* **precip\_tp\_sum\_week** · *mm* (weekly sum of ERA5 *total precipitation* “tp”)
* **precip\_mtpr\_sum\_week** · *mm* (weekly sum of the alternate precipitation field used in your pipeline; treat as a second precip source for robustness/ensemble)

> Where both exist, you can compare or use both to stabilize signal.

### Wetness and heat event indicators

* **wet\_days\_ge10\_tp / wet\_days\_ge10\_mtpr** · *days/week* (count of days with daily precip ≥ 10 mm)
* **max3d\_tp / max3d\_mtpr** · *mm* (max 3-day running precip accumulation within the week)
* **wet\_spell\_maxlen\_tp / …\_mtpr** · *days* (longest consecutive wet spell in the week, threshold consistent with builder)
* **hot\_days\_ge32** · *days/week* (count of days with Tmax ≥ 32 °C)

### Exposure-weighted precipitation proxies

* **ewap\_tp, ewap\_mtpr** · *index* (dimensionless)
Exposure/extent-weighted precipitation proxy (constructed upstream; larger ⇒ broader/intense wetness). Use descriptively and compare against plain sums.

### Climatology baselines for a subset

* **precip\_tp\_sum\_week\_clim, tmax\_mean\_clim, rh\_mean\_week\_clim** · *same units as base*
  District × week-of-year long-term baseline (multi-year mean/median as defined in your feature builder). Used for anomalies and %-of-normal.

## 6) Lagged / rolling / anomaly / standardized variants

These families appear across many ERA5 variables. Interpret using the base variable’s unit.

* **`<var>_lagK`**
  Value of `<var>` from **K weeks earlier** (K ∈ {1,…,6}).
*Example:* `rh_mean_week_lag1` = RH mean last week.

* **`<var>_roll2w_mean`, `<var>_roll2w_sum`, `<var>_roll4w_mean`, `<var>_roll4w_sum`**
  Rolling **right-aligned** windows over the **previous 2 or 4 weeks**, excluding the current week unless your builder specified otherwise.
*Example:* `precip_tp_sum_week_roll4w_mean` = mean of the last 4 weekly totals.

* **`<var>_anom`**
  Anomaly relative to district × week-of-year climatology (`<var> - <var>_clim`). Units same as `<var>`.

* **`<var>_pct_normal`**
  Percent of normal vs climatology: `100 × <var> / <var>_clim`. Unit: % (watch division by very small climatology).

* **`<var>_z`**
  Standardized anomaly (z-score) vs district × WOY climatology distribution. Unitless.

> Tip: Don’t include several highly collinear forms (e.g., raw + lag1 + roll4 + anom + z). Start compact (1–3 per concept), then expand only if they add independent signal.

## 7) Outbreak baseline & labels (derived in your labeling pipeline)

* **lepto\_prev** · *numeric/int*
  Previous week’s leptospirosis count (shifted) used to build past-only baseline.

* **medW, madW** · *numeric*
  Rolling **median** and **MAD** of `lepto_prev` over the previous `W` weeks (right-aligned). Units: cases.

* **baseline** · *numeric*
  `pmax(0, medW + h_z * madW)` — robust expected level per district/week (past-only).

* **z** · *numeric*
  Standardized exceedance vs rolling baseline: `(lepto - medW) / (madW + eps)`.

* **pct** · *numeric*
  Percent exceedance vs rolling baseline: `(lepto - medW) / (medW + eps)`.

* **n\_hist** · *int*
  Number of prior weeks of history available at each row (for warm-up logic).

* **trigger** · *logical*
  Spike trigger fired this week: `(n_hist < W & lepto ≥ min_cases) OR [lepto ≥ min_cases & (z ≥ k_z OR pct ≥ p_pct)]`.

* **state** · *int* (0/1)
Outbreak state (sticky). Once `trigger` fires, `state` stays 1 until there are **H** consecutive weeks below `baseline`.

* **onset** · *int* (0/1)
First week of each run where `state` flips 0 → 1.

* **onset\_next** · *int* (0/1)
Next-week onset indicator (lead of `onset`). This is the **binomial target** if you train a “one-week-ahead” onset classifier.

## 8) ERA5 variables present in your list (quick reference)

Base names you currently have (and therefore may also appear with suffixes from §6):
  
  * Temperature family: `tmax_mean, tmin_mean, tmean_mean, tmax_p90, tmax_p95, tmax_range, tmax_mean_*`
* Humidity/heat: `rh_mean_week, vpd_mean_week, wbgt_mean_week, ssrd_MJ_mean_week`
* Precipitation: `precip_tp_sum_week, precip_mtpr_sum_week` (+ `wet_days_ge10_*`, `max3d_*`, `wet_spell_maxlen_*`)
* Events: `hot_days_ge32`
* Proxies: `ewap_tp, ewap_mtpr`
* Climatologies & standardized: `*_clim, *_anom, *_pct_normal, *_z`
* Lags & rolls: `*_lag1…6, *_roll2w_*, *_roll4w_*`

## 9) Notes & caveats

* **Units**: Follow your feature-builder script. All temp-like vars in °C; precip in mm; radiation in MJ/m²; RH in %; VPD in kPa. If any upstream transformation changed units, keep that note next to the builder.
* **tp vs mtpr**: You built two precip flavors for robustness. Treat both as mm/week after aggregation; you can model either or both but mind collinearity.
* **A/B disease columns**: “A=weekly, B=cumulative” is the standard WER convention; if a specific issue differs, prefer the issue’s legend.
* **Rates vs counts**: Use `*_100k` for cross-district comparability. Use counts + `offset(log(poptot))` for NB/Poisson models.
* **Week 53**: Standardize to 52 for cyclic splines (`knots = list(week = c(0.5, 52.5))`).
* **Collinearity**: Lags/rolls/anoms/z of the same base variable can be highly collinear. Screen and keep a compact set.

---
  
  If you want this rendered into a CSV/Markdown table automatically, say the word and I’ll drop a small script that reads `names(model_df)` and writes out a tidy table with column, type, units, and description you can edit collaboratively.
