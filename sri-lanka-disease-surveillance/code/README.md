# code/

## Quick start (most people only need this)

1. **Open the project** by double-clicking
   `sri-lanka-disease-surveillance.rproj` (in the repo root).

2. **Open & run:** `code/analysis/sri_lanka_modeling.R`

What that script does (accurately, per the header you shared):

* Loads a **prepared, analysis-ready weekly panel** that already links
  WER district×week disease counts, population, land-cover, and **ERA5-derived weekly climate features** (including any lags/rolls/anoms/z you built upstream).
* **Discovers available ERA5 columns by regex** (no recomputation in this script) and builds predictor blocks by selecting only from those existing columns.
* Creates **per-100k outcomes** from counts if needed.
* Makes a **chronological 80/20 split** (train past → test future; no leakage).
* Provides a compact **GAM/BAM harness**:

  * annual seasonality via `s(week, bs="cc")`,
  * optional district random effects `s(district, bs="re")`,
  * smooths only the ERA5 variables you pick.
* Prints **hold-out metrics** (RMSE/MAE/R²) and **lists the exact variables used**.
* You edit the `smooth_vars` list near the bottom to try different ERA5 drivers.

> This script assumes the panel CSV exists at
> `cfg$paths$intermediate/sri_lanka-disease-landcover-climate-2015_2024.csv`.
> (Define `cfg$paths` in your `.Rprofile`.)

---

## If you need to rebuild inputs (optional)

Only touch these if you’re regenerating the panel or changing upstream features.

* `code/preprocess/initial_processing.R`
  Parse WER PDFs → clean district names (e.g., “Kalmune” → “Kalmunai”), fold **Kalmunai into Ampara** per week, derive robust `date_start/end/mid`, and write a weekly disease table.

* `code/preprocess/era5_to_weekly_features.R`
  Build **weekly district-level ERA5 features** (means/sums/etc.). If you change variables, lags, or roll windows, re-run this and regenerate the ERA5 weekly dataset.

* `code/preprocess/post_processing.R`
  Merge the disease table + land cover + **ERA5 weekly** into the final **analysis panel** (the CSV the modeling script reads). Also filters years and normalizes date fields.

* `code/initial_processing.R`
  Legacy fallback of the initial parsing step—prefer the `preprocess/` version.

---

## Expected data & config

* **Config:** your `.Rprofile` should define `cfg$paths` (e.g., `raw`, `intermediate`, `processed`, `reports`).
* **Final panel CSV (consumed by the modeling script):**
  `cfg$paths$intermediate/sri_lanka-disease-landcover-climate-2015_2024.csv`
  Includes: `district`, `date` (or `date_mid`), `week`, `year`, `poptot`, **ERA5 weekly features** (and *\_lag*/ *\_roll*/ *\_anom*/ \*\_z if you made them), and optional land-cover proportions.

---

## Typical run order (only if rebuilding)

1. `code/preprocess/initial_processing.R`  → weekly disease table
2. `code/preprocess/era5_to_weekly_features.R` → ERA5 weekly features
3. `code/preprocess/post_processing.R` → final analysis panel CSV
4. `code/analysis/sri_lanka_modeling.R` → explore & model

---

## Notes

* Week 53 is coerced to 52 for cyclic splines.
* Keep smoother sets small (3–6 ERA5 variables) to avoid collinearity across lags/rolls.
* The modeling script **does not** recompute weather features; it only selects among the ERA5 columns already present in the panel.

---

