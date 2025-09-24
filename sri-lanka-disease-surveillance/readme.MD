# Sri Lanka Disease Surveillance ▸ Climate–Health Pipeline

## Overview

This repository builds an analysis-ready dataset linking Sri Lanka’s **Weekly Epidemiological Reports (WER)** to **climate, land cover, and demographic covariates**. It is part of the DGHI Climate & Health Initiative’s Southeast Asia observatory work.

The pipeline:

1. **Indexes and parses WER PDFs** to extract district-level weekly cases.
2. **Joins mid-year population estimates** and census data to compute per-capita rates.
3. **Processes station and ERA5 reanalysis data** into district-level daily and weekly features.
4. **Adds land cover proportions** to support ecological/epidemiological modeling.
5. **Outputs integrated health + climate panels** for downstream analysis (GAMs, GBMs, etc.).

> If you only want to explore/plot the prepared data, start at `code/preprocess/post_processing.R` after cloning.

---

## Quick start

### Requirements

* **R ≥ 4.3** with key packages:
  `data.table`, `stringr`, `lubridate`, `rvest`, `xml2`, `pdftools`,
  `sf`, `terra`, `exactextractr`, `arrow`, `ggplot2`, `scales`, `ragg`.
* For optional PDF table extraction: `tabulapdf` + `tabulizerjars` (requires Java + `rJava`).
* **GDAL/GEOS/PROJ** must be available for `sf/terra`.

  * macOS: `brew install gdal proj geos`.

### Environment variables (recommended)

Set in `~/.Renviron`:

```
CHI_LOCAL_WORK="C:/Users/<you>/Desktop/srilanka"
CHI_GITHUB_ROOT="C:/Users/<you>/R_Projects/dghi-southeast-asia"
JAVA_HOME="C:/Program Files/Eclipse Adoptium/jdk-17.0.16.8-hotspot"
```

---

## Running the pipeline

**Step 1. Initial ingest & extraction**

```r
source("code/preprocess/initial_processing.R")
```

* Crawls WER site, indexes PDFs, extracts weekly district counts.
* Saves:

  * `data/intermediate/sri_lanka_WER_index_of_pdfs.csv`
  * `data/intermediate/disease_counts.txt`

**Step 2. Post-processing & enrichment**

```r
source("code/preprocess/post_processing.R")
```

* Reads case counts, joins population, weather station aggregates, land cover, and ERA5 features.
* Outputs:

  * `data/intermediate/lep_analysis_panel.csv` (cases + features)
  * `data/processed/srilanka_district_weekly_era5_areawt.csv`
  * Figures under `reports/figures/`

**Step 3. Modeling (optional)**

* GAMs, GBMs, and exploratory visualization scripts in:

  * `code/analysis/gams_and_gbms.R`
  * `code/analysis/sri_lanka_modeling.R`

---

## Inputs

Located under `data/raw/`:

* `Mid-year_population_by_district_and_sex_2024.pdf`
* `population_by_district_in_census_years.csv`
* Weather station CSVs under `data/raw/station_data/`
* Land cover (GeoTIFF in `.rar`): `SriLanka_Landcover_2018`
* ERA5-derived district files produced by `code/preprocess/era5_to_weekly_features.R`

---

## Typical outputs

* `data/intermediate/disease_counts.txt` — long-format weekly counts by district × disease.
* `data/intermediate/sri_lanka_WER_index_of_pdfs.csv` — catalog of indexed WERs.
* `data/intermediate/lep_analysis_panel.csv` — merged health + features.
* `data/processed/srilanka_district_weekly_era5_areawt.csv` — climate features (means, sums, lags).
* `reports/figures/` — time series, seasonality, top districts.

---

## Repository structure

```
sri-lanka-disease-surveillance/
├─ code/
│  ├─ analysis/
│  │  ├─ gams_and_gbms.R
│  │  ├─ gams_bulk.R
│  │  ├─ initial_viz.R
│  │  └─ sri_lanka_modeling.R
│  ├─ helpers/
│  │  ├─ helpers.R
│  │  └─ plot_helpers.R
│  └─ preprocess/
│     ├─ era5_to_weekly_features.R
│     ├─ initial_processing.R
│     └─ post_processing.R
├─ data/
│  ├─ raw/
│  │  ├─ landcover/
│  │  │  └─ SriLanka_Landcover_2018.rar
│  │  ├─ Mid-year_population_by_district_and_sex_2024.pdf
│  │  ├─ population_by_district_in_census_years.csv
│  │  └─ station_data/
│  │     ├─ SriLanka_Weather_Dataset.csv
│  │     ├─ archive.zip
│  │     └─ readme.txt
│  ├─ intermediate/
│  │  ├─ disease_counts.txt
│  │  ├─ lep_analysis_panel.csv
│  │  └─ sri_lanka_WER_index_of_pdfs.csv
│  └─ processed/
│     ├─ srilanka_district_weekly_era5_areawt.csv
│     └─ WER_leptospirosis_counts.csv
├─ reports/
│  ├─ figures/
│  │  ├─ P1_national_trend.png
│  │  ├─ P2_top15_latest_week.png
│  │  ├─ P3_top6_last52_facets.png
│  │  ├─ P4_seasonality_boxplots.png
│  │  └─ P6_yoy_change_top15.png
│  └─ tables/
│     └─ temp.txt
├─ shared/helpers.R
├─ params.yaml
├─ readme.MD
└─ sri-lanka-disease-surveillance.Rproj
```

---

## Troubleshooting

* **Java/rJava/tabulizer**: If “JAVA\_HOME cannot be determined…”, ensure a JDK is installed and `JAVA_HOME` is set.
* **GDAL/GEOS/PROJ issues**: Verify CRS definitions with `st_crs()`.
* **Large ERA5 inputs**: Use `arrow` or `disk.frame` for efficiency; ensure enough memory for joins.
