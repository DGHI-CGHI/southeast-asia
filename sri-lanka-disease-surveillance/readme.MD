# Sri Lanka Disease Surveillance ▸ Climate–Health Pipeline

## Overview

This repository creates an **analysis-ready panel** linking Sri Lanka’s **Weekly Epidemiological Reports (WER)** to climate, population, and land cover features. It supports DGHI CHI observatory analyses of leptospirosis, dengue, and other notifiable diseases.

The pipeline includes three core stages:

1. **Initial processing** — Crawl and parse WER PDFs into district-week case counts.
2. **Post-processing** — Merge disease counts with population, station weather, land cover, and ERA5 daily aggregates.
3. **ERA5 features** — Generate daily and weekly district-level climate covariates with lags, rolling windows, anomalies, and derived metrics.

Outputs are formatted for direct use in GAMs, GBMs, and other modeling frameworks.

---

## Key scripts

* `code/preprocess/initial_processing.R`

  * Crawls WER site, indexes PDFs, and extracts weekly district disease counts.
  * Produces:

    * `data/intermediate/sri_lanka_WER_index_of_pdfs.csv`
    * `data/intermediate/disease_counts.txt`

* `code/preprocess/post_processing.R`

  * Consumes WER outputs and joins:

    * Mid-year population (from PDF or CSV fallback)
    * Station weather mapped to districts
    * Land cover proportions (2018 raster)
    * ERA5 district-daily aggregates
  * Produces:

    * `data/intermediate/lep_analysis_panel.csv`
    * Figures under `reports/figures/`

* `code/preprocess/era5_to_weekly_features.R`

  * Converts ERA5 hourly Parquet (S3 or local DVC cache) → daily district means/extremes.
  * Builds epidemiology-ready **weekly covariates** aligned to WER weeks.
  * Produces:

    * `data/processed/srilanka_district_daily_era5_areawt.csv`
    * `data/processed/srilanka_district_weekly_era5_areawt.csv`

---

## Data inputs

* `data/raw/Mid-year_population_by_district_and_sex_2024.pdf`
* `data/raw/population_by_district_in_census_years.csv`
* `data/raw/station_data/SriLanka_Weather_Dataset.csv`
* `data/raw/landcover/SriLanka_Landcover_2018.rar` → contains GeoTIFF
* ERA5 hourly Parquet files:

  * `s3://dghi-chi/data/se-asia/sri-lanka-disease-surveillance/era5/{year}/*.parquet`

---

## Typical outputs

* **Case counts**

  * `data/intermediate/disease_counts.txt`
  * `data/intermediate/sri_lanka_WER_index_of_pdfs.csv`

* **Merged analysis panel**

  * `data/intermediate/lep_analysis_panel.csv`

* **ERA5 features**

  * `data/processed/srilanka_district_daily_era5_areawt.csv`
  * `data/processed/srilanka_district_weekly_era5_areawt.csv`

* **Figures** (examples in `reports/figures/`):

  * National trend (`P1_national_trend.png`)
  * Top 15 districts (`P2_top15_latest_week.png`)
  * Seasonality (`P4_seasonality_boxplots.png`)

---

## Repository structure

```
sri-lanka-disease-surveillance/
├─ code/
│  ├─ preprocess/
│  │  ├─ initial_processing.R
│  │  ├─ post_processing.R
│  │  └─ era5_to_weekly_features.R
│  ├─ analysis/
│  │  ├─ gams_and_gbms.R
│  │  ├─ gams_bulk.R
│  │  ├─ initial_viz.R
│  │  └─ sri_lanka_modeling.R
│  └─ helpers/
│     ├─ helpers.R
│     └─ plot_helpers.R
├─ data/
│  ├─ raw/ … population, stations, landcover, ERA5
│  ├─ intermediate/ … extracted PDFs + merged panels
│  └─ processed/ … final daily/weekly ERA5 features
├─ reports/
│  ├─ figures/
│  └─ tables/
├─ shared/helpers.R
├─ params.yaml
└─ sri-lanka-disease-surveillance.Rproj
```

---

## Notes

* **PDF parsing is slow** (>20–30 minutes). For exploratory work, skip `initial_processing.R` and start at `post_processing.R`.
* **ERA5 files are large.** Use Apache Arrow/DVC to pull subsets locally when possible.
* **Java is only needed** if you parse population tables directly from the PDF. Otherwise, supply a CSV.
* **Land cover raster** requires 7-Zip / p7zip to extract `.rar`.

---
