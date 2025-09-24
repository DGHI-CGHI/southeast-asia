---

# Sri Lanka Disease Surveillance ▸ Climate–Health Pipeline

Link Sri Lanka’s **Weekly Epidemiological Reports (WER)** to **population**, **land cover**, and **ERA5** climate features to support modeling (e.g., leptospirosis, dengue). Outputs are analysis-ready for GAMs and other frameworks.

---

## TL;DR — Get started fast

> If you just want to explore **final disease counts linked to ERA5** and try models, you only need the *modeling script* below.

### 1) Clone the repo

```bash
# pick any folder you like
git clone https://github.com/DGHI-CHI/southeast-asia/
cd southeast-asia/sri-lanka-disease-surveillance
```

### 2) Open the RStudio project

Double-click **`sri-lanka-disease-surveillance.Rproj`** (in the repo root).
(Opening the project sets the working directory and plays nicely with `.Rprofile`.)

### 3) Minimal config (paths)

Add this to your **`~/.Rprofile`** *or* to a **project-local** `.Rprofile` in the repo root:

```r
# Example folder layout — edit to your machine
cfg <- list(
  paths = list(
    raw          = "data/raw",
    intermediate = "data/intermediate",
    processed    = "data/processed",
    reports      = "reports"
  )
)
```

> The modeling script reads the analysis panel at
> `file.path(cfg$paths$intermediate, "sri_lanka-disease-landcover-climate-2015_2024.csv")`.

### 4) Install R packages (first run)

```r
install.packages(c(
  "data.table","lubridate","stringr","mgcv",
  "PRROC","pROC"
))
```

> (If you plan to **parse PDFs** or **rebuild** upstream data, you’ll also need:
> `tabulapdf`, `tabulizerjars`, `rJava`, `arrow`, `fixest`, `glmmTMB`, `ranger`, `glmnet`.)

### 5) Open the modeling script

Open **`code/analysis/sri_lanka_modeling.R`** and run it.
This script:

* Loads the **analysis-ready weekly panel** (WER + pop + land cover + **ERA5 weekly** features),
* Auto-discovers ERA5 columns already in the panel (no recomputation here),
* Builds a compact **GAM/BAM** with cyclic seasonality and optional district RE,
* Prints **hold-out metrics** and the exact ERA5 variables used.

> Tip: edit the `smooth_vars` list near the bottom to try different ERA5 drivers.

---

## Full pipeline (if you want to rebuild inputs)

The project has three stages. You only need these if you’re **regenerating** the panel CSV.

1. **Initial processing** – parse WER PDFs → district×week case counts
   `code/preprocess/initial_processing.R`
   Produces `data/intermediate/sri_lanka_WER_index_of_pdfs.csv` and a weekly disease table.

2. **ERA5 features** – build **weekly district-level ERA5** covariates
   `code/preprocess/era5_to_weekly_features.R`
   Produces `data/processed/srilanka_district_weekly_era5_areawt.csv`

3. **Post-processing** – merge disease + population + land cover + **ERA5 weekly** → panel
   `code/preprocess/post_processing.R`
   Produces **`data/intermediate/sri_lanka-disease-landcover-climate-2015_2024.csv`** (consumed by the modeling script).

> For quick exploration, you can skip 1–2 and just use the **ready panel** if it already exists in `data/intermediate/`.

---

## Expected inputs & outputs

### Inputs (examples)

* `data/raw/Mid-year_population_by_district_and_sex_2024.pdf` (or CSV equivalent)
* `data/raw/landcover/…` (2018 land-cover raster packaged as RAR/GeoTIFF)
* ERA5 (hourly/daily) → aggregated to weekly in step 2

### Core output for modeling

* **`data/intermediate/sri_lanka-disease-landcover-climate-2015_2024.csv`**
  Columns: `district`, `date`/`date_mid`, `week`, `year`, `poptot`, **ERA5 weekly features** (and any precomputed lags/rolls/anoms/z), optional land-cover proportions, and disease counts/rates.

---

## Repo structure (short)

```
sri-lanka-disease-surveillance/
├─ code/
│  ├─ preprocess/
│  │  ├─ initial_processing.R          # WER → weekly district counts
│  │  ├─ era5_to_weekly_features.R     # ERA5 → weekly district features
│  │  └─ post_processing.R             # merge → final analysis panel
│  ├─ analysis/
│  │  └─ sri_lanka_modeling.R          # ERA5-only modeling harness (GAM/BAM)
│  └─ helpers/
├─ data/
│  ├─ raw/           # population, land cover, ERA5, etc.
│  ├─ intermediate/  # parsed WER + merged panel (modeling reads from here)
│  └─ processed/     # ERA5 daily/weekly features
├─ reports/
└─ sri-lanka-disease-surveillance.Rproj
```

---

## Troubleshooting

* **Modeling script can’t find the panel CSV**
  Make sure your `.Rprofile` defines `cfg$paths$intermediate` **and** the file
  `sri_lanka-disease-landcover-climate-2015_2024.csv` exists there.

* **SLF4J / Java warnings** while parsing PDFs
  Benign. If `tabulapdf/tabulizerjars` need Java on Windows, set (replace path):

  ```r
  Sys.setenv(JAVA_HOME = "C:/Program Files/Eclipse Adoptium/jdk-17.x.x-hotspot")
  ```

* **Week 53**
  We coerce ISO week 53 → 52 so cyclic splines close the year.

* **Large files (ERA5)**
  Prefer running the modeling script with the prebuilt panel. Rebuilding ERA5 may require Arrow/DVC and disk space.

---

## License / attribution

Add your preferred license here and any data source acknowledgements (Epidemiology Unit, Dept. of Meteorology, ECMWF ERA5, etc.).

---
