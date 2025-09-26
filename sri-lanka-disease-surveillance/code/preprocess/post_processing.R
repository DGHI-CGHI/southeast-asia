################################################################################
# CHI  Sri Lanka WER - Post-Processing Pipeline
# File: analysis/sri_lanka/post_processing.R
#
# Purpose:
#   Consume outputs from the initial WER scrape/extraction and build an
#   analysis-ready panel by merging:
#     . Weekly district disease counts (from initial_processing.R)
#     . Mid-year population (PDF  tidy, or CSV fallback)  rates per 100k
#     . Station-based weather aggregated to district × day
#     . District land-cover proportions (exactextractr over raster)
#     . ERA5 district-daily aggregates
#   Final artifact:
#     - lep_analysis_panel.csv (district × week with covariates)
#
# Upstream Inputs (must exist):
#   - outputs/disease_counts_v4.txt (from initial_processing.R)
#   - station dataset: analysis/sri_lanka/station_data/SriLanka_Weather_Dataset.csv
#   - land cover raster: analysis/sri_lanka/SriLanka_Landcover_2018.tif
#   - ERA5 district-daily CSV: analysis/sri_lanka/srilanka_district_daily_era5_areawt.csv
#
# Population Source Options:
#   1) Default: Extract tables from PDF (requires Java + tabulizer stack)
#      midyear_pop = analysis/sri_lanka/Mid-year_population_by_district_and_sex_2024.pdf
#   2) Fallback: Provide a ready CSV via env var CHI_POP_CSV with columns:
#      district, year, poptot
#
# Environment variables (set in ~/.Renviron):
#   CHI_LOCAL_WORK   - working directory for temp/outputs
#   CHI_GITHUB_ROOT  - repo root
#   (Optional) CHI_POP_CSV - path to population CSV to bypass Java/tabulizer
#
# Java Note (only needed if using PDF extraction):
#   Configure JAVA_HOME to a valid JDK and ensure rJava works.
#
# Author: Jordan Clark (DGHI CHI)
# Last updated: 2025-09-21
################################################################################

suppressPackageStartupMessages({
  # Core / IO
  library(data.table)
  library(stringr)
  library(lubridate)
  
  # HTML/PDF (PDF *tables* only used if no CSV fallback is provided)
  library(pdftools)
  library(xml2)
  library(rvest)
  
  # Spatial + raster
  library(sf)
  library(terra)
  library(exactextractr)
  
  # Plotting helpers (not used here, but harmless if downstream scripts source this)
  library(ggplot2)
  library(scales)
  library(ragg)
})

setDTthreads(8)

###############################################################################
# 0) CONFIG & PATHS ------------------------------------------------------------
###############################################################################

# -- Configuration via environment variables (with sensible defaults) ----------
# Tip: set these in ~/.Renviron or project .Renviron
#   CHI_LOCAL_WORK="C:/Users/jordan/Desktop/srilanka"
#   CHI_GITHUB_ROOT="C:/Users/jordan/R_Projects/CHI-Data"
#   JAVA_HOME="C:/Program Files/Eclipse Adoptium/jdk-17.0.16.8-hotspot"
#################
#################
# PDF table extraction (tabulizer stack via rJava)
# MUST SET JAVA_HOME FIRST before loading library.
Sys.setenv(JAVA_HOME = Sys.getenv("JAVA_HOME", unset = Sys.getenv("JAVA_HOME")))  # no-op if already set
Sys.setenv("JAVA_HOME" = "C:/Program Files/Eclipse Adoptium/jdk-17.0.16.8-hotspot")

# To install PDF related packages:
# install.packages("rJava")
# remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")

library(tabulapdf)       # ropensci fork; requires rJava

# ------------------------------------------------------------------------------
# 0) CONFIG
# ------------------------------------------------------------------------------
# Summary: Define project-relative paths using cfg (from .Rprofile) + file.path() helper.

# 1) No setwd() and no absolute C:\ pathsa
#    Everything below is relative to the Sri Lanka subproject root.

# Define which ERA5 (ERA5 or ERA5-LAND) to utilize.
# DEFINE WHERE ERA5 (ERA5 OR ERA5-LAND) TO USE.
ERA5_TYPE = 'era5land'
# ERA5_TYPE = 'era5'


# Where to put scratch outputs during processing
paths <- list(
  temp_dir = file.path(cfg$paths$intermediate, "temp"),
  work_out = file.path(cfg$paths$intermediate, "outputs"),
  helpers  = file.path('code/helpers/helpers.R')
)


# 2) Source inputs that should live in the project repo:
#    Recommend moving PDFs/CSVs under data/raw or data/raw/external.
#    For now, keep your current filenames but make them relative.
paths$midyear_pop <- file.path(cfg$paths$raw,
                       "Mid-year_population_by_district_and_sex_2024.pdf")
paths$wx_stations <- file.path(cfg$paths$raw , "station_data/SriLanka_Weather_Dataset.csv")
# 5) Define ERA5 data root location (S3 bucket - publically accessible).
if (ERA5_TYPE == 'era5land'){
  # FOR ERA5-LAND
  paths$era5_daily  <- file.path(cfg$paths$raw, "srilanka_district_daily_era5land_areawt.csv")
} else if (ERA5_TYPE == 'era5'){
  # For ERA5.
  paths$era5_daily  <- file.path(cfg$paths$raw, "srilanka_district_daily_era5_areawt.csv")
}

# 3) Output figures directory inside project
paths$fig_dir <- file.path(cfg$paths$reports, "figures")

# 4) Additional named outputs (project-relative)
paths$outputs <- list(
  era5_weekly_aggregated = file.path(
    cfg$path$processed,
    'srilanka_district_weekly_era5_areawt.csv'
  ),
  pdf_index_csv   = file.path(cfg$path$intermediate, "sri_lanka_WER_index_of_pdfs.csv"),
  case_counts_txt = file.path(cfg$path$intermediate, "disease_counts.txt")
)

# 5) ERA5 root **local path** (hydrated by DVC), not s3://
#    DVC should import/pull into data/raw/era5/
era5_root <- file.path(cfg$paths$raw, "era5")
# era5_root <- "s3://dghi-chi/data/se-asia/sri-lanka-disease-surveillance/era5/"


# 6) Create any needed directories once (safe if they already exist)
ensure_dir <- function(...)
  dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE)
ensure_dir(cfg$paths$intermediate, "temp")
ensure_dir(cfg$paths$intermediate, "outputs")
ensure_dir(cfg$paths$reports, "figures")


# 0.3 Helper function sources ------------------------------------------------
# These must define: DISEASES, POS_MAP, norm_dist(), .norm(), .is_footer(),
# .parse_ints(), .extract_district_from_row(), .pick_n(), DIST_CANON
source(paths$helpers)



################################################################################
# 1) LOAD WEEKLY DISEASE COUNTS -------------------------------------------------
################################################################################
# SECTION GOAL: Read district-week disease counts, normalize, and precompute dates.
lep <- fread(paths$outputs$case_counts_txt)

# Normalize key fields and derive mid-week date (used for daily joins)
lep[, date_mid := as.IDate(date_start + (as.integer(date_end - date_start) / 2))]
lep[, `:=`(
  district = norm_dist(district),
  date_mid = as.IDate(date_mid),
  date_end = as.IDate(date_end)
)]
lep$year = lubridate::year(lep$date_start)

################################################################################
# 2) POPULATION: PDF  TIDY (or CSV fallback) ----------------------------------
################################################################################
# SECTION GOAL: Produce pop_dt with columns [district, year, poptot].
# Preference order:
#   (A) CHI_POP_CSV provided   read directly (no Java)
#   (B) Else parse from PDF via tabulizer (requires Java stack)
# NOTE: Requires working rJava/Tabulizer on your machine; handled with tryCatch.
pop_dt <- NULL
if (file.exists(paths$midyear_pop)) {
  tabs <- tryCatch(
    extract_tables(
      paths$midyear_pop,
      method = "lattice",
      guess = TRUE,
      output = "tibble"
    ),
    error = function(e)
      NULL
  )
  if (!is.null(tabs) && length(tabs) >= 1) {
    yrs <- 2014:2023
    pieces <- list()
    for (pg in seq_along(tabs)) {
      tb <- as.data.table(tabs[[pg]])
      yrheads <- names(tb)[names(tb) %like% paste(yrs, collapse = "|")]
      if (!length(yrheads))
        next
      # For each detected year column, grab District names + "Total" col under that year
      for (yh in yrheads) {
        district_names <- tb$District[-1]
        col_idx <- which(tb[1] == "Total")[which(names(tb) == yh)]
        if (!length(col_idx))
          next
        vals <- tb[[col_idx]][-1]
        pieces[[length(pieces) + 1L]] <- data.table(
          year    = as.integer(gsub("\\*", "", yh)),
          district = district_names,
          poptot  = as.numeric(gsub(",", "", vals)) * 1000
        )
      }
    }
    pop_dt <- rbindlist(pieces, fill = TRUE)
    pop_dt[, district := norm_dist(district)]
  }
}
stopifnot(!is.null(pop_dt) && nrow(pop_dt) > 0)
pop_dt
pop_dt = pop_dt[!is.na(district)] # remove national level totals.

# lep[!district %in% pop_dt$district][,.N,by=district]
# lep[district %in% pop_dt$district][,.N,by=district][order(district)]

# current mid year pop data starts in 2014, so for years in health data before that, simply using 2014 pop for now.
lep[, year2merge := fifelse(year >= 2014, year, 2014L)]
lep <- merge(
  lep,
  pop_dt,
  by.x = c("district", "year2merge"),
  by.y = c("district", "year"),
  all.x = TRUE
)

# -- 2.4 Compute rates per 100k -------------------------------------------------
lep[, `:=`(
  lepto_100k  = (lepto  / poptot) * 1e5,
  dengue_100k = (dengue / poptot) * 1e5
)]

lep[, year := NULL] # year no longer needed after merge





# ################################################################################
# # 3) WEATHER: MAP STATIONS  DISTRICTS & AGG DAILY -----------------------------

# Commented out given uncertainty in how these weather station data were fetched
#   how, from where, etc. ERA5 should be more robust-but keeping here for reference.

# ################################################################################
# # SECTION GOAL: Assign station observations to districts, then compute
# # areal averages per district × day for selected variables.
# 
# wx <- fread(paths$wx_stations)
# stopifnot(all(c("time", "latitude", "longitude", "city") %in% names(wx)))
# if (!inherits(wx$time, "Date"))
#   wx[, time := as.IDate(time)]
# 
# # -- Assign each unique (city, lat, lon) to a district via spatial join --------
# gadm_zip <- "https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_LKA_shp.zip"
# tdir <- tempfile("lka_gadm41_")
# dir.create(tdir)
# zipfile <- file.path(tdir, "gadm41_LKA_shp.zip")
# download.file(gadm_zip, destfile = zipfile, mode = "wb", quiet = TRUE)
# unzip(zipfile, exdir = tdir)
# adm2_shp <- list.files(tdir,
#                        pattern = "^gadm41_LKA_2\\.shp$",
#                        full.names = TRUE,
#                        recursive = TRUE)
# stopifnot(length(adm2_shp) == 1)
# 
# adm2 <- st_read(adm2_shp, quiet = TRUE)[, c("GID_2", "NAME_1", "NAME_2", "geometry")]
# names(adm2) <- c("gid2", "province", "district", "geometry")
# adm2$district <- norm_dist(adm2$district)
# adm2 <- st_make_valid(adm2)
# 
# stations_lu <- unique(wx[, .(city, latitude, longitude)])
# stations_sf <- st_as_sf(
#   stations_lu,
#   coords = c("longitude", "latitude"),
#   crs = 4326,
#   remove = FALSE
# )
# stations_sf <- st_transform(stations_sf, st_crs(adm2))
# stations_join <- st_join(stations_sf, adm2["district"], left = TRUE)
# city2dist <- as.data.table(stations_join)[, .(city, latitude, longitude, district)]
# 
# wx_d <- merge(wx,
#               city2dist,
#               by = c("city", "latitude", "longitude"),
#               all.x = TRUE)
# wx_d[, district := norm_dist(district)]
# 
# # -- District × day means (areal average across stations) ----------------------
# weather_means <- c(
#   "temperature_2m_max",
#   "temperature_2m_min",
#   "temperature_2m_mean",
#   "apparent_temperature_max",
#   "apparent_temperature_min",
#   "apparent_temperature_mean",
#   "shortwave_radiation_sum",
#   "precipitation_sum",
#   "rain_sum",
#   "precipitation_hours",
#   "windspeed_10m_max",
#   "windgusts_10m_max",
#   "et0_fao_evapotranspiration"
# )
# 
# keep_cols <- c("district", "time", weather_means)
# wx_keep  <- wx_d[, intersect(names(wx_d), keep_cols), with = FALSE]
# wx_daily <- wx_keep[, lapply(.SD, mean, na.rm = TRUE), by = .(district, date = time), .SDcols = setdiff(names(wx_keep), c("district", "time"))]




# ###############################################################################
# 4) LAND COVER  Download .RAR  Extract  Raster ------------------------------
# GOAL:
#   . Download Sri Lanka land-cover archive (.rar) from NODA
#   . Extract into analysis/sri_lanka/SriLanka_Landcover_2018/
#   . Load the extracted .tif with terra::rast()
#
# REQUIREMENTS (one-time per machine):
#   . 7-Zip or unrar installed and accessible (script auto-detects common paths)
#       - Windows: install 7-Zip (https://www.7-zip.org/)
#       - macOS:   brew install p7zip
#       - Ubuntu:  sudo apt-get install p7zip-full
#
# # Landcover citation.
# ZHONG Bo, HU Longfei, WU Junjun, YANG Aixia. 30 m Land Cover Dataset of Sri Lanka (2018)[J/DB/OL].
# Digital Journal of Global Change Data Repository, 2020. https://doi.org/10.3974/geodb.2020.09.06.V1.
###############################################################################
lc_dir = file.path(cfg$paths$raw, "landcover")
landcover_file = file.path(lc_dir, "SriLanka_Landcover_2018")
# lc_dir = gsub("\\.rar", "", landcover_file)
if (!file.exists(landcover_file)) {
  # download from first URL below is not always succesful, but can try if second URL doesn't work.
  # landcover_data_url = 'https://www.noda.ac.cn/en/knowledgehub/downloadAttachmentid?=64257e8a29d1210627d613d3&type=DATASET'
  landcover_data_url = 'https://www.geodoi.ac.cn/weben/down.aspx?fileID=5329'
  options(timeout = 180) # increase timeout for download since connection to website is slow.
  download.file(
    landcover_data_url, landcover_file,
    method = "curl",
    extra = "-L -A 'Mozilla/5.0'"
  )
}
# 
extract_rar(landcover_file, lc_dir)

# Point terra to the extracted .tif (adjust name if different inside the RAR)
tif_file <- paste0(landcover_file, ".tif")
# file.exists(tif_file)
if (is.na(tif_file))
  stofile.path("No .tif found after extraction in: ", lc_dir)

r <- rast(tif_file)

# Dissolve to district and match raster CRS
adm2_diss <- adm2 |>
  dplyr::select(province, geometry) |>
  dplyr::group_by(province) |>
  dplyr::summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_transform(crs(r))

names(adm2_diss) <- c("district", "geometry")

# Legend mapping
class_map <- data.table(
  code  = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 255),
  label = c(
    "Water",
    "BuiltUp",
    "Cropland",
    "Forest",
    "Shrub",
    "Grass",
    "Bare",
    "Wetland",
    "Paddy",
    "NoData"
  )
)

# Weighted frequency  proportions within polygon
freq_fun <- function(values, coverage_fraction) {
  dt <- data.table(code = values, w = coverage_fraction)[!is.na(code)]
  if (!nrow(dt))
    return(data.frame(code = integer(), prop = numeric()))
  dt <- dt[, .(w = sum(w, na.rm = TRUE)), by = code][, prop := w / sum(w)]
  as.data.frame(dt[, .(code, prop)])
}

lc_list <- exact_extract(
  r,
  adm2_diss,
  fun = function(values, coverage_fraction)
    list(freq_fun(values, coverage_fraction)),
  progress = TRUE
)

lc_dt <- rbindlist(lc_list, idcol = "idx")
lc_dt[, district := adm2_diss$district[idx]][, idx := NULL]
lc_dt <- lc_dt[code %in% class_map$code]
lc_dt <- merge(lc_dt, class_map, by = "code", all.x = TRUE)
lc_wide <- dcast(lc_dt, district ~ label, value.var = "prop", fill = 0)

# Re-normalize after dropping NoData to ensure rows sum ~1
prop_cols <- setdiff(names(lc_wide), "district")
row_sums <- lc_wide[, rowSums(.SD), .SDcols = prop_cols]
lc_wide[, (prop_cols) := lapply(.SD, function(z)
  ifelse(row_sums > 0, z / row_sums, 0)), .SDcols = prop_cols]

file.remove(tif_file) # remove raster, keeping .tar file.

################################################################################
# 5) ERA5 + MERGES + OUTPUT -----------------------------------------------------
################################################################################
# SECTION GOAL: Join district-daily ERA5 to weekly disease by aligning on
# date_mid; add weather & land cover; write final analysis panel.

# -- 5.1 ERA5 input -------------------------------------------------------------
era5_weekly_features <- fread(paths$outputs$era5_weekly_aggregated)

era5_weekly_features[, date := as.IDate(date_start + (as.integer(date_end - date_start) / 2))]
era5_weekly_features$date_start=NULL
era5_weekly_features$date_end=NULL


stopifnot(all(c("date", "district") %in% names(era5_weekly_features)))
era5_weekly_features[, district := norm_dist(district)]
era5_weekly_features[, date := as.IDate(date)]

# # -- 5.2 Join weather STATION (district × day) to weekly (by date_mid) -----------------
# setnames(wx_daily, "date", "wx_date")
# setkey(lep, district, date_mid)
# setkey(wx_daily, district, wx_date)
# lep <- wx_daily[lep, on = .(district, wx_date = date_mid)]  # left join onto lep

# names(lep)

# -- 5.3 Join land cover + era5_weekly_features (align era5_weekly_features date to date_mid) ------------------
lep <- merge(lep, lc_wide, by = "district", all.x = TRUE)

# Recompute date_mid defensively (ensures correct type after merges)
lep[, date_mid := as.IDate(date_start + (as.integer(date_end - date_start) / 2))]

lep <- merge(
  lep,
  era5_weekly_features,
  by.x = c("district", "date_mid"),
  by.y = c("district", "date"),
  all.x = TRUE
)

lep = lep[year(date_start) <= 2024] # climate data persists through 2024 as of now.
# start in 2014, since that is when mid year population data starts
lep = lep[year(date_start) >= 2014]

# names(lep)
outfilename = sprintf("sri_lanka-disease-landcover-climate-2015_2024-%s.csv", ERA5_TYPE)
# # -- 5.4 Persist final analysis dataset ----------------------------------------
fwrite(lep, file.path(cfg$paths$intermediate, outfilename))

message("Saved analysis panel: ",
        normalizePath(
          file.path(cfg$paths$intermediate, outfilename),
          winslash = "/"
        ))


# # End core workflow.
################################################################################
