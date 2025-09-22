################################################################################
# CHI ??? Sri Lanka ERA5 ??? District Daily + Weekly Features
# File: analysis/sri_lanka/era5_to_weekly_features.R
#
# Purpose
#   Convert ERA5 hourly gridded fields to district-level daily aggregates
#   using area-weighting, then derive epidemiology-ready weekly features
#   aligned to WER week windows (date_start ??? date_end).
#
# What this script produces
#   1) District × day (area-weighted) climate table
#      ??? srilanka_district_daily_era5_areawt.csv
#   2) District × epi-week feature table (lags, rolls, anomalies, etc.)
#      ??? srilanka_district_weekly_era5_areawt.csv
#
# Key steps
#   . Ingest hourly ERA5 parquet (per year), normalize units, localize timestamps
#   . Summarize to cell × day (means/mins/maxes; sums for precipitation)
#   . Overlay grid cells with ADM2 districts (equal-area CRS) to compute
#     intersection areas and area-weighted district × day stats
#   . Aggregate district daily ??? district weekly aligned to WER week windows,
#     and compute features (wet days, 3-day maxima, EWAP, lags, rolls, anomalies)
#
# Inputs (expected layout)
#   era5_root/???year???/*.parquet  with columns (or equivalents):
#     - validTime (seconds since epoch, UTC), ssrd (J/m² per hour),
#     - ta_scaled10_degC, td_scaled10_degC, wbgt_scaled10_degC,
#     - wind2m_scaled10_ms1 (or wind2m_scaled10_ms),
#     - tp (m of total precipitation), mtpr (mean total precip rate; units vary),
#     - lat_idx_x4, lon_idx_x4 (0.25° grid centers = index/4)
#   . Sri Lanka ADM2 boundaries (downloaded automatically from GADM 4.1)
#   . Weeks table from WER data: district, date_start, date_end
#
# Assumptions & conventions
#   . Timezone: Asia/Colombo (+5:30); hourly timestamps are localized before
#     daily aggregation so "days" reflect local health-relevant timing.
#   . Area weighting uses an equal-area projection for intersection areas.
#   . Negative precipitation artifacts are clamped to zero.
#   . Weekly coverage rule is configurable (MIN_DAYS_PER_WEEK).
#
# Author: Jordan Clark (DGHI CHI)
# Last updated: 2025-09-21
################################################################################

suppressPackageStartupMessages({
  library(arrow)        # fast parquet + lazy scanning
  library(dplyr)        # used only inside Arrow pipelines
  library(data.table)   # primary table manipulation
  library(lubridate)    # date-time helpers
  library(here)         # project-relative paths
  
  library(sf)           # vector GIS (districts, intersections)
  library(terra)        # raster, CRSs
  library(exactextractr)# not used in this file, but often paired with terra
  library(stringr)      # string cleanup for names
  library(zoo)          # rolling helpers (weekly features)
})

# ------------------------------------------------------------------------------
# 0) CONFIG
# ------------------------------------------------------------------------------

source(here("helpers", "helpers.R"))  # expected to provide norm_dist(), etc.

#################
#################

CHI_LOCAL_WORK  <- Sys.getenv("CHI_LOCAL_WORK",  unset = "C:/Users/jordan/Desktop/srilanka")
CHI_GITHUB_ROOT <- Sys.getenv("CHI_GITHUB_ROOT", unset = "C:/Users/jordan/R_Projects/CHI-Data")

setwd(CHI_GITHUB_ROOT)

paths <- list(
  temp_dir      = file.path(CHI_LOCAL_WORK, "temp"),
  work_out      = file.path(CHI_LOCAL_WORK, "outputs"),
  # Source inputs (checked into repo)
  midyear_pop   = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs/Mid-year_population_by_district_and_sex_2024.pdf"),
  wx_stations   = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs/station_data/SriLanka_Weather_Dataset.csv"),
  era5_daily    = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs/srilanka_district_daily_era5_areawt.csv"),
  fig_dir       = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/outputs/figures")
)


paths = c(paths, 
          
          outputs = list(
            
            pdf_index_csv = file.path(CHI_GITHUB_ROOT, 'analysis/sri_lanka/outputs/sri_lanka_WER_index_of_pdfs.csv'),
            case_counts_txt = file.path(CHI_GITHUB_ROOT, 'analysis/sri_lanka/outputs/disease_counts_v4.txt')
            
          )
)


# Path holding ERA5 parquet folders, one folder per year (e.g., 1980, 1981, .)
era5_root <- "E:/data/gridded/era5-srilanka/processed"



# Analysis years (adjust as needed)
years_all   <- 2006:2024

# Timezone handling:
#   ERA5 validTime is UTC; we localize to Asia/Colombo (+5:30) BEFORE daily agg
tz_local      <- "Asia/Colombo"
tz_offset_sec <- as.integer(5.5 * 3600)

# ------------------------------------------------------------------------------
# 1) COLUMN MAP (edit if parquet schema differs)
#   Many ERA5 variables are stored as "scaled by 10" integers; we unscale below.
#   Precipitation:
#     - tp   : total precipitation (m) over the accumulation window
#     - mtpr : mean total precipitation rate (unit depends on chain; treated here
#              as rate aggregated to daily sum; negative artifacts are clamped)
# ------------------------------------------------------------------------------

VAR_VALIDTIME <- "validTime"             # epoch seconds (UTC) or timestamp
VAR_SSRD      <- "ssrd"                  # surface solar radiation downwards (J/m² per hour)
VAR_TA10      <- "ta_scaled10_degC"     # air temp (°C * 10)
VAR_TD10      <- "td_scaled10_degC"     # dewpoint (°C * 10)
VAR_WBGT10    <- "wbgt_scaled10_degC"   # WBGT (°C * 10); may be absent
VAR_WIND10    <- "wind2m_scaled10_ms1"  # wind (m/s * 10); fallback handled below

VAR_MTPR      <- "mtpr"                  # precip rate (assumed mm/day after daily sum here)
VAR_TP        <- "tp"                    # total precipitation (m); sum to daily

VAR_LATIDX    <- "lat_idx_x4"            # 0.25° grid center index ??? lat = idx/4
VAR_LONIDX    <- "lon_idx_x4"            # 0.25° grid center index ??? lon = idx/4

# ------------------------------------------------------------------------------
# 2) SMALL HELPERS
# ------------------------------------------------------------------------------

scale10 <- function(x) as.numeric(x) / 10

# Simple "workability" curve from WBGT (kept for compatibility; not used here)
workability <- function(wbgt, steepness = 0.6, midpoint = 30) {
  100 * (1 - 1 / (1 + exp(-steepness * (wbgt - midpoint))))
}

# Great-circle small-distance approximation (not used in final pipeline)
dist2_m2 <- function(lon1, lat1, lon2, lat2) {
  r <- 6371000
  x <- (lon2 - lon1) * cos((lat1 + lat2) * pi / 360)
  y <- (lat2 - lat1)
  (r * pi / 180)^2 * (x^2 + y^2)
}

# ------------------------------------------------------------------------------
# 3) ADMIN BOUNDARIES (ADM2): read + normalize + dissolve to district
#    We'll use an equal-area CRS for intersections to get area in m².
# ------------------------------------------------------------------------------

get_lka_districts <- function() {
  gadm_zip <- "https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_LKA_shp.zip"
  tdir <- tempfile("lka_gadm41_"); dir.create(tdir)
  zipfile <- file.path(tdir, "gadm41_LKA_shp.zip")
  download.file(gadm_zip, destfile = zipfile, mode = "wb", quiet = TRUE)
  unzip(zipfile, exdir = tdir)
  shp <- list.files(tdir, pattern = "^gadm41_LKA_2\\.shp$", full.names = TRUE, recursive = TRUE)
  
  sf_use_s2(TRUE)
  adm2 <- st_read(shp, quiet = TRUE)[, c("GID_2","NAME_1","NAME_2","geometry")]
  names(adm2) <- c("gid2","province","district","geometry")
  
  # Normalize names to project convention
  norm_dist_local <- function(x) {
    x1 <- str_squish(str_to_title(x))
    str_replace_all(x1, "Nuwara Eliya", "Nuwara-Eliya")
  }
  adm2$district <- norm_dist_local(adm2$district)
  
  # Dissolve to unique district polygons
  adm2_diss <- adm2 |>
    dplyr::select(province, geometry) |>
    dplyr::group_by(province) |>
    dplyr::summarise(geometry = st_union(geometry), .groups = "drop")
  
  names(adm2_diss) <- c("district","geometry")
  
  sf_use_s2(FALSE)
  adm2_diss
}

districts_sf_wgs84 <- get_lka_districts()

# Use a global equal-area for weighting; UTM works too, but this is safer nationwide
EA_CRS <- "EPSG:6933"  # World Equal Area
districts_sf_m <- st_transform(districts_sf_wgs84, EA_CRS)

# ------------------------------------------------------------------------------
# 4) BUILD GRID CELL POLYGONS FROM PARQUET INDICES
#    We derive cell edges from center coords (0.25° grid) and build polygons.
# ------------------------------------------------------------------------------

# Read any year to get unique (lat_idx_x4, lon_idx_x4) pairs ??? centers
build_latlon_lookup <- function(year_for_coords = NULL) {
  yr <- if (is.null(year_for_coords)) {
    yrs <- list.dirs(era5_root, recursive = FALSE, full.names = FALSE)
    yrs <- yrs[grepl("^\\d{4}$", yrs)]
    if (!length(yrs)) stop("No year directories found under era5_root.")
    yrs[1]
  } else as.character(year_for_coords)
  
  ds <- open_dataset(file.path(era5_root, yr), format = "parquet",
                     factory_options = list(exclude_invalid_files = TRUE))
  
  # Expect index columns; we compute lat/lon = idx/4 (for 0.25° grid)
  needed <- c("lat_idx_x4","lon_idx_x4")
  if (!all(needed %in% names(ds))) {
    stop("Parquet missing lat_idx_x4 or lon_idx_x4; include them or add a lookup.")
  }
  
  coords <- ds |>
    select(dplyr::all_of(needed)) |>
    distinct() |>
    collect()
  
  setDT(coords)
  coords[, `:=`(
    lat = as.integer(lat_idx_x4)/4,
    lon = as.integer(lon_idx_x4)/4,
    lat_idx_x4 = as.integer(lat_idx_x4)/4,
    lon_idx_x4 = as.integer(lon_idx_x4)/4
  )]
  coords[]
}

coords_lookup <- build_latlon_lookup()

# Build polygons by computing edges between adjacent centers
make_cell_polygons <- function(coords_lookup) {
  ux <- sort(unique(coords_lookup$lon))
  uy <- sort(unique(coords_lookup$lat))
  
  to_edges <- function(v) {
    dv <- diff(v)
    edges <- numeric(length(v) + 1)
    edges[1] <- v[1] - dv[1]/2
    edges[2:length(v)] <- (v[-length(v)] + v[-1]) / 2
    edges[length(v)+1] <- v[length(v)] + dv[length(dv)]/2
    edges
  }
  
  ex <- to_edges(ux); ey <- to_edges(uy)
  lon_map <- data.table(lon = ux, col = seq_along(ux))
  lat_map <- data.table(lat = uy, row = seq_along(uy))
  
  coords <- merge(coords_lookup, lon_map, by = "lon")
  coords <- merge(coords,       lat_map, by = "lat")
  coords[, `:=`(
    x_min = ex[col], x_max = ex[col + 1L],
    y_min = ey[row], y_max = ey[row + 1L]
  )]
  
  polys <- lapply(seq_len(nrow(coords)), function(i) {
    with(coords[i], {
      st_polygon(list(matrix(
        c(x_min, y_min,  x_max, y_min,  x_max, y_max,  x_min, y_max,  x_min, y_min),
        ncol = 2, byrow = TRUE
      )))
    })
  })
  
  cell_sf <- st_sf(
    lat = coords$lat, lon = coords$lon,
    lat_idx_x4 = coords$lat_idx_x4, lon_idx_x4 = coords$lon_idx_x4,
    geometry = st_sfc(polys, crs = 4326)
  )
  
  cell_sf <- st_make_valid(cell_sf)
  cell_sf_m <- st_transform(cell_sf, EA_CRS)
  cell_sf_m$cell_area_m2 <- as.numeric(st_area(cell_sf_m))
  cell_sf_m
}

cell_sf_m <- make_cell_polygons(coords_lookup)

# ------------------------------------------------------------------------------
# 5) PRECOMPUTE CELL???DISTRICT AREA WEIGHTS (m² intersections)
#    This is static across years and variables ??? compute once and reuse.
# ------------------------------------------------------------------------------

compute_cell_district_weights <- function(cell_sf_m, districts_sf_m) {
  idx <- st_intersects(cell_sf_m, districts_sf_m, sparse = TRUE)
  pairs <- data.table(
    cell_row = rep(seq_along(idx), lengths(idx)),
    dist_row = unlist(idx)
  )
  if (nrow(pairs) == 0) stop("No cell-district intersections; check CRS/extent.")
  
  out <- vector("list", length = nrow(districts_sf_m))
  for (j in seq_len(nrow(districts_sf_m))) {
    pj <- pairs[dist_row == j]
    if (nrow(pj) == 0) next
    sub_cells <- cell_sf_m[pj$cell_row, c("lat_idx_x4","lon_idx_x4","cell_area_m2")]
    inter     <- suppressWarnings(st_intersection(sub_cells, districts_sf_m[j, ]))
    if (nrow(inter) == 0) next
    inter$area_m2 <- as.numeric(st_area(inter))
    inter <- inter[inter$area_m2 > 0, ]
    if (nrow(inter) == 0) next
    
    out[[j]] <- data.table(
      lat_idx_x4 = inter$lat_idx_x4,
      lon_idx_x4 = inter$lon_idx_x4,
      district   = districts_sf_wgs84$district[j],
      area_m2    = inter$area_m2
    )
  }
  
  weights <- rbindlist(out, use.names = TRUE, fill = TRUE)
  setkey(weights, lat_idx_x4, lon_idx_x4, district)
  weights
}

weights <- compute_cell_district_weights(cell_sf_m, districts_sf_m)

# ------------------------------------------------------------------------------
# 6) READ ONE YEAR OF HOURLY ERA5 AND NORMALIZE (cell × hour)
#    - Localize timestamps (Asia/Colombo)
#    - Derive daily variables and convert units where needed
# ------------------------------------------------------------------------------

read_one_year_raw <- function(yr) {
  ydir <- file.path(era5_root, yr)
  if (!dir.exists(ydir)) return(data.table())
  
  ds <- open_dataset(ydir, format = "parquet",
                     factory_options = list(exclude_invalid_files = TRUE))
  
  # Some datasets store wind as wind2m_scaled10_ms
  VAR_WIND10_LOCAL <- if ("wind2m_scaled10_ms1" %in% names(ds)) "wind2m_scaled10_ms1" else "wind2m_scaled10_ms"
  
  dt <- ds |>
    select(dplyr::all_of(c(
      VAR_LATIDX, VAR_LONIDX, VAR_VALIDTIME,
      VAR_SSRD, VAR_TA10, VAR_TD10, VAR_WBGT10, VAR_WIND10_LOCAL,
      VAR_MTPR, VAR_TP
    ))) |>
    # Localize hour to Asia/Colombo *before* forming daily aggregates
    mutate(
      validTime = cast(cast(!!sym(VAR_VALIDTIME), int64()), timestamp("s", timezone = "UTC")),
      local_ts  = cast(cast(validTime, int64()) + tz_offset_sec, timestamp("s")),
      date      = cast(local_ts, date32()),
      hour      = hour(local_ts),
      
      # Unscale °C*10 ??? °C; round to 0.1 to remove integer artifacts
      ta   = round(cast(!!sym(VAR_TA10),   float64()) / 10, 1),
      td   = round(cast(!!sym(VAR_TD10),   float64()) / 10, 1),
      wbgt = round(cast(!!sym(VAR_WBGT10), float64()) / 10, 1),
      wind = round(cast(!!sym(VAR_WIND10_LOCAL), float64()) / 10, 1),
      
      mtpr = cast(!!sym(VAR_MTPR), float64()),    # daily sum later
      tp   = cast(!!sym(VAR_TP),   float64()),    # in meters; daily sum later
      
      # Relative humidity (%) and VPD (kPa) derived from ta/td
      rh   = exp((17.62 * td) / (243.12 + td) - (17.62 * ta) / (243.12 + ta)) * 100,
      es   = 0.6108 * exp(17.27 * ta / (ta + 237.3)),
      ea   = 0.6108 * exp(17.27 * td / (td + 237.3)),
      vpd  = if_else(es - ea > 0, es - ea, 0.0),
      
      # ssrd in J/m² per hour ??? MJ/m² per hour
      ssrd_MJ = cast(!!sym(VAR_SSRD), float64()) / 1e6
    ) |>
    select(dplyr::all_of(c(VAR_LATIDX, VAR_LONIDX)),
           date, hour, ta, td, wbgt, wind, rh, vpd, ssrd_MJ, mtpr, tp) |>
    collect()
  
  setDT(dt)
  
  # Clamp negative precip artifacts
  dt[tp   < 0, tp   := 0]
  dt[mtpr < 0, mtpr := 0]
  
  # Convert index columns to true 0.25° coordinates for merging with weights
  setnames(dt, c(VAR_LATIDX, VAR_LONIDX), c("lat_idx_x4","lon_idx_x4"))
  dt[, `:=`(lat_idx_x4 = as.integer(lat_idx_x4)/4,
            lon_idx_x4 = as.integer(lon_idx_x4)/4)]
  dt[]
}

# ------------------------------------------------------------------------------
# 7) HOURLY ??? DAILY (cell × day): means/mins/maxes; precip sums
# ------------------------------------------------------------------------------

cell_day_stats <- function(dt_hourly) {
  if (nrow(dt_hourly) == 0) return(dt_hourly)
  
  vars_mmm <- c("ta","td","wbgt","wind","rh","vpd","ssrd_MJ")  # mean/min/max per day
  vars_sum <- c("tp","mtpr")                                  # sum per day
  
  mmm <- dt_hourly[, c(
    setNames(lapply(.SD, \(x) mean(x, na.rm = TRUE)), paste0(vars_mmm, "_mean")),
    setNames(lapply(.SD, \(x) min(x,  na.rm = TRUE)), paste0(vars_mmm, "_min")),
    setNames(lapply(.SD, \(x) max(x,  na.rm = TRUE)), paste0(vars_mmm, "_max"))
  ), by = .(lat_idx_x4, lon_idx_x4, date), .SDcols = vars_mmm]
  
  sums <- dt_hourly[, c(
    setNames(lapply(.SD, \(x) sum(x, na.rm = TRUE)), paste0(vars_sum, "_sum"))
  ), by = .(lat_idx_x4, lon_idx_x4, date), .SDcols = vars_sum]
  
  out <- merge(mmm, sums, by = c("lat_idx_x4","lon_idx_x4","date"), all = TRUE)
  setorder(out, lat_idx_x4, lon_idx_x4, date)
  out[]
}

# Area-weighted mean helper
aw_mean <- function(val, w) {
  i <- is.finite(val)
  if (!any(i)) return(NA_real_)
  sum(val[i] * w[i]) / sum(w[i])
}

# ------------------------------------------------------------------------------
# 8) CELL-DAILY ??? DISTRICT-DAILY using area weights (m²)
# ------------------------------------------------------------------------------

district_day_from_cell_day <- function(cell_daily, weights) {
  if (nrow(cell_daily) == 0) return(cell_daily)
  x <- merge(cell_daily, weights,
             by = c("lat_idx_x4","lon_idx_x4"),
             allow.cartesian = TRUE)
  
  vars_mmm <- c("ta","td","wbgt","wind","rh","vpd","ssrd_MJ")
  
  x[, {
    means <- setNames(
      lapply(vars_mmm, \(v) aw_mean(get(paste0(v, "_mean")), area_m2)),
      paste0(vars_mmm, "_mean")
    )
    mins  <- setNames(
      lapply(vars_mmm, \(v) suppressWarnings(min(get(paste0(v, "_min")), na.rm = TRUE))),
      paste0(vars_mmm, "_min")
    )
    maxs  <- setNames(
      lapply(vars_mmm, \(v) suppressWarnings(max(get(paste0(v, "_max")), na.rm = TRUE))),
      paste0(vars_mmm, "_max")
    )
    # Area-weighted means of daily precipitation sums ??? district daily totals
    tp_aw   <- aw_mean(tp_sum,   area_m2)
    mtpr_aw <- aw_mean(mtpr_sum, area_m2)
    
    as.list(c(means, mins, maxs,
              tp_sum = tp_aw,
              mtpr_sum = mtpr_aw))
  }, by = .(district, date)][order(district, date)]
}

# ------------------------------------------------------------------------------
# 9) DRIVER: loop over years ??? write district-daily CSV
# ------------------------------------------------------------------------------

summarize_years_area_weighted <- function(years) {
  res_list <- vector("list", length(years))
  for (i in seq_along(years)) {
    yr <- years[i]
    message("Processing ", yr, " .")
    hr <- read_one_year_raw(yr)
    if (!all(c("lat_idx_x4","lon_idx_x4","date") %in% names(hr))) {
      stop("Hourly frame missing lat_idx_x4/lon_idx_x4/date after read for year ", yr)
    }
    cd <- cell_day_stats(hr)
    dd <- district_day_from_cell_day(cd, weights)
    dd[, year := yr]
    res_list[[i]] <- dd
  }
  rbindlist(res_list, use.names = TRUE, fill = TRUE)
}

# Example: run for selected years and write to CSV
# years <- 2006:2024
# out_daily <- summarize_years_area_weighted(years)
# fwrite(out_daily, here("analysis/sri_lanka/srilanka_district_daily_era5_areawt.csv"))

# ------------------------------------------------------------------------------
# 10) WEEKLY WEATHER FEATURE PIPELINE
#     Inputs:
#       - daily_dt: district × day (from step 9)
#       - weeks_dt: district, date_start, date_end (from WER)
#     Output: district × epi week with covariates ready for modeling
# ------------------------------------------------------------------------------

# Safe helpers for aggregation
q_na    <- function(x, p) if (all(is.na(x))) NA_real_ else as.numeric(quantile(x, p, na.rm = TRUE))
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sum_na  <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)

winsor <- function(x, probs = c(0.01, 0.99)) {
  if (all(is.na(x))) return(x)
  qs <- quantile(x, probs = probs, na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

# Longest run of days meeting a precip threshold
longest_wet_spell <- function(prec, thr = 10) {
  if (all(is.na(prec))) return(NA_real_)
  wet <- as.integer(prec >= thr)
  rle_wet <- rle(wet)
  if (!length(rle_wet$lengths)) return(0)
  max(c(0, rle_wet$lengths[rle_wet$values == 1]))
}

iso_week <- function(d) as.integer(strftime(d, format = "%V"))

# Rolling with partial coverage (within a district's weekly time series)
roll_mean_partial <- function(x, k, min_obs = 1L) {
  zoo::rollapplyr(x, k, function(z) if (sum(!is.na(z)) >= min_obs) mean(z, na.rm = TRUE) else NA_real_, fill = NA)
}
roll_sum_partial <- function(x, k, min_obs = 1L) {
  zoo::rollapplyr(x, k, function(z) if (sum(!is.na(z)) >= min_obs) sum(z, na.rm = TRUE) else NA_real_, fill = NA)
}

zscore_robust <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  mu <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv <= 1e-8) {
    madv <- mad(x, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(madv) || madv <= 1e-8) return(rep(0, length(x)))
    return((x - mu) / madv)
  } else (x - mu) / sdv
}

pct_of_normal <- function(x, clim, eps = 1e-6) {
  out <- x / pmax(clim, eps)
  out[x == 0 & clim == 0] <- 1
  out
}

# Configuration for features
CFG <- list(
  WINSORIZE    = TRUE,
  WINSOR_PROBS = c(0.01, 0.99),
  HOT_TMAX     = 32,    # °C threshold for hot-day counts
  WET_MM       = 10,    # "wet day" threshold (mm)
  MAX_LAG_WEEKS= 6,
  ROLL_WINDOWS = c(2, 4),
  EWAP_ALPHA   = 0.8,   # exponentially-weighted antecedent precip
  EWAP_K       = 4,
  ANOMALY_VARS = c("precip_tp_sum_week","tmax_mean","rh_mean_week"),
  ZSCORE_VARS  = c("precip_tp_sum_week","precip_mtpr_sum_week","tmax_mean","vpd_mean_week"),
  MIN_DAYS_PER_WEEK = 5 # require ???5/7 days in week to compute metrics
)

weekly_weather_features <- function(daily_dt, weeks_dt, cfg = CFG) {
  # Required columns in inputs
  stopifnot(all(c("district","date") %in% names(daily_dt)))
  stopifnot(all(c("district","date_start","date_end") %in% names(weeks_dt)))
  
  DTd <- as.data.table(copy(daily_dt))
  DTw <- as.data.table(copy(weeks_dt))
  
  # Ensure types
  if (!inherits(DTd$date, "Date")) DTd[,  date := as.Date(date)]
  if (!inherits(DTw$date_start, "Date")) DTw[, date_start := as.Date(date_start)]
  if (!inherits(DTw$date_end, "Date"))   DTw[, date_end   := as.Date(date_end)]
  
  # Winsorize some daily inputs to reduce outlier impact (optional)
  num_cols <- intersect(
    c("ta_mean","ta_min","ta_max","wbgt_mean","rh_mean","vpd_mean",
      "ssrd_MJ_mean","tp_sum","mtpr_sum"),
    names(DTd)
  )
  if (cfg$WINSORIZE && length(num_cols)) {
    DTd[, (num_cols) := lapply(.SD, winsor, probs = cfg$WINSOR_PROBS), .SDcols = num_cols]
  }
  
  # Build a daily "calendar" for each epi week; join daily climate onto it
  idx <- DTw[, .(district, week_id = .I,
                 date_start, date_end,
                 year = year(date_end),
                 week_of_year = iso_week(date_end))]
  
  wk_days <- idx[
    , .(date = seq(date_start, date_end, by = "1 day")),
    by = .(district, week_id, date_start, date_end, year, week_of_year)
  ]
  
  setkey(DTd, district, date)
  setkey(wk_days, district, date)
  wk <- DTd[wk_days, on = .(district, date)]
  
  # Aggregate within each epi week (district × week_id)
  agg <- wk[, {
    ndays <- sum(!is.na(date))
    tmax <- ta_max; tmin <- ta_min; tmean <- ta_mean
    wbgt <- wbgt_mean; rh <- rh_mean; vpd <- vpd_mean
    ssrd <- ssrd_MJ_mean
    tp   <- tp_sum        # daily precipitation (units consistent with earlier steps)
    mtpr <- mtpr_sum
    
    # Rolling 3-day maxima within the week (requires ???3 non-NA days)
    max3_tp   <- if (sum(!is.na(tp))   >= 3) max(zoo::rollapplyr(tp,   3, sum, na.rm = TRUE), na.rm = TRUE) else NA_real_
    max3_mtpr <- if (sum(!is.na(mtpr)) >= 3) max(zoo::rollapplyr(mtpr, 3, sum, na.rm = TRUE), na.rm = TRUE) else NA_real_
    
    # Weekly metrics
    res <- list(
      n_days_week            = as.integer(ndays),
      
      # Temperature summaries (°C)
      tmax_mean              = mean_na(tmax),
      tmax_p90               = q_na(tmax, 0.90),
      tmax_p95               = q_na(tmax, 0.95),
      tmax_range             = if (all(is.na(tmax))) NA_real_ else max(tmax, na.rm = TRUE) - min(tmax, na.rm = TRUE),
      tmin_mean              = mean_na(tmin),
      tmean_mean             = mean_na(tmean),
      wbgt_mean_week         = mean_na(wbgt),     # °C
      
      # Moisture & evaporative demand
      rh_mean_week           = mean_na(rh),       # %
      vpd_mean_week          = mean_na(vpd),      # kPa
      
      # Radiation (MJ/m² per day ??? weekly mean)
      ssrd_MJ_mean_week      = mean_na(ssrd),
      
      # Weekly precipitation totals (units follow daily sums)
      precip_tp_sum_week     = sum_na(tp),
      precip_mtpr_sum_week   = sum_na(mtpr),
      
      # Wet-day counts (threshold in mm)
      wet_days_ge10_tp       = as.numeric(sum(tp   >= cfg$WET_MM, na.rm = TRUE)),
      wet_days_ge10_mtpr     = as.numeric(sum(mtpr >= cfg$WET_MM, na.rm = TRUE)),
      
      # 3-day maxima and longest wet spells
      max3d_tp               = max3_tp,
      max3d_mtpr             = max3_mtpr,
      wet_spell_maxlen_tp    = as.numeric(longest_wet_spell(tp,   thr = cfg$WET_MM)),
      wet_spell_maxlen_mtpr  = as.numeric(longest_wet_spell(mtpr, thr = cfg$WET_MM)),
      
      # Heat exceedances
      hot_days_ge32          = as.numeric(sum(ta_max >= cfg$HOT_TMAX, na.rm = TRUE))
    )
    
    # Coverage rule: mask if insufficient days (keep n_days_week visible)
    if (ndays < cfg$MIN_DAYS_PER_WEEK) {
      keep <- "n_days_week"
      res[names(res) != keep] <- lapply(res[names(res) != keep], function(.) NA_real_)
    }
    res
  }, by = .(district, week_id, date_start, date_end, year, week_of_year)]
  
  # Ordered by time for rolling ops
  setorder(agg, district, date_end)
  
  # Lags and rolling windows (past-only) for selected variables
  LAG_VARS <- intersect(
    c("precip_tp_sum_week","precip_mtpr_sum_week","rh_mean_week","tmax_mean","vpd_mean_week"),
    names(agg)
  )
  for (v in LAG_VARS) {
    for (L in seq_len(cfg$MAX_LAG_WEEKS)) {
      agg[, paste0(v, "_lag", L) := data.table::shift(get(v), L), by = district]
    }
  }
  
  # Rolling windows with partial coverage (2w mean/sum, 4w mean/sum)
  for (v in LAG_VARS) {
    agg[, paste0(v, "_roll2w_mean") := roll_mean_partial(get(v), k = 2, min_obs = 1L), by = district]
    agg[, paste0(v, "_roll2w_sum")  := roll_sum_partial (get(v), k = 2, min_obs = 1L), by = district]
    agg[, paste0(v, "_roll4w_mean") := roll_mean_partial(get(v), k = 4, min_obs = 2L), by = district]
    agg[, paste0(v, "_roll4w_sum")  := roll_sum_partial (get(v), k = 4, min_obs = 2L), by = district]
  }
  
  # Exponentially-weighted antecedent precipitation (EWAP)
  K <- cfg$EWAP_K; a <- cfg$EWAP_ALPHA
  if ("precip_tp_sum_week" %in% names(agg)) {
    agg[, ewap_tp := {
      x <- precip_tp_sum_week; ew <- x
      for (k in 1:K) ew <- ew + (a^k) * data.table::shift(x, k, fill = 0)
      ew
    }, by = district]
  }
  if ("precip_mtpr_sum_week" %in% names(agg)) {
    agg[, ewap_mtpr := {
      x <- precip_mtpr_sum_week; ew <- x
      for (k in 1:K) ew <- ew + (a^k) * data.table::shift(x, k, fill = 0)
      ew
    }, by = district]
  }
  
  # Climatology (district × ISO week) and anomalies / percent-of-normal
  CLIM_VARS <- intersect(cfg$ANOMALY_VARS, names(agg))
  if (length(CLIM_VARS)) {
    clim <- agg[, lapply(.SD, mean_na), by = .(district, week_of_year), .SDcols = CLIM_VARS]
    setnames(clim, CLIM_VARS, paste0(CLIM_VARS, "_clim"))
    agg <- clim[agg, on = .(district, week_of_year)]
    for (v in CLIM_VARS) {
      vc <- paste0(v, "_clim")
      agg[, paste0(v, "_anom")       := get(v) - get(vc)]
      agg[, paste0(v, "_pct_normal") := pct_of_normal(get(v), get(vc), eps = 1e-6)]
    }
  }
  
  # Robust within-district z-scores for select variables
  ZV <- intersect(cfg$ZSCORE_VARS, names(agg))
  if (length(ZV)) {
    for (v in ZV) {
      agg[, paste0(v, "_z") := zscore_robust(get(v)), by = district]
    }
  }
  
  setorder(agg, district, date_end)
  agg[]
}

# ------------------------------------------------------------------------------
# 11) EXAMPLE I/O (uncomment to run end-to-end)
# ------------------------------------------------------------------------------

# 1) Generate district-daily file (only once per update)
years <- 2006
out_daily <- summarize_years_area_weighted(years)
fwrite(out_daily, here("analysis/sri_lanka/srilanka_district_daily_era5_areawt.csv"))

# 2) Weekly features (requires WER week windows)
names(out_daily) <- gsub("_aw","",names(out_daily))  # if earlier suffixes existed



lepto <- fread(file.path(paths$work_out, "lep_analysis_panel.csv"))
               
               
weeks_dt <- unique(lepto[, .(district, date_start, date_end)])  # from your WER table
features_weekly <- weekly_weather_features(out_daily, weeks_dt)
fwrite(features_weekly, here("analysis/sri_lanka/srilanka_district_weekly_era5_areawt.csv"))

# 3) Join back to health data for modeling
features_weekly <- fread(here("analysis/sri_lanka/srilanka_district_weekly_era5_areawt.csv"))
lepto_feat <- features_weekly[lepto, on = .(district, date_start, date_end)]
# Fit GLM/GLMM with offset(log(poptot)) and chosen features.






