################################################################################
# CHI ??? Sri Lanka ??? ERA5 ??? District Daily & Weekly Features
# File: analysis/sri_lanka/era5_to_weekly_features.R
#
# PURPOSE
#   Build epidemiology-ready climate covariates from ERA5 hourly fields:
#   (1) aggregate to cell×day in local time (Asia/Colombo),
#   (2) area-weight to ADM2 districts, and
#   (3) summarize to epi-weeks aligned with WER (date_start/date_end), with
#       lags, rolling windows, anomalies, and derived moisture/heat metrics.
#
# DATA LOCATION & ACCESS (S3 + local)
#   . Primary store: Amazon S3 bucket
#       s3://dghi-chi/data/se-asia/sri-lanka-disease-surveillance/era5/{YEAR}/*.parquet
#   . Access path in this script: `era5_root` (can be S3 or local path).
#   . Read layer: Apache Arrow `open_dataset()` reads Parquet either from:
#       - S3 (data is public, no credentials required)

# WHAT THIS SCRIPT PRODUCES
#   1) District × day climate table (area-weighted):
#        analysis/sri_lanka/srilanka_district_daily_era5_areawt.csv
#   2) District × epi-week feature table (WER weeks):
#        analysis/sri_lanka/srilanka_district_weekly_era5_areawt.csv
#
# WHO USES THESE OUTPUTS
#   . DGHI CHI climate-health modeling (e.g., lepto & dengue early-warning).
#   . Exploratory analyses (distributed lags, thresholding, onset modeling).
#
# KEY DESIGN CHOICES & ASSUMPTIONS
#   . Time: ERA5 hours (UTC) are localized to Asia/Colombo (+05:30) *before*
#     daily aggregation so days reflect health-relevant local timing.
#   . Spatial: district values are area-weighted using cell×district
#     intersections computed in an equal-area CRS (EPSG:6933) to avoid bias.
#   . Precip artifacts: negative totals/rates are clamped to zero.n
#   . Week coverage: per-week metrics are masked if <MIN_DAYS_PER_WEEK (default 5).
#   . Units:
#       - Temperature/TD/WBGT: °C (scaled×10 in raw ??? unscaled here)
#       - Wind: m s???¹ (scaled×10 in raw ??? unscaled here)
#       - Radiation: ssrd in MJ m???² per hour ??? daily means/min/max + weekly mean
#       - Precipitation: `tp` in meters (summed to daily) and `mtpr` as rate
#         (treated as daily totals after sum; keep consistent downstream)
#
# GOVERNANCE & REPRODUCIBILITY
#   . Code: deterministic transformations; no in-place mutation of raw data.
#   . Caching: cell×district weights are static; consider persisting to RDS/Parquet.
#
# AUTHOR: Jordan Clark (DGHI CHI)
# LAST UPDATED: 2025-09-24
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
# 0) CONFIG  ??? Paths, project layout, and one-time directory setup
# WHAT:  Defines project-relative paths (scratch, outputs, figures) and points
#        `era5_root` to S3 or local DVC cache. Creates scratch dirs if missing.
# WHY:   Keeps code portable across dev laptops, EC2, and CI by avoiding setwd()
#        and absolute paths; allows switching between S3 and local seamlessly.
# WHO:   Anyone running the pipeline end-to-end; CI jobs; students onboarding.
# NOTE:  `paths$outputs$era5_weekly_aggregated` currently points to the *daily*
#        CSV. Consider renaming that key to avoid confusion.
# ------------------------------------------------------------------------------
# 1) No setwd() and no absolute C:\ paths
#    Everything below is relative to the Sri Lanka subproject root.

# Where to put scratch outputs during processing
paths <- list(
  temp_dir = file.path(cfg$paths$intermediate, "temp")
)

# 
# 2) Source inputs that should live in the project repo:
#    Recommend moving PDFs/CSVs under data/raw or data/raw/external.
#    For now, keep your current filenames but make them relative.
paths$midyear_pop <- file.path(cfg$paths$raw, "Mid-year_population_by_district_and_sex_2024.pdf")
paths$wx_stations <- file.path(cfg$paths$raw , "SriLanka_Weather_Dataset.csv")
paths$era5_daily  <- file.path(cfg$paths$raw, "srilanka_district_daily_era5_areawt.csv")


# 3) Output figures directory inside project
paths$fig_dir <- file.path(cfg$paths$reports, "figures")

# 4) Additional named outputs (project-relative)
paths$outputs <- list(
  era5_weekly_aggregated = file.path(cfg$path$intermediate, 'srilanka_district_daily_era5_areawt.csv'),
  pdf_index_csv   = file.path(cfg$path$intermediate,
                      "sri_lanka_WER_index_of_pdfs.csv"),
  case_counts_txt = file.path(cfg$path$intermediate,
                      "disease_counts.txt")
)

# DEFINE WHERE ERA5 (ERA5 OR ERA5-LAND) TO USE.
ERA5_TYPE = 'era5land'
# ERA5_TYPE = 'era5'

# 5) Define ERA5 data root location (S3 bucket - publically accessible).
if (ERA5_TYPE == 'era5land'){
  # FOR ERA5-LAND
  era5_root <- 's3://dghi-chi/data/se-asia/sri-lanka-disease-surveillance/era5land'
} else if (ERA5_TYPE == 'era5'){
  # For ERA5.
  era5_root <- 's3://dghi-chi/data/se-asia/sri-lanka-disease-surveillance/era5'
}



# 6) Create any needed directories once (safe if they already exist)
ensure_dir <- function(...) dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE)
ensure_dir(cfg$paths$intermediate, "temp")
ensure_dir(cfg$paths$intermediate, "outputs")
ensure_dir(cfg$paths$reports, "figures")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 1) COLUMN MAP  ??? ERA5 variable names and scaling conventions
# WHAT:  Centralized mapping of expected column names in Parquet and how we
#        unscale/convert them (e.g., °C×10 ??? °C). Enumerates precipitation vars.
# WHY:   ERA5 derivatives (and our pre-processing) sometimes differ in names
#        and scaling; a single map prevents drift and simplifies future changes.
# WHO:   Script maintainers; anyone aligning upstream parquet schemas.
# ------------------------------------------------------------------------------
# 
# # No longer used - simply pulling in all variables instead.
# VAR_VALIDTIME <- "validTime"             # epoch seconds (UTC) or timestamp
# VAR_SSRD      <- "ssrd"                  # surface solar radiation downwards (J/m² per hour)
# VAR_TA10      <- "ta_scaled10_degC"     # air temp (°C * 10)
# VAR_TD10      <- "td_scaled10_degC"     # dewpoint (°C * 10)
# VAR_WBGT10    <- "wbgt_scaled10_degC"   # WBGT (°C * 10); may be absent
# VAR_WIND10    <- "wind2m_scaled10_ms1"  # wind (m/s * 10); fallback handled below
# VAR_TP        <- "tp"                    # total precipitation (m); sum to daily

if (ERA5_TYPE == 'era5land'){
  # FOR ERA5-LAND
  # 0.1° grid center index  lon/lon = idx/10
  VAR_LATIDX    <- "lat_idx_x10"
  VAR_LONIDX    <- "lon_idx_x10"           
} else if (ERA5_TYPE == 'era5'){
  # FOR ERA5 
  # 0.25° grid center index  lon = idx/4
  VAR_LATIDX    <- "lat_idx_x4"
  VAR_LONIDX    <- "lon_idx_x4"
}

# ------------------------
# Timezone handling:
#   ERA5 validTime is UTC; we localize to Asia/Colombo (+5:30) BEFORE daily agg
tz_local      <- "Asia/Colombo"
tz_offset_sec <- as.integer(5.5 * 3600)


# ------------------------
# Analysis years (adjust as needed)
years_all   <- 2014:2024

# ------------------------------------------------------------------------------
# 2) SMALL HELPERS  ??? Math utilities for scaling, robust stats, and weighting
# WHAT:  scale10(), robust z-scores, great-circle area approx, area-weighted mean.
# WHY:   Encapsulate repeated numeric patterns and edge-case handling (all NA,
#        zero-variance, etc.). Keeps main pipeline readable.
# WHO:   Internal use by downstream aggregation steps.
# ------------------------------------------------------------------------------
scale10 <- function(x) as.numeric(x) / 10

# Simple "workability" curve from WBGT (kept for compatibility; not used here)
workability <- function(wbgt, steepness = 0.6, midpoint = 30) {
  100 * (1 - 1 / (1 + exp(-steepness * (wbgt - midpoint))))
}

# Great-circle small-distance approximation
dist2_m2 <- function(lon1, lat1, lon2, lat2) {
  r <- 6371000
  x <- (lon2 - lon1) * cos((lat1 + lat2) * pi / 360)
  y <- (lat2 - lat1)
  (r * pi / 180)^2 * (x^2 + y^2)
}

# Area-weighted mean helper
aw_mean <- function(val, w) {
  i <- is.finite(val)
  if (!any(i)) return(NA_real_)
  sum(val[i] * w[i]) / sum(w[i])
}

# ------------------------------------------------------------------------------
# 3) ADMIN BOUNDARIES (ADM2)  ??? Read, normalize, and dissolve to districts
# WHAT:  Downloads GADM 4.1 LKA level-2, normalizes names, and dissolves to
#        district polygons; reprojects to equal-area for intersection areas.
# WHY:   Accurate area weights require valid, clean district shapes in an EA CRS.
# WHO:   Core to area-weighting cell???district.
# NOTE:  The current dissolve groups by `province` but renames the result to
#        `district`. If the intent is *district* shapes, change:
#          dplyr::select(district, geometry) |>
#          dplyr::group_by(district) |> dplyr::summarise(...)
#        Otherwise you'll compute province-level weights with district labels.
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
# 4) GRID POLYGONS  ??? Construct ERA5 0.25° cell polygons from center indices
# WHAT:  Reads unique (lat_idx_x4, lon_idx_x4) from any year, converts centers
#        to edges, and builds polygons; computes cell area in equal-area CRS.
# WHY:   We need true polygon geometries for precise cell×district intersections
#        and correct area weights (point-in-poly is insufficient at borders).
# WHO:   Spatial core; computed once per domain.
# ------------------------------------------------------------------------------
# Read any year to get unique (lat_idx_x4, lon_idx_x4) pairs  centers
build_latlon_lookup <- function() {
  yr <- 2020
  ds <- open_dataset(file.path(era5_root, yr), format = "parquet",
                     factory_options = list(exclude_invalid_files = TRUE))
  
  # Expect index columns; we compute lat/lon = idx/4 (for 0.25° grid)
  # needed <- c("lat_idx_x4","lon_idx_x4")
  # 0.25° grid center index  lon = idx/4
  
  needed <- c(VAR_LATIDX, VAR_LONIDX)
  if (!all(needed %in% names(ds))) {
    stop("Parquet missing lat_idx_x4 or lon_idx_x4; include them or add a lookup.")
  }
  
  coords <- ds |>
    select(dplyr::all_of(needed)) |>
    distinct() |>
    collect()
  
  setDT(coords)
  
  if (all(needed %like% "x10")){
    coords[, `:=`(
      lat = as.integer(lat_idx_x10)/10,
      lon = as.integer(lon_idx_x10)/10,
      lat_idx_x10 = as.integer(lat_idx_x10)/10,
      lon_idx_x10 = as.integer(lon_idx_x10)/10
    )]
  } else {
    coords[, `:=`(
      lat = as.integer(lat_idx_x4)/4,
      lon = as.integer(lon_idx_x4)/4,
      lat_idx_x4 = as.integer(lat_idx_x4)/4,
      lon_idx_x4 = as.integer(lon_idx_x4)/4
    )]
  }
  
  coords[]
}

coords_lookup <- build_latlon_lookup()

# Build polygons from center coords, dynamically keeping index cols
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
  
  # Retain whichever index columns exist (x4 or x10)
  idx_cols <- grep("^lat_idx_|^lon_idx_", names(coords), value = TRUE)
  
  cell_sf <- st_sf(
    lat = coords$lat, lon = coords$lon,
    coords[, ..idx_cols],
    geometry = st_sfc(polys, crs = 4326)
  )
  
  cell_sf <- st_make_valid(cell_sf)
  cell_sf_m <- st_transform(cell_sf, EA_CRS)
  cell_sf_m$cell_area_m2 <- as.numeric(st_area(cell_sf_m))
  cell_sf_m
}

cell_sf_m <- make_cell_polygons(coords_lookup)

# ------------------------------------------------------------------------------
# 5) WEIGHTS  ??? Cell×district intersection areas (m²)
# WHAT:  Intersects every cell polygon with every district polygon, recording
#        the m² of overlap (per-pair weights). Produces a lookup table keyed by
#        (lat_idx_x4, lon_idx_x4, district).
# WHY:   These weights are reused across years/variables; separating the step
#        avoids repeated heavy spatial ops and guarantees consistency.
# WHO:   All downstream district aggregations rely on this table.
# PERF:  Cache to Parquet/RDS for re-runs; this is deterministic for domain/CRS.
# ------------------------------------------------------------------------------
compute_cell_district_weights <- function(cell_sf_m, districts_sf_m) {
  # Detect which index columns are present
  idx_cols <- grep("^lat_idx_|^lon_idx_", names(cell_sf_m), value = TRUE)
  if (length(idx_cols) != 2L) {
    stop("Expected exactly two index columns (lat/lon), got: ", paste(idx_cols, collapse = ", "))
  }
  
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
    
    # keep index cols + area
    sub_cells <- cell_sf_m[pj$cell_row, c(idx_cols, "cell_area_m2")]
    inter     <- suppressWarnings(st_intersection(sub_cells, districts_sf_m[j, ]))
    if (nrow(inter) == 0) next
    
    inter$area_m2 <- as.numeric(st_area(inter))
    inter <- inter[inter$area_m2 > 0, ]
    if (nrow(inter) == 0) next
    
    tmp <- as.data.table(inter)[, ..idx_cols]   # keep index cols dynamically
    tmp[, district := districts_sf_wgs84$district[j]]
    tmp[, area_m2  := inter$area_m2]
    
    out[[j]] <- tmp
  }
  
  weights <- rbindlist(out, use.names = TRUE, fill = TRUE)
  setkeyv(weights, c(idx_cols, "district"))
  weights
}

weights <- compute_cell_district_weights(cell_sf_m, districts_sf_m)

# ------------------------------------------------------------------------------
# 6) READ YEAR (cell×hour)  ??? S3/local Arrow read + unit normalization
# WHAT:  Reads `era5_root/{year}` Parquet via Arrow; localizes timestamps to
#        Asia/Colombo (+5:30); unscales temps/wind; derives RH/VPD; converts
#        ssrd to MJ m???² h???¹; clamps negative precipitation to zero.
# WHY:   Ensures all hourly inputs are in epidemiology-relevant local time and
#        consistent physical units before daily aggregation.
# WHO:   Annual batch step; used by the driver in Section 9.
# ------------------------------------------------------------------------------
read_one_year_raw <- function(yr) {
  ydir <- if (grepl("^s3://", era5_root)) sprintf("%s/%d", era5_root, yr) else file.path(era5_root, yr)
  ds <- open_dataset(ydir, format = "parquet",
                     factory_options = list(exclude_invalid_files = TRUE))
  
  # Optional inputs present in some schemas
  HAS <- function(x) x %in% names(ds)
  
  # Wind: two possible names
  VAR_WIND10_LOCAL <- if ("wind2m_scaled10_ms1" %in% names(ds)) "wind2m_scaled10_ms1" else "wind2m_scaled10_ms"
  
  # Build the base mutation (timestamp ??? local date/hour). Always available: VAR_VALIDTIME.
  tbl <- ds |>
    mutate(
      validTime = cast(cast(!!sym(VAR_VALIDTIME), int64()), timestamp("s", timezone = "UTC")),
      local_ts  = cast(cast(validTime, int64()) + tz_offset_sec, timestamp("s")),
      date      = cast(local_ts, date32()),
      hour      = hour(local_ts)
    )
  
  # Add derived variables ONLY if sources exist
  if (HAS(VAR_TA10)) {
    tbl <- tbl |> mutate(ta = round(cast(!!sym(VAR_TA10), float64()) / 10, 1))
  }
  if (HAS(VAR_TD10)) {
    tbl <- tbl |> mutate(td = round(cast(!!sym(VAR_TD10), float64()) / 10, 1))
  }
  if (HAS(VAR_WBGT10)) {
    tbl <- tbl |> mutate(wbgt = round(cast(!!sym(VAR_WBGT10), float64()) / 10, 1))
  }
  if (HAS(VAR_WIND10_LOCAL)) {
    tbl <- tbl |> mutate(wind = round(cast(!!sym(VAR_WIND10_LOCAL), float64()) / 10, 1))
  }
  if (HAS(VAR_TP)) {
    tbl <- tbl |> mutate(tp = cast(!!sym(VAR_TP), float64()))
  }
  if (HAS(VAR_SSRD)) {
    tbl <- tbl |> mutate(ssrd_MJ = cast(!!sym(VAR_SSRD), float64()) / 1e6)
  }
  # RH / VPD only if both ta and td exist
  if (HAS(VAR_TA10) && HAS(VAR_TD10)) {
    tbl <- tbl |>
      mutate(
        es  = 0.6108 * exp(17.27 * ta / (ta + 237.3)),
        ea  = 0.6108 * exp(17.27 * td / (td + 237.3)),
        vpd = if_else(es - ea > 0, es - ea, 0.0),
        rh  = exp((17.62 * td) / (243.12 + td) - (17.62 * ta) / (243.12 + ta)) * 100
      )
  }
  
  # Optional raw fields we just pass through if present
  opt_raw <- intersect(
    c("slhf","sshf","evabs","ro","ssro","swvl1","swvl2","swvl3","swvl4"),
    names(ds)
  )
  
  # Index columns must exist per your config (here using VAR_LATIDX/VAR_LONIDX);
  # if you occasionally get x4 instead of x10, swap the config OR detect dynamically.
  base_select <- c(VAR_LATIDX, VAR_LONIDX, "date", "hour")
  
  # Only keep derived columns that were actually created above
  derived_present <- intersect(
    c("ta","td","wbgt","wind","rh","vpd","ssrd_MJ","tp"),
    names(schema(tbl))
  )
  
  sel <- c(base_select, derived_present, opt_raw)
  
  dt <- tbl |>
    select(dplyr::all_of(sel)) |>
    collect()
  
  setDT(dt)
  
  # Clamp negative precip artifacts if present
  if ("tp" %in% names(dt)) dt[tp < 0, tp := 0]
  
  # Sanity check: make sure we kept the two index columns
  idx_cols <- grep("^lat_idx_|^lon_idx_", names(dt), value = TRUE)
  if (length(idx_cols) != 2L) {
    stop("Expected exactly 2 index columns; found: ", paste(idx_cols, collapse = ", "))
  }
  
  dt
}



# ------------------------------------------------------------------------------
# 7) HOURLY ??? DAILY (cell)  ??? Daily means/mins/maxes; precip sums
# WHAT:  Aggregates each cell's 24 hourly records into daily summaries for
#        temperature-related vars (mean/min/max), moisture/radiation (means),
#        and precipitation totals (sums of `tp` and `mtpr`).
# WHY:   Daily is the natural bridge to epi-weeks and many vector/pathogen lags.
# WHO:   Used immediately before area-weighting to district.
# ------------------------------------------------------------------------------
cell_day_stats <- function(dt_hourly) {
  if (nrow(dt_hourly) == 0) return(dt_hourly)
  
  # Detect grid index columns (works for x4 or x10)
  idx_cols <- grep("^lat_idx_|^lon_idx_", names(dt_hourly), value = TRUE)
  if (length(idx_cols) != 2L) stop("Expected 2 index columns in dt_hourly")
  
  # Variables to aggregate
  vars_mmm <- intersect(
    c("ta","td","wbgt","wind","rh","vpd",
      "swvl1","swvl2","swvl3","swvl4",
      "ssrd_MJ","sshf","slhf"),
    names(dt_hourly)
  )
  vars_sum <- intersect(
    c("tp","sro","ssro","ro","evabs"),
    names(dt_hourly)
  )
  
  # --- If sshf/slhf are accumulated J m^-2 per hour and you prefer daily mean W m^-2:
  #     You could pre-compute hourly W m^-2 here as dt_hourly[, sshf := sshf/3600] etc.,
  #     then treat them as vars_mmm. Leaving as-is (mmm of hourly values).
  
  # Mean/Min/Max per day (for state variables and flux means)
  mmm <- if (length(vars_mmm)) {
    dt_hourly[, c(
      setNames(lapply(.SD, \(x) mean(x, na.rm = TRUE)), paste0(vars_mmm, "_mean")),
      setNames(lapply(.SD, \(x) min(x,  na.rm = TRUE)), paste0(vars_mmm, "_min")),
      setNames(lapply(.SD, \(x) max(x,  na.rm = TRUE)), paste0(vars_mmm, "_max"))
    ), by = c(idx_cols, "date"), .SDcols = vars_mmm]
  } else {
    unique(dt_hourly[, c(idx_cols, "date"), with = FALSE])  # empty shell to merge onto
  }
  
  # Sums per day (for additive water/energy depths)
  sums <- if (length(vars_sum)) {
    dt_hourly[, c(
      setNames(lapply(.SD, \(x) sum(x, na.rm = TRUE)), paste0(vars_sum, "_sum"))
    ), by = c(idx_cols, "date"), .SDcols = vars_sum]
  } else {
    NULL
  }
  
  out <- if (!is.null(sums)) merge(mmm, sums, by = c(idx_cols, "date"), all = TRUE) else mmm
  # setorder(out, idx_cols, "date")
  out[]
}
# ------------------------------------------------------------------------------
# 8) cell???district (daily)  ??? Area-weighted aggregation using m² weights
# WHAT:  Joins cell-daily stats to weights and computes district-daily metrics:
#        . For "level" vars (e.g., ta_mean), uses area-weighted mean.
#        . For extremes (e.g., ta_max), uses min/max over member cells. 
#        . For daily precip sums, uses area-weighted mean of cell totals.
# WHY:   Mimics spatial integration over irregular district shapes without bias.
# WHO:   Produces the canonical district×day climate table.
# NOTE:  min/max are not area-weighted by design (we want spatial extremes).
# ------------------------------------------------------------------------------
district_day_from_cell_day <- function(cell_daily, weights) {
  if (nrow(cell_daily) == 0) return(cell_daily)
  
  # Detect index columns
  idx_cols <- grep("^lat_idx_|^lon_idx_", names(cell_daily), value = TRUE)
  if (length(idx_cols) != 2L) stop("Expected 2 index columns in cell_daily")
  
  # Join weights
  x <- merge(cell_daily, weights, by = idx_cols, allow.cartesian = TRUE)
  
  # Base variable names present in mmm/sum outputs
  # (strip suffixes to identify bases)
  bases_mmm <- unique(sub("_(mean|min|max)$", "", grep("_(mean|min|max)$", names(x), value = TRUE)))
  bases_sum <- unique(sub("_sum$", "", grep("_sum$", names(x), value = TRUE)))
  
  # Area-weighted aggregation by district × date
  x[, {
    # For each base in mmm: area-weighted mean of *_mean ; plain min of *_min ; plain max of *_max
    means <- if (length(bases_mmm)) {
      setNames(lapply(bases_mmm, \(v) aw_mean(get(paste0(v, "_mean")), area_m2)),
               paste0(bases_mmm, "_mean"))
    } else list()
    
    mins <- if (length(bases_mmm)) {
      setNames(lapply(bases_mmm, \(v) suppressWarnings(min(get(paste0(v, "_min")), na.rm = TRUE))),
               paste0(bases_mmm, "_min"))
    } else list()
    
    maxs <- if (length(bases_mmm)) {
      setNames(lapply(bases_mmm, \(v) suppressWarnings(max(get(paste0(v, "_max")), na.rm = TRUE))),
               paste0(bases_mmm, "_max"))
    } else list()
    
    # For sums (tp/sro/ssro/ro/evabs): area-weighted mean of cell totals ??? district daily areal mean
    sums <- if (length(bases_sum)) {
      setNames(lapply(bases_sum, \(v) aw_mean(get(paste0(v, "_sum")), area_m2)),
               paste0(bases_sum, "_sum"))
    } else list()
    
    as.list(c(means, mins, maxs, sums))
  }, by = .(district, date)][order(district, date)]
}


# ------------------------------------------------------------------------------
# 9) DRIVER  ??? Loop years ??? hourly read ??? cell-daily ??? district-daily
# WHAT:  Orchestrates the annual workflow and binds results; attaches `year`.
# WHY:   Batch execution with clear checkpoints (helpful for CI and resume).
# WHO:   Primary entrypoint when generating/refreshing the daily CSV.
# I/O:   Writes srilanka_district_daily_era5_areawt.csv to project outputs.
# ------------------------------------------------------------------------------
# summarize_years_area_weighted <- function(years) {
#   res_list <- vector("list", length(years))
#   for (i in seq_along(years)) {
#     yr <- years[i]
#     message("Processing ", yr, " .")
#     hr <- read_one_year_raw(yr)
#     if (!all(c("lat_idx_x4","lon_idx_x4","date") %in% names(hr))) {
#       stop("Hourly frame missing lat_idx_x4/lon_idx_x4/date after read for year ", yr)
#     }
#     cd <- cell_day_stats(hr)
#     dd <- district_day_from_cell_day(cd, weights)
#     dd[, year := yr]
#     res_list[[i]] <- dd
#   }
#   rbindlist(res_list, use.names = TRUE, fill = TRUE)
# }

# ------------------------------------------------------------------------------
# Driver: orchestrate year-by-year workflow
# ------------------------------------------------------------------------------
summarize_years_area_weighted <- function(years, weights_tbl = NULL) {
  if (is.null(weights_tbl)) {
    weights_tbl <- get0("weights", inherits = TRUE)
    if (is.null(weights_tbl)) stop("No 'weights' provided and none found in parent/global env.")
  }
  
  res_list <- vector("list", length(years))
  for (i in seq_along(years)) {
    yr <- years[i]
    message("Processing ", yr, " .")
    
    hr <- read_one_year_raw(yr)  # robust to missing cols
    setDT(hr)
    
    # Ensure we have a Date column (read_one_year_raw adds it)
    if (!"date" %in% names(hr)) stop("No 'date' in hourly data for year ", yr)
    
    # Index columns present?
    idx_cols <- grep("^lat_idx_|^lon_idx_", names(hr), value = TRUE)
    if (length(idx_cols) != 2L) {
      stop("Hourly frame missing grid index columns for year ", yr,
           ". Found: ", paste(names(hr), collapse = ", "))
    }
    
    # Daily cell aggregates (robust)
    cd <- cell_day_stats(hr)
    
    # --- Ensure join keys match the weights' key space (degrees vs indices) ---
    w_idx <- grep("^lat_idx_|^lon_idx_", names(weights_tbl), value = TRUE)
    if (length(w_idx) != 2L) stop("Weights table must have two index columns.")
    
    # If cd indices look like integers but weights' indices look like decimals, rescale cd
    intish <- function(x) all(is.finite(x)) && mean(abs(x - round(x)), na.rm = TRUE) < 1e-6
    decish <- function(x) any(is.finite(x)) && mean(abs(x - round(x)), na.rm = TRUE) > 1e-6
    
    need_scale <- intish(cd[[idx_cols[1]]]) && decish(weights_tbl[[w_idx[1]]])
    
    if (need_scale) {
      # Pick scale from suffix (x10 or x4)
      if (grepl("x10$", idx_cols[1])) {
        scl <- 10
      } else if (grepl("x4$", idx_cols[1])) {
        scl <- 4
      } else {
        stop("Cannot infer scale from index column names: ", paste(idx_cols, collapse = ", "))
      }
      cd[, (idx_cols) := lapply(.SD, function(z) as.numeric(z) / scl), .SDcols = idx_cols]
    }
    
    dd <- district_day_from_cell_day(cd, weights_tbl)
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
# 10) WEEKLY FEATURES  ??? District×week covariates aligned to WER windows
# WHAT:  Joins district-daily climate onto epi-week calendars (WER date_start/
#        date_end) and computes weekly summaries (means/sums, 3-day maxima,
#        wet-day counts/spells, radiation), plus:
#          . lags (1..MAX_LAG_WEEKS)
#          . rolling means/sums (2w, 4w) with partial coverage rules
#          . exponentially-weighted antecedent precip (EWAP)
#          . climatologies (district×ISO week) ??? anomalies and %-of-normal
#          . robust within-district z-scores for selected vars
# WHY:   Produces epidemiology-ready covariates for onset/spike models without
#        leaking future information (all lags/rolls are past-only).
# WHO:   Modeling and visualization teams; inputs to lepto/dengue GAM/GLM/logit.
# QA:    Masks weeks with <MIN_DAYS_PER_WEEK to avoid biased summaries.
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

# # Configuration for features
# CFG <- list(
#   WINSORIZE    = TRUE,
#   WINSOR_PROBS = c(0.01, 0.99),
#   HOT_TMAX     = 32,    # °C threshold for hot-day counts
#   WET_MM       = 10,    # "wet day" threshold (mm)
#   MAX_LAG_WEEKS= 6,
#   ROLL_WINDOWS = c(2, 4),
#   EWAP_ALPHA   = 0.8,   # exponentially-weighted antecedent precip
#   EWAP_K       = 4,
#   ANOMALY_VARS = c("precip_tp_sum_week","tmax_mean","rh_mean_week"),
#   ZSCORE_VARS  = c("precip_tp_sum_week","precip_mtpr_sum_week","tmax_mean","vpd_mean_week"),
#   MIN_DAYS_PER_WEEK = 5 # require 5/7 days in week to compute metrics
# )

CFG <- list(
  WINSORIZE           = TRUE,
  WINSOR_PROBS        = c(0.01, 0.99),
  HOT_TMAX            = 32,    # °C threshold for hot-day counts
  WET_MM              = 10,    # "wet day" threshold (mm)
  MAX_LAG_WEEKS       = 6,
  ROLL_WINDOWS        = c(2, 4),
  EWAP_ALPHA          = 0.8,   # exponentially-weighted antecedent precip
  EWAP_K              = 4,
  # Use intersection with names(agg) later, so extra entries are harmless.
  ANOMALY_VARS        = c("precip_tp_sum_week",
                          "runoff_ro_sum_week",
                          "tmax_mean",
                          "rh_mean_week","vpd_mean_week",
                          "ssrd_MJ_mean_week","sshf_mean_week","slhf_mean_week",
                          "swvl1_mean_week","swvl2_mean_week","swvl3_mean_week","swvl4_mean_week"),
  ZSCORE_VARS         = c("precip_tp_sum_week",
                          "runoff_ro_sum_week",
                          "tmax_mean","vpd_mean_week",
                          "sshf_mean_week","slhf_mean_week"),
  MIN_DAYS_PER_WEEK   = 5
)


weekly_weather_features <- function(daily_dt, weeks_dt, cfg = CFG) {
  # Required columns in inputs
  stopifnot(all(c("district","date") %in% names(daily_dt)))
  stopifnot(all(c("district","date_start","date_end") %in% names(weeks_dt)))
  
  DTd <- as.data.table(copy(daily_dt))
  DTw <- as.data.table(copy(weeks_dt))
  
  # Ensure types
  if (!inherits(DTd$date, "Date"))        DTd[,  date := as.Date(date)]
  if (!inherits(DTw$date_start, "Date"))  DTw[, date_start := as.Date(date_start)]
  if (!inherits(DTw$date_end, "Date"))    DTw[, date_end   := as.Date(date_end)]
  
  # Optional winsorization of daily inputs (only those present)
  winsor_targets <- intersect(
    c(
      # temps & humidity/flux
      "ta_mean","ta_min","ta_max","wbgt_mean","rh_mean","vpd_mean",
      "ssrd_MJ_mean","sshf_mean","slhf_mean",
      # soil moisture layers
      "swvl1_mean","swvl2_mean","swvl3_mean","swvl4_mean",
      # additive hydrology
      "tp_sum","sro_sum","ssro_sum","ro_sum","evabs_sum"
    ),
    names(DTd)
  )
  if (cfg$WINSORIZE && length(winsor_targets)) {
    DTd[, (winsor_targets) := lapply(.SD, winsor, probs = cfg$WINSOR_PROBS), .SDcols = winsor_targets]
  }
  
  # Build a daily calendar for each epi week; join daily climate onto it
  idx <- DTw[, .(district, week_id = .I,
                 date_start, date_end,
                 year = data.table::year(date_end),
                 week_of_year = iso_week(date_end))]
  
  wk_days <- idx[
    , .(date = seq(date_start, date_end, by = "1 day")),
    by = .(district, week_id, date_start, date_end, year, week_of_year)
  ]
  
  setkey(DTd, district, date)
  setkey(wk_days, district, date)
  wk <- DTd[wk_days, on = .(district, date)]
  
  # Helper to safely compute weekly mean of a daily column if present
  wmean_col <- function(col) if (col %in% names(wk)) mean_na(get(col)) else NA_real_
  wsum_col  <- function(col) if (col %in% names(wk)) sum_na(get(col))  else NA_real_
  wroll3sum_max <- function(col) {
    if (!(col %in% names(wk))) return(NA_real_)
    x <- get(col)
    if (sum(!is.na(x)) >= 3) max(zoo::rollapplyr(x, 3, sum, na.rm = TRUE), na.rm = TRUE) else NA_real_
  }
  wwet_spell <- function(col, thr = cfg$WET_MM) {
    if (!(col %in% names(wk))) return(NA_real_)
    as.numeric(longest_wet_spell(get(col), thr = thr))
  }
  wcount_ge <- function(col, thr) {
    if (!(col %in% names(wk))) return(NA_real_)
    as.numeric(sum(get(col) >= thr, na.rm = TRUE))
  }
  
  # Weekly aggregation (district × week_id)
  agg <- wk[, {
    ndays <- sum(!is.na(date))
    
    # Temperature: daily columns expected from earlier daily aggregation
    #   ta_max, ta_min, ta_mean  ??? weekly summaries below
    tmax_mean <- if ("ta_max" %in% names(wk)) mean_na(ta_max) else NA_real_
    tmax_p90  <- if ("ta_max" %in% names(wk)) q_na(ta_max, 0.90) else NA_real_
    tmax_p95  <- if ("ta_max" %in% names(wk)) q_na(ta_max, 0.95) else NA_real_
    tmax_rng  <- if ("ta_max" %in% names(wk)) {
      if (all(is.na(ta_max))) NA_real_ else max(ta_max, na.rm = TRUE) - min(ta_max, na.rm = TRUE)
    } else NA_real_
    tmin_mean <- if ("ta_min" %in% names(wk)) mean_na(ta_min) else NA_real_
    tmean_mean<- if ("ta_mean" %in% names(wk)) mean_na(ta_mean) else NA_real_
    
    # Core means (weekly means of daily means)
    wbgt_mean_week   <- wmean_col("wbgt_mean")
    rh_mean_week     <- wmean_col("rh_mean")
    vpd_mean_week    <- wmean_col("vpd_mean")
    ssrd_MJ_mean_week<- wmean_col("ssrd_MJ_mean")
    sshf_mean_week   <- wmean_col("sshf_mean")
    slhf_mean_week   <- wmean_col("slhf_mean")
    
    # Soil moisture (weekly means of daily means)
    swvl1_mean_week  <- wmean_col("swvl1_mean")
    swvl2_mean_week  <- wmean_col("swvl2_mean")
    swvl3_mean_week  <- wmean_col("swvl3_mean")
    swvl4_mean_week  <- wmean_col("swvl4_mean")
    
    # Additive hydrology (weekly sums of daily sums)
    precip_tp_sum_week <- wsum_col("tp_sum")
    runoff_sro_sum_week<- wsum_col("sro_sum")
    runoff_ssro_sum_week<- wsum_col("ssro_sum")
    runoff_ro_sum_week <- wsum_col("ro_sum")
    evabs_sum_week     <- wsum_col("evabs_sum")
    
    # Wet-day metrics and 3-day maxima (for precipitation; tp only)
    max3d_tp            <- wroll3sum_max("tp_sum")
    wet_days_ge10_tp    <- wcount_ge("tp_sum", cfg$WET_MM)
    
    # Heat exceedances (hot day count by ta_max threshold)
    hot_days_ge32       <- if ("ta_max" %in% names(wk)) wcount_ge("ta_max", cfg$HOT_TMAX) else NA_real_
    
    res <- list(
      n_days_week            = as.integer(ndays),
      # temps
      tmax_mean              = tmax_mean,
      tmax_p90               = tmax_p90,
      tmax_p95               = tmax_p95,
      tmax_range             = tmax_rng,
      tmin_mean              = tmin_mean,
      tmean_mean             = tmean_mean,
      # core
      wbgt_mean_week         = wbgt_mean_week,
      rh_mean_week           = rh_mean_week,
      vpd_mean_week          = vpd_mean_week,
      ssrd_MJ_mean_week      = ssrd_MJ_mean_week,
      sshf_mean_week         = sshf_mean_week,
      slhf_mean_week         = slhf_mean_week,
      # soil moisture
      swvl1_mean_week        = swvl1_mean_week,
      swvl2_mean_week        = swvl2_mean_week,
      swvl3_mean_week        = swvl3_mean_week,
      swvl4_mean_week        = swvl4_mean_week,
      # hydrology sums
      precip_tp_sum_week     = precip_tp_sum_week,
      runoff_sro_sum_week    = runoff_sro_sum_week,
      runoff_ssro_sum_week   = runoff_ssro_sum_week,
      runoff_ro_sum_week     = runoff_ro_sum_week,
      evabs_sum_week         = evabs_sum_week,
      # precip-derived
      max3d_tp               = max3d_tp,
      wet_days_ge10_tp       = wet_days_ge10_tp,
      # heat exceedances
      hot_days_ge32          = hot_days_ge32
    )
    
    # Coverage rule: mask if insufficient days (keep n_days_week visible)
    if (ndays < cfg$MIN_DAYS_PER_WEEK) {
      keep <- "n_days_week"
      res[names(res) != keep] <- lapply(res[names(res) != keep], function(.) NA_real_)
    }
    res
  }, by = .(district, week_id, date_start, date_end, year, week_of_year)]
  
  # Order for rolling ops
  data.table::setorder(agg, district, date_end)
  
  # Lags / rolling windows (past-only) for a curated set
  LAG_VARS <- intersect(
    c("precip_tp_sum_week",
      "runoff_ro_sum_week",
      "tmax_mean","vpd_mean_week","rh_mean_week",
      "sshf_mean_week","slhf_mean_week"),
    names(agg)
  )
  for (v in LAG_VARS) {
    for (L in seq_len(cfg$MAX_LAG_WEEKS)) {
      agg[, paste0(v, "_lag", L) := data.table::shift(get(v), L), by = district]
    }
  }
  for (v in LAG_VARS) {
    agg[, paste0(v, "_roll2w_mean") := roll_mean_partial(get(v), k = 2, min_obs = 1L), by = district]
    agg[, paste0(v, "_roll2w_sum")  := roll_sum_partial (get(v), k = 2, min_obs = 1L), by = district]
    agg[, paste0(v, "_roll4w_mean") := roll_mean_partial(get(v), k = 4, min_obs = 2L), by = district]
    agg[, paste0(v, "_roll4w_sum")  := roll_sum_partial (get(v), k = 4, min_obs = 2L), by = district]
  }
  
  # EWAP (antecedent precipitation) - tp only
  if ("precip_tp_sum_week" %in% names(agg)) {
    K <- cfg$EWAP_K; a <- cfg$EWAP_ALPHA
    agg[, ewap_tp := {
      x <- precip_tp_sum_week; ew <- x
      for (k in 1:K) ew <- ew + (a^k) * data.table::shift(x, k, fill = 0)
      ew
    }, by = district]
  }
  
  # Climatology (district × ISO week) and anomalies / percent-of-normal
  CLIM_VARS <- intersect(cfg$ANOMALY_VARS, names(agg))
  if (length(CLIM_VARS)) {
    clim <- agg[, lapply(.SD, mean_na), by = .(district, week_of_year), .SDcols = CLIM_VARS]
    data.table::setnames(clim, CLIM_VARS, paste0(CLIM_VARS, "_clim"))
    agg <- clim[agg, on = .(district, week_of_year)]
    for (v in CLIM_VARS) {
      vc <- paste0(v, "_clim")
      agg[, paste0(v, "_anom")       := get(v) - get(vc)]
      agg[, paste0(v, "_pct_normal") := pct_of_normal(get(v), get(vc), eps = 1e-6)]
    }
  }
  
  # Robust within-district z-scores
  ZV <- intersect(cfg$ZSCORE_VARS, names(agg))
  if (length(ZV)) {
    for (v in ZV) {
      agg[, paste0(v, "_z") := zscore_robust(get(v)), by = district]
    }
  }
  
  data.table::setorder(agg, district, date_end)
  agg[]
}

weekly_weather_features <- function(daily_dt, weeks_dt, cfg = CFG) {
  # Required columns in inputs
  stopifnot(all(c("district","date") %in% names(daily_dt)))
  stopifnot(all(c("district","date_start","date_end") %in% names(weeks_dt)))
  
  DTd <- data.table::as.data.table(data.table::copy(daily_dt))
  DTw <- data.table::as.data.table(data.table::copy(weeks_dt))
  
  # Ensure types
  if (!inherits(DTd$date, "Date"))        DTd[,  date := as.Date(date)]
  if (!inherits(DTw$date_start, "Date"))  DTw[, date_start := as.Date(date_start)]
  if (!inherits(DTw$date_end, "Date"))    DTw[, date_end   := as.Date(date_end)]
  
  # Winsorize daily inputs (only those present)
  winsor_targets <- intersect(
    c(
      # temps & humidity/flux
      "ta_mean","ta_min","ta_max","wbgt_mean","rh_mean","vpd_mean",
      "ssrd_MJ_mean","sshf_mean","slhf_mean",
      # soil moisture layers
      "swvl1_mean","swvl2_mean","swvl3_mean","swvl4_mean",
      # additive hydrology (daily totals)
      "tp_sum","sro_sum","ssro_sum","ro_sum","evabs_sum"
    ),
    names(DTd)
  )
  if (isTRUE(cfg$WINSORIZE) && length(winsor_targets)) {
    DTd[, (winsor_targets) := lapply(.SD, winsor, probs = cfg$WINSOR_PROBS), .SDcols = winsor_targets]
  }
  
  # Build a daily "calendar" for each epi week; join daily climate onto it
  idx <- DTw[, .(district, week_id = .I,
                 date_start, date_end,
                 year = data.table::year(date_end),
                 week_of_year = iso_week(date_end))]
  
  wk_days <- idx[
    , .(date = seq(date_start, date_end, by = "1 day")),
    by = .(district, week_id, date_start, date_end, year, week_of_year)
  ]
  
  data.table::setkey(DTd, district, date)
  data.table::setkey(wk_days, district, date)
  wk <- DTd[wk_days, on = .(district, date)]
  
  # Aggregate within each epi week (district × week_id)
  agg <- wk[, {
    ndays <- sum(!is.na(date))
    
    # ---- helpers that read from the current group (.SD) safely ----
    gmean <- function(nm)  if (nm %in% names(.SD)) mean_na(.SD[[nm]]) else NA_real_
    gsum  <- function(nm)  if (nm %in% names(.SD))  sum_na(.SD[[nm]]) else NA_real_
    gq    <- function(nm,p)if (nm %in% names(.SD))     q_na(.SD[[nm]], p) else NA_real_
    grng  <- function(nm)  if (nm %in% names(.SD)) {
      x <- .SD[[nm]]
      if (all(is.na(x))) NA_real_ else max(x, na.rm=TRUE) - min(x, na.rm=TRUE)
    } else NA_real_
    gcount_ge <- function(nm, thr) if (nm %in% names(.SD)) as.numeric(sum(.SD[[nm]] >= thr, na.rm=TRUE)) else NA_real_
    groll3sum_max <- function(nm) {
      if (!(nm %in% names(.SD))) return(NA_real_)
      x <- .SD[[nm]]
      if (sum(!is.na(x)) >= 3) max(zoo::rollapplyr(x, 3, sum, na.rm = TRUE), na.rm = TRUE) else NA_real_
    }
    gwet_spell <- function(nm, thr) {
      if (!(nm %in% names(.SD))) return(NA_real_)
      as.numeric(longest_wet_spell(.SD[[nm]], thr = thr))
    }
    # ----------------------------------------------------------------
    
    res <- list(
      n_days_week            = as.integer(ndays),
      # Temperature summaries (°C)
      tmax_mean              = gmean("ta_max"),
      tmax_p90               = gq("ta_max", 0.90),
      tmax_p95               = gq("ta_max", 0.95),
      tmax_range             = grng("ta_max"),
      tmin_mean              = gmean("ta_min"),
      tmean_mean             = gmean("ta_mean"),
      # Core means (weekly means of daily means)
      wbgt_mean_week         = gmean("wbgt_mean"),
      rh_mean_week           = gmean("rh_mean"),
      vpd_mean_week          = gmean("vpd_mean"),
      ssrd_MJ_mean_week      = gmean("ssrd_MJ_mean"),
      sshf_mean_week         = gmean("sshf_mean"),
      slhf_mean_week         = gmean("slhf_mean"),
      # Soil moisture (weekly means of daily means)
      swvl1_mean_week        = gmean("swvl1_mean"),
      swvl2_mean_week        = gmean("swvl2_mean"),
      swvl3_mean_week        = gmean("swvl3_mean"),
      swvl4_mean_week        = gmean("swvl4_mean"),
      # Hydrology (weekly sums of daily sums)
      precip_tp_sum_week     = gsum("tp_sum"),
      runoff_sro_sum_week    = gsum("sro_sum"),
      runoff_ssro_sum_week   = gsum("ssro_sum"),
      runoff_ro_sum_week     = gsum("ro_sum"),
      evabs_sum_week         = gsum("evabs_sum"),
      # Precip-derived
      max3d_tp               = groll3sum_max("tp_sum"),
      wet_days_ge10_tp       = gcount_ge("tp_sum", cfg$WET_MM),
      # Heat exceedances
      hot_days_ge32          = gcount_ge("ta_max", cfg$HOT_TMAX)
    )
    
    # Coverage rule: mask if insufficient days (keep n_days_week visible)
    if (ndays < cfg$MIN_DAYS_PER_WEEK) {
      keep <- "n_days_week"
      res[names(res) != keep] <- lapply(res[names(res) != keep], function(.) NA_real_)
    }
    res
  }, by = .(district, week_id, date_start, date_end, year, week_of_year)]
  
  # Ordered by time for rolling ops
  data.table::setorder(agg, district, date_end)
  
  # Lags / rolling windows (past-only) for variables that actually exist
  LAG_VARS <- intersect(
    c("precip_tp_sum_week",
      "runoff_ro_sum_week",
      "tmax_mean","vpd_mean_week","rh_mean_week",
      "sshf_mean_week","slhf_mean_week"),
    names(agg)
  )
  for (v in LAG_VARS) {
    for (L in seq_len(cfg$MAX_LAG_WEEKS)) {
      agg[, paste0(v, "_lag", L) := data.table::shift(get(v), L), by = district]
    }
  }
  for (v in LAG_VARS) {
    agg[, paste0(v, "_roll2w_mean") := roll_mean_partial(get(v), k = 2, min_obs = 1L), by = district]
    agg[, paste0(v, "_roll2w_sum")  := roll_sum_partial (get(v), k = 2, min_obs = 1L), by = district]
    agg[, paste0(v, "_roll4w_mean") := roll_mean_partial(get(v), k = 4, min_obs = 2L), by = district]
    agg[, paste0(v, "_roll4w_sum")  := roll_sum_partial (get(v), k = 4, min_obs = 2L), by = district]
  }
  
  # EWAP (antecedent precipitation) - tp only if present
  if ("precip_tp_sum_week" %in% names(agg)) {
    K <- cfg$EWAP_K; a <- cfg$EWAP_ALPHA
    agg[, ewap_tp := {
      x <- precip_tp_sum_week; ew <- x
      for (k in 1:K) ew <- ew + (a^k) * data.table::shift(x, k, fill = 0)
      ew
    }, by = district]
  }
  
  # Climatology (district × ISO week) and anomalies / percent-of-normal
  CLIM_VARS <- intersect(cfg$ANOMALY_VARS, names(agg))
  if (length(CLIM_VARS)) {
    clim <- agg[, lapply(.SD, mean_na), by = .(district, week_of_year), .SDcols = CLIM_VARS]
    data.table::setnames(clim, CLIM_VARS, paste0(CLIM_VARS, "_clim"))
    agg <- clim[agg, on = .(district, week_of_year)]
    for (v in CLIM_VARS) {
      vc <- paste0(v, "_clim")
      agg[, paste0(v, "_anom")       := get(v) - get(vc)]
      agg[, paste0(v, "_pct_normal") := pct_of_normal(get(v), get(vc), eps = 1e-6)]
    }
  }
  
  # Robust within-district z-scores
  ZV <- intersect(cfg$ZSCORE_VARS, names(agg))
  if (length(ZV)) {
    for (v in ZV) {
      agg[, paste0(v, "_z") := zscore_robust(get(v)), by = district]
    }
  }
  
  data.table::setorder(agg, district, date_end)
  agg[]
}

# ------------------------------------------------------------------------------
# 11) EXAMPLE I/O  ??? End-to-end run helpers
# WHAT:  Illustrates running the pipeline across all years, timing execution,
#        and writing outputs. Includes examples for joining weekly features back
#        to WER case data for modeling.
# WHY:   Copy-paste starter for ad-hoc runs
# ------------------------------------------------------------------------------
jj::timed('start')
# Takes ~1 minute per year for ERA5-Land - J Clark - 20250926 1900 UTC
out_daily <- summarize_years_area_weighted(years_all, weights)
fwrite(out_daily, paths$outputs$era5_weekly_aggregated)



lepto <- fread(paths$outputs$case_counts_txt)
weeks_dt <- unique(lepto[, .(district, date_start, date_end)])  # from your WER table

# weekly level features, including lags and rolling sums. 
features_weekly <- weekly_weather_features(out_daily, weeks_dt)
weekly_out_file <- sprintf("srilanka_district_weekly_%s_areawt.csv", ERA5_TYPE)
fwrite(features_weekly, file.path(cfg$paths$processed, weekly_out_file))


# End Script
#################!
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------





