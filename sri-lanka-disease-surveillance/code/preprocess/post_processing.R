################################################################################
# CHI ??? Sri Lanka WER - Post-Processing Pipeline
# File: analysis/sri_lanka/post_processing.R
#
# Purpose:
#   Consume outputs from the initial WER scrape/extraction and build an
#   analysis-ready panel by merging:
#     . Weekly district disease counts (from initial_processing.R)
#     . Mid-year population (PDF ??? tidy, or CSV fallback) ??? rates per 100k
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
Sys.setenv("JAVA_HOME"="C:/Program Files/Eclipse Adoptium/jdk-17.0.16.8-hotspot")

# To install PDF related packages:
# install.packages("rJava")
# remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")

library(tabulapdf)       # ropensci fork; requires rJava
library(tabulizerjars)

#################
#################


CHI_LOCAL_WORK  <- Sys.getenv("CHI_LOCAL_WORK",  unset = "C:/Users/jordan/Desktop/srilanka")
CHI_GITHUB_ROOT <- Sys.getenv("CHI_GITHUB_ROOT", unset = "C:/Users/jordan/R_Projects/CHI-Data")

paths <- list(
  temp_dir      = file.path(CHI_LOCAL_WORK, "temp"),
  work_out      = file.path(CHI_LOCAL_WORK, "outputs"),
  # Source inputs (checked into repo)
  midyear_pop   = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs/Mid-year_population_by_district_and_sex_2024.pdf"),
  wx_stations   = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs/station_data/SriLanka_Weather_Dataset.csv"),
  era5_daily    = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs/srilanka_district_daily_era5_areawt.csv"),
  fig_dir       = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/outputs/figures")
)

dir.create(paths$temp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$fig_dir, recursive = TRUE, showWarnings = FALSE)

# Outputs from initial_processing and targets from this script
paths = c(paths, 
          
          outputs = list(
            
            pdf_index_csv = file.path(CHI_GITHUB_ROOT, 'analysis/sri_lanka/outputs/sri_lanka_WER_index_of_pdfs.csv'),
            case_counts_txt = file.path(CHI_GITHUB_ROOT, 'analysis/sri_lanka/outputs/disease_counts_v4.txt')
            
          )
)

# -- 0.4 Helpers: naming + parsing utilities -----------------------------------
# Must provide: norm_dist(), .norm(), DIST_CANON, etc.
source(file.path(CHI_GITHUB_ROOT, "/helpers/helpers.R"))
source(file.path(CHI_GITHUB_ROOT, "/analysis/sri_lanka/helpers.R"))


################################################################################
# 1) LOAD WEEKLY DISEASE COUNTS -------------------------------------------------
################################################################################
# SECTION GOAL: Read district-week disease counts, normalize, and precompute dates.
lep <- fread(paths$outputs.case_counts_txt)

# Normalize key fields and derive mid-week date (used for daily joins)
lep[, date_mid := as.IDate(date_start + (as.integer(date_end - date_start) / 2))]
lep[, `:=`(district = norm_dist(district),
           date_mid = as.IDate(date_mid),
           date_end = as.IDate(date_end))]

################################################################################
# 2) POPULATION: PDF ??? TIDY (or CSV fallback) ----------------------------------
################################################################################
# SECTION GOAL: Produce pop_dt with columns [district, year, poptot].
# Preference order:
#   (A) CHI_POP_CSV provided  ??? read directly (no Java)
#   (B) Else parse from PDF via tabulizer (requires Java stack)
# NOTE: Requires working rJava/Tabulizer on your machine; handled with tryCatch.
pop_dt <- NULL
if (file.exists(paths$midyear_pop)) {
  tabs <- tryCatch(extract_tables(paths$midyear_pop, method = "lattice",
                                  guess = TRUE, output = "tibble"),
                   error = function(e) NULL)
  if (!is.null(tabs) && length(tabs) >= 1) {
    yrs <- 2014:2023
    pieces <- list()
    for (pg in seq_along(tabs)) {
      tb <- as.data.table(tabs[[pg]])
      yrheads <- names(tb)[names(tb) %like% paste(yrs, collapse = "|")]
      if (!length(yrheads)) next
      # For each detected year column, grab District names + "Total" col under that year
      for (yh in yrheads) {
        district_names <- tb$District[-1]
        col_idx <- which(tb[1] == "Total")[which(names(tb) == yh)]
        if (!length(col_idx)) next
        vals <- tb[[col_idx]][-1]
        pieces[[length(pieces)+1L]] <- data.table(
          year    = as.integer(gsub("\\*", "", yh)),
          district= district_names,
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


# -- 2.3 Keep only valid districts; fill pre-2014 with 2014 pop ----------------
lep = lep[district %in% pop_dt$district]

# current mid year pop data starts in 2014, so for years in health data before that, simply using 2014 pop for now.
lep[, year2merge := fifelse(year >= 2014, year, 2014L)]
lep <- merge(lep, pop_dt, by.x = c("district","year2merge"), by.y = c("district","year"), all.x = TRUE)

# -- 2.4 Compute rates per 100k -------------------------------------------------
lep[, `:=`(
  lepto_100k  = (lepto  / poptot) * 1e5,
  dengue_100k = (dengue / poptot) * 1e5
)]

lep[, year := NULL] # year no longer needed after merge



################################################################################
# 3) WEATHER: MAP STATIONS ??? DISTRICTS & AGG DAILY -----------------------------
################################################################################
# SECTION GOAL: Assign station observations to districts, then compute
# areal averages per district × day for selected variables.

wx <- fread(paths$wx_stations)
stopifnot(all(c("time","latitude","longitude","city") %in% names(wx)))
if (!inherits(wx$time, "Date")) wx[, time := as.IDate(time)]

# -- Assign each unique (city, lat, lon) to a district via spatial join --------
gadm_zip <- "https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_LKA_shp.zip"
tdir <- tempfile("lka_gadm41_"); dir.create(tdir)
zipfile <- file.path(tdir, "gadm41_LKA_shp.zip")
download.file(gadm_zip, destfile = zipfile, mode = "wb", quiet = TRUE)
unzip(zipfile, exdir = tdir)
adm2_shp <- list.files(tdir, pattern = "^gadm41_LKA_2\\.shp$", full.names = TRUE, recursive = TRUE)
stopifnot(length(adm2_shp) == 1)

adm2 <- st_read(adm2_shp, quiet = TRUE)[, c("GID_2","NAME_1","NAME_2","geometry")]
names(adm2) <- c("gid2","province","district","geometry")
adm2$district <- norm_dist(adm2$district)
adm2 <- st_make_valid(adm2)

stations_lu <- unique(wx[, .(city, latitude, longitude)])
stations_sf <- st_as_sf(stations_lu, coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
stations_sf <- st_transform(stations_sf, st_crs(adm2))
stations_join <- st_join(stations_sf, adm2["district"], left = TRUE)
city2dist <- as.data.table(stations_join)[, .(city, latitude, longitude, district)]

wx_d <- merge(wx, city2dist, by = c("city","latitude","longitude"), all.x = TRUE)
wx_d[, district := norm_dist(district)]

# -- District × day means (areal average across stations) ----------------------
weather_means <- c("temperature_2m_max","temperature_2m_min","temperature_2m_mean",
                   "apparent_temperature_max","apparent_temperature_min","apparent_temperature_mean",
                   "shortwave_radiation_sum","precipitation_sum","rain_sum",
                   "precipitation_hours","windspeed_10m_max","windgusts_10m_max",
                   "et0_fao_evapotranspiration")

keep_cols <- c("district","time", weather_means)
wx_keep  <- wx_d[, intersect(names(wx_d), keep_cols), with = FALSE]
wx_daily <- wx_keep[, lapply(.SD, mean, na.rm = TRUE),
                    by = .(district, date = time),
                    .SDcols = setdiff(names(wx_keep), c("district","time"))]

###############################################################################
# 4) LAND COVER ??? Download .RAR ??? Extract ??? Raster ------------------------------
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
###############################################################################

# Ensure the .rar archive is extracted into a fixed subdir
rar_file <- file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs", "SriLanka_Landcover_2018.rar")
lc_dir   <- file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs", "SriLanka_Landcover_2018")

landcover_data_url = 'https://www.noda.ac.cn/en/knowledgehub/downloadAttachment?id=64257e8a29d1210627d613d3&type=DATASET'
download.file(landcover_data_url, rar_file)


if (!dir.exists(lc_dir) || !length(list.files(lc_dir, pattern = "\\.tif(f)?$", recursive = TRUE))) {
  extract_rar(rar_file, lc_dir)
}

# Point terra to the extracted .tif (adjust name if different inside the RAR)
tif_file <- list.files(lc_dir, pattern = "\\.tif(f)?$", full.names = TRUE, recursive = TRUE)[1]
if (is.na(tif_file)) stop("No .tif found after extraction in: ", lc_dir)

r <- rast(tif_file)

# Dissolve to district and match raster CRS
adm2_diss <- adm2 |>
  dplyr::select(province, geometry) |>
  dplyr::group_by(province) |>
  dplyr::summarise(geometry = st_union(geometry), .groups = "drop")|>
  st_transform(crs(r))

names(adm2_diss) <- c("district","geometry")

# Legend mapping
class_map <- data.table(
  code  = c(10,20,30,40,50,60,70,80,90,255),
  label = c("Water","BuiltUp","Cropland","Forest","Shrub","Grass","Bare","Wetland","Paddy","NoData")
)

# Weighted frequency ??? proportions within polygon
freq_fun <- function(values, coverage_fraction) {
  dt <- data.table(code = values, w = coverage_fraction)[!is.na(code)]
  if (!nrow(dt)) return(data.frame(code = integer(), prop = numeric()))
  dt <- dt[, .(w = sum(w, na.rm = TRUE)), by = code][, prop := w / sum(w)]
  as.data.frame(dt[, .(code, prop)])
}

lc_list <- exact_extract(r, adm2_diss,
                         fun = function(values, coverage_fraction) list(freq_fun(values, coverage_fraction)),
                         progress = TRUE)

lc_dt <- rbindlist(lc_list, idcol = "idx")
lc_dt[, district := adm2_diss$district[idx]][, idx := NULL]
lc_dt <- lc_dt[code %in% class_map$code]
lc_dt <- merge(lc_dt, class_map, by = "code", all.x = TRUE)
lc_wide <- dcast(lc_dt, district ~ label, value.var = "prop", fill = 0)

# Re-normalize after dropping NoData to ensure rows sum ~1
prop_cols <- setdiff(names(lc_wide), "district")
row_sums <- lc_wide[, rowSums(.SD), .SDcols = prop_cols]
lc_wide[, (prop_cols) := lapply(.SD, function(z) ifelse(row_sums > 0, z / row_sums, 0)), .SDcols = prop_cols]


################################################################################
# 5) ERA5 + MERGES + OUTPUT -----------------------------------------------------
################################################################################
# SECTION GOAL: Join district-daily ERA5 to weekly disease by aligning on
# date_mid; add weather & land cover; write final analysis panel.

# -- 5.1 ERA5 input -------------------------------------------------------------
era5 <- fread(paths$era5_daily)
stopifnot(all(c("date","district") %in% names(era5)))
era5[, district := norm_dist(district)]
era5[, date := as.IDate(date)]

# -- 5.2 Join weather (district × day) to weekly (by date_mid) -----------------
setnames(wx_daily, "date", "wx_date")
setkey(lep, district, date_mid)
setkey(wx_daily, district, wx_date)
lep <- wx_daily[lep, on = .(district, wx_date = date_mid)]  # left join onto lep

# -- 5.3 Join land cover + ERA5 (align ERA5 date to date_mid) ------------------
lep <- merge(lep, lc_wide, by = "district", all.x = TRUE)

# Recompute date_mid defensively (ensures correct type after merges)
lep[, date_mid := as.IDate(date_start + (as.integer(date_end - date_start) / 2))]

lep <- merge(lep, era5, by.x = c("district","date_mid"), by.y = c("district","date"), all.x = TRUE)

lep = lep[year(date_start) <= 2024] # climate data persists through 2024 as of now.

# -- 5.4 Persist final analysis dataset ----------------------------------------
fwrite(lep, file.path(paths$work_out, "lep_analysis_panel.csv"))
message("Saved analysis panel: ", normalizePath(file.path(paths$work_out, "lep_analysis_panel.csv"), winslash = "/"))


# End core workflow.
################################################################################



# ###############################################################################
# ###############################################################################

# Initial investigatory plotting - optional. 

# ###############################################################################
# # 5) PANELS & PLOTS ------------------------------------------------------------
# ###############################################################################
# 
# # -- Build district-week panel, rolling metrics, and national series -----------
# lep_week <- lep[, .(cases = sum(lepto, na.rm = TRUE)), by = .(district, date_end)]
# lep_week <- lep_week[!is.na(date_end)]
# weeks_seq <- seq(min(lep_week$date_end), max(lep_week$date_end), by = "7 days")
# panel <- CJ(district = sort(unique(lep_week$district)), date_end = weeks_seq, unique = TRUE)
# panel <- lep_week[panel, on = .(district, date_end)]
# panel[is.na(cases), cases := 0L]
# panel[, `:=`(iso_week = isoweek(date_end), iso_year = isoyear(date_end),
#              month = month(date_end), year = year(date_end))]
# setorder(panel, district, date_end)
# panel[, `:=`(
#   ma4   = frollmean(cases,  4, align = "right", na.rm = TRUE),
#   ma12  = frollmean(cases, 12, align = "right", na.rm = TRUE),
#   sum12 = frollsum(cases,  12, align = "right", na.rm = TRUE)
# ), by = district]
# 
# prev <- panel[, .(district, iso_week, iso_year, cases_prev = cases)]
# prev[, iso_year := iso_year + 1L]
# panel <- merge(panel, prev, by = c("district","iso_week","iso_year"), all.x = TRUE, sort = FALSE)
# panel[, `:=`(
#   yoy_abs = cases - cases_prev,
#   yoy_pct = ifelse(!is.na(cases_prev) & cases_prev > 0, 100 * (cases - cases_prev) / cases_prev, NA_real_)
# )]
# 
# nat <- panel[, .(cases = sum(cases)), by = date_end][order(date_end)]
# nat[, `:=`(ma4 = frollmean(cases, 4, align = "right"),
#            ma12 = frollmean(cases,12, align = "right"))]
# latest <- max(panel$date_end, na.rm = TRUE)
# 
# # -- Quick plot helpers --------------------------------------------------------
# fmt_int   <- function(x) format(as.integer(x), big.mark = ",")
# safe_ylim <- function(v, pad = 0.06) {
#   v <- v[is.finite(v)]
#   if (!length(v)) return(c(0,1))
#   r <- range(v, na.rm = TRUE); d <- diff(r); if (d == 0) d <- max(1, r[2])
#   c(max(0, r[1] - d*pad), r[2] + d*pad)
# }
# set_par <- function() {
#   par(family = "sans", cex.axis = 1.05, cex.lab = 1.2, cex.main = 1.25,
#       mar = c(4.2, 4.5, 3.0, 1.0), las = 1, xaxs = "i", yaxs = "i")
# }
# 
# # -- P1: National weekly trend with 4w and 12w MAs -----------------------------
# ragg::agg_png(file.path(paths$fig_dir, "P1_national_trend.png"),
#               width = 10, height = 5.5, units = "in", res = 450)
# on.exit(dev.off(), add = TRUE)
# set_par()
# yl <- safe_ylim(c(nat$cases, nat$ma12, nat$ma4))
# plot(nat$date_end, nat$cases, type = "l", lwd = 2.2, col = "#2C3E50",
#      xlab = "Week ending", ylab = "Cases",
#      main = "Sri Lanka leptospirosis - national weekly cases", ylim = yl)
# grid(col = "grey90"); box()
# lines(nat$date_end, nat$ma12, lwd = 2.8, col = "#E74C3C")
# lines(nat$date_end, nat$ma4,  lwd = 2.0, col = "#3498DB")
# legend("topleft", legend = c("Weekly cases","12-week MA","4-week MA"),
#        col = c("#2C3E50","#E74C3C","#3498DB"), lwd = c(2.2,2.8,2.0), bty = "n")
# dev.off()
# 
# # -- P2: Top 15 districts at latest week --------------------------------------
# top15 <- panel[date_end == latest][order(-cases)][1:15]
# ragg::agg_png(file.path(paths$fig_dir, "P2_top15_latest_week.png"),
#               width = 9, height = 7.5, units = "in", res = 450)
# set_par(); par(mar = c(4.2, 10, 3.0, 1.0))
# bp <- barplot(rev(top15$cases), horiz = TRUE, col = "#2E86C1", border = NA,
#               names.arg = rev(top15$district),
#               xlab = "Cases (weekly)", main = paste0("Top districts - week ending ", latest))
# grid(nx = NA, ny = NULL, col = "grey90"); box()
# text(x = rev(top15$cases), y = bp, labels = fmt_int(rev(top15$cases)),
#      pos = 4, xpd = NA, cex = 0.9, offset = 0.4)
# dev.off()
# 
# # -- P3: Small multiples (top 6 by last-52-week burden) -----------------------
# top6 <- panel[date_end > (latest - 7*52), .(cases_52 = sum(cases)), by = district][order(-cases_52)][1:6, district]
# ragg::agg_png(file.path(paths$fig_dir, "P3_top6_last52_facets.png"),
#               width = 12, height = 8, units = "in", res = 450)
# set_par(); par(mfrow = c(2,3), mar = c(3.8, 4.5, 2.8, 1.0))
# for (d in top6) {
#   dt <- panel[district == d & date_end > (latest - 7*52)]
#   yl <- safe_ylim(c(dt$cases, dt$ma12))
#   plot(dt$date_end, dt$cases, type = "l", lwd = 2.0, col = "#2C3E50",
#        xlab = "Week", ylab = "Cases", main = d, ylim = yl)
#   grid(col = "grey92"); box()
#   lines(dt$date_end, dt$ma12, lwd = 2.6, col = "#E74C3C")
# }
# dev.off()
# 
# # -- P4: Seasonality boxplots (ISO week) for top 8 total ----------------------
# top8 <- panel[, .(tot = sum(cases)), by = district][order(-tot)][1:8, district]
# ragg::agg_png(file.path(paths$fig_dir, "P4_seasonality_boxplots.png"),
#               width = 14, height = 8, units = "in", res = 450)
# set_par(); par(mfrow = c(2,4), mar = c(4.2, 4.5, 2.8, 0.8))
# for (d in top8) {
#   dt <- panel[district == d]
#   boxplot(cases ~ iso_week, data = dt, outline = FALSE,
#           xaxt = "n", col = "#AED6F1", border = "#2E86C1",
#           xlab = "ISO week", ylab = "Weekly cases", main = d)
#   axis(1, at = seq(1, 52, by = 4), labels = seq(1, 52, by = 4), tick = TRUE)
#   grid(nx = NA, ny = NULL, col = "grey90"); box()
# }
# dev.off()
# 
# # -- P5: Heatmap of rolling 12-week burden (last 3 years) ---------------------
# start3y <- as.Date(latest) - 365*3
# sub <- panel[date_end >= start3y, .(district, date_end, sum12)]
# ord <- sub[, .(burden = sum(sum12, na.rm = TRUE)), by = district][order(-burden), district]
# sub[, district := factor(district, levels = ord)]
# dates <- sort(unique(sub$date_end))
# dists <- levels(sub$district)
# zmat <- vapply(dists, function(d) {
#   df <- sub[district == d]
#   v  <- rep(NA_real_, length(dates))
#   v[match(df$date_end, dates)] <- df$sum12
#   v
# }, numeric(length(dates)))
# pal <- colorRampPalette(c("#FFFFFF","#FF9F66","#E74C3C","#A93226"))(256)
# 
# ragg::agg_png(file.path(paths$fig_dir, "P5_heatmap_sum12_last3y.png"),
#               width = 12, height = 8.5, units = "in", res = 450)
# set_par()
# x <- seq_along(dates); y <- seq_along(dists)
# image(x, y, zmat, col = pal, xlab = "Week ending", ylab = "", xaxt = "n", yaxt = "n", useRaster = TRUE)
# box()
# xticks <- unique(round(seq(1, length(dates), length.out = 9)))
# axis(1, at = xticks, labels = format(dates[xticks], "%Y-%m"))
# axis(2, at = y, labels = dists, las = 2, cex.axis = 0.9)
# title("Rolling 12-week burden by district (last 3 years)")
# dev.off()
# 
# # -- P6: YoY change bars (latest week) ----------------------------------------
# yoy <- panel[date_end == latest & !is.na(yoy_abs), .(district, yoy_abs)][order(-abs(yoy_abs))][1:15]
# cols <- ifelse(yoy$yoy_abs >= 0, "#E74C3C", "#2E86C1")
# ragg::agg_png(file.path(paths$fig_dir, "P6_yoy_change_top15.png"),
#               width = 9, height = 7.5, units = "in", res = 450)
# set_par(); par(mar = c(4.2, 10, 3.0, 1.0))
# bp <- barplot(rev(yoy$yoy_abs), horiz = TRUE,
#               names.arg = rev(yoy$district), col = rev(cols), border = NA,
#               xlab = "?? cases vs same ISO week last year",
#               main = paste0("YoY change - week ending ", latest))
# abline(v = 0, col = "grey40", lwd = 1.5)
# grid(nx = NA, ny = NULL, col = "grey90"); box()
# text(x = rev(yoy$yoy_abs), y = bp, labels = rev(fmt_int(yoy$yoy_abs)),
#      pos = ifelse(rev(yoy$yoy_abs) >= 0, 4, 2), xpd = NA, cex = 0.9, offset = 0.4)
# dev.off()
# 
# message("Saved: ", normalizePath(paths$fig_dir, winslash = "/"))
