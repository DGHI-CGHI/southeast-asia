################################################################################
# Sri Lanka | Climate-Health Modeling Sandbox (GAM-focused, ERA5-only)
# ------------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   . Loads a prepared, analysis-ready panel that already links:
#       - WER weekly counts at district × week (e.g., leptospirosis),
#       - population (poptot),
#       - district-level land cover proportions,
#       - ERA5-derived *weekly* climate features (and their lags/rolls/anoms/z).
#   . Discovers which ERA5 weekly columns are available by regex and builds
#     predictor blocks by *selecting from those columns only* (no recomputation).
#   . Prepares "per 100k" outcomes from counts when needed.
#   . Splits the series chronologically (80/20) to mimic deployment:
#       train on the past ??? evaluate on the future (prevents leakage).
#   . Provides a compact GAM harness:
#       - annual seasonality via cyclic spline s(week, bs = "cc"),
#       - optional district random effects s(district, bs = "re"),
#       - smooths only the ERA5 variables you request (if present).
#   . Prints minimal test metrics (RMSE/MAE/R²) and the exact variables used.
#
# WHY (AND WHAT CHANGED)
#   . You requested we *do not* use any station/observational weather fields
#     included earlier (e.g., temperature_2m_*, apparent_temperature_*, rain_sum,
#     precipitation_sum, etc.). This script *excludes* those by construction.
#   . "Single source of truth" for climate drivers is the ERA5 weekly feature
#     set you produced upstream. We only select from it; we do not re-derive
#     lags/rolling windows/anomalies/z-scores here.
#
# HOW TO RUN
#   1) Ensure `.Rprofile` defines cfg$paths (raw/intermediate/processed/reports).
#   2) Ensure your modeling panel CSV exists at:
#        cfg$paths$intermediate/sri_lanka-disease-landcover-climate-2015_2024.csv
#      and includes columns like:
#        district, date or date_mid, week_of_year (or derive from date),
#        poptot, ERA5 weekly fields (tmax_mean, rh_mean_week, precip_tp_sum_week,
#        their *_lag*, *_roll2w_*, *_roll4w_*, *_anom, *_pct_normal, *_z), and
#        optional LC proportions (e.g., Paddy, BuiltUp).
#   3) Source this script; then tweak the `smooth_vars` in §6 and re-run.
#
# DESIGN NOTES
#   . Chronological split avoids look-ahead (random splits would leak seasonality).
#   . Week 53 is coerced to 52 so the cyclic spline "closes" the year.
#   . Keep the smoother block small (3-6 variables) to avoid overfitting highly
#     correlated lags/rolls. Use your earlier screening tables to choose.
################################################################################

# ============================== SETUP =========================================
suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(stringr)
  library(mgcv)   # GAM/BAM
})

setDTthreads(percent = 100)

# ------------------------------------------------------------------------------
# 0) CONFIG - project-relative paths (no setwd; no absolute paths)
# ------------------------------------------------------------------------------
if (!exists("cfg") || is.null(cfg$paths)) {
  stop("`cfg$paths` not found. Define it in .Rprofile (raw/intermediate/processed/reports).")
}

paths <- list(
  helpers     = file.path("code/helpers/helpers.R"),  # optional
  panel_csv   = file.path(cfg$paths$intermediate, "sri_lanka-disease-landcover-climate-2015_2024.csv")
)

if (file.exists(paths$helpers)) source(paths$helpers)

# ------------------------------------------------------------------------------
# Utility helpers - metrics, type safety, factor alignment, cyclic week
# ------------------------------------------------------------------------------
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
mae  <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
r2   <- function(y, yhat) 1 - sum((y - yhat)^2, na.rm=TRUE) /
  sum((y - mean(y, na.rm=TRUE))^2, na.rm=TRUE)

harmonize_factors <- function(TR, TE, facs = c("district")) {
  for (f in facs) {
    TR[[f]] <- factor(TR[[f]])
    TE[[f]] <- factor(TE[[f]], levels = levels(TR[[f]]))
  }
  list(TR = TR, TE = TE)
}

prep_cyclic_week <- function(D) {
  if (!"week" %in% names(D)) stop("`week` column required (1..52); add before modeling.")
  D[, week := as.integer(pmax(1L, pmin(52L, week)))]
  D
}

ensure_numeric <- function(DT, cols) {
  ok <- intersect(cols, names(DT))
  for (c in ok) if (!is.numeric(DT[[c]])) suppressWarnings(DT[[c]] <- as.numeric(DT[[c]]))
  ok
}

# ------------------------------------------------------------------------------
# 1) Load modeling panel (WER + pop + LC + ERA5) and standardize time fields
# ------------------------------------------------------------------------------
if (!file.exists(paths$panel_csv)) {
  stop("Modeling panel not found: ", paths$panel_csv,
       "\nMake sure your ETL/post_processing produced it.")
}
DT <- fread(paths$panel_csv)

# Use `date_mid` if `date` is absent; set Date type
if (!"date" %in% names(DT) && "date_mid" %in% names(DT)) DT[, date := date_mid]
if (!inherits(DT$date, "Date")) DT[, date := as.IDate(date)]

# Prefer existing week column if present; else derive from date
if ("week_of_year" %in% names(DT)) {
  DT[, week := as.integer(pmax(1L, pmin(52L, week_of_year)))]
} else {
  DT[, week := as.integer(isoweek(date))]
  DT[week == 53L, week := 52L]
}
DT[, `:=`(year = year(date), yday = yday(date))]

# If *_A counts exist and poptot is present, create *_rate_100k outcomes
a_cols   <- grep("_A$", names(DT), value = TRUE)
diseases <- sub("_A$", "", a_cols)
if ("poptot" %in% names(DT) && length(a_cols)) {
  for (d in diseases) {
    a   <- paste0(d, "_A")
    rpk <- paste0(d, "_rate_100k")
    if (!rpk %in% names(DT)) DT[, (rpk) := 1e5 * get(a) / poptot]
  }
}

# Outcome resolution (robust to your column variants)
resolve_outcome <- function(D) {
  cand <- c("lepto_rate_100k","lepto_100k")
  for (nm in cand) if (nm %in% names(D)) return(nm)
  # fallbacks from counts
  if ("leptospirosis_A" %in% names(D) && "poptot" %in% names(D)) {
    D[, lepto_rate_100k := fifelse(poptot > 0, 1e5 * leptospirosis_A / poptot, NA_real_)]
    return("lepto_rate_100k")
  }
  if ("lepto" %in% names(D) && "poptot" %in% names(D)) {
    D[, lepto_rate_100k := fifelse(poptot > 0, 1e5 * lepto / poptot, NA_real_)]
    return("lepto_rate_100k")
  }
  stop("Could not resolve a leptospirosis rate outcome.")
}
y <- resolve_outcome(DT)

# ------------------------------------------------------------------------------
# 2) ERA5-only predictor discovery (NO station fields)
# ------------------------------------------------------------------------------
# WHAT:   Build a clean whitelist of ERA5 weekly feature names (and their
#         lag/roll/anom/z variants) from *existing columns only*.
# WHY:    You asked to drop station-series fields entirely.
# HOW:    Regex anchored on ERA5 base names, plus allowed suffix families.

# ERA5 weekly base names you've produced upstream
era5_base <- c(
  # temperature / humidity / moisture / radiation / heat stress
  "tmax_mean","tmin_mean","tmean_mean",
  "rh_mean_week","vpd_mean_week","ssrd_MJ_mean_week","wbgt_mean_week",
  # precipitation
  "precip_tp_sum_week","precip_mtpr_sum_week",
  # event counts / spells
  "wet_days_ge10_tp","wet_days_ge10_mtpr","max3d_tp","max3d_mtpr",
  "wet_spell_maxlen_tp","wet_spell_maxlen_mtpr","hot_days_ge32",
  # exposure-weighted precipitation proxies
  "ewap_tp","ewap_mtpr",
  # book-keeping
  "n_days_week"
)

# Allowed suffix families (lags/rolls/anoms/z)
era5_suffix <- "(?:_lag[1-6]|_roll(?:2w|4w)_(?:mean|sum)|_anom|_pct_normal|_z)?"
era5_regex  <- paste0("^(", paste(era5_base, collapse = "|"), ")", era5_suffix, "$")

era5_vars <- grep(era5_regex, names(DT), value = TRUE)

# (Optional) Land cover (linear terms only)
lc_vars <- intersect(c("BuiltUp","Cropland","Grass","Paddy","Shrub","Water","Wetland"), names(DT))

# Final candidate pool (ERA5 only + optional LC), pruned for NA and variance
PREDICTORS <- unique(c(era5_vars, lc_vars))

if (length(PREDICTORS)) {
  na_share <- DT[, lapply(.SD, function(z) mean(is.na(z))), .SDcols = PREDICTORS]
  keep <- names(na_share)[na_share < 0.30]  # keep vars with <30% missing
  if (length(keep)) {
    zv <- DT[, sapply(.SD, \(z) data.table(z)[, uniqueN(z, na.rm = TRUE)]), .SDcols = keep]
    keep <- keep[zv > 1]
  }
  PREDICTORS <- keep
}

# ------------------------------------------------------------------------------
# 3) Temporal split (chronological 80/20 by `date`)
# ------------------------------------------------------------------------------
# WHAT:  Split the panel into an "earlier" training subset and a "later"
#        evaluation subset at the 80th percentile of time.
# WHY:   This simulates real early-warning deployment (fit on past, score future)
#        and avoids look-ahead leakage that inflates performance.
# HOW:   Sort by `date`, cut at floor(0.8*N), and align factor levels so
#        predict() can handle district factors cleanly.

setorder(DT, date)
cut_date <- DT$date[floor(0.8 * nrow(DT))]
TR <- DT[date <= cut_date]
TE <- DT[date  >  cut_date]

hz <- harmonize_factors(TR, TE, facs = c("district"))
TR <- prep_cyclic_week(hz$TR)
TE <- prep_cyclic_week(hz$TE)

# ------------------------------------------------------------------------------
# 4) GAM harness (cyclic week, optional district RE, ERA5 smooths only)
# ------------------------------------------------------------------------------
# WHAT:  mgcv::gam with:
#         . s(week, bs='cc') to encode annual cycle,
#         . optional s(district, bs='re') when >1 district,
#         . smooths limited to ERA5 variables you pass in (if present),
#         . optional linear LC terms.
# WHY:   Parsimonious + interpretable. Avoids dumping dozens of collinear lags.
# HOW:   Build formula from intersected columns; use k=6 cubic regression splines.

fit_predict_gam_safe <- function(TR, TE, outcome,
                                 smooth_vars,
                                 linear_vars = character(0),
                                 use_bam = FALSE,
                                 verbose = TRUE) {
  # Intersect with columns that actually exist in TR
  smooth_vars <- intersect(smooth_vars, names(TR))
  linear_vars <- intersect(linear_vars, names(TR))
  
  # Keep only what we need on each side
  keep <- unique(c(outcome, "district", "week", smooth_vars, linear_vars))
  TR2 <- copy(TR)[, ..keep]
  TE2 <- copy(TE)[, ..keep]
  
  ensure_numeric(TR2, setdiff(names(TR2), c("district","week", outcome)))
  ensure_numeric(TE2, setdiff(names(TE2), c("district","week", outcome)))
  
  have_RE <- "district" %in% names(TR2) && nlevels(factor(TR2$district)) >= 2
  re_term <- if (have_RE) " + s(district, bs='re')" else ""
  
  s_terms <- if (length(smooth_vars)) sprintf("s(%s, bs='cr', k=6)", smooth_vars) else character(0)
  s_part  <- paste(s_terms, collapse = " + ")
  lin_part<- if (length(linear_vars)) paste(linear_vars, collapse = " + ") else NULL
  
  rhs <- paste(c("s(week, bs='cc', k=20)", s_part, lin_part), collapse = " + ")
  rhs <- paste0(rhs, re_term)
  
  fml <- as.formula(paste0(outcome, " ~ ", rhs))
  if (verbose) message("GAM formula: ", deparse(fml))
  
  fit_fun <- if (use_bam) mgcv::bam else mgcv::gam
  mod <- fit_fun(fml, data = TR2, method = "REML",
                 knots = list(week = c(0.5, 52.5)))
  
  pred <- as.numeric(predict(mod, newdata = TE2))
  list(model = mod, pred = pred,
       used = list(smooth_vars = smooth_vars, linear_vars = linear_vars))
}

# ------------------------------------------------------------------------------
# 5) Example run - ERA5-only compact smoother set
# ------------------------------------------------------------------------------
# Tip: keep 3-6 smooths and pick from your screening table. Below is a reasonable
# starting set: temperature + moisture + precip + one lag + heat/radiation proxy.

smooth_pref <- c(
  "tmax_mean", "tmax_mean_lag2",
  "rh_mean_week",
  "precip_tp_sum_week", "precip_tp_sum_week_lag2",
  "wbgt_mean_week",
  "ssrd_MJ_mean_week",
  "vpd_mean_week",
  "ewap_tp"
)

smooth_vars <- intersect(smooth_pref, PREDICTORS)  # ERA5-only
linear_vars <- intersect(lc_vars,      PREDICTORS)  # optional land cover (linear)

y <- "lepto_100k"
res <- fit_predict_gam_safe(
  TR, TE, outcome = y,
  smooth_vars = smooth_vars,
  linear_vars = linear_vars
)

cat("\nERA5-only variables USED:\n  smooth: ",
    if (length(res$used$smooth_vars)) paste(res$used$smooth_vars, collapse = ", ") else "(none)",
    "\n  linear: ",
    if (length(res$used$linear_vars)) paste(res$used$linear_vars, collapse = ", ") else "(none)",
    "\n", sep = "")

cat("\nHold-out metrics for ", y, "\n",
    "  RMSE: ", round(rmse(TE[[y]], res$pred), 3),
    "  MAE: ",  round(mae (TE[[y]], res$pred), 3),
    "  R²: ",   round(r2  (TE[[y]], res$pred), 3), "\n", sep = "")
