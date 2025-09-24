################################################################################
# Sri Lanka | Climate-Health Modeling Sandbox (GAM-focused)
# ------------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   . Loads a prepared, analysis-ready panel that already links:
#       - WER weekly counts (e.g., leptospirosis) at district × week,
#       - population (poptot),
#       - district-level land cover proportions,
#       - ERA5-derived *weekly* climate features (incl. lags, rolling windows,
#         anomalies, and z-scores) that were produced upstream.
#   . Discovers which climate columns are available and defines thematic
#     predictor blocks (temperature, moisture, radiation, heat, etc.) by
#     *selecting from existing columns* - no recomputation of lags/rolls here.
#   . Prepares outcomes in "per 100k" units if you loaded counts with poptot.
#   . Splits the time series *chronologically* (80/20) to simulate deployment:
#     training on the "past", evaluating on the "future" (prevents leakage).
#   . Provides a compact, safe GAM harness (`fit_predict_gam_safe`) that:
#       - encodes *annual seasonality* with a cyclic spline s(week, bs='cc'),
#       - adds *district random effects* when >1 district is present,
#       - *smooths only the climate variables you request* that are present,
#       - leaves room to include linear terms if you want them.
#   . Outputs simple regression metrics (RMSE/MAE/R²) for quick model checks.
#
# WHY THIS EXISTS
#   . To give you a reproducible "playground" where you can quickly try
#     combinations of climate covariates (picked from your upstream features)
#     without re-deriving features or risking data leakage.
#   . To keep model code minimal, explicit, and aligned with how an *early
#     warning* system would actually be trained (past) and used (future).
#
# HOW TO RUN
#   1) Ensure `.Rprofile` defines `cfg$paths` used below (raw/intermediate/
#      processed/reports). This avoids setwd() and hard-coded paths.
#   2) Ensure your modeling panel (referenced below) contains at minimum:
#         district (chr), date or date_mid (Date/IDate), week (1..52),
#         poptot (num), climate features (weekly ERA5 with lags/rolls/anoms/z),
#         and optional land cover (e.g., Paddy, BuiltUp).
#   3) Source this script in an interactive session. After running, explore:
#         - objects: DT/TR/TE, predictor blocks, and `fit_predict_gam_safe()`.
#         - try alternative `smooth_vars` sets and re-fit.
#
# KEY DESIGN CHOICES & CAVEATS
#   . *Single source of truth for features*: weekly ERA5 features (lags, rolls,
#     anomalies, z-scores) are computed upstream; this script only *selects*
#     them. That keeps provenance clean and consistent with your prior work.
#   . *Chronological split*: avoids look-ahead (e.g., random splits contaminate
#     "future" with "past"). This matters for early warning realism.
#   . *Cyclic seasonality*: `s(week, bs='cc')` closes the yearly cycle. Week 53
#     is coerced to 52 to maintain the cycle.
#   . *District effects*: add `s(district, bs='re')` when appropriate; this
#     absorbs static between-district offsets so climate smooths model the
#     remaining temporal signal.
#   . *Parsimony*: use a *small* set of smooths (e.g., 3-6) guided by your
#     screening (??PR-AUC / ??NLL), to avoid overfitting highly correlated lags.
#
# OUTPUTS
#   . No files are written by default. The example shows in-memory RMSE/MAE/R².
#     If you want artifacts, add `fwrite()` or `saveRDS()` at the bottom.
################################################################################

# ============================== SETUP =========================================
suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(stringr)
  
  # Model families; handy to swap among in this sandbox
  library(mgcv)         # GAM/BAM with smooths and random effects
  library(fixest)       # feols/fepois (sometimes useful)
  library(glmmTMB)      # NB GLMMs (optional)
  library(ranger)       # Random Forest (optional)
  library(glmnet)       # LASSO/Elastic Net (optional)
  library(car)          # VIF (optional diagnostics)
})

setDTthreads(percent = 100)

# ------------------------------------------------------------------------------
# 0) CONFIG - project-relative paths (no setwd; no absolute local paths)
# ------------------------------------------------------------------------------
# WHAT: Define where inputs live and where optional outputs go using `cfg$paths`.
# WHY:  Keeps code portable across machines/EC2/CI; zero hard-coded drive roots.
# HOW:  `cfg` is expected in .Rprofile; guard with informative warnings.

if (!exists("cfg") || is.null(cfg$paths)) {
  stop("`cfg$paths` not found. Define it in .Rprofile (raw/intermediate/processed/reports).")
}

paths <- list(
  temp_dir  = file.path(cfg$paths$intermediate, "temp"),
  work_out  = file.path(cfg$paths$intermediate, "outputs"),
  helpers   = file.path("code/helpers/helpers.R"),
  fig_dir   = file.path(cfg$paths$reports, "figures"),
  era5_daily = file.path(cfg$paths$raw, "srilanka_district_daily_era5_areawt.csv"),
  # Header-only discovery of weekly ERA5 table (optional, see §1)
  era5_weekly = file.path(cfg$paths$processed, "srilanka_district_weekly_era5_areawt.csv"),
  # Main modeling panel (WER + pop + LC + ERA5) produced in your ETL
  panel_csv  = file.path(cfg$paths$intermediate, "sri_lanka-disease-landcover-climate-2015_2024.csv")
)

ensure_dir <- function(...) dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE)
invisible(list(
  ensure_dir(cfg$paths$intermediate, "temp"),
  ensure_dir(cfg$paths$intermediate, "outputs"),
  ensure_dir(cfg$paths$reports, "figures")
))

# Helper source is harmless here (kept for completeness if others reuse script)
if (file.exists(paths$helpers)) source(paths$helpers)

# ------------------------------------------------------------------------------
# Utility helpers - metrics, type safety, factor alignment, cyclic week
# ------------------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
mae  <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
r2   <- function(y, yhat) 1 - sum((y - yhat)^2, na.rm=TRUE) /
  sum((y - mean(y, na.rm=TRUE))^2, na.rm=TRUE)

# Add back dropped regressors for fixest predict(); handy if you use FE models.
safe_newdata_fixest <- function(mod, TE, fe_keys = c("district","year")) {
  used_beta <- names(coef(mod))
  dropped   <- mod$collin.var %||% character(0)
  toks      <- all.vars(mod$fml_all$linear)
  keep      <- unique(c(fe_keys, used_beta, dropped, toks))
  miss <- setdiff(keep, names(TE))
  if (length(miss)) for (m in miss) TE[[m]] <- NA_real_
  TE[, ..keep]
}

# Keep factor levels aligned across splits (critical for predict()).
harmonize_factors <- function(TR, TE, facs = c("district")) {
  for (f in facs) {
    TR[[f]] <- factor(TR[[f]])
    TE[[f]] <- factor(TE[[f]], levels = levels(TR[[f]]))
  }
  list(TR = TR, TE = TE)
}

# Force week into [1..52]; coerce ISO week 53 to 52 so cyclic spline "closes".
prep_cyclic_week <- function(D) {
  if (!"week" %in% names(D)) stop("`week` column required (1..52); add before modeling.")
  D[, week := as.integer(pmax(1L, pmin(52L, week)))]
  D
}

# Ensure listed columns are numeric; safe coercion for mixed types.
ensure_numeric <- function(DT, cols) {
  ok <- intersect(cols, names(DT))
  for (c in ok) if (!is.numeric(DT[[c]])) suppressWarnings(DT[[c]] <- as.numeric(DT[[c]]))
  ok
}

# If using all LC proportions in linear models, drop one to avoid simplex collinearity.
drop_one_landcover <- function(PV, drop = c("NoData","Water","BuiltUp")) {
  lc <- c("BuiltUp","Cropland","Grass","Paddy","Shrub","Water","Wetland","NoData")
  present <- intersect(lc, PV)
  if (length(present) <= 1) return(PV)
  d <- intersect(drop, present)
  if (length(d)) return(setdiff(PV, d[1L]))
  setdiff(PV, present[1L])
}

# ------------------------------------------------------------------------------
# 1) Discover available ERA5 weekly feature names (optional, header only)
# ------------------------------------------------------------------------------
# WHAT:  Reads only the column names from the weekly ERA5 table, so you can
#        confirm the schema and align regex pickers later.
# WHY:   This is a cheap sanity check; the model section relies on columns
#        actually present in the main panel you load next.
# HOW:   `fread(..., nrows=1)` to fetch header without loading full data.

if (file.exists(paths$era5_weekly)) {
  era5_var_names <- names(fread(paths$era5_weekly, nrows = 1L))
} else {
  warning("ERA5 weekly feature file not found: ", paths$era5_weekly)
  era5_var_names <- character(0)
}

# ------------------------------------------------------------------------------
# 2) Load modeling panel (WER + pop + LC + ERA5) and standardize time fields
# ------------------------------------------------------------------------------
# WHAT:  Ingest the unified panel produced by your ETL. Ensure it has
#        district, date/date_mid, week (or enough info to derive it), poptot,
#        and the ERA5 weekly features you want to use.
# WHY:   Keeps modeling isolated from scraping/ETL logic; this script assumes
#        the panel is "ready to model".
# HOW:   Read CSV, set Date, and compute week/yday if needed.

if (!file.exists(paths$panel_csv)) {
  stop("Modeling panel not found: ", paths$panel_csv,
       "\nMake sure your ETL/post_processing produced it.")
}
DT <- fread(paths$panel_csv)

# Alias date_mid ??? date if that's your canonical weekly stamp
if (!"date" %in% names(DT) && "date_mid" %in% names(DT)) DT[, date := date_mid]
if (!inherits(DT$date, "Date")) DT[, date := as.IDate(date)]

# Derive time parts if missing
if (!"week" %in% names(DT)) DT[, week := as.integer(isoweek(date))]
DT[week == 53L, week := 52L]
DT[, `:=`(year = year(date), yday = yday(date))]

# If counts are present, generate rate-per-100k outcomes (keeps your naming)
a_cols   <- grep("_A$", names(DT), value = TRUE)
diseases <- sub("_A$", "", a_cols)
if ("poptot" %in% names(DT) && length(a_cols)) {
  for (d in diseases) {
    a   <- paste0(d, "_A")
    rpk <- paste0(d, "_rate_100k")
    if (!rpk %in% names(DT)) DT[, (rpk) := 1e5 * get(a) / poptot]
  }
}
outcome_cols <- grep("_rate_100k$", names(DT), value = TRUE)
if (!length(outcome_cols)) outcome_cols <- intersect(names(DT), c("lepto_100k","dengue_100k"))
stopifnot(length(outcome_cols) > 0)

# ------------------------------------------------------------------------------
# 3) Select predictor families from *existing* weekly features (no recompute)
# ------------------------------------------------------------------------------
# WHAT:  Define climate blocks by regex directly from columns you already made
#        upstream (lags/rolls/anoms/z). No new lags/rolls here.
# WHY:   Prevents feature drift/duplication and keeps provenance clear.
# HOW:   Narrow regexes pick variables, then we prune for NA share and variance.

has <- function(x) x[x %in% names(DT)]  # quick helper

# 3.1 Base weekly summaries (non-lagged)
base_weekly <- has(c(
  "tmax_mean","tmin_mean","tmean_mean","wbgt_mean_week",
  "rh_mean_week","vpd_mean_week","ssrd_MJ_mean_week",
  "precip_tp_sum_week","precip_mtpr_sum_week",
  "wet_days_ge10_tp","wet_days_ge10_mtpr",
  "max3d_tp","max3d_mtpr","wet_spell_maxlen_tp","wet_spell_maxlen_mtpr",
  "hot_days_ge32","n_days_week"
))






# 3.2 Fixed lags (1..6)
lag_vars  <- grep("_(lag[1-6])$", names(DT), value = TRUE)

# 3.3 Rolling windows (2w/4w, mean/sum)
roll_vars <- grep("_roll(2w|4w)_(mean|sum)$", names(DT), value = TRUE)

# 3.4 Anomalies and % of normal (climatology-based)
anom_vars <- grep("(_anom|_pct_normal)$", names(DT), value = TRUE)

# 3.5 Z-scores (within-district standardization)
z_vars    <- grep("_z$", names(DT), value = TRUE)

# 3.6 Land cover (proportions)
lc_vars <- intersect(c("BuiltUp","Cropland","Grass","Paddy","Shrub","Water","Wetland"),
                     names(DT))

# Optional mechanistic interactions (only if both parents exist)
DT[, `:=`(
  precipX_paddy = if ("precip_tp_sum_week" %in% names(DT) && "Paddy" %in% names(DT))
    precip_tp_sum_week * Paddy else NA_real_,
  tmaxX_built   = if ("tmax_mean" %in% names(DT) && "BuiltUp" %in% names(DT))
    tmax_mean * BuiltUp else NA_real_
)]
lc_int_vars <- c("precipX_paddy","tmaxX_built")[c("precipX_paddy","tmaxX_built") %in% names(DT)]

# Kitchen-sink candidate list (we'll prune before modeling)
PREDICTORS <- unique(c(base_weekly, lag_vars, roll_vars, anom_vars, z_vars,
                       lc_vars, lc_int_vars))

# Prune high-NA and zero-variance columns (keeps script robust to partial panels)
if (length(PREDICTORS)) {
  na_share <- DT[, lapply(.SD, function(z) mean(is.na(z))), .SDcols = PREDICTORS]
  keep <- names(na_share)[na_share < 0.30]  # allow up to 30% missing
  if (length(keep)) {
    zv <- DT[, sapply(.SD, \(z) data.table(z)[, uniqueN(z, na.rm = TRUE)]), .SDcols = keep]
    keep <- keep[zv > 1]
  }
  PREDICTORS <- keep
}
P_SMALL <- PREDICTORS  # keep name used elsewhere


# ------------------------------------------------------------------------------
# 4) Temporal split (chronological 80/20 by `date`)
# ------------------------------------------------------------------------------
# WHAT:  Split the series into an "earlier" training subset and a "later"
#        evaluation subset using a *time-ordered* cut at the 80th percentile.
# WHY:   Mimics *real deployment*: you will fit on the past and score the
#        future. Random splits leak future information (seasonal structure,
#        trend) into training and inflate performance estimates.
# HOW:   Sort by `date`, take a cut index at floor(0.8*N), and split. Then
#        harmonize factor levels (e.g., `district`) so predict() won't fail.

setorder(DT, date)
cut_date <- DT$date[floor(0.8 * nrow(DT))]
TR <- DT[date <= cut_date]
TE <- DT[date  >  cut_date]

hz <- harmonize_factors(TR, TE, facs = c("district"))
TR <- prep_cyclic_week(hz$TR)
TE <- prep_cyclic_week(hz$TE)

TE$lepto_100k
# ------------------------------------------------------------------------------
# 5) GAM harness (cyclic week, optional district RE, slim climate smooths)
# ------------------------------------------------------------------------------
# WHAT:  Builds a mgcv GAM that (1) captures annual seasonality with a cyclic
#        spline, (2) optionally absorbs district offsets via a random effect,
#        and (3) *smooths only* the climate variables you pass in (that exist).
# WHY:   Keeps the model parsimonious and interpretable; avoids overfitting
#        dozens of correlated lags/rolls; aligns with your prior screening.
# HOW:   Construct formula from available columns; `k=6` cubic regression
#        splines for each climate smoother; `REML` for stable smoothing.

fit_predict_gam_safe <- function(TR, TE, outcome,
                                 smooth_vars = c(
                                   # sensible compact defaults; change per outcome
                                   "tmax_mean", "tmax_mean_lag2",
                                   "rh_mean_week", "vpd_mean_week",
                                   "precip_tp_sum_week", "precip_tp_sum_week_lag2",
                                   "ssrd_MJ_mean_week", "wbgt_mean_week",
                                   "ewap_tp"
                                 ),
                                 linear_vars  = character(0),
                                 use_bam = FALSE, verbose = TRUE) {
  
  # Keep only columns that exist in TR
  smooth_vars <- smooth_vars[smooth_vars %in% names(TR)]
  linear_vars <- linear_vars[linear_vars %in% names(TR)]
  keep <- unique(c(outcome, "district", "week", smooth_vars, linear_vars))
  TR2 <- copy(TR)[, ..keep]
  TE2 <- copy(TE)[, ..keep]
  
  # Defensive numeric coercion for non-factor columns
  ensure_numeric(TR2, setdiff(names(TR2), c("district","week", outcome)))
  ensure_numeric(TE2, setdiff(names(TE2), c("district","week", outcome)))
  
  # Random effect for district if multiple levels
  have_RE <- "district" %in% names(TR2) && nlevels(factor(TR2$district)) >= 2
  re_term <- if (have_RE) " + s(district, bs='re')" else ""
  
  # Build smooth and (optional) linear parts
  s_terms <- if (length(smooth_vars)) sprintf("s(%s, bs='cr', k=6)", smooth_vars) else character(0)
  s_part  <- paste(s_terms, collapse = " + ")
  lin_part <- if (length(linear_vars)) paste(linear_vars, collapse = " + ") else NULL
  
  # Cyclic annual seasonality (explicit cc knots close the loop)
  week_term <- "s(week, bs='cc', k=20)"
  
  rhs <- paste(c(week_term, s_part, lin_part), collapse = " + ")
  rhs <- paste0(rhs, re_term)
  fml <- as.formula(paste0(outcome, " ~ ", rhs))
  if (verbose) message("GAM formula: ", deparse(fml))
  
  fit_fun <- if (use_bam) mgcv::bam else mgcv::gam
  mod <- fit_fun(
    fml, data = TR2,
    method = "REML",
    knots = list(week = c(0.5, 52.5))
  )
  
  pred <- as.numeric(predict(mod, newdata = TE2))
  list(model = mod, pred = pred)
}

# ------------------------------------------------------------------------------
# 6) Example: fit a compact climate GAM and print quick metrics
# ------------------------------------------------------------------------------
# WHAT:  Minimal end-to-end sanity check: pick an outcome, choose 3-6 smooths,
#        fit on TR, score on TE, and print RMSE/MAE/R². Adjust `smooth_vars`
#        using your screening table (??VAL PR-AUC / ??NLL) for each outcome.
# WHY:   Quick feedback loop without touching upstream ETL or recomputing lags.
# HOW:   Uncomment and run.

y <- "lepto_rate_100k"
res <- fit_predict_gam_safe(
  TR, TE, outcome = y,
  smooth_vars = c("tmax_mean_lag2",
                  "rh_mean_week",
                  "precip_tp_sum_week",
                  "precip_tp_sum_week_lag2",
                  "ewap_tp")
)
cat("Hold-out metrics for", y, "\n",
    "  RMSE:", round(rmse(TE[[y]], res$pred), 3),
    " MAE:",  round(mae (TE[[y]], res$pred), 3),
    "  R²:",  round(r2  (TE[[y]], res$pred), 3), "\n")

# OPTIONAL: If you want district-adjusted partial dependence or smooth summaries:
# summary(res$model)
# plot(res$model, pages = 1, shade = TRUE)

################################################################################
# Notes on redundancy cleanup
# ------------------------------------------------------------------------------
# . We removed all code that recomputed lags/rolls - your ERA5 weekly features
#   already include: *_lag1..6, *_roll2w_*, *_roll4w_*, *_anom, *_pct_normal,
#   and *_z. This script solely *selects* from them.
# . We kept helper functions only once (no duplicate `harmonize_factors`,
#   `ensure_numeric`, etc.). If you paste into an existing file, delete any
#   duplicate definitions to avoid masking warnings and confusion.
# . If you later want EWS classification scoring (PR-AUC vs. MAD/ABS/PCT
#   baselines), tie this modeling harness into your existing alert evaluation
#   functions; the chronological split and seasonality handling here are
#   already aligned with that goal.
################################################################################






















































































# 
# for (y in outcome_cols) {
#   message("=== OUTCOME: ", y, " ===")
#   
#   # ---------- (A) FE-OLS on rate --------------------------------------------
#   f_rhs <- paste(P_SMALL, collapse = " + ")
#   f_ols <- as.formula(paste0(y, " ~ ", f_rhs, " | district + year"))
#   # mod_ols <- feols(f_ols, data = TR, vcov = ~district)
#   
#   mod_ols  <- feols(f_ols, data = TR, vcov = ~district)
#   new_ols  <- safe_newdata_fixest(mod_ols, TE)   # <-- use patched helper
#   yhat_ols <- as.numeric(predict(mod_ols, newdata = new_ols))
#   
#   
#   # new_ols <- safe_newdata_fixest(mod_ols, TE)
#   # yhat_ols <- as.numeric(predict(mod_ols, newdata = new_ols))
#   perf_ols <- data.table(model = "FE-OLS", outcome = y,
#                          RMSE = rmse(TE[[y]], yhat_ols),
#                          MAE  = mae(TE[[y]],  yhat_ols),
#                          R2   = r2(TE[[y]],   yhat_ols))
#   
#   # ---------- (B) FE-Poisson (counts + offset) ------------------------------
#   cnt <- sub("_rate_100k$", "_A", y)
#   perf_pois <- perf_nb <- NULL
#   
#   
#   if (cnt %in% names(TR)) {
#     # Put the offset inside the formula (critical!)
#     f_pois <- as.formula(paste0(
#       cnt, " ~ ", f_rhs, " + offset(log(poptot)) | district + year"
#     ))
#     
#     mod_pois <- fepois(f_pois, data = TR, vcov = ~district)
#     
#     # Build newdata that includes FE keys, used slopes, dropped vars, and poptot
#     new_pois <- safe_newdata_fixest(mod_pois, TE)
#     # predict() will evaluate offset(log(poptot)) on new_pois
#     mu_pois  <- as.numeric(predict(mod_pois, newdata = new_pois, type = "response"))
#     
#     perf_pois <- data.table(model = "FE-Poisson", outcome = y,
#                             RMSE = rmse(TE[[cnt]], mu_pois),
#                             MAE  = mae(TE[[cnt]],  mu_pois),
#                             R2   = r2(TE[[cnt]],   mu_pois))
#     
#     # ---------- (C) NegBin (glmmTMB FE via factors) ----------
#     f_nb <- as.formula(paste0(
#       cnt, " ~ ", f_rhs,
#       " + factor(district) + factor(year) + offset(log(poptot))"
#     ))
#     mod_nb <- try(glmmTMB(family = nbinom2(), data = TR, formula = f_nb), silent = TRUE)
#     if (!inherits(mod_nb, "try-error")) {
#       # TE must have poptot too
#       stopifnot("poptot" %in% names(TE))
#       mu_nb <- as.numeric(predict(mod_nb, newdata = TE, type = "response"))
#       perf_nb <- data.table(model = "NB (glmmTMB FE)", outcome = y,
#                             RMSE = rmse(TE[[cnt]], mu_nb),
#                             MAE  = mae(TE[[cnt]],  mu_nb),
#                             R2   = r2(TE[[cnt]],   mu_nb))
#     }
#   }
#   
#   
#   # ---------- (D) GAM (rate) ------------------------------------------------
#   # smooths for a few key ERA5 vars + cyclic week + RE for district
#   f_gam <- as.formula(
#     paste0(y, " ~ s(temperature_2m_mean) + s(precipitation_sum) + ",
#            "s(temperature_2m_mean_lag1) + s(precipitation_sum_lag1) + ",
#            "s(week, bs='cc', k=20) + s(district, bs='re') + ",
#            paste(setdiff(P_SMALL, c("temperature_2m_mean","precipitation_sum",
#                                     "temperature_2m_mean_lag1","precipitation_sum_lag1",
#                                     "sin52","cos52")), collapse = " + "))
#   )
#   
#   # TR$week
#   
#   gp <- fit_predict_gam_safe(TR, TE, outcome = y, P_SMALL = P_SMALL, use_bam = FALSE)
#   mod_gam  <- gp$model
#   yhat_gam <- gp$pred
#   
#   
#   perf_gam <- data.table(model = "GAM (rate)", outcome = y,
#                          RMSE = rmse(TE[[y]], yhat_gam),
#                          MAE  = mae(TE[[y]],  yhat_gam),
#                          R2   = r2(TE[[y]],   yhat_gam))
#   
#   
#   
#   
#   # ---------- (E) Random Forest (rate) --------------------------------------
#   Xcols <- unique(c("district","year", P_SMALL))
#   X_tr  <- TR[, ..Xcols]; y_tr <- TR[[y]]
#   X_te  <- TE[, ..Xcols]; y_te <- TE[[y]]
#   for (cc in c("district","year")) {
#     if (cc %in% names(X_tr)) { X_tr[[cc]] <- as.factor(X_tr[[cc]]); X_te[[cc]] <- as.factor(X_te[[cc]]) }
#   }
#   X_tr$y_resp <- y_tr
#   mod_rf <- ranger(
#     y_resp ~ ., data = X_tr,
#     num.trees = 800, mtry = max(3, floor(sqrt(ncol(X_tr)))),
#     min.node.size = 10, importance = "impurity", seed = 42
#   )
#   yhat_rf <- predict(mod_rf, data = X_te)$predictions
#   perf_rf <- data.table(model = "RandomForest (rate)", outcome = y,
#                         RMSE = rmse(y_te, yhat_rf),
#                         MAE  = mae(y_te,  yhat_rf),
#                         R2   = r2(y_te,   yhat_rf))
#   
#   # ---------- (F) LASSO / Elastic-Net (rate) --------------------------------
#   # Build a numeric model matrix (drop district/year here; add as dummies if needed)
#   mm_form <- as.formula(paste0("~ 0 + ", paste(P_SMALL, collapse = " + ")))
#   Xmm_tr  <- model.matrix(mm_form, data = TR)
#   Xmm_te  <- model.matrix(mm_form, data = TE)
#   y_tr_mm <- TR[[y]]; y_te_mm <- TE[[y]]
#   
#   cvfit <- cv.glmnet(Xmm_tr, y_tr_mm, alpha = 1, family = "gaussian", nfolds = 5)  # LASSO
#   yhat_lasso <- as.numeric(predict(cvfit, newx = Xmm_te, s = "lambda.min"))
#   perf_lasso <- data.table(model = "LASSO (rate)", outcome = y,
#                            RMSE = rmse(y_te_mm, yhat_lasso),
#                            MAE  = mae(y_te_mm,  yhat_lasso),
#                            R2   = r2(y_te_mm,   yhat_lasso))
#   
#   quiet_fixest <- function(expr) {
#     # Suppress fixest's collinearity chatter within a block
#     suppressMessages(suppressWarnings(capture.output(res <- force(expr))))
#     invisible(res); res
#   }
#   
#   # ---------- Rolling-Origin CV (brief) -------------------------------------
#   # cv_gauss <- rolling_origin_eval(DT, y, P_SMALL, k = 5, family = "gaussian")
#   # if (cnt %in% names(DT)) {
#   #   cv_pois  <- rolling_origin_eval(DT, y, P_SMALL, k = 5, family = "poisson")
#   #   cv_negb  <- rolling_origin_eval(DT, y, P_SMALL, k = 5, family = "negbin")
#   # } else {
#   #   cv_pois <- cv_negb <- data.table()
#   # }
#   
#   # Map outcome (rate) -> count column
#   resolve_count_name <- function(outcome, count_name = NULL, count_map = NULL) {
#     # 1) explicit
#     if (!is.null(count_name)) return(count_name)
#     # 2) named list map
#     if (!is.null(count_map) && outcome %in% names(count_map)) return(count_map[[outcome]])
#     # 3) common heuristics
#     #    *_rate_100k -> *_A
#     if (grepl("_rate_100k$", outcome)) return(sub("_rate_100k$", "_A", outcome))
#     #    dengue_100k -> dengue_A
#     if (grepl("_100k$", outcome))       return(sub("_100k$", "_A", outcome))
#     #    lepto_100k -> leptospirosis_A (special case)
#     if (outcome == "lepto_100k")        return("leptospirosis_A")
#     stop("Can't infer count column for outcome '", outcome,
#          "'. Provide count_name= or count_map=.")
#   }
#   
#   
#   # You have these columns per your names():
#   #   rates: dengue_100k, lepto_100k
#   #   counts: dengue_A,   leptospirosis_A
#   count_map <- list(
#     dengue_100k = "dengue_A",
#     lepto_100k  = "leptospirosis_A"
#   )
#   
#   # Example calls:
#   cv_gauss <- rolling_origin_eval(DT, "dengue_100k", P_SMALL, k = 5, family = "gaussian")
#   
#   cv_pois  <- rolling_origin_eval(DT, "dengue_100k", P_SMALL, k = 5,
#                                   family = "poisson",  count_map = count_map)
#   
#   # cv_negb  <- rolling_origin_eval(DT, "lepto_100k",  P_SMALL, k = 5,
#   #                                 family = "negbin",  count_map = count_map)
#   # 
#   # ---------- Collect --------------------------------------------------------
#   perf_all <- rbindlist(list(perf_ols, perf_pois, perf_nb, perf_gam, perf_rf, perf_lasso), fill = TRUE)
#   results_all[[y]] <- list(
#     perf_test = perf_all,
#     cv_gauss  = cv_gauss,
#     cv_pois   = cv_pois,
#     # cv_negb   = cv_negb,
#     mod_ols   = mod_ols,
#     mod_gam   = mod_gam,
#     mod_rf    = mod_rf,
#     lasso     = cvfit
#   )
#   
#   cat(y, '\n')
# }
# 
# # Test performance across outcomes
# perf_summary <- rbindlist(lapply(results_all, `[[`, "perf_test"), use.names = TRUE, fill = TRUE)
# print(perf_summary[order(outcome, -R2)])
# 
# # ====================== OPTIONAL DIAGNOSTICS ==================================
# # VIF (simple OLS on TR for one or two outcomes)
# for (y in head(outcome_cols, 2L)) {
#   cols_needed <- unique(c(y, P_SMALL))
#   Z <- TR[, ..cols_needed]
#   Z <- Z[complete.cases(Z)]
#   if (nrow(Z) > 200) {
#     m_tmp <- lm(as.formula(paste0(y, " ~ ", paste(P_SMALL, collapse = " + "))), data = Z)
#     cat("\nVIF for", y, "\n"); print(try(vif(m_tmp), silent = TRUE))
#   }
# }
# 
# # ====================== EXAMPLES: USE RESULTS =================================
# # Coef table for one outcome (FE-OLS)
# if ("lepto_rate_100k" %in% names(results_all)) print(etable(results_all$lepto_rate_100k$mod_ols))
# 
# # RF importance for one outcome
# if ("lepto_rate_100k" %in% names(results_all)) {
#   imp <- as.data.table(results_all$lepto_rate_100k$mod_rf$variable.importance, keep.rownames = TRUE)
#   setnames(imp, c("rn","V1"), c("variable","importance"))
#   print(imp[order(-importance)][1:20])
# }
# 
# 
# 
# 
# 
