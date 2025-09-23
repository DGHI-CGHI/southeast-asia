# ============================== SETUP =========================================
library(data.table)
library(lubridate)
library(stringr)

library(fixest)       # feols/fepois
library(glmmTMB)      # negbin
library(mgcv)         # GAM
library(ranger)       # RF
library(glmnet)       # LASSO/EN
library(car)          # vif

setDTthreads(percent = 100)

`%||%` <- function(a,b) if (!is.null(a)) a else b
rmse   <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
mae    <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
r2     <- function(y, yhat) 1 - sum((y - yhat)^2, na.rm=TRUE) / sum((y - mean(y, na.rm=TRUE))^2, na.rm=TRUE)

# PATCH: safe newdata for fixest when some regressors were dropped (collin.var)
# Helper: ensure newdata has everything predict() expects
# Keep your existing safe_newdata_fixest; repeating here for clarity
safe_newdata_fixest <- function(mod, TE, fe_keys = c("district","year")) {
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  used_beta <- names(coef(mod))
  dropped   <- mod$collin.var %||% character(0)
  lin_fml   <- mod$fml_all$linear
  toks      <- all.vars(lin_fml)               # includes offset vars
  keep      <- unique(c(fe_keys, used_beta, dropped, toks))
  miss <- setdiff(keep, names(TE))
  if (length(miss)) for (m in miss) TE[[m]] <- NA_real_
  TE[, ..keep]
}

# Ensure TR/TE share factor levels for district & year
harmonize_factors <- function(TR, TE, facs = c("district","year")) {
  for (f in facs) {
    if (!is.factor(TR[[f]])) TR[[f]] <- factor(TR[[f]])
    TE[[f]] <- factor(TE[[f]], levels = levels(TR[[f]]))
  }
  list(TR = TR, TE = TE)
}


ensure_numeric <- function(DT, cols) {
  ok <- intersect(cols, names(DT))
  for (c in ok) {
    if (!is.numeric(DT[[c]])) {
      suppressWarnings({ DT[[c]] <- as.numeric(DT[[c]]) })
    }
  }
  ok
}

harmonize_factors <- function(TR, TE, facs = c("district")) {
  for (f in facs) {
    TR[[f]] <- factor(TR[[f]])
    TE[[f]] <- factor(TE[[f]], levels = levels(TR[[f]]))
  }
  list(TR=TR, TE=TE)
}

# Set cyclic week (1..52) and knots for bs='cc'

# Cyclic week (1..52)
prep_week <- function(D) {
  if (!"week" %in% names(D)) stop("week column required")
  D[, week := as.integer(pmax(1L, pmin(52L, week)))]
  D
}
# Optionally drop one LC to avoid simplex collinearity
drop_one_landcover <- function(PV, drop = c("NoData","Water","BuiltUp")) {
  lc <- c("BuiltUp","Cropland","Grass","Paddy","Shrub","Water","Wetland","NoData")
  present <- intersect(lc, PV)
  if (length(present) <= 1) return(PV)
  d <- intersect(drop, present)
  if (length(d)) return(setdiff(PV, d[1L]))
  setdiff(PV, present[1L]) # fallback: drop the first LC found
}


# ---- main: fit/predict GAM safely ----

# ---- main: fit/predict GAM safely ----
fit_predict_gam_safe <- function(TR, TE, outcome, P_SMALL,
                                 use_bam = FALSE, verbose = TRUE) {
  # factors & week
  hz <- harmonize_factors(TR, TE, facs = c("district"))
  TR <- hz$TR; TE <- hz$TE
  TR <- prep_cyclic_week(as.data.table(TR))
  TE <- prep_cyclic_week(as.data.table(TE))
  
  # keep only columns we'll need
  keep <- unique(c(outcome, "district","week", P_SMALL))
  TR <- TR[, ..keep]; TE <- TE[, ..keep]
  
  # numeric linear terms
  lin_terms <- setdiff(P_SMALL, c("sin52","cos52"))  # you're already using week cyclic
  lin_terms <- ensure_numeric(TR, lin_terms)
  ensure_numeric(TE, lin_terms)
  
  # build RHS pieces that actually exist & are not all NA in TR
  exists_nonNA <- function(cols) cols[cols %in% names(TR) & sapply(cols, \(z) any(is.finite(TR[[z]])))]
  lin_terms <- exists_nonNA(lin_terms)
  
  # district RE only if ???2 levels in TR
  have_RE <- nlevels(TR$district) >= 2
  re_term <- if (have_RE) " + s(district, bs='re')" else ""
  
  # smooths (include only if present)
  s_terms <- c(
    if ("temperature_2m_mean"       %in% lin_terms) "s(temperature_2m_mean)",
    if ("precipitation_sum"         %in% lin_terms) "s(precipitation_sum)",
    if ("temperature_2m_mean_lag1"  %in% lin_terms) "s(temperature_2m_mean_lag1)",
    if ("precipitation_sum_lag1"    %in% lin_terms) "s(precipitation_sum_lag1)"
  )
  s_terms <- paste(s_terms, collapse = " + ")
  
  # linear add-ons besides the smooths (drop the ones already smoothed)
  drop_from_lin <- c("temperature_2m_mean","precipitation_sum",
                     "temperature_2m_mean_lag1","precipitation_sum_lag1")
  lin_keep <- setdiff(lin_terms, drop_from_lin)
  lin_keep <- paste(lin_keep, collapse = " + ")
  
  # cyclic seasonality: bs='cc' requires knots
  week_term <- "s(week, bs='cc', k=20)"
  
  rhs <- paste(
    c(week_term,
      if (nzchar(s_terms)) s_terms else NULL,
      if (nzchar(lin_keep)) lin_keep else NULL),
    collapse = " + "
  )
  rhs <- paste0(rhs, re_term)
  
  f_gam <- as.formula(paste0(outcome, " ~ ", rhs))
  if (verbose) message("GAM formula: ", deparse(f_gam))
  
  # fit
  knots_list <- list(week = c(0.5, 52.5))  # close the circle
  fit_fun <- if (use_bam) mgcv::bam else mgcv::gam
  mod <- fit_fun(f_gam, data = TR, method = "REML",
                 knots = knots_list)
  
  # predict
  yhat <- as.numeric(predict(mod, newdata = TE))
  list(model = mod, pred = yhat)
}

# ============================== DATA ==========================================
DT <- as.data.table(lepto)
if (!inherits(DT$date, "Date")) DT[, date := as.IDate(date)]
DT[, `:=`(year = year(date), week = isoweek(date), yday = yday(date))]

# Build weekly rate per 100k for every *_A count (if not already present)
a_cols    <- grep("_A$", names(DT), value = TRUE)
diseases  <- sub("_A$", "", a_cols)
if ("poptot" %in% names(DT)) {
  for (d in diseases) {
    a   <- paste0(d, "_A")
    rpk <- paste0(d, "_rate_100k")
    if (!rpk %in% names(DT)) DT[, (rpk) := 1e5 * get(a) / poptot]
  }
}

# Choose outcomes (rates). You can restrict to c("lepto_rate_100k", .) if desired.
outcome_cols <- grep("_rate_100k$", names(DT), value = TRUE)
if (!length(outcome_cols)) outcome_cols <- intersect(names(DT), c("lepto_100k","dengue_100k"))
stopifnot(length(outcome_cols) > 0)

# ======================= FEATURE ENGINEERING ==================================
# All ERA5/weather vars you currently have
wx_all <- intersect(c(
  "temperature_2m_max","temperature_2m_min","temperature_2m_mean",
  "apparent_temperature_max","apparent_temperature_min","apparent_temperature_mean",
  "shortwave_radiation_sum","precipitation_sum","rain_sum","precipitation_hours",
  "windspeed_10m_max","windgusts_10m_max","et0_fao_evapotranspiration"
), names(DT))

# Make lags/rolling means by district
setkey(DT, district, date)
make_lags_rolls <- function(DT, v, lags = c(1,2,3,4), k_roll = c(2,4,8)) {
  for (L in lags) DT[, paste0(v, "_lag", L) := shift(get(v), L), by = district]
  for (K in k_roll) DT[, paste0(v, "_roll", K) := frollmean(get(v), K, align="right"), by = district]
  invisible(NULL)
}
for (v in wx_all) make_lags_rolls(DT, v)

# Seasonality
DT[, `:=`(sin52 = sin(2*pi*week/52), cos52 = cos(2*pi*week/52))]

# Land cover + interactions
lc_vars <- intersect(c("BuiltUp","Cropland","Grass","Paddy","Shrub","Water","Wetland"), names(DT))
DT[, `:=`(
  precipX_paddy = (precipitation_sum %||% NA_real_) * (Paddy %||% 0),
  tempX_built   = (temperature_2m_mean %||% NA_real_) * (BuiltUp %||% 0)
)]

# Candidate predictors (kitchen sink; we'll screen below)
PREDICTORS <- unique(c(
  wx_all,
  grep("_(lag|roll)\\d+$", names(DT), value = TRUE),
  lc_vars, "precipX_paddy", "tempX_built",
  "sin52","cos52"
))

# Drop highly NA features / zero variance
na_share <- DT[, lapply(.SD, function(z) mean(is.na(z))), .SDcols = PREDICTORS]
ok_cols  <- names(na_share)[na_share < 0.30]  # <30% missing allowed
zv       <- DT[, sapply(.SD, function(z) data.table(z)[, uniqueN(z, na.rm = TRUE)]), .SDcols = ok_cols]
ok_cols  <- ok_cols[zv > 1]
PREDICTORS <- ok_cols

# Compact set for stability (you can switch to PREDICTORS)
P_SMALL <- unique(c(
  "temperature_2m_mean","precipitation_sum","rain_sum","precipitation_hours",
  "temperature_2m_mean_lag1","temperature_2m_mean_lag2",
  "precipitation_sum_lag1","precipitation_sum_lag2",
  "apparent_temperature_mean","apparent_temperature_mean_lag1",
  "shortwave_radiation_sum","shortwave_radiation_sum_lag1",
  "windspeed_10m_max","windgusts_10m_max",
  "temperature_2m_mean_roll4","precipitation_sum_roll4",
  lc_vars, "precipX_paddy","tempX_built","sin52","cos52"
))
P_SMALL <- intersect(P_SMALL, PREDICTORS)
stopifnot(length(P_SMALL) > 0)

# ========================== SPLITS ============================================
setorder(DT, date)
cut_date <- DT[, date][floor(.N * 0.8)]
TR <- DT[date <= cut_date]
TE <- DT[date  > cut_date]

hz <- harmonize_factors(TR, TE, facs = c("district","year"))
TR <- hz$TR; TE <- hz$TE

# Rolling-origin CV
rolling_origin_eval <- function(DT, outcome, preds, k = 5, family = c("gaussian","poisson","negbin")) {
  family <- match.arg(family)
  dates <- sort(unique(DT$date))
  slices <- split(dates, cut(seq_along(dates), breaks = k, labels = FALSE))
  out <- vector("list", k)
  
  for (i in seq_len(k)) {
    test_dates  <- slices[[i]]
    train_dates <- dates[dates < min(test_dates)]
    tr <- DT[date %in% train_dates]
    te <- DT[date %in% test_dates]
    if (nrow(tr) < 100 || nrow(te) < 20) { out[[i]] <- data.table(); next }
    
    f_rhs <- paste(preds, collapse = " + ")
    
    if (family == "gaussian") {
      m <- feols(as.formula(paste0(outcome, " ~ ", f_rhs, " | district + year")),
                 data = tr, vcov = ~district)
      yhat <- as.numeric(predict(m, newdata = safe_newdata_fixest(m, te)))
      y    <- te[[outcome]]
      
    } else if (family == "poisson") {
      cnt  <- sub("_rate_100k$", "_A", outcome)
      if (!cnt %in% names(te)) { out[[i]] <- data.table(); next }
      m <- fepois(as.formula(paste0(cnt, " ~ ", f_rhs, " | district + year")),
                  data = tr, offset = log(tr$poptot), vcov = ~district)
      yhat <- as.numeric(predict(m, newdata = safe_newdata_fixest(m, te), type = "response"))
      y    <- te[[cnt]]
      
    } else {
      cnt  <- sub("_rate_100k$", "_A", outcome)
      if (!cnt %in% names(te)) { out[[i]] <- data.table(); next }
      f <- as.formula(paste0(cnt, " ~ ", f_rhs, " + factor(district) + factor(year) + offset(log(poptot))"))
      m <- try(glmmTMB(family = nbinom2(), data = tr, formula = f), silent = TRUE)
      if (inherits(m, "try-error")) { out[[i]] <- data.table(); next }
      yhat <- as.numeric(predict(m, newdata = te, type = "response"))
      y    <- te[[cnt]]
    }
    
    out[[i]] <- data.table(fold = i, outcome = outcome,
                           RMSE = rmse(y, yhat), MAE = mae(y, yhat), R2 = r2(y, yhat))
  }
  rbindlist(out, fill = TRUE)
}

rolling_origin_eval <- function(
    DT, outcome, P_SMALL, k = 5, family = "gaussian",
    count_name = NULL, count_map = NULL
) {
  stopifnot(all(c("district","year","week","poptot") %in% names(DT)))
  
  # resolve count column if needed
  cnt <- NA_character_
  if (family %in% c("poisson","negbin")) {
    cnt <- resolve_count_name(outcome, count_name, count_map)
  }
  
  # drop one LC to reduce simplex collinearity noise (optional)
  drop_one_landcover <- function(PV, drop = c("NoData","Water","BuiltUp")) {
    lc <- c("BuiltUp","Cropland","Grass","Paddy","Shrub","Water","Wetland","NoData")
    present <- intersect(lc, PV)
    if (length(present) <= 1) return(PV)
    d <- intersect(drop, present)
    if (length(d)) return(setdiff(PV, d[1L]))
    setdiff(PV, present[1L])
  }
  P_SMALL <- drop_one_landcover(P_SMALL)
  
  keep <- unique(c("district","year","week","poptot", outcome, P_SMALL))
  if (family %in% c("poisson","negbin")) keep <- unique(c(keep, cnt))
  D <- data.table::copy(DT[, ..keep])
  
  # week to 1..52
  if (!"week" %in% names(D)) stop("week column required in DT")
  D[, week := as.integer(pmax(1L, pmin(52L, week)))]
  
  # ensure numeric for linear terms
  ensure_numeric <- function(DT, cols) {
    ok <- intersect(cols, names(DT))
    for (c in ok) if (!is.numeric(DT[[c]]))
      suppressWarnings(DT[[c]] <- as.numeric(DT[[c]]))
    invisible(NULL)
  }
  ensure_numeric(D, setdiff(P_SMALL, c("district","year","week")))
  
  data.table::setorder(D, year, week, district)
  n <- nrow(D)
  splits <- floor(seq(0, n, length.out = k + 1L))
  out <- vector("list", k)
  
  harmonize_factors <- function(TR, TE, facs = c("district","year")) {
    for (f in facs) {
      TR[[f]] <- factor(TR[[f]])
      TE[[f]] <- factor(TE[[f]], levels = levels(TR[[f]]))
    }
    list(TR=TR, TE=TE)
  }
  quiet_fixest <- function(expr) {
    suppressMessages(suppressWarnings(capture.output(res <- force(expr))))
    res
  }
  rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))
  mae  <- function(a, b) mean(abs(a - b), na.rm = TRUE)
  r2   <- function(y, yhat) 1 - sum((y - yhat)^2, na.rm = TRUE) /
    sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  
  for (i in seq_len(k)) {
    lo <- splits[i] + 1L; hi <- splits[i + 1L]
    if (hi <= lo) next
    te_idx <- lo:hi; tr_idx <- setdiff(seq_len(n), te_idx)
    tr <- D[tr_idx]; te <- D[te_idx]
    hz <- harmonize_factors(tr, te); tr <- hz$TR; te <- hz$TE
    
    f_rhs <- paste(P_SMALL, collapse = " + ")
    
    if (family == "gaussian") {
      f_ols <- as.formula(paste0(outcome, " ~ ", f_rhs, " | district + year"))
      m <- quiet_fixest(fixest::feols(f_ols, data = tr, vcov = ~district))
      yhat <- as.numeric(predict(m, newdata = safe_newdata_fixest(m, te)))
      obs  <- te[[outcome]]
      
    } else if (family == "poisson") {
      stopifnot(cnt %in% names(tr), cnt %in% names(te))
      f_pois <- as.formula(paste0(
        cnt, " ~ ", f_rhs, " + offset(log(poptot)) | district + year"
      ))
      m <- quiet_fixest(fixest::fepois(f_pois, data = tr, vcov = ~district))
      yhat <- as.numeric(predict(m, newdata = safe_newdata_fixest(m, te), type = "response"))
      obs  <- te[[cnt]]
      
    } else if (family == "negbin") {
      stopifnot(cnt %in% names(tr), cnt %in% names(te))
      f_nb <- as.formula(paste0(
        cnt, " ~ ", f_rhs, " + (1|district) + (1|year) + offset(log(poptot))"
      ))
      m <- try(glmmTMB::glmmTMB(formula = f_nb, data = tr, family = glmmTMB::nbinom2(),
                                control = glmmTMB::glmmTMBControl(parallel = 1L,
                                                                  optimizer = optim,
                                                                  optArgs = list(method = "BFGS"))),
               silent = TRUE)
      if (inherits(m, "try-error")) {
        yhat <- rep(NA_real_, nrow(te))
      } else {
        yhat <- as.numeric(predict(m, newdata = te, type = "response", allow.new.levels = TRUE))
      }
      obs <- te[[cnt]]
      
    } else stop("family must be 'gaussian', 'poisson', or 'negbin'.")
    
    out[[i]] <- data.table::data.table(
      fold = i, family = family, outcome = outcome,
      RMSE = rmse(obs, yhat), MAE = mae(obs, yhat), R2 = r2(obs, yhat)
    )
  }
  data.table::rbindlist(out, fill = TRUE)
}


# ========================== MODEL LOOP ========================================
results_all <- list()

library(jj)
outcome_cols = outcome_cols[outcome_cols %likeany% c("dengue","lepto")]

outcome_cols = c("lepto","dengue")


outcome_cols = 'lepto'



for (y in outcome_cols) {
  message("=== OUTCOME: ", y, " ===")
  
  # ---------- (A) FE-OLS on rate --------------------------------------------
  f_rhs <- paste(P_SMALL, collapse = " + ")
  f_ols <- as.formula(paste0(y, " ~ ", f_rhs, " | district + year"))
  # mod_ols <- feols(f_ols, data = TR, vcov = ~district)
  
  mod_ols  <- feols(f_ols, data = TR, vcov = ~district)
  new_ols  <- safe_newdata_fixest(mod_ols, TE)   # <-- use patched helper
  yhat_ols <- as.numeric(predict(mod_ols, newdata = new_ols))
  
  
  # new_ols <- safe_newdata_fixest(mod_ols, TE)
  # yhat_ols <- as.numeric(predict(mod_ols, newdata = new_ols))
  perf_ols <- data.table(model = "FE-OLS", outcome = y,
                         RMSE = rmse(TE[[y]], yhat_ols),
                         MAE  = mae(TE[[y]],  yhat_ols),
                         R2   = r2(TE[[y]],   yhat_ols))
  
  # ---------- (B) FE-Poisson (counts + offset) ------------------------------
  cnt <- sub("_rate_100k$", "_A", y)
  perf_pois <- perf_nb <- NULL
  
  
  if (cnt %in% names(TR)) {
    # Put the offset inside the formula (critical!)
    f_pois <- as.formula(paste0(
      cnt, " ~ ", f_rhs, " + offset(log(poptot)) | district + year"
    ))
    
    mod_pois <- fepois(f_pois, data = TR, vcov = ~district)
    
    # Build newdata that includes FE keys, used slopes, dropped vars, and poptot
    new_pois <- safe_newdata_fixest(mod_pois, TE)
    # predict() will evaluate offset(log(poptot)) on new_pois
    mu_pois  <- as.numeric(predict(mod_pois, newdata = new_pois, type = "response"))
    
    perf_pois <- data.table(model = "FE-Poisson", outcome = y,
                            RMSE = rmse(TE[[cnt]], mu_pois),
                            MAE  = mae(TE[[cnt]],  mu_pois),
                            R2   = r2(TE[[cnt]],   mu_pois))
    
    # ---------- (C) NegBin (glmmTMB FE via factors) ----------
    f_nb <- as.formula(paste0(
      cnt, " ~ ", f_rhs,
      " + factor(district) + factor(year) + offset(log(poptot))"
    ))
    mod_nb <- try(glmmTMB(family = nbinom2(), data = TR, formula = f_nb), silent = TRUE)
    if (!inherits(mod_nb, "try-error")) {
      # TE must have poptot too
      stopifnot("poptot" %in% names(TE))
      mu_nb <- as.numeric(predict(mod_nb, newdata = TE, type = "response"))
      perf_nb <- data.table(model = "NB (glmmTMB FE)", outcome = y,
                            RMSE = rmse(TE[[cnt]], mu_nb),
                            MAE  = mae(TE[[cnt]],  mu_nb),
                            R2   = r2(TE[[cnt]],   mu_nb))
    }
  }
  
  
  # ---------- (D) GAM (rate) ------------------------------------------------
  # smooths for a few key ERA5 vars + cyclic week + RE for district
  f_gam <- as.formula(
    paste0(y, " ~ s(temperature_2m_mean) + s(precipitation_sum) + ",
           "s(temperature_2m_mean_lag1) + s(precipitation_sum_lag1) + ",
           "s(week, bs='cc', k=20) + s(district, bs='re') + ",
           paste(setdiff(P_SMALL, c("temperature_2m_mean","precipitation_sum",
                                    "temperature_2m_mean_lag1","precipitation_sum_lag1",
                                    "sin52","cos52")), collapse = " + "))
  )
  
  # TR$week
  
  gp <- fit_predict_gam_safe(TR, TE, outcome = y, P_SMALL = P_SMALL, use_bam = FALSE)
  mod_gam  <- gp$model
  yhat_gam <- gp$pred
  
  
  perf_gam <- data.table(model = "GAM (rate)", outcome = y,
                         RMSE = rmse(TE[[y]], yhat_gam),
                         MAE  = mae(TE[[y]],  yhat_gam),
                         R2   = r2(TE[[y]],   yhat_gam))
  
  
  
  
  # ---------- (E) Random Forest (rate) --------------------------------------
  Xcols <- unique(c("district","year", P_SMALL))
  X_tr  <- TR[, ..Xcols]; y_tr <- TR[[y]]
  X_te  <- TE[, ..Xcols]; y_te <- TE[[y]]
  for (cc in c("district","year")) {
    if (cc %in% names(X_tr)) { X_tr[[cc]] <- as.factor(X_tr[[cc]]); X_te[[cc]] <- as.factor(X_te[[cc]]) }
  }
  X_tr$y_resp <- y_tr
  mod_rf <- ranger(
    y_resp ~ ., data = X_tr,
    num.trees = 800, mtry = max(3, floor(sqrt(ncol(X_tr)))),
    min.node.size = 10, importance = "impurity", seed = 42
  )
  yhat_rf <- predict(mod_rf, data = X_te)$predictions
  perf_rf <- data.table(model = "RandomForest (rate)", outcome = y,
                        RMSE = rmse(y_te, yhat_rf),
                        MAE  = mae(y_te,  yhat_rf),
                        R2   = r2(y_te,   yhat_rf))
  
  # ---------- (F) LASSO / Elastic-Net (rate) --------------------------------
  # Build a numeric model matrix (drop district/year here; add as dummies if needed)
  mm_form <- as.formula(paste0("~ 0 + ", paste(P_SMALL, collapse = " + ")))
  Xmm_tr  <- model.matrix(mm_form, data = TR)
  Xmm_te  <- model.matrix(mm_form, data = TE)
  y_tr_mm <- TR[[y]]; y_te_mm <- TE[[y]]
  
  cvfit <- cv.glmnet(Xmm_tr, y_tr_mm, alpha = 1, family = "gaussian", nfolds = 5)  # LASSO
  yhat_lasso <- as.numeric(predict(cvfit, newx = Xmm_te, s = "lambda.min"))
  perf_lasso <- data.table(model = "LASSO (rate)", outcome = y,
                           RMSE = rmse(y_te_mm, yhat_lasso),
                           MAE  = mae(y_te_mm,  yhat_lasso),
                           R2   = r2(y_te_mm,   yhat_lasso))
  
  quiet_fixest <- function(expr) {
    # Suppress fixest's collinearity chatter within a block
    suppressMessages(suppressWarnings(capture.output(res <- force(expr))))
    invisible(res); res
  }
  
  # ---------- Rolling-Origin CV (brief) -------------------------------------
  # cv_gauss <- rolling_origin_eval(DT, y, P_SMALL, k = 5, family = "gaussian")
  # if (cnt %in% names(DT)) {
  #   cv_pois  <- rolling_origin_eval(DT, y, P_SMALL, k = 5, family = "poisson")
  #   cv_negb  <- rolling_origin_eval(DT, y, P_SMALL, k = 5, family = "negbin")
  # } else {
  #   cv_pois <- cv_negb <- data.table()
  # }
  
  # Map outcome (rate) -> count column
  resolve_count_name <- function(outcome, count_name = NULL, count_map = NULL) {
    # 1) explicit
    if (!is.null(count_name)) return(count_name)
    # 2) named list map
    if (!is.null(count_map) && outcome %in% names(count_map)) return(count_map[[outcome]])
    # 3) common heuristics
    #    *_rate_100k -> *_A
    if (grepl("_rate_100k$", outcome)) return(sub("_rate_100k$", "_A", outcome))
    #    dengue_100k -> dengue_A
    if (grepl("_100k$", outcome))       return(sub("_100k$", "_A", outcome))
    #    lepto_100k -> leptospirosis_A (special case)
    if (outcome == "lepto_100k")        return("leptospirosis_A")
    stop("Can't infer count column for outcome '", outcome,
         "'. Provide count_name= or count_map=.")
  }
  
  
  # You have these columns per your names():
  #   rates: dengue_100k, lepto_100k
  #   counts: dengue_A,   leptospirosis_A
  count_map <- list(
    dengue_100k = "dengue_A",
    lepto_100k  = "leptospirosis_A"
  )
  
  # Example calls:
  cv_gauss <- rolling_origin_eval(DT, "dengue_100k", P_SMALL, k = 5, family = "gaussian")
  
  cv_pois  <- rolling_origin_eval(DT, "dengue_100k", P_SMALL, k = 5,
                                  family = "poisson",  count_map = count_map)
  
  # cv_negb  <- rolling_origin_eval(DT, "lepto_100k",  P_SMALL, k = 5,
  #                                 family = "negbin",  count_map = count_map)
  # 
  # ---------- Collect --------------------------------------------------------
  perf_all <- rbindlist(list(perf_ols, perf_pois, perf_nb, perf_gam, perf_rf, perf_lasso), fill = TRUE)
  results_all[[y]] <- list(
    perf_test = perf_all,
    cv_gauss  = cv_gauss,
    cv_pois   = cv_pois,
    # cv_negb   = cv_negb,
    mod_ols   = mod_ols,
    mod_gam   = mod_gam,
    mod_rf    = mod_rf,
    lasso     = cvfit
  )
  
  cat(y, '\n')
}

# Test performance across outcomes
perf_summary <- rbindlist(lapply(results_all, `[[`, "perf_test"), use.names = TRUE, fill = TRUE)
print(perf_summary[order(outcome, -R2)])

# ====================== OPTIONAL DIAGNOSTICS ==================================
# VIF (simple OLS on TR for one or two outcomes)
for (y in head(outcome_cols, 2L)) {
  cols_needed <- unique(c(y, P_SMALL))
  Z <- TR[, ..cols_needed]
  Z <- Z[complete.cases(Z)]
  if (nrow(Z) > 200) {
    m_tmp <- lm(as.formula(paste0(y, " ~ ", paste(P_SMALL, collapse = " + "))), data = Z)
    cat("\nVIF for", y, "\n"); print(try(vif(m_tmp), silent = TRUE))
  }
}

# ====================== EXAMPLES: USE RESULTS =================================
# Coef table for one outcome (FE-OLS)
if ("lepto_rate_100k" %in% names(results_all)) print(etable(results_all$lepto_rate_100k$mod_ols))

# RF importance for one outcome
if ("lepto_rate_100k" %in% names(results_all)) {
  imp <- as.data.table(results_all$lepto_rate_100k$mod_rf$variable.importance, keep.rownames = TRUE)
  setnames(imp, c("rn","V1"), c("variable","importance"))
  print(imp[order(-importance)][1:20])
}





