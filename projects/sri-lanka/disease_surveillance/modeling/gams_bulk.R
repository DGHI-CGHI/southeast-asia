# ?????? Dependencies ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
library(data.table)

# ?????? Utilities ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
`%||%` <- function(a,b) if (!is.null(a)) a else b

drop_constant_cols <- function(DT) {
  stopifnot(inherits(DT, "data.table"))
  keep <- names(DT)[vapply(DT, function(x) length(unique(x)) > 1L, logical(1))]
  DT[, ..keep]
}

# Safer predict() newdata builder for fixest: keeps only used slope vars and FE keys if present
safe_newdata_fixest <- function(mod, newdata) {
  out <- copy(as.data.table(newdata))
  # slope names actually in the model (exclude intercept)
  used <- names(mod$coefficients)
  used <- setdiff(used, "(Intercept)")
  # add any missing slope vars as 0 so predict() won't choke
  miss <- setdiff(used, names(out))
  if (length(miss)) for (m in miss) out[[m]] <- 0
  out
}

# ?????? Seasonality helpers ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# Adds ISO week if missing, then sin/cos seasonal terms with given period (default 52)
add_seasonal_terms <- function(DT, date_col = "date", week_col = "week", period = 52L) {
  stopifnot(inherits(DT, "data.table"))
  if (!inherits(DT[[date_col]], "Date")) DT[, (date_col) := as.IDate(get(date_col))]
  if (!week_col %in% names(DT)) {
    DT[, (week_col) := as.integer(strftime(get(date_col), "%V"))]
    DT[get(week_col) < 1L | get(week_col) > period, (week_col) := pmin(period, pmax(1L, get(week_col)))]
  }
  DT[, sin52 := sin(2*pi*get(week_col)/period)]
  DT[, cos52 := cos(2*pi*get(week_col)/period)]
  DT[]
}

# ?????? Rolling means (optional) ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# Fast within-group trailing rolling means; aligns to the right (past only)
make_rollmeans <- function(DT, cols, k = 3L, by = "district", order_by = "date", suffix = NULL) {
  stopifnot(inherits(DT, "data.table"))
  suff <- suffix %||% paste0("_roll", k)
  setorderv(DT, c(by, order_by))
  for (v in cols) {
    newv <- paste0(v, suff)
    DT[, (newv) := frollmean(get(v), n = k, align = "right", na.rm = TRUE), by = by]
  }
  DT[]
}

# ?????? Generic lag maker ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# Creates group-wise lags for multiple variables and multiple lag orders.
# Example: make_lags(DT, c("precipitation_sum","temperature_2m_mean"), lags = 1:4)
make_lags <- function(DT, cols, lags = 1L, by = "district", order_by = "date",
                      suffix = "lag", inplace = FALSE) {
  stopifnot(inherits(DT, "data.table"))
  X <- if (isTRUE(inplace)) DT else copy(DT)
  setorderv(X, c(by, order_by))
  for (v in cols) {
    for (k in lags) {
      newv <- paste0(v, "_", suffix, k)
      X[, (newv) := shift(get(v), n = k, type = "lag"), by = by]
    }
  }
  X[]
}

# ?????? Autoregressive lags for the outcome ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
# Adds AR lags of the outcome you plan to model (either a *rate* column or a *count* column).
# Example (rate): add_ar_start(DT, outcome_col = "lepto_100k", lags = 1:4)
# Example (count): add_ar_start(DT, outcome_col = "leptospirosis_A", lags = 1:4)
add_ar_start <- function(DT, outcome_col, lags = 1L, by = "district", order_by = "date",
                         suffix = "lag", keep_first_nas = TRUE) {
  stopifnot(inherits(DT, "data.table"))
  stopifnot(outcome_col %in% names(DT))
  X <- copy(DT)
  setorderv(X, c(by, order_by))
  for (k in lags) {
    newv <- paste0(outcome_col, "_", suffix, k)
    X[, (newv) := shift(get(outcome_col), n = k, type = "lag"), by = by]
  }
  if (!keep_first_nas) {
    # Optionally drop rows that have incomplete AR history at series starts
    needed <- paste0(outcome_col, "_", suffix, lags)
    X <- X[complete.cases(X[, ..needed])]
  }
  X[]
}

# ?????? Tiny scoring helpers (for completeness) ??????????????????????????????????????????????????????????????????????????????????????????????????????
rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))
mae  <- function(a, b) mean(abs(a - b), na.rm = TRUE)
r2   <- function(y, yhat) 1 - sum((y - yhat)^2, na.rm = TRUE) / sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
to_rate100k <- function(mu_counts, pop) 100000 * (mu_counts / pmax(pop, 1))

# Make weekly seasonality if absent
ensure_week_terms <- function(DT, period = 52L) {
  if (!("week" %in% names(DT))) {
    DT[, week := as.integer(strftime(as.IDate(date), "%V"))]
    DT[is.na(week) | week < 1 | week > period, week := pmin(period, pmax(1L, week))]
  }
  DT[]
}

# ---------- formula builders ----------
smooth_term <- function(var, bs = c("tp","cr","cc","re"), k = NULL) {
  bs <- match.arg(bs)
  if (bs == "re")   return(sprintf("s(%s, bs='re')", var))
  if (bs == "cc")   return(sprintf("s(%s, bs='cc'%s)", var, if (!is.null(k)) paste0(", k=", k) else ""))
  # 'tp' or 'cr'
  sprintf("s(%s, bs='%s'%s)", var, bs, if (!is.null(k)) paste0(", k=", k) else "")
}

linear_term <- function(var) var

# Build term vector from available columns safely
avail_terms <- function(vars, data, fun = linear_term, ...) {
  vars <- intersect(vars, names(data))
  if (!length(vars)) return(character(0))
  vapply(vars, function(v) fun(v, ...), character(1))
}




# Join back onto your epi table (district×week outcomes):
lepto_feat <- features_weekly[lepto, on = .(district, date_start, date_end)]
# Now model with Poisson/NB using offset(log(poptot)) and the features you like.

station_vars = c("temperature_2m_max","temperature_2m_min",'temperature_2m_mean',
                 'apparent_temperature_mean',
                 "apparent_temperature_max","apparent_temperature_min",
                 "shortwave_radiation_sum","precipitation_sum",
                 "rain_sum","precipitation_hours",
                 "windspeed_10m_max","windgusts_10m_max",
                 "et0_fao_evapotranspiration")
lepto_feat = lepto_feat[,!station_vars,with=FALSE]




as.data.frame(
  countna(lepto_feat))


lepto_feat = lepto_feat[!is.na(vpd_mean_week_z)]
# Use your table
DT <- as.data.table(lepto_feat)

# Outcome choices (rates)
OUTCOME <- "lepto"   # or "dengue_100k"

# # Weather/ERA-5 pool (tailor as you like)
# WX_BASE <- c(
#   "temperature_2m_mean", "precipitation_sum",
#   "apparent_temperature_mean", "shortwave_radiation_sum",
#   "windspeed_10m_max", "windgusts_10m_max", "et0_fao_evapotranspiration"
# )

# WX_BASE <- c('ta_max','td_mean_aw','rh_mean_aw')



WX_BASE = names(DT) [names(DT) %likeany% c("ta_","precip","td","rh","wind",'vpd','wet_days','tmax')]
WX_BASE = (WX_BASE)[!(WX_BASE) %like% "lag"]



# Land cover (drop one later)
LC_BASE <- c("BuiltUp","Cropland","Grass","Paddy","Shrub","Water","Wetland","NoData")

# Ensure existence
# require_cols(DT, c("district","date", OUTCOME, WX_BASE, LC_BASE))

# Add derived columns

DT[, week := as.integer(strftime(date, "%V"))]
DT[week < 1 | week > 52, week := pmin(52L, pmax(1L, week))]

DT = make_lags(DT, cols = WX_BASE, lags = c(1,2,3,4))
# DT = add_ar_start(DT)
DT <- add_ar_start(DT, outcome_col = "lepto", lags = 1:4, by = "district", order_by = "date")

# DT$lepto

# Optional: scale continuous predictors (mgcv handles scale, but scaling can help)
scale_cols <- function(D, cols) {
  for (c in cols) if (is.numeric(D[[c]])) {
    mu <- mean(D[[c]], na.rm = TRUE); sdv <- sd(D[[c]], na.rm = TRUE)
    if (is.finite(sdv) && sdv > 0) D[[c]] <- (D[[c]] - mu)/sdv
  }
  D
}
CONT <- intersect(c(WX_BASE, grep("^temperature_2m_mean_lag|^precipitation_sum_lag", names(DT), value = TRUE), LC_BASE), names(DT))
DT <- scale_cols(DT, CONT)

# Keep complete rows for modeling
MODEL_COLS <- unique(c("district","date","week", OUTCOME, WX_BASE, LC_BASE,
                       grep("_lag[1-4]$", names(DT), value = TRUE)))
DTm <- DT[complete.cases(DT[, ..MODEL_COLS])]




# Which land-cover to drop (avoid perfect collinearity)
LC_USE <- drop_one_landcover(LC_BASE, drop = c("NoData","Water","BuiltUp"))

CFG_GRID <- CJ(
  use_lags      = c("none","lag1","lag12","lag14"),
  use_lc        = c(FALSE, TRUE),
  use_tensor    = c(FALSE, TRUE),
  bs_kind       = c("cs","tp"),
  use_ar1       = c(FALSE, TRUE)
)





# Construct a GAM formula string given the outcome and cfg
build_gam_formula <- function(outcome,
                              cfg,
                              data_train = NULL,
                              pop_var = "poptot") {
  # cfg fields (all optional; defaults set here)
  wx_vars     <- cfg$wx_vars     %||% character(0)      # e.g., c("ta_max","td_mean_aw","rh_mean_aw")
  wx_lags     <- cfg$wx_lags     %||% integer(0)        # e.g., 1:4
  lc_vars     <- cfg$lc_vars     %||% character(0)      # land cover terms (linear)
  ar_lags     <- cfg$ar_lags     %||% integer(0)        # AR lags for outcome
  season_type <- cfg$season      %||% "cyclic"          # "cyclic"|"fourier"|"none"
  k_cyclic    <- cfg$k_cyclic    %||% 20
  wx_bs       <- cfg$wx_bs       %||% "tp"              # "tp","cr"
  district_re <- isTRUE(cfg$district_re %||% TRUE)      # add s(district, bs='re')
  family      <- cfg$family      %||% "poisson"         # "poisson" for counts(+offset) or "gaussian" for rates
  use_offset  <- isTRUE(cfg$use_offset %||% (family=="poisson"))
  include_pop <- isTRUE(cfg$include_pop %||% TRUE)
  
  # If we know training data, only keep columns that exist
  if (!is.null(data_train)) {
    wx_vars <- intersect(wx_vars, names(data_train))
    lc_vars <- intersect(lc_vars, names(data_train))
  }
  
  # Weather (contemporaneous + lags)
  wx_all <- wx_vars
  if (length(wx_lags)) {
    lag_names <- as.vector(outer(wx_vars, paste0("_lag", wx_lags), paste0))
    wx_all <- unique(c(wx_all, lag_names))
  }
  wx_s <- if (length(wx_all)) paste(avail_terms(wx_all, data_train %||% data.frame(), smooth_term, bs = wx_bs), collapse = " + ") else ""
  
  # Seasonality
  seas <- switch(season_type,
                 "none"    = "",
                 "fourier" = "sin52 + cos52",
                 "cyclic"  = smooth_term("week", bs = "cc", k = k_cyclic),
                 smooth_term("week", bs = "cc", k = k_cyclic)
  )
  
  # District RE
  dist_re <- if (district_re) smooth_term("district", bs = "re") else ""
  
  # Land cover (linear)
  lc <- if (length(lc_vars)) paste(avail_terms(lc_vars, data_train %||% data.frame(), linear_term), collapse = " + ") else ""
  
  # AR lags for outcome (linear)
  ar <- if (length(ar_lags)) paste(paste0(outcome, "_lag", ar_lags), collapse = " + ") else ""
  
  rhs_parts <- c(wx_s, seas, dist_re, lc, ar)
  rhs <- paste(rhs_parts[nzchar(rhs_parts)], collapse = " + ")
  if (!nzchar(rhs)) rhs <- "1"
  
  # Offset for counts
  if (family == "poisson" && use_offset) {
    rhs <- paste("offset(log(", pop_var, ")) + ", rhs)
  }
  
  as.formula(paste0(outcome, " ~ ", rhs))
}





# ---------- grid evaluation ----------
evaluate_gam_grid <- function(TR, TE, outcome, CFG_GRID,
                              pop_var = "poptot",
                              eval_on = c("auto","rate","count")) {
  eval_on <- match.arg(eval_on)
  res <- vector("list", nrow(CFG_GRID))
  for (i in seq_len(nrow(CFG_GRID))) {
    cfg <- as.list(CFG_GRID[i])
    gp  <- fit_predict_gam_cfg(TR, TE, outcome, cfg, pop_var = pop_var, eval_on = eval_on)
    res[[i]] <- cbind(data.table(cfg_id = i), as.data.table(gp$metrics),
                      data.table(family = cfg$family %||% "poisson",
                                 season  = cfg$season %||% "cyclic",
                                 district_re = cfg$district_re %||% TRUE,
                                 wx_bs   = cfg$wx_bs %||% "tp",
                                 wx_vars = paste(cfg$wx_vars %||% character(0), collapse=","), 
                                 wx_lags = paste(cfg$wx_lags %||% integer(0), collapse=","), 
                                 lc_vars = paste(cfg$lc_vars %||% character(0), collapse=","),
                                 ar_lags = paste(cfg$ar_lags %||% integer(0), collapse=",")),
                      data.table(formula = deparse(gp$formula)))
  }
  rbindlist(res, fill = TRUE)
}

# ---------- convenience: time split ----------
time_split <- function(DT, prop = 0.8) {
  X <- copy(DT)
  if (!inherits(X$date, "Date")) X[, date := as.IDate(date)]
  setorder(X, date, district)
  u <- sort(unique(X$date))
  cut <- u[floor(length(u)*prop)]
  list(TR = X[date <= cut], TE = X[date > cut])
}



# Synchronize factor levels in TR and TE so prediction won't fail on unseen levels
factor_sync <- function(TR, TE, facs = c("district","year")) {
  TR <- as.data.table(TR); TE <- as.data.table(TE)
  for (f in facs) {
    if (!f %in% names(TR) || !f %in% names(TE)) next
    TR[[f]] <- as.factor(TR[[f]])
    # Restrict TE levels to those seen in TR; unknowns become NA
    TE[[f]] <- factor(as.character(TE[[f]]), levels = levels(TR[[f]]))
  }
  list(TR = TR, TE = TE)
}

# Utility: extract LHS variable name from a formula
lhs_var <- function(fml) {
  # as.character(~) returns c("~", "lhs", "rhs"); index 2 is lhs
  as.character(fml)[2]
}

# Detect whether an offset(log(poptot)) term is present in the formula
has_pop_offset <- function(fml, pop_var = "poptot") {
  grepl(sprintf("offset\\s*\\(\\s*log\\s*\\(\\s*%s\\s*\\)\\s*\\)", pop_var), deparse(fml))
}

# Guess family from the formula (Poisson if offset is present; otherwise Gaussian), unless provided explicitly
guess_family <- function(fml, family = NULL, pop_var = "poptot") {
  if (!is.null(family)) return(family)
  if (has_pop_offset(fml, pop_var = pop_var)) "poisson" else "gaussian"
}

# Optionally estimate a rough AR(1) rho from training residuals of a quick Gaussian fit
# If that seems too heuristic, set a fixed rho (e.g., 0.3).
estimate_rho <- function(y, group, time) {
  DT <- data.table(y = y, g = group, t = time)
  setorder(DT, g, t)
  # simple pooled lag-1 correlation within groups
  rhos <- DT[, {
    z <- y
    if (length(z) < 3) return(NA_real_)
    stats::cor(z[-1], z[-length(z)], use = "pairwise.complete.obs")
  }, by = g]$V1
  rhos <- rhos[is.finite(rhos)]
  rho <- median(rhos, na.rm = TRUE)
  if (!is.finite(rho)) rho <- 0.3
  pmin(pmax(rho, 0.0), 0.95)
}


fit_predict_score <- function(TR, TE, fml, use_ar1 = FALSE) {
  # harmonize factors
  sync <- factor_sync(TR, TE, facs = "district")
  TR <- sync$TR; TE <- sync$TE
  
  # AR(1) setup
  if (use_ar1) {
    stopifnot("AR_start" %in% names(TR), "AR_start" %in% names(TE))
    # estimate rho on train via residual autocorrelation of a simple model
    # heuristic: rho = 0.3 if unknown
    rho_est <- 0.3
    mod <- bam(fml, data = TR, family = gaussian(), method = "fREML",
               discrete = TRUE, select = TRUE, rho = rho_est, AR.start = TR$AR_start)
    yhat <- as.numeric(predict(mod, newdata = TE, discrete = TRUE,
                               AR.start = TE$AR_start))
  } else {
    mod <- bam(fml, data = TR, family = gaussian(), method = "fREML",
               discrete = TRUE, select = TRUE)
    yhat <- as.numeric(predict(mod, newdata = TE, discrete = TRUE))
  }
  
  list(
    model = mod,
    metrics = data.table(
      RMSE = rmse(TE[[all.vars(fml)[1]]], yhat),
      MAE  = mae(TE[[all.vars(fml)[1]]],  yhat),
      R2   = r2(TE[[all.vars(fml)[1]]],   yhat),
      AIC  = AIC(mod),
      GCV  = mod$gcv.ubre %||% NA_real_,
      DevExpl = summary(mod)$dev.expl %||% NA_real_
    )
  )
}
`%||%` <- function(a,b) if (!is.null(a)) a else b





safe_fit <- function(expr) {
  out <- try(eval.parent(substitute(expr)), silent = TRUE)
  if (inherits(out, "try-error")) {
    list(ok = FALSE, error = attr(out, "condition")$message)
  } else list(ok = TRUE, value = out)
}



library(future)
library(future.apply)

plan(multisession, workers = 6)

split <- time_split(DTm)
TR <- split$TR; TE <- split$TE

OUTCOME <- "lepto"  # change to dengue_100k, etc.

# Evaluate all configurations
results <- future_lapply(seq_len(nrow(CFG_GRID)), function(i) {
  cfg <- as.list(CFG_GRID[i])
  fml <- build_gam_formula(OUTCOME, cfg)
  
  fit <- try(fit_predict_score(TR, TE, fml, use_ar1 = cfg$use_ar1), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(data.table(
      i = i, use_lags = cfg$use_lags, use_lc = cfg$use_lc, use_tensor = cfg$use_tensor,
      bs_kind = cfg$bs_kind, use_ar1 = cfg$use_ar1,
      RMSE = NA_real_, MAE = NA_real_, R2 = NA_real_, AIC = NA_real_, GCV = NA_real_,
      DevExpl = NA_real_, ok = FALSE
    ))
  } else {
    cbind(
      data.table(i = i, use_lags = cfg$use_lags, use_lc = cfg$use_lc, use_tensor = cfg$use_tensor,
                 bs_kind = cfg$bs_kind, use_ar1 = cfg$use_ar1, ok = TRUE),
      fit$metrics
    )
  }
})

plan(sequential)

GRID_SCORES <- rbindlist(results, fill = TRUE)
setorder(GRID_SCORES, RMSE)
GRID_SCORES[]

WX_BASE


lepto$doy = yday(lepto$date)
lepto$year = lepto$year.x
lepto$district_factor = as.factor(lepto$district)


agam = gam(lepto ~ s(ta_max, bs = 'cr', k =10)
           + s(rh_mean_aw, bs = 'cr', k =10, fx = TRUE)
           # + s(rh_mean_aw, by = district_factor, bs = 'cr', k =5, fx = FALSE)
           + s(precip_tp_sum_week, k = 10, bs = 'cc') 
           + s(td_max, bs = 'cr', k =10, fx = TRUE)
           + s(year, bs = 'cc')
           + s(doy, bs = 'cc')
           + district
           , data = lepto_feat, offset = log(poptot))


summary(agam)



plot(agam)
abline(h=0)



quantile(lepto$rh_mean_aw)





# lepto_feat$precip_tp_sum_week_clim
names(lepto_feat)


names(lepto_feat)[names(lepto_feat) == 'precip_tp_sum_week_clim']



range(lepto_feat$precip_tp_sum_week_clim, na.rm=TRUE)
range(lepto_feat$precip_tp_sum_week, na.rm=TRUE)



# 
range(lepto_feat$precip_tp_sum_week_anom, na.rm=TRUE)





















