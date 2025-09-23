# ===========================  SETUP  ==========================================
suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)     # gam/bam
})

# Use your merged table
DT <- as.data.table(lepto_feat)


# How many unique, finite values?
n_unique <- function(x) length(unique(x[is.finite(x)]))

# Decide how to treat a predictor given the data
# - If very few unique values -> linear term
# - Else -> smooth with k capped by n_unique-1
decide_term <- function(var, data, bs = "cr", k_target = 10L,
                        min_k_smooth = 3L, linear_if_le = 3L) {
  if (!var %in% names(data)) return(NULL)
  nu <- n_unique(data[[var]])
  if (nu <= linear_if_le) {
    # treat as linear (or factor if integer 0/1 and you prefer)
    return(var)
  }
  k_use <- max(min_k_smooth, min(k_target, nu - 1L))
  sprintf("s(%s, bs='%s', k=%d)", var, bs, k_use)
}
# Build a RHS given a set of candidate continuous vars (auto k),
# plus optional pure-linear vars (e.g., land cover).
build_rhs_auto <- function(cont_vars, data, bs = "cr", k_target = 10L,
                           linear_vars = character(0)) {
  # smooths with auto k
  s_terms <- unlist(lapply(cont_vars, decide_term, data = data,
                           bs = bs, k_target = k_target))
  # linear add-ons
  lin_terms <- intersect(linear_vars, names(data))
  rhs <- c(s_terms, lin_terms)
  rhs <- rhs[nzchar(rhs)]
  if (!length(rhs)) "1" else paste(rhs, collapse = " + ")
}

# Example usage on your feature table









# --- small safety: drop constants (no variance -> useless/infinite edf) ----
drop_constant_cols <- function(DT) {
  keep <- names(DT)[vapply(DT, function(x) length(unique(x)) > 1L, logical(1))]
  DT[, ..keep]
}
DT <- drop_constant_cols(DT)

# ---- housekeeping (unchanged but with safe offset) ----
stopifnot(all(c("district","date","poptot","lepto") %in% names(DT)))
if (!inherits(DT$date, "Date")) DT[, date := as.IDate(date)]
DT[, district := as.factor(district)]
DT[, week := as.integer(strftime(date, "%V"))]
DT[week < 1 | week > 53, week := pmin(53L, pmax(1L, week))]
DT[, year := as.integer(format(date, "%Y"))]
DT[, AR_start := (shift(district) != district) | (shift(week) != week-1), by = district]
DT[is.na(AR_start), AR_start := TRUE]

# ---- split ----
setorder(DT, date, district)
u <- sort(unique(DT$date))
cut_date <- u[floor(0.8*length(u))]
TR <- DT[date <= cut_date]
TE <- DT[date >  cut_date]

# ====================== CANDIDATE VARIABLE GROUPS =============================
keep_if_exist <- function(vars, D=TR) intersect(vars, names(D))

WX_CONTEMP <- keep_if_exist(c(
  "tmax_mean","tmin_mean","tmean_mean","rh_mean_week","vpd_mean_week",
  "ssrd_MJ_mean_week","precip_tp_sum_week","precip_mtpr_sum_week",
  "wet_days_ge10_tp","wet_days_ge10_mtpr","max3d_tp","max3d_mtpr",
  "wet_spell_maxlen_tp","wet_spell_maxlen_mtpr","hot_days_ge32",
  "ta_max","td_max","wbgt_mean_aw","rh_mean_aw","vpd_mean_aw"
))
WX_LAGS <- keep_if_exist(c(
  paste0("tmax_mean_lag", 1:3),
  paste0("rh_mean_week_lag", 1:3),
  paste0("vpd_mean_week_lag", 1:3),
  paste0("precip_tp_sum_week_lag", 1:3),
  paste0("precip_mtpr_sum_week_lag", 1:3)
))
WX_ROLLS <- keep_if_exist(c(
  "precip_tp_sum_week_roll2w_sum","precip_tp_sum_week_roll4w_sum",
  "vpd_mean_week_roll2w_mean","vpd_mean_week_roll4w_mean",
  "tmax_mean_roll2w_mean","tmax_mean_roll4w_mean",
  "rh_mean_week_roll2w_mean","rh_mean_week_roll4w_mean"
))
WX_ANOM <- keep_if_exist(c("precip_tp_sum_week_anom","tmax_mean_anom","rh_mean_week_anom"))
WX_CLIM <- keep_if_exist(c("precip_tp_sum_week_clim","tmax_mean_clim","rh_mean_week_clim"))
WX_PCT  <- keep_if_exist("precip_tp_sum_week_pct_normal")
WX_Z    <- keep_if_exist(c("precip_tp_sum_week_z","precip_mtpr_sum_week_z","tmax_mean_z","vpd_mean_week_z"))
WX_EWAP <- keep_if_exist(c("ewap_tp","ewap_mtpr"))

# AR lags (if not present, create 1:4)
AR_LAGS <- keep_if_exist(paste0("lepto_lag", 1:4))
if (!length(AR_LAGS)) {
  for (k in 1:4) TR[, paste0("lepto_lag", k) := shift(lepto, k), by=district]
  for (k in 1:4) TE[, paste0("lepto_lag", k) := shift(lepto, k), by=district]
  AR_LAGS <- paste0("lepto_lag", 1:4)
}

# ==================== ADAPTIVE TERM BUILDERS ==================================
# Use your decide_term()/build_rhs_auto defined earlier

SEASON <- "s(week, bs='cc', k=20)"
DISTRE <- "s(district, bs='re')"
OFFSET <- "offset(log(pmax(poptot, 1)))"   # <-- protect offset

rhs_auto <- function(vars, data, bs="cr", k=10, linear_vars = character(0)) {
  build_rhs_auto(cont_vars = vars, data = data, bs = bs, k_target = k, linear_vars = linear_vars)
}

# tensors (only if both vars present)
te_term <- function(v1, v2, bs=c("tp","tp"), k=NULL) {
  if (!all(c(v1,v2) %in% names(TR))) return(NULL)
  if (is.null(k)) sprintf("te(%s,%s, bs=c('%s','%s'))", v1, v2, bs[1], bs[2])
  else            sprintf("te(%s,%s, bs=c('%s','%s'), k=%d)", v1, v2, bs[1], bs[2], k)
}
TE_ta_rh     <- te_term("ta_max","rh_mean_aw", bs=c("tp","tp"), k=10)
TE_tmax_vpd  <- te_term("tmax_mean","vpd_mean_week", bs=c("tp","tp"), k=10)
TE_prec_tmax <- te_term("precip_tp_sum_week","tmax_mean", bs=c("tp","tp"), k=10)

# convenience
rhs_c <- function(...) {
  parts <- unlist(list(...))
  parts <- parts[nzchar(parts)]
  if (!length(parts)) "1" else paste(parts, collapse=" + ")
}

# ==================== DEFINE 18 FORMULA CONFIGS (adaptive) ====================
FMLS <- list(
  f1  = paste("lepto ~", rhs_c(OFFSET, rhs_auto(WX_CONTEMP, TR, "cr", 10), SEASON, DISTRE)),
  f2  = paste("lepto ~", rhs_c(OFFSET, rhs_auto(c(WX_CONTEMP, WX_LAGS), TR, "cr", 8), SEASON, DISTRE)),
  f3  = paste("lepto ~", rhs_c(OFFSET, rhs_auto(c(WX_CONTEMP, WX_ROLLS), TR, "tp", 8), SEASON, DISTRE)),
  f4  = paste("lepto ~", rhs_c(OFFSET, rhs_auto(WX_ANOM, TR, "tp", 8), SEASON, DISTRE)),
  f5  = paste("lepto ~", rhs_c(OFFSET, rhs_auto(c(WX_CLIM, WX_ANOM), TR, "tp", 8), SEASON, DISTRE)),
  f6  = paste("lepto ~", rhs_c(OFFSET, rhs_auto(c(WX_Z, WX_PCT), TR, "tp", 8), SEASON, DISTRE)),
  f7  = paste("lepto ~", rhs_c(OFFSET, rhs_auto(c("precip_tp_sum_week", WX_EWAP), TR, "tp", 10), SEASON, DISTRE)),
  f8  = paste("lepto ~", rhs_c(OFFSET, TE_ta_rh, SEASON, DISTRE)),
  f9  = paste("lepto ~", rhs_c(OFFSET, TE_tmax_vpd, rhs_auto(WX_LAGS, TR, "cr", 8), SEASON, DISTRE)),
  f10 = paste("lepto ~", rhs_c(OFFSET, TE_prec_tmax, rhs_auto(WX_ANOM, TR, "tp", 8), SEASON, DISTRE)),
  f11 = paste("lepto ~", rhs_c(OFFSET, rhs_auto(WX_CONTEMP, TR, "cr", 10), paste(AR_LAGS, collapse="+"), SEASON, DISTRE)),
  f12 = paste("lepto ~", rhs_c(OFFSET, rhs_auto(WX_LAGS, TR, "cr", 8), SEASON, DISTRE)),
  f13 = paste("lepto ~", rhs_c(OFFSET, rhs_auto(WX_ROLLS, TR, "tp", 8), SEASON, DISTRE)),
  f14 = paste("lepto ~", rhs_c(OFFSET, rhs_auto(c("precip_tp_sum_week", paste0("precip_tp_sum_week_lag",1:3)), TR, "cr", 8), SEASON, DISTRE)),
  f15 = paste("lepto ~", rhs_c(OFFSET, rhs_auto(c("ssrd_MJ_mean_week","hot_days_ge32","rh_mean_week"), TR, "tp", 5), SEASON, DISTRE)),
  f16 = paste("lepto ~", rhs_c(OFFSET, rhs_auto(c(WX_CONTEMP, WX_LAGS, WX_ROLLS, WX_ANOM, WX_EWAP), TR, "tp", 8), SEASON, DISTRE)),
  f17 = paste("lepto ~", rhs_c(OFFSET, rhs_auto(WX_CLIM, TR, "tp", 8), SEASON, DISTRE)),
  f18 = paste("lepto ~", rhs_c(OFFSET, rhs_auto(c("ta_max","td_max","wbgt_mean_aw","rh_mean_aw","vpd_mean_aw"), TR, "cr", 10),
                               TE_ta_rh, SEASON, DISTRE))
)

# Drop any empty RHS (just in case)
FMLS <- FMLS[vapply(FMLS, function(ff) grepl("~\\s*.+", ff), logical(1))]

# ======================= FIT / SCORE / COMPARE ================================
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
mae  <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
r2   <- function(y, yhat) {
  den <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(den) || den == 0) return(NA_real_)
  1 - sum((y - yhat)^2, na.rm = TRUE) / den
}

fit_one <- function(fml_chr) {
  fml <- as.formula(fml_chr)
  m <- try(bam(fml, data = TR, family = poisson(link="log"),
               method = "fREML", discrete = TRUE, select = TRUE,
               na.action = na.exclude),
           silent = TRUE)
  if (inherits(m, "try-error")) {
    return(data.table(formula = fml_chr, ok = FALSE,
                      RMSE_cnt = NA_real_, MAE_cnt = NA_real_, R2_cnt = NA_real_,
                      RMSE_rpk = NA_real_, MAE_rpk = NA_real_, R2_rpk = NA_real_,
                      AIC = NA_real_, GCV = NA_real_, DevExpl = NA_real_))
  }
  mu <- as.numeric(predict(m, newdata = TE, type = "response", discrete = TRUE,
                           na.action = na.exclude))
  data.table(
    formula = fml_chr, ok = TRUE,
    edf     = sum(m$edf),
    DevExpl = tryCatch(summary(m)$dev.expl, error = function(e) NA_real_),
    GCV     = tryCatch(m$gcv.ubre,        error = function(e) NA_real_),
    AIC     = tryCatch(AIC(m),            error = function(e) NA_real_),
    RMSE_cnt = rmse(TE$lepto, mu),
    MAE_cnt  = mae (TE$lepto, mu),
    R2_cnt   = r2  (TE$lepto, mu),
    RMSE_rpk = rmse(1e5*TE$lepto/pmax(TE$poptot,1), 1e5*mu/pmax(TE$poptot,1)),
    MAE_rpk  = mae (1e5*TE$lepto/pmax(TE$poptot,1), 1e5*mu/pmax(TE$poptot,1)),
    R2_rpk   = r2  (1e5*TE$lepto/pmax(TE$poptot,1), 1e5*mu/pmax(TE$poptot,1))
  )
}

SCORES <- rbindlist(lapply(FMLS, fit_one), fill = TRUE)
setorder(SCORES, RMSE_cnt)
SCORES[]






#   fml <- as.formula(fml_chr)
m <- try(bam(as.formula(SCORES$formula[1]), data = TR, family = poisson(link="log"),
             method = "fREML", discrete = TRUE, select = TRUE, offset = log(poptot),
             na.action = na.exclude),
         silent = TRUE)

plot(m)
summary(m)
        
        
        # agam = gam(lepto ~ offset(log(pmax(poptot, 1))) + s(tmax_mean, bs='cr', k=10) + s(tmin_mean, bs='cr', k=10) + s(tmean_mean, bs='cr', k=10) 
#            + s(rh_mean_week, bs='cr', k=10) + s(vpd_mean_week, bs='cr', k=10) + s(ssrd_MJ_mean_week, bs='cr', k=10) + s(precip_tp_sum_week, bs='cr', k=10) 
#            + s(precip_mtpr_sum_week, bs='cr', k=10) + s(wet_days_ge10_tp, bs='cr', k=7) + s(wet_days_ge10_mtpr, bs='cr', k=7) + s(max3d_tp, bs='cr', k=10) 
#            + s(max3d_mtpr, bs='cr', k=10) + wet_spell_maxlen_tp + wet_spell_maxlen_mtpr + s(hot_days_ge32, bs='cr', k=7) + s(ta_max, bs='cr', k=10) 
#            + s(td_max, bs='cr', k=10) + s(wbgt_mean_aw, bs='cr', k=10) + s(rh_mean_aw, bs='cr', k=10) + s(vpd_mean_aw, bs='cr', k=10) + 
#              lepto_lag1+lepto_lag2+lepto_lag3+lepto_lag4 + s(week, bs='cc', k=20) + s(district, bs='re'), data = lepto_feat, offset = log(poptot))
# 
# 


m2 <- bam(
  lepto ~ offset(log(pmax(poptot,1))) +
    s(tmax_mean, bs="cr", k=9) +
    s(rh_mean_week, bs="cr", k=9) +
    s(precip_tp_sum_week, bs="cr", k=9) +
    s(precip_tp_sum_week_lag2, bs="cr", k=7) +   # one plausible lag
    s(precip_tp_sum_week_lag3, bs="cr", k=7) +   # one plausible lag
    s(max3d_tp, bs="cr", k=7) +                  # one intensity metric
    s(week, bs="cc", k=6) +
    s(district, bs="re") +
    lepto_lag1 + lepto_lag2,
  data = TR, family = poisson(), method = "fREML", discrete = TRUE, select = TRUE
)
mu2 <- predict(m2, newdata = TE, type = "response", discrete = TRUE)

summary(m2)





plot_gam_smooths_base(m2, pages = c(2,3), show_rug = TRUE)






mA <- bam(
  lepto ~ offset(log(poptot)) +
    s(week, bs="cc", k=20) +
    s(district, bs="re") +
    precip_tp_sum_week + s(district, by = precip_tp_sum_week, bs="re") +
    s(tmax_mean, bs="cr", k=9) + s(rh_mean_week, bs="cr", k=9) +
    lepto_lag1 + lepto_lag2,
  data=TR, family=poisson(), method="fREML", discrete=TRUE, select=TRUE
)

mB <- bam(
  lepto ~ offset(log(poptot)) +
    te(week, district, bs=c("cc","re"), k=c(20, NA)) +
    s(precip_tp_sum_week, district, bs="fs", k=6) +     # FS interaction
    s(tmax_mean, bs="cr", k=9) + s(rh_mean_week, bs="cr", k=9) +
    lepto_lag1 + lepto_lag2,
  data=TR, family=poisson(), method="fREML", discrete=TRUE, select=TRUE
)

m2 <- bam(
  lepto ~ offset(log(poptot)) +
    # main weather
    s(tmax_mean, bs="cr", k=9) +
    s(rh_mean_week, bs="cr", k=9) +
    s(precip_tp_sum_week_clim, bs="cr", k=7) +
    s(precip_tp_sum_week_anom, bs="cr", k=9) +
    # interaction
    ti(tmax_mean, rh_mean_week, bs=c("cr","cr"), k=c(7,7)) +
    # short DLNM (lags as separate smooths)
    s(precip_tp_sum_week_lag1, bs="cr", k=7) +
    s(precip_tp_sum_week_lag2, bs="cr", k=7) +
    # season: global + district deviation
    s(week, bs="cc", k=20) +
    s(week, district, bs="fs", k=8) +
    # spatial pooling (pick one of these)
    s(district, bs="re"),
  # or: s(district, bs="mrf", xt=list(nb = nb_mgcv)),
  data = TR, family = poisson(), method="fREML", discrete=TRUE, select=TRUE
)




summary(m2)




summary(m_vc)





plot_gam_smooths_base(mA, pages = c(2,3), show_rug = TRUE)

plot_gam_smooths_base(m2, pages = c(2,3), show_rug = TRUE)

# 
# 
# 
# lepto$doy = yday(lepto$date)
# lepto$year = lepto$year.x
# lepto$district_factor = as.factor(lepto$district)
# 
# 
# agam = gam(lepto ~ s(ta_max, bs = 'cr', k =10)
#            + s(rh_mean_aw, bs = 'cr', k =10, fx = TRUE)
#            # + s(rh_mean_aw, by = district_factor, bs = 'cr', k =5, fx = FALSE)
#            + s(precip_tp_sum_week, k = 10, bs = 'cc') 
#            + s(td_max, bs = 'cr', k =10, fx = TRUE)
#            + s(year, bs = 'cc')
#            + s(doy, bs = 'cc')
#            + district
#            , data = lepto_feat, offset = log(poptot))
# 
# 
# summary(agam)
# 
# 
# 
# plot(agam)
# abline(h=0)
# 





# install.packages("gbm")



# ===================== GBM: Build, Fit, Predict, Visualize ====================

suppressPackageStartupMessages({
  library(data.table)
  library(gbm)
})

# ---- 0) Start from TR/TE already split and cleaned ----
# Assumes you already prepared TR, TE with columns:
#   district (factor), date (Date), lepto (count), poptot (offset)
#   + all engineered features (lags, rolls, anomalies, sin52/cos52, week, etc.)

suppressPackageStartupMessages({
  library(data.table)
  library(gbm)
})

# ---- Assume TR/TE already exist and contain: lepto, poptot, district (factor), date, and features ----

# 1) Choose feature set (no targets/offsets here)
FEATS0 <- c(
  # contemporaneous
  "tmax_mean","tmin_mean","tmean_mean","rh_mean_week","vpd_mean_week",
  "ssrd_MJ_mean_week","precip_tp_sum_week","precip_mtpr_sum_week",
  "wet_days_ge10_tp","wet_days_ge10_mtpr","max3d_tp","max3d_mtpr","hot_days_ge32",
  # lags/rolls
  paste0("tmax_mean_lag",1:3),
  paste0("rh_mean_week_lag",1:3),
  paste0("vpd_mean_week_lag",1:3),
  paste0("precip_tp_sum_week_lag",1:3),
  "precip_tp_sum_week_roll2w_sum","tmax_mean_roll2w_mean","rh_mean_week_roll2w_mean",
  # anomalies / ewap / z
  "tmax_mean_anom","rh_mean_week_anom","precip_tp_sum_week_anom","ewap_tp",
  "tmax_mean_z","vpd_mean_week_z","precip_tp_sum_week_z",
  # AR lags
  paste0("lepto_lag",1:4),
  # seasonality
  "week","sin52","cos52"
)
FEATS0 <- intersect(FEATS0, names(TR))
FEATS0 <- FEATS0[vapply(TR[, ..FEATS0], is.numeric, logical(1))]

# 2) District one-hot, aligned
TR[, district := factor(district)]
TE[, district := factor(district, levels = levels(TR$district))]
MM_TR <- model.matrix(~ district - 1, data = TR)
MM_TE <- model.matrix(~ district - 1, data = TE)
miss_cols <- setdiff(colnames(MM_TR), colnames(MM_TE))
if (length(miss_cols)) {
  add <- matrix(0, nrow = nrow(MM_TE), ncol = length(miss_cols))
  colnames(add) <- miss_cols
  MM_TE <- cbind(MM_TE, add)
}
MM_TE <- MM_TE[, colnames(MM_TR), drop = FALSE]

# 3) Final numeric design matrices
XTR <- cbind(as.matrix(TR[, ..FEATS0]), MM_TR)
XTE <- cbind(as.matrix(TE[, ..FEATS0]), MM_TE)

# Keep column names stable (no mangling)
colnames(XTR) <- make.names(colnames(XTR), unique = TRUE)
colnames(XTE) <- colnames(XTR) # must match

FEATS <- colnames(XTR)

# 4) Targets & exposure
y_rate <- TR$lepto / pmax(TR$poptot, 1)
w_expo <- pmax(TR$poptot, 1)

# 5) Row sanitation (AFTER the above builds)
is_finite_row <- function(M) apply(M, 1L, function(r) all(is.finite(r)))
mask <- is.finite(y_rate) & (y_rate >= 0) &
  is.finite(w_expo) & (w_expo >= 0) &
  is_finite_row(XTR)

if (!any(mask)) stop("No valid rows after cleaning train.")
if (any(!mask)) message("Dropping ", sum(!mask), " invalid train rows.")
XTRc   <- XTR[mask, , drop = FALSE]
y_ratec <- as.numeric(y_rate[mask])
w_expoc <- as.numeric(w_expo[mask])

# Zero-variance feature drop (optional but helpful)
nzv <- apply(XTRc, 2L, function(x) length(unique(x)) > 1L)
if (any(!nzv)) {
  drop_cols <- FEATS[!nzv]
  message("Dropping constant features: ", paste(drop_cols, collapse = ", "))
  XTRc <- XTRc[, nzv, drop = FALSE]
  XTE  <- XTE[,  nzv, drop = FALSE]
  FEATS <- colnames(XTRc)
}

# Final belt-and-suspenders: clamp any tiny negatives (shouldn't exist, but be safe)
y_ratec[!is.finite(y_ratec)] <- 0
y_ratec[y_ratec < 0] <- 0

# 6) Fit with gbm.fit (avoid formula/model.frame pitfalls)
set.seed(42)
gbm_fit <- gbm.fit(
  x                  = XTRc,
  y                  = y_ratec,
  distribution       = "gaussian",     # Poisson on non-negative rate
  w                  = w_expoc,       # exposure as weights
  n.trees            = 100,
  interaction.depth  = 5,
  shrinkage          = 0.03,
  n.minobsinnode     = 20,
  bag.fraction       = 0.8,
  train.fraction     = 1.0,           # use CV instead of internal holdout
  cv.folds           = 5,
  keep.data          = FALSE,
  verbose            = TRUE,
  n.cores            = parallel::detectCores()  # << works in gbm.fit
)



best_iter <- gbm.perf(gbm_fit, method = "cv", plot.it = FALSE)

# 7) Predict on held-out test
mu_rate_te <- predict(gbm_fit, newdata = XTE, n.trees = best_iter, type = "response")
mu_cnt_te  <- mu_rate_te * pmax(TE$poptot, 1)

# 8) Simple scores
rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))
mae  <- function(a, b) mean(abs(a - b), na.rm = TRUE)
r2   <- function(y, yhat) 1 - sum((y - yhat)^2, na.rm = TRUE) / sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)

test_cnt <- TE$lepto
cat(sprintf("GBM @ %d trees | RMSE_cnt=%.3f  MAE_cnt=%.3f  R2_cnt=%.3f\n",
            best_iter, rmse(test_cnt, mu_cnt_te), mae(test_cnt, mu_cnt_te), r2(test_cnt, mu_cnt_te)))

# 9) Variable importance & partial effects (base graphics)
# Importance
rel_inf <- summary(gbm_fit, n.trees = best_iter, plotit = FALSE)
print(head(rel_inf, 20))

# Partial dependence for top variables
op <- par(mfrow = c(2,2), mar = c(4,4,2,1))
on.exit(par(op), add = TRUE)

top_vars <- rel_inf$var[seq_len(min(4, nrow(rel_inf)))]
for (v in top_vars) {
  # gbm::plot will handle transforming the single feature grid internally
  plot(gbm_fit, i.var = v, n.trees = best_iter, main = paste("Partial:", v))
}

# ===================== NOTES =============================
# 1) We used one-hot districts-safer than integer codes.
# 2) The offset is passed as base_margin (XGB) and as 'offset' (GBM).
# 3) Split is strictly time-based to avoid leakage.
# 4) We impute NA???0 for tree models (simple & robust).
# 5) Expand FEATS with your roll4w, pct_normal, z-scores, etc., as desired.





suppressPackageStartupMessages({
  library(data.table)
  library(gbm)
})

## -------- 0) Inputs assumed ----------
## TR, TE already exist (your time split); contain:
##   lepto (counts), poptot (exposure), district (factor), date (Date),
##   engineered weather features (lags/rolls/anoms), week, sin52, cos52, etc.

stopifnot(all(c("district","date","lepto","poptot") %in% names(TR)))
stopifnot(all(c("district","date","lepto","poptot") %in% names(TE)))

## -------- 1) Choose features (no target/offset here) ----------
FEATS0 <- c(
  # contemporaneous
  "tmax_mean","tmin_mean","tmean_mean","rh_mean_week","vpd_mean_week",
  "ssrd_MJ_mean_week","precip_tp_sum_week","precip_mtpr_sum_week",
  "wet_days_ge10_tp","wet_days_ge10_mtpr","max3d_tp","max3d_mtpr","hot_days_ge32",
  # lags/rolls
  paste0("tmax_mean_lag",1:3),
  paste0("rh_mean_week_lag",1:3),
  paste0("vpd_mean_week_lag",1:3),
  paste0("precip_tp_sum_week_lag",1:3),
  "precip_tp_sum_week_roll2w_sum","tmax_mean_roll2w_mean","rh_mean_week_roll2w_mean",
  # anomalies / ewap / z
  "tmax_mean_anom","rh_mean_week_anom","precip_tp_sum_week_anom","ewap_tp",
  "tmax_mean_z","vpd_mean_week_z","precip_tp_sum_week_z",
  # AR lags for outcome
  paste0("lepto_lag",1:4),
  # seasonality
  "week","sin52","cos52"
)

FEATS0 <- intersect(FEATS0, names(TR))
FEATS0 <- FEATS0[vapply(TR[, ..FEATS0], is.numeric, logical(1))]

## -------- 2) One-hot district and align columns ----------
TR[, district := factor(district)]
TE[, district := factor(district, levels = levels(TR$district))]

MM_TR <- model.matrix(~ district - 1, data = TR)
MM_TE <- model.matrix(~ district - 1, data = TE)

miss_cols <- setdiff(colnames(MM_TR), colnames(MM_TE))
if (length(miss_cols)) {
  add <- matrix(0, nrow = nrow(MM_TE), ncol = length(miss_cols))
  colnames(add) <- miss_cols
  MM_TE <- cbind(MM_TE, add)
}
MM_TE <- MM_TE[, colnames(MM_TR), drop = FALSE]

## -------- 3) Build design matrices ----------
XTR <- cbind(as.matrix(TR[, ..FEATS0]), MM_TR)
XTE <- cbind(as.matrix(TE[, ..FEATS0]), MM_TE)

FEATS <- colnames(XTR) # final list

## -------- 4) Targets, weights, and hard sanitation ----------
y_rate <- TR$lepto / pmax(TR$poptot, 1)      # non-negative
w_expo <- pmax(TR$poptot, 1)                  # non-negative

# Row-wise finite mask on everything we send to the learner
is_finite_row <- function(M) apply(as.matrix(M), 1L, function(r) all(is.finite(r)))
mask <- is.finite(y_rate) & (y_rate >= 0) &
  is.finite(w_expo) & (w_expo >= 0) &
  is_finite_row(XTR)

if (!any(mask)) stop("No valid training rows after sanitation.")
if (any(!mask)) message("Dropping ", sum(!mask), " invalid training rows.")

XTRc    <- XTR[mask, , drop = FALSE]
y_ratec <- pmax(y_rate[mask], 0)             # floor any machine eps < 0
w_expoc <- pmax(w_expo[mask], 0)

# Drop zero-variance features (after masking)
nzv <- vapply(as.data.frame(XTRc), function(x) length(unique(x)) > 1L, logical(1))
if (any(!nzv)) {
  message("Dropping constant features: ", paste(colnames(XTRc)[!nzv], collapse = ", "))
  XTRc <- XTRc[, nzv, drop = FALSE]
  FEATS <- colnames(XTRc)
}

## -------- 5) External (time-blocked) validation split ----------
# Use the *kept* rows' dates to split 85/15 by time
dates_trc <- TR$date[mask]
ord <- order(dates_trc)
n   <- length(dates_trc)
cut <- floor(0.85 * n)
idx_tr <- ord[1:cut]
idx_vl <- ord[(cut+1):n]

X_tr <- XTRc[idx_tr, , drop = FALSE]
y_tr <- y_ratec[idx_tr]
w_tr <- w_expoc[idx_tr]

X_vl <- XTRc[idx_vl, , drop = FALSE]
y_vl <- y_ratec[idx_vl]
w_vl <- w_expoc[idx_vl]

# Belt-and-suspenders prints (so if it ever fails again, you'll see it)
cat("y_tr range:", range(y_tr), "  any <0? ", any(y_tr < 0), "\n")
cat("y_vl range:", range(y_vl), "  any <0? ", any(y_vl < 0), "\n")

## -------- 6) Fit GBM via gbm.fit() (no formula, no hidden CV) ----------
set.seed(42)
gbm_fit <- gbm.fit(
  x                 = X_tr,
  y                 = y_tr,
  distribution      = "gaussian",   # Poisson on rate with exposure in weights
  w                 = w_tr,
  n.trees           = 1000,
  interaction.depth = 5,
  shrinkage         = 0.03,
  n.minobsinnode    = 20,
  bag.fraction      = 0.8,
  nTrain            = nrow(X_tr),  # replaces train.fraction
  keep.data         = FALSE,
  verbose           = TRUE
)

## -------- 7) Pick best n.trees on the validation block ----------
poisson_dev <- function(obs_rate, pred_rate, w) {
  eps <- 1e-12
  mu  <- pmax(pred_rate, eps)
  y   <- pmax(obs_rate, 0)
  # mean weighted Poisson deviance
  term <- ifelse(y > 0, y*log(y/mu) - (y - mu), mu)
  2 * sum(w * term, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

grid_trees <- seq(100, 6000, by = 100)
vl_dev <- sapply(grid_trees, function(nt) {
  pr <- predict(gbm_fit, X_vl, n.trees = nt, type = "response")  # rate
  poisson_dev(y_vl, pr, w_vl)
})
best_nt <- grid_trees[which.min(vl_dev)]
cat("Chosen n.trees (external val):", best_nt, "\n")

## -------- 8) Final test predictions (held-out TE) ----------
# Align TE feature columns to the kept FEATS
XTE_use <- XTE[, FEATS, drop = FALSE]
# Any missing columns (shouldn't happen, but safe):
miss_te <- setdiff(FEATS, colnames(XTE_use))
if (length(miss_te)) {
  add <- matrix(0, nrow = nrow(XTE_use), ncol = length(miss_te))
  colnames(add) <- miss_te
  XTE_use <- cbind(XTE_use, add)[, FEATS, drop = FALSE]
}

rate_hat_te <- predict(gbm_fit, XTE_use, n.trees = best_nt, type = "response")
count_hat_te <- pmax(TE$poptot, 1) * rate_hat_te

rmse <- function(a,b) sqrt(mean((a-b)^2, na.rm = TRUE))
mae  <- function(a,b) mean(abs(a-b), na.rm = TRUE)
r2   <- function(y,yhat){
  den <- sum((y-mean(y, na.rm=TRUE))^2, na.rm=TRUE)
  if (!is.finite(den) || den==0) return(NA_real_)
  1 - sum((y-yhat)^2, na.rm=TRUE)/den
}

cat(
  "Test RMSE (counts): ", rmse(TE$lepto, count_hat_te), "\n",
  "Test MAE  (counts): ", mae (TE$lepto, count_hat_te), "\n",
  "Test R2   (counts): ", r2  (TE$lepto, count_hat_te), "\n"
)

## -------- 9) Quick variable importance & partial dependence ----------
rel_inf <- summary.gbm(gbm_fit, n.trees = best_nt, cBars = 0, plotit = FALSE)
print(head(rel_inf, 20))  # top 20 features

# Base-R partial dependence for a few core variables
pd_plot <- function(fit, X, var, nt = best_nt) {
  if (!var %in% colnames(X)) { message("Missing: ", var); return(invisible()) }
  pd <- plot.gbm(fit, i.var = var, n.trees = nt, return.grid = TRUE)
  plot(pd[[var]], pd$y, type = "l", xlab = var, ylab = "f(var)", main = paste("PD:", var))
  rug(X[, var])
}

par(mfrow = c(2,3))
for (v in intersect(c("tmax_mean","rh_mean_week","precip_tp_sum_week","vpd_mean_week","week","sin52"), FEATS)) {
  pd_plot(gbm_fit, XTRc, v)
}
par(mfrow = c(1,1))





























suppressPackageStartupMessages({
  library(data.table)
  library(gbm)
})

## -------- Assumptions --------
## You already have TR, TE (time split) with:
##   district (factor), date (Date), lepto (counts >=0), poptot (>=0),
##   engineered features (lags/rolls/anoms), week, sin52, cos52, etc.

stopifnot(all(c("district","date","lepto","poptot") %in% names(TR)))
stopifnot(all(c("district","date","lepto","poptot") %in% names(TE)))
if (!inherits(TR$date, "Date")) TR[, date := as.IDate(date)]
if (!inherits(TE$date, "Date")) TE[, date := as.IDate(date)]

## -------- 1) Feature set (no target/offset here) --------
FEATS0 <- c(
  "tmax_mean","tmin_mean","tmean_mean","rh_mean_week","vpd_mean_week",
  "ssrd_MJ_mean_week","precip_tp_sum_week","precip_mtpr_sum_week",
  "wet_days_ge10_tp","wet_days_ge10_mtpr","max3d_tp","max3d_mtpr","hot_days_ge32",
  paste0("tmax_mean_lag",1:3),
  paste0("rh_mean_week_lag",1:3),
  paste0("vpd_mean_week_lag",1:3),
  paste0("precip_tp_sum_week_lag",1:3),
  "precip_tp_sum_week_roll2w_sum","tmax_mean_roll2w_mean","rh_mean_week_roll2w_mean",
  "tmax_mean_anom","rh_mean_week_anom","precip_tp_sum_week_anom","ewap_tp",
  "tmax_mean_z","vpd_mean_week_z","precip_tp_sum_week_z",
  paste0("lepto_lag",1:4),
  "week","sin52","cos52"
)
FEATS0 <- intersect(FEATS0, names(TR))
FEATS0 <- FEATS0[vapply(TR[, ..FEATS0], is.numeric, logical(1))]

# FEATS0 = FEATS0[!FEATS0 %like% "lepto"]

## -------- 2) One-hot for district; align TR/TE --------
TR[, district := factor(district)]
TE[, district := factor(district, levels = levels(TR$district))]

MM_TR <- model.matrix(~ district - 1, data = TR)
MM_TE <- model.matrix(~ district - 1, data = TE)

if (length(setdiff(colnames(MM_TR), colnames(MM_TE)))) {
  miss <- setdiff(colnames(MM_TR), colnames(MM_TE))
  MM_TE <- cbind(MM_TE, matrix(0, nrow(MM_TE), length(miss), dimnames = list(NULL, miss)))
}
MM_TE <- MM_TE[, colnames(MM_TR), drop = FALSE]

## -------- 3) Final design matrices --------
XTR <- cbind(as.matrix(TR[, ..FEATS0]), MM_TR)
XTE <- cbind(as.matrix(TE[, ..FEATS0]), MM_TE)
FEATS <- colnames(XTR)

## -------- 4) Targets (COUNTS) and OFFSETS (log pop) --------
y_cnt   <- pmax(TR$lepto, 0)
off_tr  <- log(pmax(TR$poptot, 1))
off_te  <- log(pmax(TE$poptot, 1))

## Row-wise finite mask on X, y, offset
is_finite_row <- function(M) apply(as.matrix(M), 1L, function(r) all(is.finite(r)))
mask <- is_finite_row(XTR) & is.finite(y_cnt) & (y_cnt >= 0) & is.finite(off_tr)

if (!any(mask)) stop("No valid training rows after sanitation.")
if (any(!mask)) message("Dropping ", sum(!mask), " invalid training rows.")

XTRc   <- XTR[mask, , drop = FALSE]
y_cntc <- y_cnt[mask]
off_tc <- off_tr[mask]

## Drop zero-variance features (after masking)
nzv <- vapply(as.data.frame(XTRc), function(x) length(unique(x)) > 1L, logical(1))
if (any(!nzv)) {
  message("Dropping constant features: ", paste(colnames(XTRc)[!nzv], collapse = ", "))
  XTRc <- XTRc[, nzv, drop = FALSE]
  FEATS <- colnames(XTRc)
}

## -------- 5) External time-blocked validation (85/15) --------
dates_trc <- TR$date[mask]
ord <- order(dates_trc)
n   <- length(dates_trc)
cut <- floor(0.85 * n)
idx_tr <- ord[1:cut]
idx_vl <- ord[(cut+1):n]

X_tr  <- XTRc[idx_tr, , drop = FALSE]
y_tr  <- y_cntc[idx_tr]
off_tr_sub <- off_tc[idx_tr]

X_vl  <- XTRc[idx_vl, , drop = FALSE]
y_vl  <- y_cntc[idx_vl]
off_vl_sub <- off_tc[idx_vl]

cat("y_tr (counts) range:", range(y_tr), "  any <0? ", any(y_tr < 0), "\n")
cat("offset_tr range:", range(off_tr_sub), "\n")


summary(TE$lepto)
quantile(TE$lepto, c(.1,.25,.5,.75,.9))

# naive persistence baseline
naive_mu <- pmax(0, TE$lepto_lag1)  # or clamp at 0
cat(
  "Naive RMSE:", sqrt(mean((TE$lepto - naive_mu)^2, na.rm=TRUE)), "\n",
  "Naive MAE :", mean(abs(TE$lepto - naive_mu), na.rm=TRUE), "\n"
)



num_of_tress = 6000

## -------- 6) Fit GBM (Poisson with log offset) --------
## IMPORTANT: use counts + offset; do NOT put exposure in weights.
set.seed(42)
gbm_fit <- gbm.fit(
  x                 = X_tr,
  y                 = y_tr,
  distribution      = "poisson",
  offset            = off_tr_sub,        # <<< exposure enters HERE
  n.trees           = num_of_tress,
  interaction.depth = 5,
  shrinkage         = 0.6,
  n.minobsinnode    = 20,
  bag.fraction      = 0.8,
  nTrain            = nrow(X_tr),
  keep.data         = FALSE,
  verbose           = TRUE
)

## -------- 7) Pick n.trees via validation deviance --------
## For gbm.fit, predict(type="link") gives f(x); mean = exp( offset + f(x) )
poisson_dev <- function(y, mu) {
  ## y: counts (>=0); mu: mean (>0)
  eps <- 1e-12
  mu  <- pmax(mu, eps)
  y   <- pmax(y, 0)
  2 * sum(ifelse(y > 0, y*log(y/mu) - (y - mu), mu), na.rm = TRUE) / length(y)
}

grid_trees <- seq(100, num_of_tress, by = 100)
vl_dev <- sapply(grid_trees, function(nt) {
  fx   <- predict(gbm_fit, X_vl, n.trees = nt, type = "link")
  mu   <- exp(off_vl_sub + fx)   # add offset on link scale
  poisson_dev(y_vl, mu)
})
best_nt <- grid_trees[which.min(vl_dev)]
cat("Chosen n.trees:", best_nt, "\n")

## -------- 8) Final test scoring --------
## Align TE columns to FEATS kept in train:
XTE_use <- XTE[, FEATS, drop = FALSE]
miss <- setdiff(FEATS, colnames(XTE_use))
if (length(miss)) {
  add <- matrix(0, nrow = nrow(XTE_use), ncol = length(miss))
  colnames(add) <- miss
  XTE_use <- cbind(XTE_use, add)[, FEATS, drop = FALSE]
}

fx_te    <- predict(gbm_fit, XTE_use, n.trees = best_nt, type = "link")
mu_count <- exp(off_te + fx_te)           # expected counts

rmse <- function(a,b) sqrt(mean((a-b)^2, na.rm = TRUE))
mae  <- function(a,b) mean(abs(a-b), na.rm = TRUE)
r2   <- function(y,yhat){
  den <- sum((y-mean(y, na.rm=TRUE))^2, na.rm=TRUE)
  if (!is.finite(den) || den==0) return(NA_real_)
  1 - sum((y-yhat)^2, na.rm=TRUE)/den
}

cat(
  "Test RMSE (counts): ", rmse(TE$lepto, mu_count), "\n",
  "Test MAE  (counts): ", mae (TE$lepto, mu_count), "\n",
  "Test R2   (counts): ", r2  (TE$lepto, mu_count), "\n"
)

## -------- 9) Variable importance & partial dependence (base R) --------
rel_inf <- summary.gbm(gbm_fit, n.trees = best_nt, cBars = 0, plotit = FALSE)
print(head(rel_inf, 20))

pd_plot <- function(fit, X, var, nt = best_nt, offset_val = median(off_tr_sub, na.rm=TRUE)) {
  if (!var %in% colnames(X)) return(invisible())
  ## For gbm.fit with matrices, we'll do a simple grid along var, hold others at their medians,
  ## and add a constant offset so the PD is on the response (count) scale.
  x0 <- as.data.frame(X)
  z  <- x0
  z[] <- lapply(z, function(col) if (is.numeric(col)) median(col, na.rm=TRUE) else col)
  grid <- quantile(x0[[var]], probs = seq(0.02, 0.98, length.out = 60), na.rm = TRUE, names = FALSE)
  z_rep <- z[rep(1, length(grid)), , drop = FALSE]
  z_rep[[var]] <- grid
  fx <- predict(fit, as.matrix(z_rep), n.trees = nt, type = "link")
  mu <- exp(offset_val + fx)
  plot(grid, mu, type = "l", xlab = var, ylab = "partial mean count", main = paste("PD:", var))
  rug(x0[[var]])
}

par(mfrow=c(2,3))
for (v in intersect(c("tmax_mean","rh_mean_week","precip_tp_sum_week","vpd_mean_week","week","sin52"), FEATS)) {
  pd_plot(gbm_fit, XTRc, v)
}
par(mfrow=c(1,1))

