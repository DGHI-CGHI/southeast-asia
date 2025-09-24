# ================================
# Climate feature screening script
# ================================
library(data.table)
library(mgcv)
library(PRROC)

# ---------- Prep & splits ----------
DT <- as.data.table(copy(model_df))
setkey(DT, district, date_start)
DT[, district := factor(district)]
if (!"week" %in% names(DT)) DT[, week := as.integer(strftime(date_start, "%V"))]
DT[week == 53L, week := 52L]
DT <- DT[poptot > 0]

# names(dt)[names(dt) %likeany% c("ta_","rh_","tp_","mtpr_","wbgt_","wind_","ssrd")]

base_terms <- c("ta_mean_lag1","rh_mean_lag1","tp_sum_lag2")  # sensible defaults
need_base  <- c("lepto","poptot","week","district","date_start")

ALL0 <- DT[complete.cases(DT[, ..need_base])]
setorder(ALL0, date_start)

cut1 <- ALL0$date_start[floor(0.8 * nrow(ALL0))]
train_all <- ALL0[date_start <= cut1]
test      <- ALL0[date_start  >  cut1]

cut2 <- train_all$date_start[floor(0.8 * nrow(train_all))]
train_fit <- train_all[date_start <= cut2]
val       <- train_all[date_start  >  cut2]

# Safe cyclic knots for 'cc'
# k_week  <- min(18L, max(4L, uniqueN(train_fit$week) - 1L))
k_week = 6
cc_knots <- list(week = c(0.5, 52.5))

# ---------- Baselines for spike labelling ----------
week52 <- function(d) { w <- as.integer(strftime(d, "%V")); ifelse(w == 53L, 52L, w) }

make_baseline_clim <- function(train_dat, target_dat, type = c("mad","q"), h = 1.0, q = 0.8) {
  type <- match.arg(type)
  T1 <- copy(train_dat)[, .(district, date_start, lepto)]
  T1[, woy := week52(date_start)]
  CLIM <- T1[, {
    x <- as.numeric(lepto)
    if (type == "mad") {
      med <- median(x, na.rm = TRUE); md <- mad(x, na.rm = TRUE)
      .(baseline = pmax(0, med + h * md))
    } else {
      .(baseline = pmax(0, as.numeric(quantile(x, q, na.rm = TRUE))))
    }
  }, by = .(district, woy)]
  tgt <- copy(target_dat)[, .(district, date_start)]
  tgt[, woy := week52(date_start)]
  CLIM[tgt, on = .(district, woy)][, .(district, date_start, baseline)]
}

# Climatology baselines (train-derived) for VAL and TEST
BL_val  <- make_baseline_clim(train_fit, val,  type = "mad", h = 1.0)
BL_test <- make_baseline_clim(train_all, test, type = "mad", h = 1.0)

# ---------- Base model (no climate) ----------
fit_base <- function(train_dat) {
  bam(
    lepto ~ s(week, bs="cc", k = k_week) + s(district, bs="re") + offset(log(poptot)),
    data     = train_dat,
    family   = nb(link="log"),
    method   = "fREML",
    select   = TRUE, gamma = 1.4,
    discrete = TRUE, nthreads = 8,
    knots    = cc_knots
  )
}
m_base <- fit_base(train_fit)

# ---------- Scoring helpers ----------
score_counts <- function(model, newdata) {
  mu <- as.numeric(predict(model, newdata = newdata, type = "response"))
  theta <- model$family$getTheta(TRUE)
  list(mu = mu, theta = theta)
}

# Safe PR-AUC (handles degenerate cases)
safe_pr_auc <- function(scores_pos, scores_neg) {
  if (length(scores_pos) == 0L || length(scores_neg) == 0L ||
      all(!is.finite(scores_pos)) || all(!is.finite(scores_neg))) return(NA_real_)
  pr <- PRROC::pr.curve(scores.class0 = scores_pos, scores.class1 = scores_neg)
  pr$auc.davis.goadrich
}

# Compute alert PR-AUC and NB log score on a set with a given baseline table
eval_set <- function(model, set_dt, BL) {
  X <- BL[copy(set_dt), on = .(district, date_start)]
  sc <- score_counts(model, X)
  mu <- sc$mu; th <- sc$theta
  # alert score
  base_int <- pmax(0L, floor(pmax(0, X$baseline)))
  spike_prob <- 1 - pnbinom(pmax(base_int - 1L, 0L), size = th, mu = mu)
  obs_spike  <- as.integer(X$lepto >= X$baseline)
  pr_auc <- safe_pr_auc(spike_prob[obs_spike == 1L], spike_prob[obs_spike == 0L])
  # count log score (negative log-likelihood; lower is better)
  nll <- -sum(dnbinom(X$lepto, size = th, mu = mu, log = TRUE), na.rm = TRUE)
  data.table(pr_auc = pr_auc, mean_nll = nll / nrow(X))
}

# Base scores (for deltas)
val_base  <- eval_set(m_base,  val,  BL_val)
test_base <- eval_set(m_base,  test, BL_test)

# ---------- Candidate variable discovery ----------
# Keep *_lag1..lag4 and *_roll2/4/8, drop district_* and obviously non-climate
all_names <- names(DT)
is_candidate <- grepl("_(lag[1-4]|roll[248])$", all_names) &
  !grepl("^district_", all_names) &
  !all_names %in% c("district_lag1","district_lag2","district_lag3","district_lag4")
cands <- setdiff(all_names[is_candidate], c("week","district","date_start","lepto","poptot"))



# cands = jj::jsample(data.table(cands), n = 10, random = TRUE)$cands

# Filter to numeric and not (nearly) constant in training
is_ok <- function(v) {
  x <- train_fit[[v]]
  is.numeric(x) && sum(is.finite(x)) >= 0.7 * nrow(train_fit) && sd(x, na.rm = TRUE) > 1e-8
}
cands <- Filter(is_ok, cands)

# ---------- Univariate screening (base + s(feature)) ----------
fit_with_feature <- function(train_dat, feat) {
  fm <- as.formula(paste0(
    "lepto ~ s(week, bs='cc', k=", k_week, ") + s(district, bs='re') + ",
    "s(", feat, ", bs='cr', k=20) + offset(log(poptot))"
  ))
  bam(fm, data = train_dat, family = nb(link="log"),
      method="fREML", select=TRUE, gamma=1.4,
      discrete=TRUE, nthreads=8, knots=cc_knots)
}

rank_rows <- vector("list", length(cands))
for (i in seq_along(cands)) {
  v <- cands[[i]]
  need_v <- c(need_base, v)
  trv <- train_fit[complete.cases(train_fit[, ..need_v])]
  valv <- val[complete.cases(val[, ..need_v])]
  testv <- test[complete.cases(test[, ..need_v])]
  # skip if too small after NA filtering
  if (nrow(trv) < 0.5 * nrow(train_fit) || nrow(valv) < 200) next
  m <- fit_with_feature(trv, v)
  val_sc  <- eval_set(m, valv,  BL_val)
  test_sc <- eval_set(m, testv, BL_test)
  # edf of the feature (how much wiggle it used)
  sm <- summary(m)$s.table
  row_idx <- grep(paste0("^s\\(", v, "\\)"), rownames(sm))
  edf_v <- if (length(row_idx)) sm[row_idx, "edf"] else NA_real_
  rank_rows[[i]] <- data.table(
    feature   = v,
    val_pr_auc   = val_sc$pr_auc,
    test_pr_auc  = test_sc$pr_auc,
    d_val_pr_auc = val_sc$pr_auc - val_base$pr_auc,
    val_mean_nll   = val_sc$mean_nll,
    test_mean_nll  = test_sc$mean_nll,
    d_val_mean_nll = val_sc$mean_nll - val_base$mean_nll,
    edf = edf_v
  )
  cat(i, '\n')
}
feat_rank <- rbindlist(rank_rows, use.names = TRUE, fill = TRUE)
setorder(feat_rank, -d_val_pr_auc, d_val_mean_nll)  # higher PR-AUC gain, lower NLL is better

# ---------- Optional: simple forward selection ----------
# Greedy add features while PR-AUC gain on VAL > min_gain and avoiding high correlation
min_gain <- 0.01   # require at least +0.01 PR-AUC on validation to add
max_k    <- 5      # cap number of features
picked <- character(0)
current_model <- m_base
current_need  <- need_base
current_val_pr <- val_base$pr_auc

# precompute Spearman correlation among candidates to avoid redundancy
corr_cap <- 0.8
cand_mat <- train_fit[, ..cands]
cand_cor <- suppressWarnings(cor(cand_mat, use="pairwise.complete.obs", method="spearman"))

forward_rows <- list()
repeat {
  best <- NULL
  best_gain <- 0
  for (v in cands) {
    if (v %in% picked) next
    # avoid high correlation with already picked
    if (length(picked) && any(abs(cand_cor[v, picked]) > corr_cap, na.rm=TRUE)) next
    need_v <- c(current_need, v)
    trv <- train_fit[complete.cases(train_fit[, ..need_v])]
    valv <- val[complete.cases(val[, ..need_v])]
    if (nrow(trv) < 0.5 * nrow(train_fit) || nrow(valv) < 200) next
    # fit base plus all picked plus candidate
    rhs <- c(
      sprintf("s(week, bs='cc', k=%d)", k_week),
      "s(district, bs='re')",
      sprintf("s(%s, bs='cr', k=6)", c(picked, v)),
      "offset(log(poptot))"
    )
    fm <- as.formula(paste("lepto ~", paste(rhs, collapse=" + ")))
    m_try <- bam(fm, data=trv, family=nb(link="log"),
                 method="fREML", select=TRUE, gamma=1.4,
                 discrete=TRUE, nthreads=8, knots=cc_knots)
    val_sc <- eval_set(m_try, valv, BL_val)
    gain <- val_sc$pr_auc - current_val_pr
    if (is.finite(gain) && gain > best_gain) {
      best_gain <- gain; best <- list(v=v, m=m_try, val_sc=val_sc)
    }
    cat(v, '\n')
  }
  if (is.null(best) || best_gain < min_gain || length(picked) >= max_k) break
  picked <- c(picked, best$v)
  current_model <- best$m
  current_val_pr <- best$val_sc$pr_auc
  forward_rows[[length(picked)]] <- data.table(
    step = length(picked), added = best$v,
    val_pr_auc = best$val_sc$pr_auc,
    d_val_pr_auc = best_gain
  )
}
fwd_summary <- rbindlist(forward_rows, use.names=TRUE, fill=TRUE)

# ---------- Outputs to inspect ----------
cat("\nBase model (no climate) VAL PR-AUC:", round(val_base$pr_auc, 3),
    "  TEST PR-AUC:", round(test_base$pr_auc, 3), "\n")

cat("\nTop 15 single features by ??VAL PR-AUC:\n")
print(feat_rank[order(-d_val_pr_auc)][1:15])

if (nrow(fwd_summary)) {
  cat("\nForward selection path:\n")
  print(fwd_summary)
  # Evaluate final forward model on TEST for reference
  final_test <- eval_set(current_model, test[complete.cases(test[, c(need_base, picked), with=FALSE])], BL_test)
  cat("\nFinal forward model TEST PR-AUC:", round(final_test$pr_auc, 3),
      "  mean NLL:", round(final_test$mean_nll, 3), "\n")
} else {
  cat("\nForward selection did not add features above min_gain =", min_gain, "\n")
}
