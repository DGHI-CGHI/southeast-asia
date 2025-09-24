################################################################################
# Sri Lanka WER ??? Outbreak Labels + Weekly Modeling (Onset OR Counts)
# ------------------------------------------------------------------------------
# WHAT this script does
#   1) Builds a robust, past-only baseline per district using rolling
#      median + MAD (no leakage - current week is excluded).
#   2) Converts weekly counts into:
#        . state  = 1 while "in outbreak" (sticky with hysteresis)
#        . onset  = 1 for the first in-outbreak week of each run
#        . onset_next = 1 on the week before an onset (what we can predict)
#   3) Fits ONE of two modeling tracks (pick the one you want):
#        (A) Onset classifier (binomial) ??? P(onset next week | not in outbreak)
#        (B) Count model (negative binomial with offset) ??? expected counts
#
# WHY these choices
#   . Baseline = median + h*MAD is robust to outliers and quiet weeks.
#   . Excluding the current week from the baseline window avoids leakage that
#     can make triggers artificially easy to hit.
#   . Hysteresis (needing H consecutive quiet weeks to end) prevents flicker.
#   . Binomial for onset makes the target "rare events" explicit.
#   . Negative binomial for counts handles overdispersion (Poisson underfits var).
#
# Inputs you MUST have in dt:
#   district (chr/fct), date_start (Date), week (1..52), lepto (int),
#   poptot (numeric > 0), and any climate features you plan to use.
#
# Common gotchas we've run into (now handled here)
#   . Week 53 ??? coerced to 52 so cyclic splines "close the year".
#   . Baseline leakage ??? we lag counts before rolling stats.
#   . Early-window NA ??? warm-up rule handles first W weeks.
#   . Offset(log(pop)) ??? drop/guard rows with poptot <= 0.
#   . Data.table chaining: compute aggregates first, then add derived metrics.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
})

# =============================================================================
# 0) CONFIG - outbreak labeling + modeling switches + climate inputs
# =============================================================================
# This section controls TWO distinct things:
#   A) How we turn raw weekly counts into "outbreak state" and "onset" labels.
#   B) Which modeling track(s) to run: a binomial onset model and/or a count model.
#
# Why you/your colleague care: these knobs move precision/recall in opposite
# directions. If you raise thresholds, you'll get fewer false alarms (???precision)
# but miss more true spikes (???recall). We'll tune with a validation split, but
# the *label definition* below is still the single biggest lever.

# ---- A) Outbreak labeling hyperparameters -----------------------------------
# Baseline is "typical recent level" per district, computed from the *past only*:
#   baseline_t = median(last W weeks) + h_z * MAD(last W weeks).
# We fire a trigger when the current week meaningfully exceeds that baseline.

min_cases <- 4L     # Floor on weekly cases to even consider a trigger.
#   Why: avoids calling an "outbreak" on pure noise in very quiet districts.
#   Tune: 1-5. ???min_cases => fewer small blips counted as spikes (???precision, ???recall).

W         <- 8L     # Rolling window length for the baseline (in weeks).
#   Common choices: 8, 12, 26, 52. Larger W = smoother baseline (fewer spurious spikes),
#   but slower to adapt after a regime shift.

h_z       <- 1.0    # Strength of the baseline threshold: median + h_z * MAD.
#   MAD is robust scale; h_z is your "how unusual" dial.
#   Tune: 1.0-2.0. ???h_z => stricter baseline => fewer spikes (???precision, ???recall).

k_z       <- 1.0    # Z-threshold using the same MAD baseline: z = (y - median) / (MAD + eps).
#   We treat a week as "hot" if z >= k_z OR it's >= p_pct above baseline (see below).
#   Tune: 0.5-2.0. ???k_z => require stronger deviations to trigger.

p_pct     <- 0.25   # Percent above baseline alternative to z (0.25 = +25%).
#   Helpful where MAD is tiny (very stable districts). If either z>=k_z OR pct>=p_pct
#   is true (and >= min_cases), we trigger.

H         <- 2L     # Hysteresis: need H consecutive *below-baseline* weeks to end outbreak.
#   Prevents rapid on/off flicker. Tune: 1-3. ???H => longer outbreaks (stickier state).

# NOTE on interactions:
#   . (W, h_z, k_z, p_pct) jointly shape what "unusual" is. If you raise h_z and
#     k_z together, also consider lengthening W.
#   . In quiet districts, MAD can be ~0. The p_pct rule is your safety net there.

# ---- B) Modeling tracks ------------------------------------------------------
# Choose which model(s) to run after labeling.
#   . RUN_ONSET_CLASSIFIER: binary model to predict "OUTBREAK STARTS next week?"
#   . RUN_COUNT_MODEL:     count model to explain/predict weekly counts now/next.

RUN_ONSET_CLASSIFIER <- TRUE   # Binomial: outcome = onset_next (0/1)
RUN_COUNT_MODEL      <- TRUE   # Negative binomial: outcome = lepto (counts)

# Important discipline for RUN_ONSET_CLASSIFIER:
#   Use *only features available at decision time* (e.g., lagged climate).
#   Don't include same-week climate if you're calling this a "next-week" alert.

# ---- C) Climate terms to feed the models ------------------------------------
# This is the *intent* list; we'll intersect with the columns that actually exist.
# Keep this small (3-6 smooths) to avoid collinearity among lags/rolls.


# Climate terms you intend to use (will be intersected with available cols)
climate_terms <- c(
  "tmax_mean",
  "rh_mean_week", "rh_mean_week_lag1",
  "precip_tp_sum_week", "precip_tp_sum_week_lag1",
  "precip_tp_sum_week_roll4w_mean"
)

# ---------------------------------------------------------------------
# 1) Data: expect dt present; do minimal preflight checks / cleaning
# ---------------------------------------------------------------------

# TR is created in script, sri_lanka_modeling.R ******
dt <- copy(TR)
dt = dt[!is.na(tmax_mean)]

req <- c("district","date_start","lepto","poptot")
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing required columns in dt: ", paste(miss, collapse=", "))

setDT(dt)
setkey(dt, district, date_start)

# Week (1..52) - collapse 53 to 52 for cyclic seasonality
if (!"week" %in% names(dt)) {
  dt[, week := as.integer(strftime(date_start, "%V"))]
}
dt[week == 53L, week := 52L]
dt[, week := as.integer(pmax(1L, pmin(52L, week)))]

# Guard offset rows
dt <- dt[is.finite(poptot) & poptot > 0]

# ---------------------------------------------------------------------
# 2) Past-only rolling baseline (NO leakage)
#    baseline_t = median_{t-W..t-1} + h_z * MAD_{t-W..t-1}
# ---------------------------------------------------------------------
setorder(dt, district, date_start)
dt[, lepto_prev := shift(lepto, 1L, type = "lag"), by = district]  # exclude current week

dt[, medW := frollapply(lepto_prev, W, median, align = "right", na.rm = TRUE), by = district]
dt[, madW := frollapply(lepto_prev, W, mad,    align = "right", na.rm = TRUE), by = district]

dt[, baseline := pmax(0, medW + h_z * madW)]

# ---------------------------------------------------------------------
# 3) Build triggers, state (with hysteresis), onset, onset_next
# ---------------------------------------------------------------------
# Safe z / percent exceedance
eps <- 0.5  # stabilizer to avoid huge ratios near zero
dt[, z   := fifelse(!is.na(madW) & madW > 0, (lepto - medW) / (madW + eps), NA_real_)]
dt[, pct := fifelse(!is.na(medW) & medW > 0, (lepto - medW) / (medW + eps), NA_real_)]

# Warm-up history length (# of prior weeks)
dt[, n_hist := pmax(0L, seq_len(.N) - 1L), by = district]

# Trigger fires if:
#   . we're still in warm-up and cases >= min_cases, OR
#   . cases >= min_cases AND (z >= k_z OR pct >= p_pct)
dt[, trigger := (n_hist < W & lepto >= max(min_cases, 1L)) |
     (lepto >= min_cases & (
       (!is.na(z)   & z   >= k_z) |
         (!is.na(pct) & pct >= p_pct)
     ))]

# State machine with hysteresis: once ON, require H quiet weeks (below baseline)
dt[, state := {
  s <- integer(.N)
  on <- FALSE
  below <- 0L
  for (i in seq_len(.N)) {
    if (!on && isTRUE(trigger[i])) { on <- TRUE; below <- 0L }
    if (on) {
      s[i] <- 1L
      thr <- (ifelse(is.na(medW[i]), 0, medW[i])) + h_z * (ifelse(is.na(madW[i]), 0, madW[i]))
      if (lepto[i] < thr) {
        below <- below + 1L
        if (below >= H) { on <- FALSE; below <- 0L }
      } else {
        below <- 0L
      }
    }
  }
  s
}, by = district]

# Onset = first 1 after a 0 (per district)
dt[, onset := fifelse(state == 1L & shift(state, 1L, 0L) == 0L, 1L, 0L), by = district]

# Onset-next (what we can predict one week ahead)
dt[, onset_next := shift(onset, type = "lead"), by = district]

cat("Sanity checks:\n",
    "  In-outbreak weeks  :", sum(dt$state == 1L, na.rm = TRUE), "\n",
    "  Onset weeks        :", sum(dt$onset == 1L,  na.rm = TRUE), "\n",
    "  Onset_next (global):", sum(dt$onset_next == 1L, na.rm = TRUE), "\n", sep="")

# Final modeling frame
model_df <- copy(dt)
model_df[, district := factor(district)]

# ---------------------------------------------------------------------
# 4A) Onset classifier (binomial): predict onset_next among non-outbreak weeks
# ---------------------------------------------------------------------
if (RUN_ONSET_CLASSIFIER) {
  model_df_sub <- model_df[state == 0L & !is.na(onset_next)]
  if (nrow(model_df_sub) == 0L) stop("No rows available for onset classifier (check state/onset labels).")
  
  # Handle class imbalance with simple inverse-prevalence weights
  p_pos <- mean(model_df_sub$onset_next == 1, na.rm = TRUE)
  wts   <- ifelse(model_df_sub$onset_next == 1, (1 - p_pos) / p_pos, 1)
  
  # Keep only terms that exist
  use_terms <- intersect(climate_terms, names(model_df_sub))
  s_terms   <- if (length(use_terms)) sprintf("s(%s, bs='cr', k=6)", use_terms) else NULL
  rhs <- paste(c(s_terms, "s(week, bs='cc', k=20)", "s(district, bs='re')"), collapse = " + ")
  fml <- as.formula(paste0("onset_next ~ ", rhs))
  
  m_onset <- bam(
    fml,
    data   = model_df_sub,
    family = binomial(link = "logit"),
    weights = wts,
    method = "fREML",
    discrete = TRUE
  )
  
  cat("\n--- Onset classifier (binomial) ---\n")
  print(summary(m_onset))
}

# ---------------------------------------------------------------------
# 4B) COUNT MODEL - Weekly leptospirosis counts with population offset
# ---------------------------------------------------------------------
# GOAL
#   Model the expected *count* per district-week (??_t) while accounting for:
#     . population exposure via an offset:  log(??_t) = ... + log(pop_t)
#     . smooth climate effects (additive on the log scale)
#     . annual seasonality (cyclic spline on week)
#     . between-district heterogeneity (random effect)
# CURRENTLY NO LANDCOVER INFORMATION INCLUDED.
if (RUN_COUNT_MODEL) {
  use_terms <- intersect(climate_terms, names(model_df))
  s_terms   <- if (length(use_terms)) sprintf("s(%s, bs='cr', k=6)", use_terms) else NULL
  rhs <- paste(c(s_terms, "s(week, bs='cc', k=20)", "s(district, bs='re')"),
               collapse = " + ")
  fml <- as.formula(paste0("lepto ~ ", rhs))
  
  
  m_counts <- mgcv::bam(
    lepto ~
      s(tmax_mean,                      bs = "cr", k = 6) +
      s(rh_mean_week,                   bs = "cr", k = 6) +
      s(rh_mean_week_lag1,              bs = "cr", k = 6) +
      s(precip_tp_sum_week,             bs = "cr", k = 6) +
      s(precip_tp_sum_week_lag1,        bs = "cr", k = 6) +
      s(precip_tp_sum_week_roll4w_mean, bs = "cr", k = 6) +
      s(week, bs = "cc", k = 20) +
      s(district, bs = "re"),
    data = model_df,
    
    # --------------
    # NEGATIVE BINOMIAL
    # count data with overdispersion (variance > mean).
    # family = nb(link = "log"),   # << overdispersed counts (one strategy? - check for overdispersion)
    # --------------
    
    # --------------
    # POISSON
    family = poisson(link="log"),  # relative difference to quasipoisson -- unsure.
    # seems to simply tighten the CI relative to quasipoisson... ? 
    # --------------
    
    # --------------
    # QUASIPOISSON
    # family = quasipoisson(link="log"), # seems to be better than nb(link = "log")
    #    # can't recall, but certainly can't infer anything from deviance explained when using quasipoisson
    # --------------
    
    method = "fREML",
    # select  = TRUE,                 # shrink unhelpful smooths to ~0
    # gamma   = 1.4,                  # mild extra penalty for stability
    offset   = log(poptot),
    # discrete = TRUE,
    knots    = list(week = c(0.5, 52.5))
  )
  
  
  cat("\n--- Count model (NB with offset) ---\n")
  print(summary(m_counts))
}

plot_gam_smooths_base(m_counts, pages = c(3, 3))
BuiltUp
tmax_mean




as.data.frame(jj::countna(model_df))


################################################################################
# NOTES / Troubleshooting
# . If you see almost-all-positives from the classifier: your baseline + trigger
#   may be too permissive (raise h_z, raise k_z/p_pct, or increase W). ?
# . If NB model errors with knots/week: ensure week ??? {1..52} and pass
#   knots = list(week = c(0.5, 52.5)) if you switch to gam() (bam() handles fine).
# . If AR(1) residual structure matters, you can re-fit bam(..., rho=..., AR.start=...),
#   but only after confirming residual ACF(1) > ~0.05 and building AR.start per
#   district series. Keep that as a later refinement.
################################################################################












