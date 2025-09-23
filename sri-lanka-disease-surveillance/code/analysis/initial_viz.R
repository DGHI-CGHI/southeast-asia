
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
