# ---------- Plot helpers (styling + safe open) ----------
open_file <- function(path) {
  if (requireNamespace("berryFunctions", quietly = TRUE)) {
    berryFunctions::openFile(path)
  } else {
    message("berryFunctions not installed; opening with default browser.")
    utils::browseURL(paste0("file:///", normalizePath(path)))
  }
}

# One place to set your visual defaults
set_plot_theme <- function() {
  par(mar = c(5.2, 6.0, 4.2, 1.6),   # bottom, left, top, right
      mgp = c(3.1, 1.0, 0),          # label-to-axis distance (bigger first value = farther)
      tcl = -0.35,                   # tick length (negative = inward)
      las = 1,                       # y labels horizontal
      cex.axis = 1.25,
      cex.lab  = 1.5,
      cex.main = 1.7,
      xaxs = "i", yaxs = "i")        # flush to axis range
}

draw_hgrid <- function(y) abline(h = pretty(y), col = "gray90", lwd = 1)
draw_vgrid <- function(x) abline(v = pretty(x), col = "gray95", lwd = 1)

# village colors (tweak to taste)
col_vil <- c("Mandena" = "#2F3E46", "Sarahandrano" = "#7A7D7F")
