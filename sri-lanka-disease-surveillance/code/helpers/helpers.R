# helpers.R  (add at top)

# ---- Portable paths ----
safeload <- function(pkg) { if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg); TRUE }
safeload("here"); safeload("fs")

data_root <- function() {
  dr <- Sys.getenv("CHI_DATA_ROOT", unset = NA_character_)
  if (is.na(dr) || dr == "") stop("CHI_DATA_ROOT is not set. Add it to .Renviron.", call. = FALSE)
  fs::path_norm(dr)
}

# Build paths from CHI_DATA_ROOT and/or project root
path_data   <- function(...) fs::path(data_root(), ...)
path_proj   <- function(...) fs::path(here::here(), ...)
dir_ensure  <- function(...) { fs::dir_create(fs::path(...)); fs::path(...) }

# ---- Project-specific convenience ----
# ERA5 root on disk (under CHI_DATA_ROOT)
era5_root <- function() path_data("gridded", "era5-africa", "processed")

# Central location for small "publishable" outputs (in-repo)
outputs_root <- function() path_proj("analytics", "madagascar", "outputs")

# Optional streaming directory (also under CHI_DATA_ROOT if it's large)
stream_root <- function() path_data("gridded", "era5land", "madagascar", "madagascar_hourly_daily_streamed")



plot_gam_smooths_base <- function(mod,
                                  data       = NULL,
                                  which      = NULL,
                                  pages      = c(2, 3),
                                  se_mult    = 2,
                                  show_ci    = TRUE,
                                  show_rug   = TRUE,
                                  col_line   = "black",
                                  col_ci     = "grey85",
                                  main_prefix= NULL,
                                  ref_rule   = c("median","mean","mode"),
                                  ...) {
  stopifnot(inherits(mod, "gam"))
  ref_rule <- match.arg(ref_rule)
  
  # 1) decide the reference dataset we'll pull values from
  if (is.null(data)) {
    # try to retrieve the original data if provided in the call
    data <- try({
      dname <- as.character(mod$call$data)
      if (length(dname) && nzchar(dname)) get(dname, envir = parent.frame())
    }, silent = TRUE)
    if (inherits(data, "try-error") || is.null(data)) {
      # fallback to model.frame (may miss offset base vars!)
      data <- model.frame(mod)
    }
  }
  data <- as.data.frame(data)  # we only read from it
  
  # 2) collect all variables the formula needs (covers offset/log/pmax, etc.)
  req_vars <- unique(all.vars(formula(mod)))
  
  # 3) make a single-row reference with sensible values for every required var
  ref_value <- function(x) {
    if (is.numeric(x)) {
      stats::median(x[is.finite(x)], na.rm = TRUE)
    } else if (is.factor(x)) {
      lv <- levels(x); lv[which.max(tabulate(x))]
    } else if (is.character(x)) {
      ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]
    } else {
      x[1]
    }
  }
  base_ref_df <- setNames(
    replicate(length(req_vars), NA, simplify = FALSE),
    req_vars
  )
  for (v in req_vars) {
    if (v %in% names(data)) {
      x <- data[[v]]
      # keep factor levels from data
      if (is.factor(x)) {
        base_ref_df[[v]] <- factor(ref_value(x), levels = levels(x))
      } else {
        base_ref_df[[v]] <- ref_value(x)
      }
    } else {
      # if a required symbol isn't in `data` (rare), set a neutral value
      base_ref_df[[v]] <- 0
    }
  }
  base_ref_df <- as.data.frame(base_ref_df, stringsAsFactors = FALSE)
  
  sm <- mod$smooth
  if (length(sm) == 0L) { message("No smooth terms."); return(invisible(NULL)) }
  
  # filter 'which'
  term_labels <- vapply(sm, function(s) s$label, character(1))
  if (!is.null(which)) {
    keep <- term_labels %in% which
    sm <- sm[keep]; term_labels <- term_labels[keep]
  }
  
  nr <- pages[1]; nc <- pages[2]
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
  par(mfrow = c(nr, nc), mar = c(4,4,2,1), oma = c(0,0,2,0))
  
  draw_1d <- function(sobj, term_label) {
    xname <- sobj$term
    x <- data[[xname]]
    if (!is.numeric(x)) {
      lev <- if (is.factor(x)) levels(x) else sort(unique(x))
      nd <- base_ref_df[rep(1, length(lev)), , drop = FALSE]
      nd[[xname]] <- if (is.factor(x)) factor(lev, levels = levels(x)) else lev
      pr <- predict(mod, newdata = nd, type = "terms", terms = term_label, se.fit = TRUE)
      eff <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
      lo <- eff - se_mult*se; hi <- eff + se_mult*se
      ylim <- range(c(lo, hi), na.rm = TRUE); pad <- 0.05*diff(ylim); ylim <- ylim + c(-pad, pad)
      plot(seq_along(lev), eff, xaxt="n", xlab=xname, ylab="partial effect",
           main = main_prefix %||% term_label, ylim = ylim, pch = 16, col = col_line)
      axis(1, at = seq_along(lev), labels = lev, las = 2, cex.axis = 0.8)
      if (show_ci) segments(x0 = seq_along(lev), y0 = lo, x1 = seq_along(lev), y1 = hi, col = "grey50")
      return(invisible(NULL))
    }
    xr <- range(x, finite = TRUE); xg <- seq(xr[1], xr[2], length.out = 200)
    nd <- base_ref_df[rep(1, length(xg)), , drop = FALSE]; nd[[xname]] <- xg
    pr <- predict(mod, newdata = nd, type = "terms", terms = term_label, se.fit = TRUE)
    eff <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
    lo <- eff - se_mult*se; hi <- eff + se_mult*se
    ylim <- range(c(lo, hi), na.rm = TRUE); pad <- 0.05*diff(ylim); ylim <- ylim + c(-pad, pad)
    plot(xg, eff, type="n", xlab=xname, ylab="partial effect",
         main = main_prefix %||% term_label, ylim = ylim)
    if (show_ci) polygon(c(xg, rev(xg)), c(lo, rev(hi)), border = NA, col = col_ci)
    lines(xg, eff, lwd = 2, col = col_line)
    abline(h = 0, lty = 3, col = "grey40")
    if (show_rug) rug(x[is.finite(x)], ticksize = 0.02, col = "grey40")
  }
  
  draw_2d <- function(sobj, term_label) {
    xname <- sobj$term[1]; yname <- sobj$term[2]
    x <- data[[xname]]; y <- data[[yname]]
    if (!is.numeric(x) || !is.numeric(y)) { plot.new(); title(main=term_label); return(invisible(NULL)) }
    xr <- range(x, finite = TRUE); yr <- range(y, finite = TRUE)
    xg <- seq(xr[1], xr[2], length.out = 80); yg <- seq(yr[1], yr[2], length.out = 80)
    grid <- expand.grid(xg, yg); names(grid) <- c(xname, yname)
    nd <- base_ref_df[rep(1, nrow(grid)), , drop = FALSE]
    nd[[xname]] <- grid[[xname]]; nd[[yname]] <- grid[[yname]]
    z <- matrix(as.numeric(predict(mod, newdata = nd, type = "terms", terms = term_label)),
                nrow = length(xg), ncol = length(yg))
    contour(xg, yg, z, nlevels = 10, xlab = xname, ylab = yname,
            main = main_prefix %||% paste0(term_label, " (contours)"))
  }
  
  for (j in seq_along(sm)) {
    s <- sm[[j]]; lbl <- term_labels[j]; d <- length(s$term)
    if (d == 1L) draw_1d(s, lbl)
    else if (d == 2L) draw_2d(s, lbl)
    else { plot.new(); title(main = paste0(lbl, " (", d, "D not plotted)")) }
  }
  invisible(NULL)
}




###############################################################################
# 1) HELPER FUNCTIONS ----------------------------------------------------------
###############################################################################

# -- String utils --------------------------------------------------------------
`%||%` <- function(a,b) if (!is.null(a)) a else b
.norm   <- function(x) { x <- gsub("[\u00A0]", " ", x, perl = TRUE); trimws(x) }

# -- District normalization (single source of truth) ---------------------------
DIST_CANON <- c("Colombo","Gampaha","Kalutara","Kandy","Matale","Nuwara-Eliya",
                "Galle","Matara","Hambantota","Jaffna","Kilinochchi","Mannar",
                "Vavuniya","Mullaitivu","Batticaloa","Ampara","Trincomalee",
                "Kurunegala","Puttalam","Anuradhapura","Polonnaruwa","Badulla",
                "Monaragala","Ratnapura","Kegalle","Kalmunai")

alias_map <- data.table(
  raw   = c("Sri Lanka","Nuwara Eliya","Nuwaraeliya","Nuwara-eliya","Moneragala","Rathnapura",
            "Kalmunai","Kalmuniya","Galle District","Colombo District","Ampara District",
            "Puttlam","Vavniya"),
  canon = c(NA_character_,"Nuwara-Eliya","Nuwara-Eliya","Nuwara-Eliya","Monaragala","Ratnapura",
            "Ampara","Ampara","Galle","Colombo","Ampara","Puttalam","Vavuniya")
)

norm_dist <- function(x) {
  x1 <- str_squish(x)
  x1 <- ifelse(is.na(x1) | x1 == "", NA_character_, x1)
  x1 <- str_replace_all(str_to_title(x1), "Nuwara Eliya", "Nuwara-Eliya")
  m  <- match(x1, alias_map$raw)
  x2 <- ifelse(!is.na(m), alias_map$canon[m], x1)
  i  <- match(tolower(x2), tolower(DIST_CANON))
  ifelse(!is.na(i), DIST_CANON[i], x2)
}

# -- Disease column layout map (A/B pairs by position) -------------------------
DISEASES <- c(
  "dengue","dysentery","encephalitis","enteric_fever",
  "food_poisoning","leptospirosis","typhus_f","viral_hep",
  "rabies","chickenpox","meningitis","leishmania","tuberculosis","wrcd"
)

.make_pos_map <- function(diseases) {
  setNames(lapply(seq_along(diseases), function(i) list(A = 2L*i - 1L, B = 2L*i)), diseases)
}
POS_MAP <- .make_pos_map(DISEASES)

# -- Row parsing helpers for PDF tables ----------------------------------------
.is_footer  <- function(x) grepl("(?i)^(total|source|key to table|page|wer\\s+sri\\s+lanka)", .norm(x %||% ""))
.parse_ints <- function(x) {
  xs <- str_extract_all(.norm(x %||% ""), "\\b\\d{1,7}\\b")[[1]]
  if (!length(xs)) integer(0) else as.integer(xs)
}
.extract_dist_from_row <- function(s) {
  m <- str_match(s, "^(.*?)(?=\\b\\d)")[,2]
  d <- .norm(m %||% "")
  if (.is_footer(d)) return("")
  gsub("(?i)^sri\\s*lanka\\s*$", "Sri Lanka", d, perl = TRUE)
}
.pick_n <- function(ints, n) if (length(ints) >= n) ints[n] else NA_integer_





`%||%` <- function(a,b) if (!is.null(a)) a else b
.norm <- function(x) { x <- gsub("[\u00A0]", " ", x, perl=TRUE); trimws(x) }
.is_footer <- function(x) grepl("(?i)^(total|source|key to table|page|wer\\s+sri\\s+lanka)", .norm(x %||% ""))

# Extract all integers (in order) from a string
.parse_ints <- function(x) {
  xs <- str_extract_all(.norm(x %||% ""), "\\b\\d{1,7}\\b")[[1]]
  if (!length(xs)) integer(0) else as.integer(xs)
}

# District name = text before first number (after collapsing the row)
.extract_district_from_row <- function(s_row) {
  m <- str_match(s_row, "^(.*?)(?=\\b\\d)")[,2]
  d <- .norm(m %||% "")
  if (.is_footer(d)) return("")
  d <- gsub("(?i)^sri\\s*lanka\\s*$", "Sri Lanka", d, perl=TRUE)
  d
}

# Build a 2*k position map from a disease vector
.make_pos_map <- function(diseases) {
  setNames(
    lapply(seq_along(diseases), function(i) list(A = 2L*i - 1L, B = 2L*i)),
    diseases
  )
}

# Default order from your screenshot (14 diseases -> 28 positions)
DISEASES <- c(
  "dengue", "dysentery", "encephalitis", "enteric_fever",
  "food_poisoning", "leptospirosis", "typhus_f", "viral_hep",
  "rabies", "chickenpox", "meningitis", "leishmania",
  "tuberculosis", "wrcd"
)
POS_MAP <- .make_pos_map(DISEASES)  # e.g., dengue A=1,B=2; dysentery A=3,B=4; .; wrcd A=27,B=28

# Helper to safely pick nth number
.pick_n <- function(ints, n) if (length(ints) >= n) ints[n] else NA_integer_



alias_map <- data.table(
  raw = c("Sri Lanka","Nuwara Eliya","Nuwaraeliya","Nuwara-eliya",
          "Moneragala","Rathnapura","Kalmunai","Kalmuniya",
          "Galle District","Colombo District","Ampara District",
          "Puttlam","Vavniya"),
  canon = c(NA_character_,"Nuwara-Eliya","Nuwara-Eliya","Nuwara-Eliya",
            "Monaragala","Ratnapura","Ampara","Ampara",
            "Galle","Colombo","Ampara",
            "Puttalam","Vavuniya")
)

norm_dist <- function(x) {
  x1 <- str_squish(x)
  x1 <- ifelse(is.na(x1) | x1 == "", NA_character_, x1)
  # title case but keep hyphenated Eliya
  x1 <- str_replace_all(str_to_title(x1), "Nuwara Eliya", "Nuwara-Eliya")
  # apply alias map
  m <- match(x1, alias_map$raw)
  x2 <- ifelse(!is.na(m), alias_map$canon[m], x1)
  x2
}

# ---- RAR extractor (tries unrar, then 7z; cross-platform) --------------------
extract_rar <- function(rar_file, out_dir) {
  if (!file.exists(rar_file)) stop("RAR file not found: ", rar_file)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Normalize to absolute paths
  rar_file <- normalizePath(rar_file, winslash = "/", mustWork = TRUE)
  out_dir  <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  
  cand <- c(
    Sys.which("unrar"),
    Sys.which("7z"),
    "C:/Program Files/7-Zip/7z.exe",
    "C:/Program Files (x86)/7-Zip/7z.exe",
    "/usr/local/bin/7z",
    "/opt/homebrew/bin/7z",
    "/usr/bin/7z"
  )
  cand <- cand[nzchar(cand)]
  
  if (!length(cand)) {
    stop(paste0(
      "No extractor found. Install one of:\n",
      "  . Windows: 7-Zip (https://www.7-zip.org/) and ensure 7z.exe is on PATH\n",
      "  . macOS:  brew install p7zip\n",
      "  . Linux:  sudo apt-get install p7zip-full\n",
      "Then re-run."
    ))
  }
  
  for (exe in cand) {
    if (grepl("unrar", basename(exe), ignore.case = TRUE)) {
      args <- c("x", "-o+", shQuote(rar_file), shQuote(out_dir))
    } else {
      args <- c("x", shQuote(rar_file), paste0("-o", shQuote(out_dir)), "-y")
    }
    msg <- tryCatch(
      system2(exe, args, stdout = TRUE, stderr = TRUE),
      error = function(e) e$message
    )
    tif_found <- length(list.files(out_dir, pattern = "\\.tif(f)?$", recursive = TRUE)) > 0
    if (tif_found) return(invisible(TRUE))
  }
  
  stop("Extraction failed; tried: ", paste(basename(cand), collapse = ", "),
       ". Check permissions or use a different extractor.")
}




