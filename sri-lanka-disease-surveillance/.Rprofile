Sys.setenv(ARROW_S3_NO_SIGN_REQUEST = "1")  # public, no creds

# Activate renv if present (scoped to this subproject)
if (file.exists("renv.lock")) {
  try(renv::activate(), silent = TRUE)
}

# Anchor relative paths to this subproject (works with or without 'here')
if (requireNamespace("here", quietly = TRUE)) {
  here::i_am("sri-lanka-disease-surveillance.Rproj")
}

# Minimal config loader + path helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

load_cfg <- function() {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Please install 'yaml'")
  yaml::read_yaml("params.yaml")
}

p <- function(...) normalizePath(file.path(...), winslash = "/", mustWork = FALSE)

cfg <- try(load_cfg(), silent = TRUE)
if (inherits(cfg, "try-error") || is.null(cfg$paths)) {
  cfg <- list(paths = list(raw = "data/raw",
                           processed = "data/processed",
                           intermediate = "data/intermediate",
                           reports = "reports"))
}

# --- Java / rJava bootstrap ----------------------------------------------------
chi_java_ok <- FALSE

chi_try_load_rJava <- function() {
  suppressWarnings(suppressMessages(requireNamespace("rJava", quietly = TRUE)))
}

chi_try_common_windows_jdk <- function() {
  if (Sys.info()[["sysname"]] != "Windows") return(FALSE)
  candidates <- c(
    "C:/Program Files/Eclipse Adoptium/jdk-17",                 # Temurin default
    "C:/Program Files/Eclipse Adoptium",                        # parent (scan)
    "C:/Program Files/Java/jdk-17",                             # Oracle default
    "C:/Program Files/Java"                                     # parent (scan)
  )
  # Expand any parent dirs by scanning for jdk-17* folders
  expanded <- unlist(lapply(candidates, function(p) {
    if (dir.exists(p)) {
      subs <- list.dirs(p, full.names = TRUE, recursive = FALSE)
      c(p, subs)
    } else character(0)
  }))
  expanded <- unique(expanded[grepl("jdk-17", expanded, ignore.case = TRUE)])
  for (d in expanded) {
    Sys.setenv(JAVA_HOME = d)
    if (chi_try_load_rJava()) return(TRUE)
  }
  FALSE
}

# 1) If rJava loads, we’re done
if (chi_try_load_rJava()) {
  chi_java_ok <- TRUE
} else {
  # 2) Try to auto-point JAVA_HOME on Windows
  if (chi_try_common_windows_jdk()) {
    chi_java_ok <- TRUE
  } else {
    # 3) Still not OK: instruct the user exactly what to do
    msg <- paste(
      "\n--- Java setup needed for PDF extraction (tabulizer) ---",
      "What to do:",
      "  1) Install a 64-bit JDK 17 (e.g., Adoptium Temurin 17).",
      "  2) Set JAVA_HOME to that JDK folder (the one containing 'bin', e.g.,",
      "     C:/Program Files/Eclipse Adoptium/jdk-17.x.x), then restart R.",
      "",
      "Windows (recommended): System Properties → Advanced → Environment Variables →",
      "  New… User variable  Name=JAVA_HOME  Value=C:\\Program Files\\Eclipse Adoptium\\jdk-17.x.x",
      "  Then add %JAVA_HOME%\\bin to your PATH, or just restart RStudio.",
      "",
      "R-only fallback (not permanent):",
      '  In this project’s .Rprofile you can temporarily set:',
      '    Sys.setenv(JAVA_HOME = "C:/Program Files/Eclipse Adoptium/jdk-17.x.x")',
      "",
      "After Java is set, install PDF packages once:",
      '  install.packages("rJava")',
      '  if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")',
      '  remotes::install_github(c("ropensci/tabulizerjars","ropensci/tabulizer"), INSTALL_opts="--no-multiarch")',
      sep = "\n"
    )
    if (interactive()) {
      message(msg)
    } else {
      stop(msg, call. = FALSE)
    }
  }
}

rm(chi_try_load_rJava, chi_try_common_windows_jdk)
# --- end Java / rJava bootstrap ------------------------------------------------

strip_trailing_slash <- function(x) sub("/+$", "", x)
# choose ONE root (local-first fallback shown here)
get_era5_root <- function() {
  local <- p(cfg$paths$raw, "era5")           # after `dvc pull`
  if (dir.exists(local)) return(local)
  "s3://dghi-chi/data/se-asia/sri-lanka-disease-surveillance/era5"
}
era5_root <- strip_trailing_slash(get_era5_root())
