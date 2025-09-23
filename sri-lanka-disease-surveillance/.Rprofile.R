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