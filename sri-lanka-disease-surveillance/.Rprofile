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




## .Rprofile - put this in your project root
.local_open_welcome <- function() {
  # Prefer Quarto ??? then R Markdown ??? else fall back to README.md in source pane
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    proj <- rprojroot::find_root_file(".", criterion = rprojroot::is_rstudio_project)
    qmd  <- file.path(proj, "WELCOME.qmd")
    rmd  <- file.path(proj, "WELCOME.Rmd")
    readme_md <- file.path(proj, "README.md")
    
    # Render Quarto if present
    if (file.exists(qmd) && requireNamespace("quarto", quietly = TRUE)) {
      out <- try(quarto::quarto_render(qmd, quiet = TRUE), silent = TRUE)
      if (!inherits(out, "try-error")) {
        # Quarto returns the output path; open in Viewer
        rstudioapi::viewer(out)
        return(invisible())
      }
    }
    
    # Render Rmd if present
    if (file.exists(rmd) && requireNamespace("rmarkdown", quietly = TRUE)) {
      html <- file.path(tempdir(), "WELCOME.html")
      ok <- try(rmarkdown::render(rmd, output_file = html, quiet = TRUE), silent = TRUE)
      if (!inherits(ok, "try-error") && file.exists(html)) {
        rstudioapi::viewer(html)
        return(invisible())
      }
    }
    
    # Fallback: open README.md in a source tab
    if (file.exists(readme_md)) {
      rstudioapi::navigateToFile(readme_md)
    }
  } else {
    # Not in RStudio ??? open README.md in the default editor if available
    if (file.exists("README.md")) try(file.edit("README.md"), silent = TRUE)
  }
}

if (interactive()) {
  setHook("rstudio.sessionInit", function(isNewSession) {
    # Run both: open your working scripts and the welcome page
    if (!isTRUE(getOption("project_docs_opened"))) {
      .local_open_files()       # your existing function that opens scripts
      .local_open_welcome()     # NEW: opens the tidy README/welcome
      options(project_docs_opened = TRUE)
    }
  }, action = "append")
}



.local_open_files <- function() {
  wanted <- c(
    
    "code/analysis/sri_lanka_modeling.R",
    "code/preprocess/era5_to_weekly_features.R",
    "code/preprocess/initial_processing.R",  # preferred location
    "code/preprocess/post_processing.R",
    "code/initial_processing.R"              # fallback if it's here
  )
  
  files <- normalizePath(wanted, winslash = "/", mustWork = FALSE)
  files <- files[file.exists(files)]
  if (!length(files)) {
    message("Auto-open: no target files found (skipping).")
    return(invisible())
  }
  
  # Prefer RStudio API when available, else use file.edit()
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    message("Auto-open via rstudioapi: ", paste(basename(files), collapse = ", "))
    for (fp in files) try(rstudioapi::navigateToFile(fp), silent = TRUE)
    try(rstudioapi::executeCommand("activateConsole"), silent = TRUE)
  } else {
    message("Auto-open via file.edit(): ", paste(basename(files), collapse = ", "))
    for (fp in files) try(file.edit(fp), silent = TRUE)
  }
  
  # Don't re-open repeatedly in this session
  options(project_docs_opened = TRUE)
}

if (interactive()) {
  # Run after RStudio finishes initializing the session
  setHook("rstudio.sessionInit", function(isNewSession) {
    if (!isTRUE(getOption("project_docs_opened"))) .local_open_files()
  }, action = "append")
  
  # Fallback: if we're NOT in RStudio (or the hook doesn't fire), do it now.
  if (Sys.getenv("RSTUDIO") != "1" && !isTRUE(getOption("project_docs_opened"))) {
    .local_open_files()
  }
}

.local_open_readme_html <- function() {
  if (!requireNamespace("rstudioapi", quietly=TRUE) || !rstudioapi::isAvailable()) return(invisible())
  rd <- "README.md"
  if (!file.exists(rd)) return(invisible())
  html <- file.path(tempdir(), "README.html")
  # use rmarkdown if present; otherwise quick pandoc call
  if (requireNamespace("rmarkdown", quietly=TRUE)) {
    ok <- try(rmarkdown::render(rd, output_file = html, quiet = TRUE), silent = TRUE)
    if (!inherits(ok, "try-error") && file.exists(html)) rstudioapi::viewer(html)
  } else if (nzchar(Sys.which("pandoc"))) {
    cmd <- sprintf('"%s" "%s" -o "%s"', Sys.which("pandoc"), rd, html)
    system(cmd)
    if (file.exists(html)) rstudioapi::viewer(html)
  } else {
    rstudioapi::navigateToFile(rd)  # last resort: open as plain text
  }
}


