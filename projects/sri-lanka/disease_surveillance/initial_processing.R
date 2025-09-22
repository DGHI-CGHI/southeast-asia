################################################################################
# CHI ??? Sri Lanka WER - Initial Processing Pipeline
# File: analysis/sri_lanka/initial_processing.R
#
# **
# Parsing of PDFs takes >20-30 minutes. Start with post_processing.R script if
#     wanting to skip this for investigation.
# **
#
#
# Purpose:
#   Index Weekly Epidemiological Report (WER) PDFs from the Sri Lanka Ministry
#   of Health website and extract weekly district-level disease counts from
#   tabular content (position-mapped A/B columns). Saves two artifacts:
#     1) CSV index of WER PDFs with (vol, issue, week, date range)
#     2) Disease count wide-format table per district-week
#
# Upstream / Assumptions:
#   . Internet access to crawl the WER index and download PDFs
#   . Helper files provide the following objects/functions:
#       - DIST_CANON, DISEASES, POS_MAP
#       - norm_dist(), .norm(), .is_footer(), .parse_ints(),
#         .extract_district_from_row(), .pick_n()
#     (Sourced below from: helpers/helpers.R, analysis/sri_lanka/helpers.R)
#
# Downstream:
#   . Outputs used by later scripts to join population, weather, land cover,
#     and ERA5 aggregates, and to build analysis panels/plots.
#
# Inputs (paths are auto-derived from env vars with sensible defaults):
#   CHI_LOCAL_WORK   - working directory (temp + outputs)
#   CHI_GITHUB_ROOT  - repo root containing analysis/ and helpers/
#   JAVA_HOME        - JDK path (needed by tabulizer/tabulapdf via rJava)
#
# Outputs:
#   - analysis/sri_lanka/outputs/sri_lanka_WER_index_of_pdfs.csv
#   - analysis/sri_lanka/outputs/disease_counts_v4.txt
#
# Java Notes:
#   . Windows: set JAVA_HOME to a valid JDK (not just JRE), e.g.
#       "C:/Program Files/Eclipse Adoptium/jdk-17.0.16.8-hotspot"
#   . Linux: typically /usr/lib/jvm/java-17-openjdk-amd64
#   . If rJava complains about registry/JAVA_HOME on Windows, try UNSETTING
#     JAVA_HOME (let R/Java discover it), or ensure PATH contains the JDK bin.
#
# Repro Tips:
#   . Put CHI_LOCAL_WORK, CHI_GITHUB_ROOT, JAVA_HOME in ~/.Renviron
#   . Keep helper functions centralized; this script should only orchestrate.
#
# Author: Jordan Clark (DGHI CHI)
# Last updated: 2025-09-21
################################################################################


suppressPackageStartupMessages({
  # Core
  library(data.table)
  library(stringr)
  library(lubridate)
  
  # HTML crawl / PDFs
  library(pdftools)
  library(xml2)
  library(rvest)

  # Spatial (used in later stages; harmless to load here)
  library(sf)
  library(terra)
  library(exactextractr)
  
  # Plotting helpers (used downstream)
  library(ggplot2)
  library(scales)
  library(ragg)
})

setDTthreads(8)

###############################################################################
# 0) CONFIG & PATHS ------------------------------------------------------------
###############################################################################

# -- Configuration via environment variables (with sensible defaults) ----------
# Tip: set these in ~/.Renviron or project .Renviron
#   CHI_LOCAL_WORK="C:/Users/jordan/Desktop/srilanka"
#   CHI_GITHUB_ROOT="C:/Users/jordan/R_Projects/CHI-Data"
#   JAVA_HOME="C:/Program Files/Eclipse Adoptium/jdk-17.0.16.8-hotspot"
#################
#################  
# PDF table extraction (tabulizer stack via rJava)
# MUST SET JAVA_HOME FIRST before loading library.
Sys.setenv(JAVA_HOME = Sys.getenv("JAVA_HOME", unset = Sys.getenv("JAVA_HOME")))  # no-op if already set
Sys.setenv("JAVA_HOME"="C:/Program Files/Eclipse Adoptium/jdk-17.0.16.8-hotspot")

# To install PDF related packages:
# install.packages("rJava")
# remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")

library(tabulapdf)       # ropensci fork; requires rJava
library(tabulizerjars)

#################
#################

CHI_LOCAL_WORK  <- Sys.getenv("CHI_LOCAL_WORK",  unset = "C:/Users/jordan/Desktop/srilanka")
CHI_GITHUB_ROOT <- Sys.getenv("CHI_GITHUB_ROOT", unset = "C:/Users/jordan/R_Projects/CHI-Data")

setwd(CHI_GITHUB_ROOT)

paths <- list(
  temp_dir      = file.path(CHI_LOCAL_WORK, "temp"),
  work_out      = file.path(CHI_LOCAL_WORK, "outputs"),
  # Source inputs (checked into repo)
  midyear_pop   = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs/Mid-year_population_by_district_and_sex_2024.pdf"),
  wx_stations   = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs/station_data/SriLanka_Weather_Dataset.csv"),
  era5_daily    = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/inputs/srilanka_district_daily_era5_areawt.csv"),
  fig_dir       = file.path(CHI_GITHUB_ROOT, "analysis/sri_lanka/outputs/figures")
)


paths = c(paths, 
          
          outputs = list(
            
            pdf_index_csv = file.path(CHI_GITHUB_ROOT, 'analysis/sri_lanka/outputs/sri_lanka_WER_index_of_pdfs.csv'),
            case_counts_txt = file.path(CHI_GITHUB_ROOT, 'analysis/sri_lanka/outputs/disease_counts_v4.txt')
            
          )
)


dir.create(paths$temp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$fig_dir, recursive = TRUE, showWarnings = FALSE)

# ?????? 0.3 Helper function sources ------------------------------------------------
# These must define: DISEASES, POS_MAP, norm_dist(), .norm(), .is_footer(),
# .parse_ints(), .extract_district_from_row(), .pick_n(), DIST_CANON
source(file.path(CHI_GITHUB_ROOT, "/helpers/helpers.R"))
source(file.path(CHI_GITHUB_ROOT, "/analysis/sri_lanka/helpers.R"))

################################################################################
# 1) INDEX WER PDFS -------------------------------------------------------------
################################################################################
# Goal: Crawl the WER landing page and construct an absolute list of .pdf URLs,
#       then parse volume/issue to derive a weekly date window where possible.

# ?????? 1.1 Crawl & collect absolute PDF URLs -------------------------------------
wer_base <- "https://www.epid.gov.lk"
wer_url  <- paste0(wer_base, "/weekly-epidemiological-report")

idx_html <- read_html(wer_url)
hrefs    <- html_attr(html_elements(idx_html, "a"), "href")
pdfs     <- unique(grep("\\.pdf$", hrefs, value = TRUE))
pdfs     <- ifelse(startsWith(pdfs, "http"), pdfs, paste0(wer_base, pdfs))


# ?????? 1.2 Parse vol/issue ??? weekly date range (Vol ??? 34) ------------------------
# Anchor: Vol 34 No 1 -> week ending 2007-01-05 (Sat-Fri window).
parse_issue_from_filename_v34plus <- function(u) {
  f <- basename(u)
  vol   <- suppressWarnings(as.integer(str_match(f, "(?i)vol[._ -]*(\\d+)")[,2]))
  issue <- suppressWarnings(as.integer(str_match(f, "(?i)no[._ -]*(\\d+)")[,2]))
  
  # Anchor: Vol 34 No 1 -> Week ending 2007-01-05 (Sat-Fri window)
  base_start <- as.Date("2006-12-30")
  base_end   <- as.Date("2007-01-05")
  
  offset <- if (!is.na(vol) && !is.na(issue)) (vol - 34) * 52 + (issue - 1) else NA_integer_
  date_start <- if (!is.na(offset)) base_start + 7*offset else as.Date(NA)
  date_end   <- if (!is.na(offset)) base_end   + 7*offset else as.Date(NA)
  
  data.table(
    url        = u,
    file       = f,
    volume     = vol,
    issue      = issue,
    year       = if (!is.na(date_end)) year(date_end) else NA_integer_,
    week       = if (!is.na(date_end)) isoweek(date_end) else NA_integer_,
    date_start = date_start,
    date_end   = date_end,
    language   = fifelse(grepl("eng|english", f, TRUE), "English",
                         fifelse(grepl("tam|tamil", f, TRUE), "Tamil",
                                 fifelse(grepl("sin|sinhala", f, TRUE), "Sinhala", NA_character_)))
  )
}

idx <- rbindlist(lapply(pdfs, parse_issue_from_filename_v34plus), fill = TRUE)
idx <- idx[!is.na(url)]

# ?????? 1.3 Persist the index for reproducibility ---------------------------------
idx_path <- paths$outputs.pdf_index_csv
fwrite(idx, idx_path)


################################################################################
# 2) EXTRACT DISEASE TABLES FROM PDFs ------------------------------------------
################################################################################

# Goal: For each WER PDF, parse the district-level disease table using a
#       position-based mapping (two columns per disease: *_A, *_B).
#       Returns a wide format with one row per district per PDF table row.

# ?????? 2.1 Core extractor (position-based A/B pairs across all pages) ------------
extract_all_diseases_by_position <- function(
    pdf_path,
    diseases = DISEASES,
    pos_map  = POS_MAP,
    keep_total = FALSE,
    debug = FALSE
) {
  tabs <- tryCatch(
    extract_tables(pdf_path, guess = TRUE, method = "stream", output = "tibble"),
    error = function(e) NULL
  )
  if (is.null(tabs) || !length(tabs)) return(NULL)
  
  max_idx <- max(unlist(lapply(pos_map, unlist)), na.rm = TRUE)
  
  out_all <- list()
  
  for (ti in seq_along(tabs)) {
    tab <- as.data.table(tabs[[ti]])  # <- as requested
    if (!is.data.table(tab) || nrow(tab) < 2 || ncol(tab) < 1) next
    
    # Force character; normalize whitespace
    tab_chr <- as.data.frame(lapply(tab, as.character), stringsAsFactors = FALSE)
    tab_chr[] <- lapply(tab_chr, .norm)
    
    rows_out <- vector("list", nrow(tab_chr)); n_out <- 0L
    
    for (r in seq_len(nrow(tab_chr))) {
      s_row <- .norm(paste(tab_chr[r, ], collapse = " "))
      if (!nzchar(s_row) || .is_footer(s_row)) next
      
      district <- .extract_district_from_row(s_row)
      if (!nzchar(district)) next
      if (!keep_total && grepl("(?i)^sri\\s*lanka$", district)) next
      
      ints <- .parse_ints(s_row)
      
      # Build a named list of values like dengue_A, dengue_B, .
      vals <- list(table_id = ti, district = district, n_numbers_in_row = length(ints))
      for (d in diseases) {
        idxA <- pos_map[[d]]$A
        idxB <- pos_map[[d]]$B
        vals[[paste0(d, "_A")]] <- .pick_n(ints, idxA)
        vals[[paste0(d, "_B")]] <- .pick_n(ints, idxB)
      }
      
      n_out <- n_out + 1L
      rows_out[[n_out]] <- as.data.table(vals)
    }
    
    if (n_out) out_all[[length(out_all)+1L]] <- rbindlist(rows_out[seq_len(n_out)], use.names = TRUE, fill = TRUE)
  }
  
  for (ti in 1:length(out_all)){
    if (sum(out_all[[ti]]$n_numbers_in_row, na.rm=TRUE) < 300){
      out_all[[ti]] <- data.table()
    }
  }
  
  if (!length(out_all)) return(NULL)
  out <- rbindlist(out_all, use.names = TRUE, fill = TRUE)
  
  out = out[district != 'Table']
  
  # Optional diagnostics
  if (debug) {
    if (any(out$n_numbers_in_row < max_idx, na.rm = TRUE)) {
      message(sprintf("Some rows have fewer than %d numbers; corresponding *_B (and possibly later diseases) will be NA.", max_idx))
      print(out[n_numbers_in_row < max_idx])
    }
  }
  
  # Drop diagnostic column if you prefer
  out[, n_numbers_in_row := NULL]
  
  out[]
}

# ?????? 2.2 Batch extractor over indexed PDFs -------------------------------------
# Strategy:
#   . If index row contains a URL, download to tempfile; else use local path
#   . Extract tables; attach normalized district names + WER week dates
allresults <- vector("list", nrow(idx))
for (i in seq_len(nrow(idx))) {
  cur <- idx[i]
  file <- cur$file
  if (grepl("^https?://", cur$url, ignore.case = TRUE)) {
    tf <- tempfile(fileext = ".pdf")
    ok <- try(utils::download.file(cur$url, tf, mode = "wb", quiet = TRUE), silent = TRUE)
    if (inherits(ok, "try-error")) { allresults[[i]] <- data.table(); next }
    file <- tf
  }
  # Extract and attach week dates
  res <- extract_all_diseases_by_position(file, debug = FALSE)
  if (!is.null(res) && nrow(res)) {
    res[, `:=`(
      district   = norm_dist(district),
      date_start = cur$date_start,
      date_end   = cur$date_end
    )]
    allresults[[i]] <- res
  } else {
    allresults[[i]] <- data.table()
  }
  if (i %% 25 == 0) message(".processed ", i, " PDFs")
}

# ?????? 2.3 Post-process + persist for downstream scripts -------------------------
lepto_dt <- rbindlist(allresults, fill = TRUE)

lepto_dt[, `:=`(lepto = as.integer(leptospirosis_A),
                dengue = as.integer(dengue_A),
                year   = year(date_start))]

fwrite(lepto_dt, paths$outputs.case_counts_txt)

# End of script.
################################################################################
