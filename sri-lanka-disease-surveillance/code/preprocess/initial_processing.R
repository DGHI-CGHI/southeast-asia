################################################################################
# CHI  Sri Lanka WER - Initial Processing Pipeline
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





setwd("C:/Users/jordan/R_Projects/dghi-southeast-asia/sri-lanka-disease-surveillance")

source("C:/Users/jordan/R_Projects/dghi-southeast-asia/sri-lanka-disease-surveillance/.Rprofile")
pdf_dir <- "C:/Users/jordan/R_Projects/sri-lanka-disease-surveillance/pdf_cache"
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)




# To install PDF related packages:
# install.packages("rJava")
# remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")

library(tabulapdf)       # ropensci fork; requires rJava

# ------------------------------------------------------------------------------
# 0) CONFIG
# ------------------------------------------------------------------------------
# Summary: Define project-relative paths using cfg (from .Rprofile) + file.path() helper.

# 1) No setwd() and no absolute C:\ paths
#    Everything below is relative to the Sri Lanka subproject root.

# Where to put scratch outputs during processing
paths <- list(
  temp_dir = file.path(cfg$paths$intermediate, "temp"),
  work_out = file.path(cfg$paths$intermediate, "outputs"),
  helpers  = file.path('code/helpers/helpers.R')
)


# 2) Source inputs that should live in the project repo:
#    Recommend moving PDFs/CSVs under data/raw or data/raw/external.
#    For now, keep your current filenames but make them relative.
paths$midyear_pop <- file.path(cfg$paths$raw, "Mid-year_population_by_district_and_sex_2024.pdf")
paths$wx_stations <- file.path(cfg$paths$raw , "SriLanka_Weather_Dataset.csv")
paths$era5_daily  <- file.path(cfg$paths$raw, "srilanka_district_daily_era5_areawt.csv")


# 3) Output figures directory inside project
paths$fig_dir <- file.path(cfg$paths$reports, "figures")

# 4) Additional named outputs (project-relative)
paths$outputs <- list(
  era5_weekly_aggregated = file.path(cfg$path$intermediate, 'srilanka_district_daily_era5_areawt.csv'),
  pdf_index_csv   = file.path(cfg$path$intermediate, 
                      "sri_lanka_WER_index_of_pdfs.csv"),
  case_counts_txt = file.path(cfg$path$intermediate, 
                      "disease_counts.txt")
)

# 5) ERA5 root **local path** (hydrated by DVC), not s3://
#    DVC should import/pull into data/raw/era5/
era5_root <- file.path(cfg$paths$raw, "era5")
# era5_root <- "s3://dghi-chi/data/se-asia/sri-lanka-disease-surveillance/era5/"


# 6) Create any needed directories once (safe if they already exist)
ensure_dir <- function(...) dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE)
ensure_dir(cfg$paths$intermediate, "temp")
ensure_dir(cfg$paths$intermediate, "outputs")
ensure_dir(cfg$paths$reports, "figures")


# 0.3 Helper function sources ------------------------------------------------
# These must define: DISEASES, POS_MAP, norm_dist(), .norm(), .is_footer(),
# .parse_ints(), .extract_district_from_row(), .pick_n(), DIST_CANON
source(paths$helpers)
# 
# ################################################################################
# # 1) INDEX WER PDFS -------------------------------------------------------------
# ################################################################################
# # Goal: Crawl the WER landing page and construct an absolute list of .pdf URLs,
# #       then parse volume/issue to derive a weekly date window where possible.
# 
# #  1.1 Crawl & collect absolute PDF URLs -------------------------------------
# wer_base <- "https://www.epid.gov.lk"
# wer_url  <- paste0(wer_base, "/weekly-epidemiological-report")
# 
# idx_html <- read_html(wer_url)
# hrefs    <- html_attr(html_elements(idx_html, "a"), "href")
# pdfs     <- unique(file.path("\\.pdf$", hrefs, value = TRUE))
# pdfs     <- ifelse(startsWith(pdfs, "http"), pdfs, paste0(wer_base, pdfs))
# pdfs = pdfs[pdfs != 'https://www.epid.gov.lk/storage/post/pdfs/en_63ec87e11b430_WER.pdf']
# 
# 
# # # https://www.epid.gov.lk/storage/post/pdfs/vol_34_no_01_english.pdf 12/30/06 -> 01/08/07
# # # https://www.epid.gov.lk/storage/post/pdfs/vol_34_no_02_english.pdf 2007-01-06 -> 2007-01-12
# # 
# # idx = 
# # idx[order(date_start)]
# # 
# # idx[45:55]
# 
# # tail(idx, 30)
# 
# #  1.2 Parse vol/issue  weekly date range (Vol  34) ------------------------
# # Anchor: Vol 34 No 1 -> week ending 2007-01-05 (Sat-Fri window).
# parse_issue_from_filename_v34plus <- function(u) {
#   f <- basename(u)
#   # vol   <- suppressWarnings(as.integer(str_match(f, "(i)vol[._ -]*(\\d+)")[,2]))
#   # issue <- suppressWarnings(as.integer(str_match(f, "(i)no[._ -]*(\\d+)")[,2]))
#   vol   <- str_match(f, "(?i)vol[._ -]*(\\d+)")[,2]
#   issue <- str_match(f, "(?i)no[._ -]*(\\d+)")[,2]
#   
#   vol   <- as.integer(vol)
#   issue <- as.integer(issue)
#   # 1
#   # 49
#   # 53
#   # 48
#   # 
#   # 
#   # Anchor: Vol 34 No 1 -> Week ending 2007-01-05 (Sat-Fri window)
#   base_start <- as.Date("2006-12-30")
#   base_end   <- as.Date("2007-01-05")
#   
#   offset <- if (!is.na(vol) && !is.na(issue)) (vol - 34) * 52 + (issue - 1) else NA_integer_
#   date_start <- if (!is.na(offset)) base_start + 7*offset else as.Date(NA)
#   date_end   <- if (!is.na(offset)) base_end   + 7*offset else as.Date(NA)
#   
#   data.table(
#     url        = u,
#     file       = f,
#     volume     = vol,
#     issue      = issue,
#     year       = if (!is.na(date_end)) year(date_end) else NA_integer_,
#     week       = if (!is.na(date_end)) isoweek(date_end) else NA_integer_,
#     date_start = date_start,
#     date_end   = date_end,
#     language   = fifelse(grepl("eng|english", f, TRUE), "English",
#                          fifelse(grepl("tam|tamil", f, TRUE), "Tamil",
#                                  fifelse(grepl("sin|sinhala", f, TRUE), "Sinhala", NA_character_)))
#   )
# }
# library(data.table)
# library(stringr)
# library(ISOweek)
# 
# # Helper: number of ISO weeks in a calendar year
# iso_weeks_in_year <- function(y) {
#   # ISO week-year rule: the week containing Jan 4 is week 1
#   # Using ISOweek: last ISO week of y is the week of Dec 28
#   iw <- ISOweek(sprintf("%d-12-28", y))  # e.g., "2021-W52-2"
#   as.integer(str_match(iw, "W(\\d{2})")[,2])
# }
# # 
# 
# parse_issue_from_filename_wer <- function(u, lang = "English") {
#   f <- basename(u)
#   m <- str_match(f, "(?i)vol[._ -]*(\\d+)[^\\d]+no[._ -]*(\\d+)")
#   if (any(is.na(m))) file.path("Could not parse volume/issue from: ", f)
#   
#   vol   <- as.integer(m[, 2])
#   issue <- as.integer(m[, 3])
#   year_cal <- 1973L + vol                 # WER convention: volume = year - 1973
#   
#   jan1  <- as.Date(sprintf("%04d-01-01", year_cal))
#   w0    <- as.POSIXlt(jan1)$wday          # 0=Sun,...,5=Fri,6=Sat
#   d2fri <- (5L - w0 + 7L) %% 7L           # days from Jan 1 to first Friday
#   first_friday <- jan1 + d2fri
#   
#   date_end   <- first_friday + 7L * (issue - 1L)  # Nth Friday
#   date_start <- date_end - 6L                     # Sat..Fri window
#   
#   # Week-in-volume index (1..52/53). Safe integer arithmetic for Dates:
#   week_in_volume <- as.integer((as.integer(date_end - first_friday)) / 7L) + 1L
#   
#   data.table(
#     url        = u,
#     file       = f,
#     volume     = vol,
#     issue      = issue,
#     year       = as.integer(format(date_end, "%Y")),  # calendar year of that Friday
#     week       = week_in_volume,                      # WER's own week numbering
#     date_start = date_start,
#     date_end   = date_end,
#     language   = lang
#   )
# }
# 
# 
# # s
# # 
# # # --- quick sanity checks on your examples ---
# # parse_issue_from_filename_wer("https://.../vol_48_no_53-english.pdf")
# # parse_issue_from_filename_wer("https://.../vol_49_no_01-english.pdf")
# # 
# # 
# # # vol 48 no 01 is 12/26/20 - 01/01/21
# # # last pdf for 2022 is vol 49 no 42 22-28 oct 2022
# # 
# # 
# # parse_issue_from_filename_wer('https://www.epid.gov.lk/storage/post/pdfs/vol_48_no_01-english.pdf')
# # parse_issue_from_filename_wer('https://www.epid.gov.lk/storage/post/pdfs/vol_48_no_53-english.pdf')
# # parse_issue_from_filename_wer('https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_01-english.pdf')
# # parse_issue_from_filename_wer('https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_02-english.pdf')
# # 
# 
# for (aa in 1:length(pdfs)){
#   print(parse_issue_from_filename_wer(pdfs[aa]))
# }
# ]
# # idx <- rbindlist(lapply(pdfs, parse_issue_from_filename_v34plus), fill = TRUE)
# idx <- rbindlist(lapply(pdfs, parse_issue_from_filename_wer), fill = TRUE)
# 
# # 
# # 
# # 
# # idx[,c("file","date_start","date_end")]
# # idx2[,c("file","date_start","date_end")]
# # 
# # ff = merge(idx[,c("file","date_start","date_end",'url')], idx2[,c("file","date_start","date_end")], by = c("file"))
# # ff[date_start.x === date_start.y]
# # 
# 
# 
# 
# 
# idx <- idx[!is.na(url)]
# 
# idx = idx[!is.na(year)]
# idx = idx[year >= 2014]
# 
# # which(idx$url == 'https://www.epid.gov.lk/storage/post/pdfs/vol_43_no_02-english.pdf')
# # 
# # which(idx$url == 'https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_01-english.pdf')
# # 
# # which(idx$url == "https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_02-english_2.pdf")
# # 
# # u = 'https://www.epid.gov.lk/storage/post/pdfs/vol_48_no_53-english.pdf'
# # u = 'https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_01-english.pdf'
# # lepto_dt[date_start == "2021-12-18"]$url[1]
# # 
# # 
# # parse_issue_from_filename_v34plus('https://www.epid.gov.lk/storage/post/pdfs/vol_48_no_53-english.pdf')
# # parse_issue_from_filename_v34plus('https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_01-english.pdf')
# 
# #  1.3 Persist the index for reproducibility ---------------------------------
# fwrite(idx, paths$outputs$pdf_index_csv)
# 
# ################################################################################
# # 2) EXTRACT DISEASE TABLES FROM PDFs ------------------------------------------
# ################################################################################
# 
# # Goal: For each WER PDF, parse the district-level disease table using a
# #       position-based mapping (two columns per disease: *_A, *_B).
# #       Returns a wide format with one row per district per PDF table row.
# # Try to read "Week ending ..." from the first page text.
# extract_week_ending <- function(pdf_path) {
#   txt <- extract_text(pdf_path, pages = 1)
#   # Examples found in headers: "WEEK ENDING 17th December 2021" or "Week ending 17.12.2021"
#   pat <- "(?i)week\\s*ending[^0-9]*(\\d{1,2})[\\.\\-/\\s]*([A-Za-z]{3,9}|\\d{1,2})[\\.\\-/\\s]*(\\d{4})"
#   m <- str_match(txt, pat)
#   if (any(is.na(m))) return(NA)
#   d <- as.integer(m[,2]); mon <- m[,3]; y <- as.integer(m[,4])
#   # Month can be a word or a number
#   if (grepl("^[A-Za-z]+$", mon)) {
#     mon <- match(tolower(substr(mon,1,3)), tolower(month.abb))
#   } else mon <- as.integer(mon)
#   as.Date(sprintf("%04d-%02d-%02d", y, mon, d))
# }
# 
# extract_week_ending(pdf_path)
# 
# # 
# 

# # from helpers.R????
# POS_MAP = readRDS("C:/Users/jordan/R_Projects/sri-lanka-disease-surveillance/posmap.Rds")
# 

#  2.1 Core extractor (position-based A/B pairs across all pages) ------------
extract_all_diseases_by_position <- function(
    pdf_path,
    diseases = DISEASES,
    pos_map  = POS_MAP,
    keep_total = FALSE,
    debug = FALSE
) {
  
  tabs <- tryCatch(
    suppressWarnings(
      suppressMessages(
        extract_tables(pdf_path, guess = TRUE, method = "stream", output = "tibble")
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(tabs) || !length(tabs)) return(NULL)
  
  max_idx <- max(unlist(lapply(pos_map, unlist)), na.rm = TRUE)
  
  out_all <- list()
  # # berryFunctions::openFile(pdf_path)
  # for (ti in seq_along(tabs)) {
  # 
  #   tab <- as.data.table(tabs[[ti]])  # <- as requested
  #   if (!is.data.table(tab) || nrow(tab) < 2 || ncol(tab) < 1) next
  #   
  #   # Force character; normalize whitespace
  #   tab_chr <- as.data.frame(lapply(tab, as.character), stringsAsFactors = FALSE)
  #   tab_chr[] <- lapply(tab_chr, .norm)
  #   
  #   rows_out <- vector("list", nrow(tab_chr)); n_out <- 0L
  #   
  #   # if (names(tabs[[ti]])[1] %like% "Colombo"){
  #   #   tab_chr = rbind(tab_chr, names(tabs[[ti]]))
  #   # }
  #   #   
  #     
  #   for (r in seq_len(nrow(tab_chr))) {
  #     s_row <- .norm(paste(tab_chr[r, ], collapse = " "))
  #     if (!nzchar(s_row) || .is_footer(s_row)) next
  #     
  #     district <- .extract_district_from_row(s_row)
  #     if (!nzchar(district)) next
  #     if (!keep_total && grepl("(i)^sri\\s*lanka$", district)) next
  #     
  #     ints <- .parse_ints(s_row)
  #     
  #     # Build a named list of values like dengue_A, dengue_B, .
  #     vals <- list(table_id = ti, district = district, n_numbers_in_row = length(ints))
  #     for (d in diseases) {
  #       idxA <- pos_map[[d]]$A
  #       idxB <- pos_map[[d]]$B
  #       vals[[paste0(d, "_A")]] <- .pick_n(ints, idxA)
  #       vals[[paste0(d, "_B")]] <- .pick_n(ints, idxB)
  #     }
  #     
  #     n_out <- n_out + 1L
  #     rows_out[[n_out]] <- as.data.table(vals)
  #   }
  #   
  #   if (n_out) out_all[[length(out_all)+1L]] <- rbindlist(rows_out[seq_len(n_out)], use.names = TRUE, fill = TRUE)
  # }
  
  
  out_all <- list()
  
  for (ti in seq_along(tabs)) {
    tab <- as.data.table(tabs[[ti]])
    if (!is.data.table(tab) || nrow(tab) < 2 || ncol(tab) < 1) next
    
    # Force character; normalize whitespace
    tab_chr <- as.data.frame(lapply(tab, as.character), stringsAsFactors = FALSE)
    tab_chr[] <- lapply(tab_chr, .norm)
    
    rows_out <- vector("list", nrow(tab_chr) + 1L)  # +1 in case we add a header-as-row
    n_out <- 0L
    had_colombo <- FALSE
    
    ## ---------------------------------------------------------------------------
    ## 1) Parse each *row* in the table as before
    ## ---------------------------------------------------------------------------
    for (r in seq_len(nrow(tab_chr))) {
      s_row <- .norm(paste(tab_chr[r, ], collapse = " "))
      if (!nzchar(s_row) || .is_footer(s_row)) next
      
      district <- .extract_district_from_row(s_row)
      
      # Robust NA / empty check
      if (is.na(district) || !nzchar(district)) next
      
      # Skip "Sri Lanka" total row if requested
      if (!keep_total &&
          grepl("^(sri\\s*lanka)$", district, ignore.case = TRUE)) {
        next
      }
      
      ints <- .parse_ints(s_row)
      
      vals <- list(
        table_id         = ti,
        district         = district,
        n_numbers_in_row = length(ints)
      )
      for (d in diseases) {
        idxA <- pos_map[[d]]$A
        idxB <- pos_map[[d]]$B
        vals[[paste0(d, "_A")]] <- .pick_n(ints, idxA)
        vals[[paste0(d, "_B")]] <- .pick_n(ints, idxB)
      }
      
      # Guard against NA here too
      if (!is.na(district) && nchar(district) > 18L) {
        # Instead of stop(), log and skip:
        # write(
        #   paste("Long district name in table", ti, ":", district),
        #   file = "pdf_parse_issues.log",
        #   append = TRUE
        # )
        next
      }
      
      if (!is.na(district) && identical(tolower(district), "colombo")) {
        had_colombo <- TRUE
      }
      
      n_out <- n_out + 1L
      rows_out[[n_out]] <- as.data.table(vals)
    }
    
    ## ---------------------------------------------------------------------------
    ## 2) If we *never* saw a Colombo row, but the first column name
    ##    starts with "Colombo", treat the header line as a synthetic row
    ## ---------------------------------------------------------------------------
    first_name <- names(tab)[1]
    if (!had_colombo &&
        !is.null(first_name) &&
        grepl("^\\s*colombo\\b", first_name, ignore.case = TRUE)) {
      
      # Build a fake "row" string from the column names themselves
      header_str <- .norm(paste(names(tab), collapse = " "))
      
      district_hdr <- .extract_district_from_row(header_str)
      ints_hdr     <- .parse_ints(header_str)
      
      if (!is.na(district_hdr) &&
          nzchar(district_hdr) &&
          identical(tolower(district_hdr), "colombo") &&
          length(ints_hdr) > 0L) {
        
        vals_hdr <- list(
          table_id         = ti,
          district         = district_hdr,
          n_numbers_in_row = length(ints_hdr)
        )
        for (d in diseases) {
          idxA <- pos_map[[d]]$A
          idxB <- pos_map[[d]]$B
          vals_hdr[[paste0(d, "_A")]] <- .pick_n(ints_hdr, idxA)
          vals_hdr[[paste0(d, "_B")]] <- .pick_n(ints_hdr, idxB)
        }
        
        if (!is.na(district_hdr) && nchar(district_hdr) <= 18L) {
          n_out <- n_out + 1L
          rows_out[[n_out]] <- as.data.table(vals_hdr)
        } else {
          # Optional logging instead of failing:
          # write(
          #   paste("Header-based Colombo name too long in table", ti, ":", district_hdr),
          #   file = "pdf_parse_issues.log",
          #   append = TRUE
          # )
        }
      }
    }
    
    ## ---------------------------------------------------------------------------
    ## 3) Bind rows for this table
    ## ---------------------------------------------------------------------------
    if (n_out > 0L) {
      out_all[[length(out_all) + 1L]] <-
        rbindlist(rows_out[seq_len(n_out)], use.names = TRUE, fill = TRUE)
    }
  }
  
  for (ti in 1:length(out_all)){
    if (sum(out_all[[ti]]$n_numbers_in_row, na.rm=TRUE) < 300){
      out_all[[ti]] <- data.table()
    }
  }
  
  if (!length(out_all)) return(NULL)
  out <- rbindlist(out_all, use.names = TRUE, fill = TRUE)
  
  out = out[district != 'Table']
  
  
  # Remove stand-alone "NA" tokens (case-insensitive just in case)
  out[, district := gsub("\\bNA\\b", "", district, perl = TRUE, ignore.case = FALSE)]
  
  # Trim any leftover extra spaces
  out[, district := trimws(gsub("\\s{2,}", " ", district))]  
  
  
  if ( out[,.N,by=district][district == 'Ampara']$N > 1){
    stop("More than 1 Ampara")
  }
  
  
  if ('Page' %in% out$district){
    stop("Page found")
  }
  
  # if ( out[,.N,by=district][district == 'Kalmune']$N > 1){
  # #   stop("Kalmune found")
  # # }
  # if ('Kalmune' %in% out$district){
  #   stop("Kalmune found")
  # }
  # if ('Kalmunei' %in% out$district){
  #   stop("Kalmunei found")
  # }
  # 
  # if ( out[,.N,by=district][district == 'Kalmunei']$N > 1){
  #   stop("Kalmunei found")
  # }
  
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
  
  out[order(district)]
}


# pdf_path = 'C:/Users/jordan/R_Projects/sri-lanka-disease-surveillance/pdf_cache/en_68f5d334eb339_Vol_52_no_36-english.pdf'

# pdf_path = 'C:/Users/jordan/R_Projects/sri-lanka-disease-surveillance/pdf_cache/vol_47_no_32-english_1.pdf'
## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------



################################################################################
# 1) INDEX WER PDFS -------------------------------------------------------------
################################################################################
# Goal: Crawl the WER landing page and construct an absolute list of .pdf URLs,
#       then parse volume/issue to derive a weekly date window where possible.

#  1.1 Crawl & collect absolute PDF URLs -------------------------------------
wer_base <- "https://www.epid.gov.lk"
wer_url  <- paste0(wer_base, "/weekly-epidemiological-report")

idx_html <- read_html(wer_url)
hrefs    <- html_attr(html_elements(idx_html, "a"), "href")
pdfs     <- unique(grep("\\.pdf$", hrefs, value = TRUE))
pdfs     <- ifelse(startsWith(pdfs, "http"), pdfs, paste0(wer_base, pdfs))


#  1.2 Parse vol/issue  weekly date range (Vol  34) ------------------------
# Anchor: Vol 34 No 1 -> week ending 2007-01-05 (Sat-Fri window).
parse_issue_from_filename_v34plus <- function(u) {
  f <- basename(u)
  # vol   <- suppressWarnings(as.integer(str_match(f, "(i)vol[._ -]*(\\d+)")[,2]))
  # issue <- suppressWarnings(as.integer(str_match(f, "(i)no[._ -]*(\\d+)")[,2]))
  vol   <- str_match(f, "(?i)vol[._ -]*(\\d+)")[,2]
  issue <- str_match(f, "(?i)no[._ -]*(\\d+)")[,2]
  
  vol   <- as.integer(vol)
  issue <- as.integer(issue)
  
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

idx = idx[!is.na(year)]
idx = idx[year >= 2014]

#  1.3 Persist the index for reproducibility ---------------------------------
fwrite(idx, paths$outputs$pdf_index_csv)

# #  2.2 Batch extractor over indexed PDFs -------------------------------------
library(data.table)
library(future)
library(future.apply)

## ------------------------------------------------------------------
## Logging setup
## ------------------------------------------------------------------
log_path <- "C:/Users/jordan/R_Projects/dghi-southeast-asia/sri-lanka-disease-surveillance/log_extract_issues.csv"

log_issue <- function(cur, file_path, stage, msg, extra = NA_character_) {
  # cur is a 1-row data.table / data.frame from idx
  entry <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stage     = stage,
    file_id   = if ("file" %in% names(cur)) as.character(cur$file) else NA_character_,
    url       = if ("url"  %in% names(cur)) as.character(cur$url)  else NA_character_,
    file_path = file_path,
    message   = msg,
    extra     = extra
  )
  
  # append = TRUE, write header only if file does not yet exist
  data.table::fwrite(
    entry,
    file      = log_path,
    append    = file.exists(log_path),
    col.names = !file.exists(log_path)
  )
}

## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------

library(future)
library(future.apply)

## ------------------------------------------------------------------
## 0) Parallel setup: 5 workers
## ------------------------------------------------------------------
## Make sure idx is a data.table (or data.frame) with at least:
##   file, url, date_start, date_end
setDT(idx)

## ------------------------------------------------------------------
## 1) PDF cache directory (constant across runs)
## ------------------------------------------------------------------
pdf_dir <- "C:/Users/jordan/R_Projects/dghi-southeast-asia/sri-lanka-disease-surveillance/pdf_cache"
# list.files(pdf_dir)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

choose_local_pdf_path <- function(cur, pdf_dir) {
  base <- basename(cur$file)
  
  if (!nzchar(base) || !grepl("\\.pdf$", base, ignore.case = TRUE)) {
    if (!is.null(cur$url) && !is.na(cur$url) &&
        grepl("^https?://", cur$url, ignore.case = TRUE)) {
      base <- basename(cur$url)
    }
  }
  if (!nzchar(base)) {
    base <- paste0("doc_", digest::digest(paste(cur$file, cur$url)), ".pdf")
  }
  if (!grepl("\\.pdf$", base, ignore.case = TRUE)) {
    base <- paste0(base, ".pdf")
  }
  
  file.path(pdf_dir, base)
}

process_one_pdf <- function(cur, pdf_dir, log_path) {
  # cur is 1-row data.table or data.frame
  has_url <- !is.null(cur$url) && !is.na(cur$url) &&
    grepl("^https?://", cur$url, ignore.case = TRUE)
  
  file <- cur$file
  
  if (has_url) {
    local_path <- choose_local_pdf_path(cur, pdf_dir)
    
    if (!file.exists(local_path)) {
      ok <- try(
        utils::download.file(cur$url, local_path, mode = "wb", quiet = TRUE),
        silent = TRUE
      )
      if (inherits(ok, "try-error")) {
        log_issue(cur, local_path, stage = "download",
                  msg   = conditionMessage(ok),
                  extra = "download failed")
        return(data.table())
      }
    }
    file <- local_path
  } else {
    # local file; if missing, see if it exists in cache
    if (!file.exists(file)) {
      candidate <- file.path(pdf_dir, basename(file))
      if (file.exists(candidate)) {
        file <- candidate
      } else {
        log_issue(cur, file, stage = "file_missing",
                  msg   = "Local file not found",
                  extra = NA_character_)
        return(data.table())
      }
    }
  }
  
  # ------------------ run extractor with error capture -----------------------
  res <- tryCatch(
    extract_all_diseases_by_position(file, debug = FALSE),
    error = function(e) {
      log_issue(cur, file, stage = "extract",
                msg   = conditionMessage(e),
                extra = "extract_all_diseases_by_position() error")
      return(NULL)
    }
  )
  
  if (is.null(res) || !nrow(res)) {
    return(data.table())
  }
  
  setDT(res)
  
  # normalize & attach metadata
  res[, `:=`(
    district   = district, # norm_dist(district),
    date_start = cur$date_start,
    date_end   = cur$date_end
  )]
  
  res <- res[!is.na(district)]
  
  # ------------------ check for overly long district names -------------------
  max_len <- max(nchar(res$district), na.rm = TRUE)
  if (max_len > 18L) {
    long_names <- unique(res[nchar(district) > 18L]$district)
    
    log_issue(
      cur,
      file,
      stage = "district_length",
      msg   = sprintf("Found district name(s) with length > 18 (max = %d)", max_len),
      extra = paste(long_names, collapse = "; ")
    )
    
    # OPTION A: drop the offending rows and keep going
    res <- res[nchar(district) <= 18L | is.na(district)]
    
    # OPTION B (if you want to ignore the entire file instead):
    # return(data.table())
  }
  
  # Attach identifiers for traceability
  res[, `:=`(
    file_id   = cur$file,
    url       = cur$url,
    file_path = file
  )]
  
  res
}


# # 1 minute for 50 files with 10 workers.
# # 4 minute for 10 files with 1 worker.
# idx = idx[1:10]
# Script finished at 2025-12-10 15:34:01 after 9.17 mins # 10 workers.
# Script finished at 2025-12-10 16:12:48 after 7.22 mins # 15 workers
## ------------------------------------------------------------------
## 3) Parallel loop over all rows of idx
## ------------------------------------------------------------------
# plan(sequential)

plan(multisession, workers = 15)


jj::timed('start')

setDT(idx)

allresults <- future_lapply(
  X = seq_len(nrow(idx)),
  FUN = function(i) {
    cur <- idx[i]
    out <- process_one_pdf(cur, pdf_dir = pdf_dir, log_path = log_path)
    
    if (i %% 25L == 0L) {
      message(".processed ", i, " PDFs")
    }
    
    out
  },
  future.seed = TRUE
)

allresults_dt <- rbindlist(allresults, use.names = TRUE, fill = TRUE)

jj::timed('end')

plan(sequential)
gc()




## ------------------------------------------------------------------
## 4) (Optional) Bind everything into one big data.table
## ------------------------------------------------------------------
lepto_dt <- rbindlist(allresults, use.names = TRUE, fill = TRUE)


lepto_dt[, `:=`(lepto = as.integer(leptospirosis_A),
                dengue = as.integer(dengue_A),
                year   = year(date_start))]


fwrite(lepto_dt, "data/intermediate/disease_counts_v5_TMP.txt")
fwrite(lepto_dt, "xdata/intermediate/disease_counts_v5_TMP_backup.txt")











################################################################################
# 1) LOAD WEEKLY DISEASE COUNTS -------------------------------------------------
################################################################################
# SECTION GOAL: Read district-week disease counts, normalize, and precompute dates.
lep <- fread( "data/intermediate/disease_counts_v5_TMP.txt")

# # 
# # # Normalize key fields and derive mid-week date (used for daily joins)
# lep[, date_mid := as.IDate(date_start + (as.integer(date_end - date_start) / 2))]
# lep[, `:=`(
#   # district = norm_dist(district),
#   date_mid = as.IDate(date_mid),
#   date_end = as.IDate(date_end)
# )]
# 
# lep = lep[year >= 2014]
# # 
# # lep[district == 'Ampara']




################################################################################
# 2) POPULATION: PDF  TIDY (or CSV fallback) ----------------------------------
################################################################################
# SECTION GOAL: Produce pop_dt with columns [district, year, poptot].
# Preference order:
#   (A) CHI_POP_CSV provided   read directly (no Java)
#   (B) Else parse from PDF via tabulizer (requires Java stack)
# NOTE: Requires working rJava/Tabulizer on your machine; handled with tryCatch.
pop_dt <- NULL
if (file.exists(paths$midyear_pop)) {
  tabs <- tryCatch(
    extract_tables(
      paths$midyear_pop,
      method = "lattice",
      guess = TRUE,
      output = "tibble"
    ),
    error = function(e)
      NULL
  )
  if (!is.null(tabs) && length(tabs) >= 1) {
    yrs <- 2014:2024
    pieces <- list()
    for (pg in seq_along(tabs)) {
      tb <- as.data.table(tabs[[pg]])
      yrheads <- names(tb)[names(tb) %like% paste(yrs, collapse = "|")]
      if (!length(yrheads))
        next
      # For each detected year column, grab District names + "Total" col under that year
      for (yhnum in seq_along(yrheads)) {
        yh = yrheads[[yhnum]]
        district_names <- tb$District[-1]
        # col_idx <- which(tb[1] == "Total")[which(names(tb) == yh)]
        col_idx <- which(tb[1] == "Total")[yhnum]
        if (!length(col_idx))
          next
        vals <- tb[[col_idx]][-1]
        pieces[[length(pieces) + 1L]] <- data.table(
          year    = as.integer(gsub("\\*", "", yh)),
          district = district_names,
          poptot  = as.numeric(gsub(",", "", vals)) * 1000
        )
      }
    }
    pop_dt <- rbindlist(pieces, fill = TRUE)
    pop_dt[, district := norm_dist(district)]
  }
}
stopifnot(!is.null(pop_dt) && nrow(pop_dt) > 0)
pop_dt
pop_dt = pop_dt[!is.na(district)] # remove national level totals.
# 
# pop_dt[,.N,by=year]
# 
# pop_dt1 =copy(pop_dt)
# 


pop_dt[district == 'Ampara']
# pop_dt[!district %in% lep$district]
# 
# lep[!district %in% pop_dt$district][,.N,by=district]

# lep[,.N,by=district]


lep = lep[district != 'Na']
lep[, district := gsub(" Na", "", district)]
lep = lep[!district %like% c("Timeliness")]
lep = lep[!district %like% c("Timely")]
lep = lep[!district %like% c("Selected Notifiable Diseases")]
lep = lep[!district %like% c("Tab Le")]
lep = lep[district != 'Na']





lep[,.N,by=district]
lep[,.N,by=district][order(N)]







# 'Anuradha' := 'Anuradhapu'
# 'Kili-'='Kilinochchi'
# Kilinoch-='Kilinochchi'
# 
# Nuwara-Eliya
# Nuwara


districtnames = unique(lep$district)
# districtnames[districtnames %like% "Na"]

lep = lep[!tolower(district) %like% "srilanka"]

lep = lep[!is.na(leptospirosis_A) & !is.na(leptospirosis_B)] 


# -- 2.3 Keep only valid districts; fill pre-2014 with 2014 pop ----------------
lep[district == 'Nuwaraeliy', district := 'Nuwara-Eliya']
lep[district == 'Anuradhapur', district := 'Anuradhapura']
lep[district == 'Polonnaruw', district := 'Polonnaruwa']
lep[district == 'Paha', district := 'Gampaha']
lep[district == 'Trincomale', district := 'Trincomalee']
lep[district == 'Anuradhapu', district := 'Anuradhapura']
lep[district == 'M31atale', district := 'Matale']
lep[district == 'Hambantot', district := 'Hambantota']

# > alias_map
# raw        canon
# <char>       <char>
#   1:        Sri Lanka         <NA>
#   2:     Nuwara Eliya Nuwara-Eliya
# 3:      Nuwaraeliya Nuwara-Eliya
# 4:     Nuwara-eliya Nuwara-Eliya
# 5:       Moneragala   Monaragala
# 6:       Rathnapura    Ratnapura
# 7:         Kalmunai       Ampara
# 8:        Kalmuniya       Ampara
# 9:   Galle District        Galle
# 10: Colombo District      Colombo
# 11:  Ampara District       Ampara
# 12:          Puttlam     Puttalam
# 13:          Vavniya     Vavuniya

lep[district == 'Nuwara Eliya', district := 'Nuwara-Eliya']
lep[district == 'Nuwaraeliya', district := 'Nuwara-Eliya']
lep[district == 'Nuwara-eliya', district := 'Nuwara-Eliya']
lep[district == 'Moneragala', district := 'Monaragala']
lep[district == 'Rathnapura', district := 'Ratnapura']
# lep[district == 'Kalmunai', district := 'Ampara']
# lep[district == 'Kalmuniya', district := 'Ampara']
lep[district == 'Galle District', district := 'Galle']
lep[district == 'Colombo District', district := 'Colombo']
lep[district == 'Ampara District', district := 'Ampara']
lep[district == 'Puttlam', district := 'Puttalam']
lep[district == 'Vavniya', district := 'Vavuniya']

lep[district == 'NuwaraEliy', district := 'Nuwara-Eliya']
lep[district == 'NuwaraEliya', district := 'Nuwara-Eliya']


lep[district == 'Anuradhap', district := 'Anuradhapura']

lep[district == 'paha', district := 'Gampaha']



# 
# lep shows:
#   
#   Ampara: N = 604
# 
# Kalmunai + Kalmune + Kalmunei: 85 + 488 + 31 = 604
# 
# So there are exactly 604 Kal* rows and exactly 604 Ampara rows, i.e. one full "district-length" set 
# of weeks for each. That is strong evidence that Kalmunai is being reported separately but is within Ampara District. 
# That aligns with reality: Kalmunai is an area / municipal council within Ampara, not a separate district.
# 
# So:
lep[district %in% c("Kalmune", "Kalmunei"), district := "Kalmunai"]


notin = lep[!district %in% pop_dt$district] # [,.N,by=district][order(district)]
lep[!district %in% pop_dt$district][,.N,by=district][order(district)]

yesin = lep[district %in% pop_dt$district] # [,.N,by=district][order(district)]
lep[district %in% pop_dt$district][,.N,by=district][order(district)]

# For consistency with your district-level population table (25 districts), the cleanest solution is:
#   Map all three (Kalmunai, Kalmune, Kalmunei) ??? Ampara.


# 
# lep[district == 'Ampara']
# 
# lep[district == 'Page']
# 
# berryFunctions::openFile('C:/Users/jordan/R_Projects/sri-lanka-disease-surveillance/pdf_cache/vol_44_no_27-english.pdf')
# 
# lep[file_path == 'C:/Users/jordan/R_Projects/sri-lanka-disease-surveillance/pdf_cache/vol_44_no_27-english.pdf' & district == 'Page']
# 
# 
# lep[file_path == 'C:/Users/jordan/R_Projects/sri-lanka-disease-surveillance/pdf_cache/vol_44_no_27-english.pdf' & district == 'Colombo']
# 
# lep[file_path == 'C:/Users/jordan/R_Projects/sri-lanka-disease-surveillance/pdf_cache/vol_44_no_27-english.pdf'][,.N,by=district]
# a# Anuradhapur
# # Kalmune   
# 
# lep[file_path == 'C:/Users/jordan/R_Projects/sri-lanka-disease-surveillance/pdf_cache/en_68f5d334eb339_Vol_52_no_36-english.pdf'][order(district)]
# # Na Page
# # Polonnaruw    
# notin[,.N,by=district][,.N]
# 
# notin[,.N,by=district][order(district)]
# # 4:   Kalmunai     2
# 5:    Kalmune   489
# 6:   Kalmunei    31

# yesin[,.N,by=district][order(district)]
# 
# Hambantot
# 
# 
# notin[district == 'Hambantot'][,.N,by=year]
# yesin[district == 'Hambantota'][,.N,by=year]
# 




lep[,.N,by=district][order(district)]





lep[!district %in% pop_dt$district][,.N, by = district][order(district)]


# current mid year pop data starts in 2014, so for years in health data before that, simply using 2014 pop for now.
lep[, year2merge := fifelse(year >= 2014, year, 2014L)]
lep <- merge(
  lep,
  pop_dt,
  by.x = c("district", "year2merge"),
  by.y = c("district", "year"),
  all.x = TRUE
)

# -- 2.4 Compute rates per 100k -------------------------------------------------
lep[, `:=`(
  lepto_100k  = (lepto  / poptot) * 1e5,
  dengue_100k = (dengue / poptot) * 1e5
)]

lep[, year := NULL] # year no longer needed after merge


fwrite(lep, "data/intermediate/disease_counts_v5.txt")


################################################################################




