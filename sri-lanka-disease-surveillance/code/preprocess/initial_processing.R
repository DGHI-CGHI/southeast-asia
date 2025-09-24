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

# To install PDF related packages:
# install.packages("rJava")
# remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")

library(tabulapdf)       # ropensci fork; requires rJava

# ------------------------------------------------------------------------------
# 0) CONFIG
# ------------------------------------------------------------------------------
# Summary: Define project-relative paths using cfg (from .Rprofile) + p() helper.

# 1) No setwd() and no absolute C:\ paths
#    Everything below is relative to the Sri Lanka subproject root.

# Where to put scratch outputs during processing
paths <- list(
  temp_dir = p(cfg$paths$intermediate, "temp"),
  work_out = p(cfg$paths$intermediate, "outputs"),
  helpers  = p('code/helpers/helpers.R')
)


# 2) Source inputs that should live in the project repo:
#    Recommend moving PDFs/CSVs under data/raw or data/raw/external.
#    For now, keep your current filenames but make them relative.
paths$midyear_pop <- p(cfg$paths$raw, "Mid-year_population_by_district_and_sex_2024.pdf")
paths$wx_stations <- p(cfg$paths$raw , "SriLanka_Weather_Dataset.csv")
paths$era5_daily  <- p(cfg$paths$raw, "srilanka_district_daily_era5_areawt.csv")


# 3) Output figures directory inside project
paths$fig_dir <- p(cfg$paths$reports, "figures")

# 4) Additional named outputs (project-relative)
paths$outputs <- list(
  era5_weekly_aggregated = p(cfg$path$intermediate, 'srilanka_district_daily_era5_areawt.csv'),
  pdf_index_csv   = p(cfg$path$intermediate, 
                      "sri_lanka_WER_index_of_pdfs.csv"),
  case_counts_txt = p(cfg$path$intermediate, 
                      "disease_counts_v4.txt")
)

# 5) ERA5 root **local path** (hydrated by DVC), not s3://
#    DVC should import/pull into data/raw/era5/
era5_root <- p(cfg$paths$raw, "era5")
# era5_root <- "s3://dghi-chi/data/se-asia/sri-lanka-disease-surveillance/era5/"


# 6) Create any needed directories once (safe if they already exist)
ensure_dir <- function(...) dir.create(p(...), recursive = TRUE, showWarnings = FALSE)
ensure_dir(cfg$paths$intermediate, "temp")
ensure_dir(cfg$paths$intermediate, "outputs")
ensure_dir(cfg$paths$reports, "figures")


# 0.3 Helper function sources ------------------------------------------------
# These must define: DISEASES, POS_MAP, norm_dist(), .norm(), .is_footer(),
# .parse_ints(), .extract_district_from_row(), .pick_n(), DIST_CANON
source(paths$helpers)

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
pdfs = pdfs[pdfs != 'https://www.epid.gov.lk/storage/post/pdfs/en_63ec87e11b430_WER.pdf']


# # https://www.epid.gov.lk/storage/post/pdfs/vol_34_no_01_english.pdf 12/30/06 -> 01/08/07
# # https://www.epid.gov.lk/storage/post/pdfs/vol_34_no_02_english.pdf 2007-01-06 -> 2007-01-12
# 
# idx = 
# idx[order(date_start)]
# 
# idx[45:55]

# tail(idx, 30)

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
  # 1
  # 49
  # 53
  # 48
  # 
  # 
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
library(data.table)
library(stringr)
library(ISOweek)

# Helper: number of ISO weeks in a calendar year
iso_weeks_in_year <- function(y) {
  # ISO week-year rule: the week containing Jan 4 is week 1
  # Using ISOweek: last ISO week of y is the week of Dec 28
  iw <- ISOweek(sprintf("%d-12-28", y))  # e.g., "2021-W52-2"
  as.integer(str_match(iw, "W(\\d{2})")[,2])
}

parse_issue_from_filename_wer <- function(u, lang = "English") {
  f <- basename(u)
  m <- str_match(f, "(?i)vol[._ -]*(\\d+)[^\\d]+no[._ -]*(\\d+)")
  if (any(is.na(m))) stop("Could not parse volume/issue from: ", f)
  
  vol   <- as.integer(m[, 2])
  issue <- as.integer(m[, 3])
  year_cal <- 1973L + vol                 # WER convention: volume = year - 1973
  
  jan1  <- as.Date(sprintf("%04d-01-01", year_cal))
  w0    <- as.POSIXlt(jan1)$wday          # 0=Sun,...,5=Fri,6=Sat
  d2fri <- (5L - w0 + 7L) %% 7L           # days from Jan 1 to first Friday
  first_friday <- jan1 + d2fri
  
  date_end   <- first_friday + 7L * (issue - 1L)  # Nth Friday
  date_start <- date_end - 6L                     # Sat..Fri window
  
  # Week-in-volume index (1..52/53). Safe integer arithmetic for Dates:
  week_in_volume <- as.integer((as.integer(date_end - first_friday)) / 7L) + 1L
  
  data.table(
    url        = u,
    file       = f,
    volume     = vol,
    issue      = issue,
    year       = as.integer(format(date_end, "%Y")),  # calendar year of that Friday
    week       = week_in_volume,                      # WER's own week numbering
    date_start = date_start,
    date_end   = date_end,
    language   = lang
  )
}


# 
# 
# # --- quick sanity checks on your examples ---
# parse_issue_from_filename_wer("https://.../vol_48_no_53-english.pdf")
# parse_issue_from_filename_wer("https://.../vol_49_no_01-english.pdf")
# 
# 
# # vol 48 no 01 is 12/26/20 - 01/01/21
# # last pdf for 2022 is vol 49 no 42 22-28 oct 2022
# 
# 
# parse_issue_from_filename_wer('https://www.epid.gov.lk/storage/post/pdfs/vol_48_no_01-english.pdf')
# parse_issue_from_filename_wer('https://www.epid.gov.lk/storage/post/pdfs/vol_48_no_53-english.pdf')
# parse_issue_from_filename_wer('https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_01-english.pdf')
# parse_issue_from_filename_wer('https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_02-english.pdf')
# 


# idx <- rbindlist(lapply(pdfs, parse_issue_from_filename_v34plus), fill = TRUE)
idx <- rbindlist(lapply(pdfs, parse_issue_from_filename_wer), fill = TRUE)

# 
# 
# 
# idx[,c("file","date_start","date_end")]
# idx2[,c("file","date_start","date_end")]
# 
# ff = merge(idx[,c("file","date_start","date_end",'url')], idx2[,c("file","date_start","date_end")], by = c("file"))
# ff[date_start.x === date_start.y]
# 




idx <- idx[!is.na(url)]

idx = idx[!is.na(year)]
idx = idx[year >= 2014]

# which(idx$url == 'https://www.epid.gov.lk/storage/post/pdfs/vol_43_no_02-english.pdf')
# 
# which(idx$url == 'https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_01-english.pdf')
# 
# which(idx$url == "https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_02-english_2.pdf")
# 
# u = 'https://www.epid.gov.lk/storage/post/pdfs/vol_48_no_53-english.pdf'
# u = 'https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_01-english.pdf'
# lepto_dt[date_start == "2021-12-18"]$url[1]
# 
# 
# parse_issue_from_filename_v34plus('https://www.epid.gov.lk/storage/post/pdfs/vol_48_no_53-english.pdf')
# parse_issue_from_filename_v34plus('https://www.epid.gov.lk/storage/post/pdfs/vol_49_no_01-english.pdf')

#  1.3 Persist the index for reproducibility ---------------------------------
fwrite(idx, paths$outputs$pdf_index_csv)

################################################################################
# 2) EXTRACT DISEASE TABLES FROM PDFs ------------------------------------------
################################################################################

# Goal: For each WER PDF, parse the district-level disease table using a
#       position-based mapping (two columns per disease: *_A, *_B).
#       Returns a wide format with one row per district per PDF table row.
# Try to read "Week ending ..." from the first page text.
extract_week_ending <- function(pdf_path) {
  txt <- extract_text(pdf_path, pages = 1)
  # Examples found in headers: "WEEK ENDING 17th December 2021" or "Week ending 17.12.2021"
  pat <- "(?i)week\\s*ending[^0-9]*(\\d{1,2})[\\.\\-/\\s]*([A-Za-z]{3,9}|\\d{1,2})[\\.\\-/\\s]*(\\d{4})"
  m <- str_match(txt, pat)
  if (any(is.na(m))) return(NA)
  d <- as.integer(m[,2]); mon <- m[,3]; y <- as.integer(m[,4])
  # Month can be a word or a number
  if (grepl("^[A-Za-z]+$", mon)) {
    mon <- match(tolower(substr(mon,1,3)), tolower(month.abb))
  } else mon <- as.integer(mon)
  as.Date(sprintf("%04d-%02d-%02d", y, mon, d))
}

extract_week_ending(pdf_path)

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
  
  for (ti in seq_along(tabs)) {
    
    tab <- as.data.table(tabs[[ti]])  # <- as requested
    if (!is.data.table(tab) || nrow(tab) < 2 || ncol(tab) < 1) next
    
    if (any(names(tab) %like% "Colombo")){
      colmbovals = as.data.table(t(names(tab)))
      colmbovals[,] <- lapply(
        colmbovals,
        function(x) ifelse(grepl("\\.\\.\\.", x), NA_character_, x)
      )      
      names(tab) <- names(colmbovals)
      tab = rbindlist(list(colmbovals, tab))
    }

    # Force character; normalize whitespace
    tab_chr <- as.data.frame(lapply(tab, as.character), stringsAsFactors = FALSE)
    tab_chr[] <- lapply(tab_chr, .norm)
    
    rows_out <- vector("list", nrow(tab_chr)); n_out <- 0L
    
    for (r in seq_len(nrow(tab_chr))) {
      s_row <- .norm(paste(tab_chr[r, ], collapse = " "))
      if (!nzchar(s_row) || .is_footer(s_row)) next
      
      district <- .extract_district_from_row(s_row)

      if (!nzchar(district)) next
      if (!keep_total && grepl("(i)^sri\\s*lanka$", district)) next
      
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

#  2.2 Batch extractor over indexed PDFs -------------------------------------
# Strategy:
#   . If index row contains a URL, download to tempfile; else use local path
#   . Extract tables; attach normalized district names + WER week dates
allresults <- vector("list", nrow(idx))
for (i in seq_len(nrow(idx))) {
  cur <- idx[i]
  file <- cur$file
  if (grepl("^https://", cur$url, ignore.case = TRUE)) {
    tf <- tempfile(fileext = ".pdf")
    ok <- try(utils::download.file(cur$url, tf, mode = "wb", quiet = TRUE), silent = TRUE)
    if (inherits(ok, "try-error")) { allresults[[i]] <- data.table(); next }
    file <- tf
  }
  
  allresults[[i]] <- data.table()
  
  # Extract and attach week dates
  res <- extract_all_diseases_by_position(file, debug = FALSE)
  if (!is.null(res) && nrow(res)) {
    res[, `:=`(
      district   = norm_dist(district),
      date_start = cur$date_start,
      date_end   = cur$date_end
    )]
    
    res = res[!is.na(district)]

    res$file = cur$file
    res$url = cur$url
    
    if (max(nchar(res$district)) > 18) stop()
        
    
    allresults[[i]] <- res
  }
  if (i %% 25 == 0) message(".processed ", i, " PDFs")
}

#  2.3 Post-process + persist for downstream scripts -------------------------
lepto_dt <- allresults # rbindlist(allresults, fill = TRUE)

lepto_dt[, `:=`(lepto = as.integer(leptospirosis_A),
                dengue = as.integer(dengue_A),
                year   = year(date_start))]


# Remove erroneous rows included due to parsing challenges
lepto_dt = lepto_dt[!district %like% c("Timeliness")]
lepto_dt = lepto_dt[!district %like% c("Timely")]
lepto_dt = lepto_dt[!district %like% c("Selected Notifiable Diseases")]
lepto_dt = lepto_dt[!district %like% c("Tab Le")]
lepto_dt = lepto_dt[!district %like% c("Na Na Na Na")]




# fix names of district in few files resuting in appended " Na" to district names.
lepto_dt[substr(district, nchar(district)-2, nchar(district)) == ' Na', district := trimws(substr(district, 1, nchar(district)-2))]
# 
# Anuradhapur   245
# 29:         490
lepto_dt = lepto_dt[district != "Srilanka"]

lepto_dt[,.N,by=district][order(district)]
lepto_dt[,.N,by=district][order(N)]

# Rename districts where records had issues with correct parsing of names. 
lepto_dt[district == 'Paha', district := 'Gampaha']



# lepto_dt[district == 'Ampara']


dupes <- lepto_dt[district == "Ampara", .N, by = date_start][N > 1]$date_start
# Relabel the *second and beyond* records for those dates
lepto_dt[
  district == "Ampara" & date_start %in% dupes,
  district := ifelse(seq_len(.N) > 1, "aaaa Ampara", district),
  by = date_start
]


dupes2 <- lepto_dt[district == "Kalmune", .N, by = date_start][N > 1]$date_start



lepto_dt[date_start == "2021-12-11" & dengue_A == 1 & dengue_B == 67, district := 'Ampara']


ff = lepto_dt[date_start == "2021-12-11"][order(district)]
ff = unique(ff, by = c("district","dengue_A","dengue_B"))


lepto_dt[district == 'aaaa Ampara', district := 'Kalmune']
lepto_dt[district == 'Kalmunai', district := 'Kalmune']
lepto_dt[district == 'Kalmunei', district := 'Kalmune']


lepto_dt[district == 'Hambantot', district := 'Hambantota']

lepto_dt[district == 'Anuradhapur', district := 'Anuradhapura']
lepto_dt[district == 'Anuradhap', district := 'Anuradhapura']
lepto_dt[district == 'Anuradhapu', district := 'Anuradhapura']

lepto_dt[district == 'Nuwaraeliy', district := 'Nuwara-Eliya']

lepto_dt[district == 'Trincomale', district := 'Trincomalee']
lepto_dt[district == 'Polonnaruw', district := 'Polonnaruwa']

lepto_dt[district == 'M31atale', district := 'Matale']


# uniquedistricts = unique(lepto_dt$district)[order(unique(lepto_dt$district))]
# uniquedistricts[substr(uniquedistricts, nchar(uniquedistricts)-3, nchar(uniquedistricts)) == 'tale']


fwrite(lepto_dt, paths$outputs$case_counts_txt)

berryFunctions::openFile(pdf_path)
# End of script.
#####################!



################################################################################
