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
  
  # Candidates to try, in order
  cand <- c(
    Sys.which("unrar"),
    Sys.which("7z"),
    "C:/Program Files/7-Zip/7z.exe",
    "C:/Program Files (x86)/7-Zip/7z.exe",
    "/usr/local/bin/7z",           # macOS Intel (Homebrew)
    "/opt/homebrew/bin/7z",        # macOS Apple Silicon
    "/usr/bin/7z"                  # Linux
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
  
  # Try each candidate until one works
  for (exe in cand) {
    # Build args for unrar vs 7z
    if (grepl("unrar", basename(exe), ignore.case = TRUE)) {
      args <- c("x", "-o+", rar_file, out_dir)
    } else {
      args <- c("x", rar_file, paste0("-o", out_dir), "-y")
    }
    msg <- tryCatch(
      system2(exe, args, stdout = TRUE, stderr = TRUE),
      error = function(e) e$message
    )
    # Heuristic: consider success if a .tif shows up
    tif_found <- length(list.files(out_dir, pattern = "\\.tif(f)?$", recursive = TRUE)) > 0
    if (tif_found) return(invisible(TRUE))
  }
  
  stop("Extraction failed; tried: ", paste(basename(cand), collapse = ", "),
       ". Check permissions or use a different extractor.")
}



