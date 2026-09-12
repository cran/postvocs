#' Parse a single GC-MS exported TXT file
#'
#' Extracts the \code{[MC Peak Table]} and
#' \code{[MS Similarity Search Results for Spectrum Process Table]} from
#' a Shimadzu GC-MS workstation text export. The function returns the two
#' tables as data frames in a list, with element names derived from the input
#' file name. No files are written to disk.
#'
#' @param txt_file Character. Path to a single GC-MS exported TXT file
#'   (e.g., \code{"data-raw/03.qgd.txt"}). The file name (without extension)
#'   is used to name the output list elements.
#' @param debug Logical. If \code{TRUE}, prints debugging information such as
#'   total lines, start positions, and column names. Default is \code{FALSE}.
#' @param encoding Character. Encoding of the input TXT file. Default is
#'   \code{"UTF-8"}. Change to \code{"GBK"} or \code{"latin1"} if the file
#'   contains non-UTF-8 characters (e.g., in NIST library paths).
#'
#' @return A list with two data frames, named as
#'   \code{"<basename>_PeakTable"} and \code{"<basename>_SearchResults"},
#'   where \code{<basename>} is the input file name without extension.
#'   \itemize{
#'     \item \strong{PeakTable}: The extracted peak table from
#'           \code{[MC Peak Table]}, containing columns like \code{Peak#},
#'           \code{Ret.Time}, \code{Area}, etc.
#'     \item \strong{SearchResults}: The extracted similarity search results
#'           from \code{[MS Similarity Search Results ...]}, containing
#'           columns like \code{Spectrum#}, \code{CAS#}, \code{Name}, etc.
#'           CAS numbers are cleaned by removing all spaces.
#'   }
#'
#' @keywords postvocs
#'
#' @details
#' The function expects the TXT file to have the following structure:
#' \itemize{
#'   \item A header line starting with \code{"Data File Name"} (used only for
#'         debug and not for naming).
#'   \item After \code{[MC Peak Table]}, a line like \code{# of Peaks\t38} or
#'         \code{# of Peaks,38} giving the number of peaks.
#'   \item The column name line of the peak table must start with \code{"Peak#"}
#'         (e.g., \code{Peak#\tRet.Time\tArea...} or \code{Peak#,Ret.Time,Area...}).
#'   \item The mass spectral search results section starts with
#'         \code{[MS Similarity Search Results ...]}, and its column name line
#'         must start with \code{"Spectrum#"}.
#' }
#'
#' The delimiter (tab or comma) is auto-detected from the column name line of
#' each table and used for parsing all rows of that table. The two tables may
#' use different delimiters.
#'
#' @importFrom utils read.csv
#' @export
#'
#' @examples
#' \donttest{
#' # Get path to the example data directory
#' txt_dir <- system.file("extdata/txt", package = "postvocs")
#'
#' # Construct the full path to a sample TXT file
#' txt_file <- file.path(txt_dir, "sample1.txt")
#'
#' # Process the file (default encoding UTF-8)
#' res <- process_gcms_txt(txt_file, debug = TRUE)
#'
#' # If encountering encoding warnings, try GBK
#' res <- process_gcms_txt(txt_file, encoding = "GBK")
#'
#' # Access the two tables using dynamically generated names
#' peak_table <- res[["sample1_PeakTable"]]
#' search_results <- res[["sample1_SearchResults"]]
#'
#' # See all names
#' names(res)
#' }
process_gcms_txt <- function(txt_file, debug = FALSE, encoding = "UTF-8") {

  # Helper: detect delimiter (tab or comma) from a line
  detect_separator <- function(line) {
    if (is.na(line) || !is.character(line) || nchar(trimws(line)) == 0) {
      stop("Column name line is empty or missing, cannot detect delimiter")
    }
    if (grepl("\t", line)) {
      return("\t")
    } else if (grepl(",", line)) {
      return(",")
    } else {
      stop("Cannot detect delimiter in column name line (neither tab nor comma)")
    }
  }

  # Read the entire file with user-specified encoding
  lines <- readLines(txt_file, warn = FALSE, skipNul = TRUE, encoding = encoding)
  if (debug) cat("Total lines in file:", length(lines), "\n")

  # ----- 1. Extract sample name from file name -----
  sample_name <- tools::file_path_sans_ext(basename(txt_file))
  if (debug) cat("Sample name derived from file:", sample_name, "\n")

  # ----- 2. Extract Peak Table -----
  peak_start <- grep("^\\[MC Peak Table\\]", lines)
  if (length(peak_start) == 0) stop("[MC Peak Table] not found")
  if (debug) cat("Peak table start line:", peak_start, "\n")

  n_peaks_line <- lines[peak_start + 1]
  # Auto-detect delimiter and extract the number (supports tab or comma)
  if (grepl("\t", n_peaks_line)) {
    n_peaks_str <- gsub(".*\\t", "", n_peaks_line)
  } else if (grepl(",", n_peaks_line)) {
    parts <- strsplit(n_peaks_line, ",")[[1]]
    n_peaks_str <- parts[length(parts)]
  } else {
    n_peaks_str <- n_peaks_line
  }
  n_peaks <- as.numeric(trimws(n_peaks_str))
  if (is.na(n_peaks)) {
    stop("Cannot parse number of peaks from line: ", n_peaks_line)
  }
  if (debug) cat("Declared number of peaks:", n_peaks, "\n")

  # Find column name line (starts with "Peak#")
  header_idx <- peak_start + 2
  while (header_idx <= length(lines)) {
    current_line <- lines[header_idx]
    if (is.na(current_line)) {
      header_idx <- header_idx + 1
      next
    }
    if (grepl("^Peak#", current_line)) break
    header_idx <- header_idx + 1
  }
  if (header_idx > length(lines)) stop("Peak table column name line not found")

  # Detect delimiter and obtain column names
  header_line <- lines[header_idx]
  sep <- detect_separator(header_line)
  col_names <- strsplit(header_line, sep, fixed = TRUE)[[1]]
  if (debug) cat("Peak table column names:", paste(col_names, collapse = " | "), "\n")

  # Data row range
  data_start <- header_idx + 1
  data_end <- data_start + n_peaks - 1
  if (data_end > length(lines)) {
    warning("Insufficient data rows, reading until end of file")
    data_end <- length(lines)
  }
  data_lines <- lines[data_start:data_end]
  data_lines <- data_lines[!grepl("^\\s*$", data_lines)]
  if (length(data_lines) == 0) stop("No valid data rows in Peak Table")

  data_list <- strsplit(data_lines, sep, fixed = TRUE)
  max_cols <- length(col_names)
  peak_data <- do.call(rbind, lapply(data_list, function(x) {
    if (length(x) < max_cols) x <- c(x, rep(NA, max_cols - length(x)))
    x[1:max_cols]
  }))
  peak_df <- as.data.frame(peak_data, stringsAsFactors = FALSE)
  names(peak_df) <- col_names

  # Convert numeric columns (non-numeric become NA)
  numeric_cols <- c("Peak#", "Ret.Time", "Area", "Height", "A/H", "Conc.", "Ret. Index")
  for (col in intersect(numeric_cols, names(peak_df))) {
    peak_df[[col]] <- suppressWarnings(as.numeric(peak_df[[col]]))
  }

  # ----- 3. Extract Search Results -----
  search_start <- grep("^\\[MS Similarity Search Results", lines)
  if (length(search_start) == 0) stop("[MS Similarity Search Results] not found")
  if (debug) cat("Search results start line:", search_start, "\n")

  # Find column name line (starts with "Spectrum#")
  search_header_idx <- search_start + 2
  while (search_header_idx <= length(lines)) {
    current_line <- lines[search_header_idx]
    if (is.na(current_line)) {
      search_header_idx <- search_header_idx + 1
      next
    }
    if (grepl("^Spectrum#", current_line)) break
    search_header_idx <- search_header_idx + 1
  }
  if (search_header_idx > length(lines)) stop("Search results column name line not found")

  search_header_line <- lines[search_header_idx]
  sep_search <- detect_separator(search_header_line)
  search_col_names <- strsplit(search_header_line, sep_search, fixed = TRUE)[[1]]
  if (debug) cat("Search results column names:", paste(search_col_names, collapse = " | "), "\n")

  # Skip empty lines after header
  search_data_start <- search_header_idx + 1
  while (search_data_start <= length(lines) && grepl("^\\s*$", lines[search_data_start])) {
    search_data_start <- search_data_start + 1
  }
  if (search_data_start > length(lines)) stop("No data rows in Search Results")

  # Read until end of file, filter empty lines
  search_lines <- lines[search_data_start:length(lines)]
  search_lines <- search_lines[!grepl("^\\s*$", search_lines)]
  if (length(search_lines) == 0) stop("No valid data rows in Search Results")

  search_list <- strsplit(search_lines, sep_search, fixed = TRUE)
  max_cols_search <- length(search_col_names)
  search_data <- do.call(rbind, lapply(search_list, function(x) {
    if (length(x) < max_cols_search) x <- c(x, rep(NA, max_cols_search - length(x)))
    x[1:max_cols_search]
  }))
  search_df <- as.data.frame(search_data, stringsAsFactors = FALSE)
  names(search_df) <- search_col_names

  # Convert Spectrum# to numeric
  if ("Spectrum#" %in% names(search_df)) {
    search_df$`Spectrum#` <- suppressWarnings(as.numeric(search_df$`Spectrum#`))
  }

  # ----- 4. Clean CAS numbers: remove all spaces -----
  # Find the CAS column (may be named "CAS..", "CAS.", or "CAS")
  cas_col <- grep("^CAS", names(search_df), value = TRUE)[1]
  if (!is.na(cas_col) && length(cas_col) > 0) {
    search_df[[cas_col]] <- gsub(" ", "", search_df[[cas_col]])
    if (debug) cat("Cleaned CAS column:", cas_col, "\n")
  } else {
    warning("CAS column not found in SearchResults, skipping space removal")
  }

  # ----- 5. Create dynamically named list and return -----
  peak_name <- paste0(sample_name, "_PeakTable")
  search_name <- paste0(sample_name, "_SearchResults")
  result <- list()
  result[[peak_name]] <- peak_df
  result[[search_name]] <- search_df

  if (debug) {
    cat("Extracted peak table with", nrow(peak_df), "rows and",
        ncol(peak_df), "columns.\n")
    cat("Extracted search results with", nrow(search_df), "rows and",
        ncol(search_df), "columns.\n")
  }

  message("Extracted data for sample: ", sample_name)
  return(result)
}
