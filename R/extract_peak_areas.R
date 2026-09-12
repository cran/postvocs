#' Extract CAS and total peak areas from GC-MS results
#'
#' This function extracts CAS numbers and their corresponding total peak areas
#' from GC-MS PeakTable and SearchResults data. It supports three input types:
#' \itemize{
#'   \item A list from \code{process_gcms_txt()} (single sample)
#'   \item A list from \code{batch_process_gcms()} (multiple samples)
#'   \item A folder path containing PeakTable and SearchResults files (CSV or Excel)
#' }
#'
#' @param x Either:
#'   \itemize{
#'     \item A list from \code{process_gcms_txt()} containing \code{PeakTable} and \code{SearchResults}
#'     \item A list from \code{batch_process_gcms()} containing \code{PeakTables} and \code{SearchResults}
#'     \item A character path to a folder containing PeakTable and SearchResults files
#'   }
#' @param ... Additional arguments passed to \code{read.csv} or \code{readxl::read_excel}
#'   (e.g., \code{sheet = 1} for Excel files).
#'
#' @return A named list. Each element is a named numeric vector with CAS as names
#'   and total peak area as values. The list names are sample names.
#'
#' @keywords postvocs
#'
#' @details
#' The function automatically detects the input type:
#' \itemize{
#'   \item If \code{x} is a character path to an existing folder, it scans for
#'         files containing "PeakTable" and "SearchResults" in their names,
#'         pairs them by common prefix, and processes all samples.
#'   \item If \code{x} is a list containing \code{PeakTable} and \code{SearchResults},
#'         it processes as a single sample and returns a list with one element.
#'   \item If \code{x} is a list containing \code{PeakTables} and \code{SearchResults},
#'         it processes all samples in batch.
#' }
#'
#' @importFrom dplyr %>% filter mutate select distinct group_by summarise
#' @importFrom rlang sym
#' @export
#'
#' @examples
#' \donttest{
#' # Get paths to example data
#' txt_file <- system.file("extdata/txt/sample1.txt", package = "postvocs")
#' txt_dir <- system.file("extdata/txt", package = "postvocs")
#' sample_file <- system.file("extdata/SampleID.xlsx", package = "postvocs")
#'
#' # 1. From process_gcms_txt result (single sample)
#' res <- process_gcms_txt(txt_file)
#' result_single <- extract_peak_areas(res)
#' # result_single is a list with one element named "sample1"
#'
#' # 2. From batch_process_gcms result (multiple samples)
#' batch <- batch_process_gcms(txt_dir, sample_file)
#' results_batch <- extract_peak_areas(batch)
#' # results_batch is a list with elements named by mapped sample names
#' # e.g., "Site1_TreatA_1", "Site1_TreatA_2", ...
#'
#' # 3. From a folder containing PeakTable and SearchResults CSV files
#' # (This requires you to have saved the results from save_gcms_results first)
#' # For example:
#' save_gcms_results(res, tempdir(), format = "csv")
#' results_folder <- extract_peak_areas(tempdir())
#' }
extract_peak_areas <- function(x, ...) {

  # ---- Case 1: x is a character path (folder mode) ----
  if (is.character(x) && length(x) == 1) {
    if (dir.exists(x)) {
      return(.extract_from_folder(x, ...))
    } else {
      stop("'", x, "' is not an existing directory. Please provide a valid folder path.")
    }
  }

  # ---- Case 2: x is a list ----
  if (is.list(x)) {
    # Check if it's a batch_process_gcms result (has PeakTables and SearchResults)
    if ("PeakTables" %in% names(x) && "SearchResults" %in% names(x)) {
      return(.extract_batch(x))
    }
    # Check if it's a process_gcms_txt result (has PeakTable and SearchResults)
    if (any(grepl("PeakTable", names(x), ignore.case = TRUE)) &&
        any(grepl("SearchResults", names(x), ignore.case = TRUE))) {
      return(.extract_from_df(x))
    }
  }

  stop("Invalid input. Please provide:\n",
       "  1. A list from process_gcms_txt() (single sample)\n",
       "  2. A list from batch_process_gcms() (batch processing)\n",
       "  3. A folder path containing PeakTable and SearchResults files")
}

# ---- Helper: Extract from folder (auto-detect) ----
.extract_from_folder <- function(folder_path, ...) {
  files <- list.files(folder_path, full.names = TRUE)
  if (length(files) == 0) {
    stop("No files found in folder: ", folder_path)
  }

  peak_files <- grep("PeakTable", files, ignore.case = TRUE, value = TRUE)
  search_files <- grep("SearchResults", files, ignore.case = TRUE, value = TRUE)

  if (length(peak_files) == 0 || length(search_files) == 0) {
    stop("No PeakTable or SearchResults files found in folder.\n",
         "Files must contain 'PeakTable' or 'SearchResults' in their names.")
  }

  peak_names <- gsub("_PeakTable\\.(csv|xlsx|xls)$", "", basename(peak_files), ignore.case = TRUE)
  search_names <- gsub("_SearchResults\\.(csv|xlsx|xls)$", "", basename(search_files), ignore.case = TRUE)

  common_names <- intersect(peak_names, search_names)

  if (length(common_names) == 0) {
    stop("No matching PeakTable/SearchResults pairs found in folder.")
  }

  message("Found ", length(common_names), " samples to process.")

  results <- list()

  for (sname in common_names) {
    peak_file <- grep(paste0(sname, "_PeakTable"), peak_files, ignore.case = TRUE, value = TRUE)[1]
    search_file <- grep(paste0(sname, "_SearchResults"), search_files, ignore.case = TRUE, value = TRUE)[1]

    peak_df <- .read_file(peak_file, ...)
    search_df <- .read_file(search_file, ...)

    vec <- .extract_core(peak_df, search_df)

    if (!is.null(vec) && length(vec) > 0) {
      results[[sname]] <- vec
    }
  }

  if (length(results) == 0) {
    warning("No samples were successfully processed.")
    return(NULL)
  }

  return(results)
}

# ---- Helper: Read file (auto-detect CSV or Excel) ----
.read_file <- function(file_path, ...) {
  ext <- tolower(tools::file_ext(file_path))
  if (ext == "csv") {
    return(utils::read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE, ...))
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read Excel files. Please install it.")
    }
    return(readxl::read_excel(file_path, ...))
  } else {
    stop("Unsupported file extension: ", ext, ". Please use .csv, .xlsx, or .xls")
  }
}

# ---- Helper: Extract from data frames (single sample from process_gcms_txt) ----
.extract_from_df <- function(input_list) {
  peak_idx <- grep("PeakTable", names(input_list), ignore.case = TRUE, value = TRUE)[1]
  search_idx <- grep("SearchResults", names(input_list), ignore.case = TRUE, value = TRUE)[1]

  if (is.na(peak_idx) || is.na(search_idx)) {
    stop("Could not find PeakTable and SearchResults in the provided list.")
  }

  # Extract sample name from the list element name
  sample_name <- gsub("_PeakTable$|_SearchResults$", "", peak_idx)

  peak_df <- input_list[[peak_idx]]
  search_df <- input_list[[search_idx]]

  vec <- .extract_core(peak_df, search_df)

  if (is.null(vec) || length(vec) == 0) {
    warning("No data extracted for sample: ", sample_name)
    return(NULL)
  }

  # Return as a named list with sample name
  result <- list(vec)
  names(result) <- sample_name
  return(result)
}

# ---- Helper: Batch extraction from batch_process_gcms result ----
.extract_batch <- function(batch_result) {
  peak_tables <- batch_result$PeakTables
  search_results <- batch_result$SearchResults

  if (length(peak_tables) != length(search_results)) {
    warning("Number of PeakTables and SearchResults do not match.")
  }

  result_list <- list()
  sample_names <- names(peak_tables)

  for (sname in sample_names) {
    core_name <- gsub("_PeakTable$", "", sname)
    search_name <- paste0(core_name, "_SearchResults")

    if (!search_name %in% names(search_results)) {
      warning("SearchResults not found for sample: ", core_name)
      next
    }

    vec <- .extract_core(peak_tables[[sname]], search_results[[search_name]])

    if (!is.null(vec) && length(vec) > 0) {
      result_list[[core_name]] <- vec
    }
  }

  if (length(result_list) == 0) {
    warning("No samples were successfully processed.")
    return(NULL)
  }

  return(result_list)
}

# ---- Core extraction logic (shared by all paths) ----
.extract_core <- function(peak_df, search_df) {

  # ---- Detect column names ----
  peak_col <- grep("^Peak", names(peak_df), ignore.case = TRUE, value = TRUE)[1]
  area_col <- grep("^Area", names(peak_df), ignore.case = TRUE, value = TRUE)[1]

  if (is.na(peak_col) || is.na(area_col)) {
    stop("Could not find Peak and Area columns in PeakTable")
  }

  spec_col <- grep("^Spectrum", names(search_df), ignore.case = TRUE, value = TRUE)[1]
  cas_col <- grep("^CAS", names(search_df), ignore.case = TRUE, value = TRUE)[1]

  if (is.na(spec_col) || is.na(cas_col)) {
    stop("Could not find Spectrum and CAS columns in SearchResults")
  }

  # ---- Clean SearchResults ----
  search_df[[cas_col]] <- as.character(search_df[[cas_col]])

  # Remove invalid CAS
  search_clean <- search_df %>%
    dplyr::filter(!.data[[cas_col]] %in% c("0-00-0", "00-00-00")) %>%
    dplyr::filter(.data[[cas_col]] != "" & !is.na(.data[[cas_col]]))

  if (nrow(search_clean) == 0) {
    return(structure(numeric(0), names = character(0)))
  }

  # Extract Spectrum-CAS mapping
  cas_spec <- search_clean %>%
    dplyr::select(Spectrum = !!rlang::sym(spec_col), CAS = !!rlang::sym(cas_col)) %>%
    dplyr::distinct(Spectrum, CAS, .keep_all = TRUE)

  # ---- Extract Peak-Area mapping ----
  peak_area <- peak_df %>%
    dplyr::select(Peak = !!rlang::sym(peak_col), Area = !!rlang::sym(area_col)) %>%
    dplyr::mutate(
      Peak = as.numeric(.data$Peak),
      Area = as.numeric(.data$Area)
    ) %>%
    dplyr::filter(!is.na(.data$Area) & .data$Area > 0)

  if (nrow(peak_area) == 0) {
    return(structure(numeric(0), names = character(0)))
  }

  # ---- Merge and summarize ----
  merged <- cas_spec %>%
    dplyr::left_join(peak_area, by = c("Spectrum" = "Peak")) %>%
    dplyr::filter(!is.na(.data$Area))

  if (nrow(merged) == 0) {
    return(structure(numeric(0), names = character(0)))
  }

  cas_total <- merged %>%
    dplyr::group_by(.data$CAS) %>%
    dplyr::summarise(TotalArea = sum(.data$Area, na.rm = TRUE), .groups = "drop")

  # ---- Return named vector ----
  result_vec <- cas_total$TotalArea
  names(result_vec) <- cas_total$CAS
  return(result_vec)
}
