#' Save GC-MS parsing results to CSV or Excel files
#'
#' Saves the data frames from \code{process_gcms_txt()} or
#' \code{batch_process_gcms()} to individual files. For SearchResults data
#' frames, the CAS column is prefixed with a single quote to prevent Excel
#' from interpreting CAS numbers as dates.
#'
#' @param results A list returned by \code{process_gcms_txt()} or
#'   \code{batch_process_gcms()}. For \code{process_gcms_txt}, the list
#'   contains two data frames (PeakTable and SearchResults). For
#'   \code{batch_process_gcms}, the list contains components \code{PeakTables}
#'   and \code{SearchResults}, each a named list of data frames.
#' @param output_dir Character. Path to the directory where files will be saved.
#'   Created if it does not exist.
#' @param format Character. Output format: \code{"csv"} or \code{"xlsx"}.
#'   Default is \code{"csv"}.
#' @param ... Additional arguments passed to \code{write.csv} (if format = "csv")
#'   or \code{openxlsx::writeData} (if format = "xlsx"). For CSV, common
#'   arguments include \code{row.names = FALSE}.
#'
#' @return Invisibly, a character vector of saved file paths.
#'
#' @keywords postvocs
#'
#' @details
#' The function detects data frames that contain CAS-related columns (based on
#' the data frame name containing "SearchResults") and adds a single quote
#' prefix to the CAS column to avoid Excel date conversion. For PeakTable
#' data frames, no modification is applied.
#'
#' If \code{format = "xlsx"}, the \code{openxlsx} package is required.
#' If not installed, the function will prompt to install it.
#'
#' @importFrom utils write.csv
#' @export
#'
#' @examples
#' \donttest{
#' # Get paths to example data
#' txt_dir <- system.file("extdata/txt", package = "postvocs")
#' txt_file <- file.path(txt_dir, "sample1.txt")
#' sample_file <- system.file("extdata/SampleID.xlsx", package = "postvocs")
#'
#' # (1) Process a single file and save as CSV
#' res <- process_gcms_txt(txt_file)
#' save_gcms_results(res, output_dir = tempdir(), format = "csv")
#'
#' # (2) Batch process and save as Excel (all in one file)
#' batch <- batch_process_gcms(txt_dir, sample_file)
#' save_gcms_results(batch, output_dir = tempdir(), format = "xlsx")
#'
#' # (3) Save only PeakTables (custom extraction)
#' save_gcms_results(batch$PeakTables, output_dir = tempdir(), format = "csv")
#' }
save_gcms_results <- function(results, output_dir, format = c("csv", "xlsx"), ...) {

  format <- match.arg(format)

  # ---- 1. Create output directory ----
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # ---- 2. Extract all data frames from results ----
  all_dfs <- list()

  if (is.list(results) && !is.null(names(results))) {
    if ("PeakTables" %in% names(results) && "SearchResults" %in% names(results)) {
      peak_list <- results[["PeakTables"]]
      search_list <- results[["SearchResults"]]
      if (is.list(peak_list)) all_dfs <- c(all_dfs, peak_list)
      if (is.list(search_list)) all_dfs <- c(all_dfs, search_list)
    } else {
      all_dfs <- results
    }
  } else if (is.list(results) && is.null(names(results))) {
    all_dfs <- results
  } else {
    stop("results must be a list of data frames or a batch_process_gcms result list.")
  }

  is_df <- vapply(all_dfs, is.data.frame, logical(1))
  if (!all(is_df)) {
    warning("Some elements in results are not data frames and will be skipped.")
    all_dfs <- all_dfs[is_df]
  }

  if (length(all_dfs) == 0) {
    stop("No data frames found to save.")
  }

  # ---- 3. Helper to add single quote to CAS column ----
  add_quote_to_cas <- function(df, df_name) {
    if (!grepl("SearchResults", df_name, ignore.case = TRUE)) {
      return(df)
    }
    cas_col <- grep("^CAS", names(df), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(cas_col)) {
      warning("No CAS column found in ", df_name, ", skipping quote addition.")
      return(df)
    }
    df[[cas_col]] <- as.character(df[[cas_col]])
    df[[cas_col]] <- paste0("'", df[[cas_col]])
    return(df)
  }

  # ---- 4. Save each data frame ----
  saved_files <- character()

  for (df_name in names(all_dfs)) {
    df <- all_dfs[[df_name]]
    df <- add_quote_to_cas(df, df_name)

    ext <- ifelse(format == "csv", ".csv", ".xlsx")
    file_path <- file.path(output_dir, paste0(df_name, ext))

    if (format == "csv") {
      utils::write.csv(df, file_path, row.names = FALSE, ...)
    } else {
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop("Package 'openxlsx' is required for xlsx output. Please install it.")
      }
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Sheet1")
      # Write data without any formatting (default Excel style)
      openxlsx::writeData(wb, "Sheet1", df, rowNames = FALSE, ...)
      openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
    }

    saved_files <- c(saved_files, file_path)
    message("Saved: ", file_path)
  }

  invisible(saved_files)
}
