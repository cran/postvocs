#' Save postvocs results to CSV or Excel files
#'
#' This function saves the results from `build_cas_abundance()`,
#' `annotate_compounds()`, or `filter_by_frequency()` to CSV or Excel files.
#' It automatically adds a single quote prefix to CAS columns to prevent
#' Excel date conversion. It checks whether the CAS column already contains
#' a leading single quote; if so, it does not add another one.
#'
#' For Excel output (`format = "xlsx"`), if the input is a list containing
#' multiple data frames (e.g., from `filter_by_frequency()` or
#' `annotate_compounds()`), all tables are saved as separate worksheets
#' in a single Excel file. For CSV output, each table is saved as an
#' independent file.
#'
#' @param x Either:
#'   \itemize{
#'     \item A data.frame (e.g., from `build_cas_abundance()`)
#'     \item A list from `annotate_compounds()` (contains `annotation` and
#'           `abundance_updated`)
#'     \item A list from `filter_by_frequency()` (contains `summary`,
#'           `abundance`, `round1`, `round2`, `summary_table`)
#'   }
#' @param output_dir Character. Directory where files will be saved.
#'   Created if it does not exist. **Must be specified** (no default).
#' @param prefix Character. Optional prefix for output file names.
#'   If `NULL`, a prefix is automatically derived from the input:
#'   \itemize{
#'     \item For data.frame: `"abundance"`
#'     \item For annotate_compounds: `"annotated"`
#'     \item For filter_by_frequency: `"filtered"`
#'   }
#' @param format Character. Output format: `"csv"` or `"xlsx"`.
#'   Default is `"csv"`.
#' @param what Character vector. Which components to save when `x` is a list.
#'   For `annotate_compounds` results, can be `"annotation"` and/or
#'   `"abundance_updated"` (default both). For `filter_by_frequency` results,
#'   can be any of `"summary"`, `"abundance"`, `"round1"`, `"round2"`,
#'   `"summary_table"` (default all except `summary` because it's just a
#'   small stats table). Set to `"all"` to save all components.
#' @param ... Additional arguments passed to `write.csv` (if `format = "csv"`)
#'   or `openxlsx::writeData` (if `format = "xlsx"`). Common arguments include
#'   `row.names = FALSE`.
#'
#' @return Invisibly, a character vector of saved file paths.
#'
#' @keywords postvocs
#'
#' @details
#' For data.frames that contain a column named "CAS", the function adds a
#' single quote prefix to each CAS value to prevent Excel from interpreting
#' them as dates. It does not add a second quote if the value already starts
#' with a single quote.
#'
#' If `format = "xlsx"`, the `openxlsx` package is required.
#'
#' @importFrom utils write.csv
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @export
#'
#' @examples
#' \donttest{
#' # Load pre-annotated example data
#' annotated <- readRDS(system.file("extdata/annotated_example.rds", package = "postvocs"))
#'
#' # Path to sample grouping file
#' sample_file <- system.file("extdata/SampleID.xlsx", package = "postvocs")
#'
#' # Perform frequency-based screening (if needed)
#' result <- filter_by_frequency(
#'   abundance_data = annotated,
#'   sample_group_file = sample_file,
#'   group_col = "Combined_Factor",
#'   blank_indicators = c("Factor1", "Factor2", "Combined_Factor")
#' )
#'
#' # Save results
#' save_postvocs_results(result, output_dir = tempdir(), format = "csv")
#' }
save_postvocs_results <- function(
    x,
    output_dir,
    prefix = NULL,
    format = c("csv", "xlsx"),
    what = NULL,
    ...
) {

  format <- match.arg(format)

  # Check that output_dir is provided
  if (missing(output_dir) || is.null(output_dir)) {
    stop("`output_dir` must be specified. Use `tempdir()` for temporary files.")
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # ---- Helper: add single quote to CAS column (if not already present) ----
  protect_cas <- function(df) {
    if (!is.data.frame(df)) return(df)
    if ("CAS" %in% names(df)) {
      df$CAS <- as.character(df$CAS)
      has_quote <- grepl("^'", df$CAS)
      df$CAS[!has_quote] <- paste0("'", df$CAS[!has_quote])
    }
    return(df)
  }

  # ---- Helper: transpose summary table ----
  transpose_summary <- function(df) {
    if (!is.data.frame(df)) return(df)
    if (nrow(df) == 1) {
      df_t <- as.data.frame(t(df))
      colnames(df_t) <- "Value"
      df_t$Metric <- rownames(df_t)
      rownames(df_t) <- NULL
      df_t <- df_t[, c("Metric", "Value")]
      return(df_t)
    } else {
      return(df)
    }
  }

  # ---- Helper: save a single data frame to CSV ----
  save_csv <- function(df, base_name, ...) {
    df <- protect_cas(df)
    file_path <- file.path(output_dir, paste0(base_name, ".csv"))
    utils::write.csv(df, file_path, row.names = FALSE, ...)
    return(file_path)
  }

  # ---- Helper: save multiple data frames to a single Excel file ----
  save_xlsx_multi <- function(df_list, base_name, sheet_names = NULL, ...) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for xlsx output.")
    }

    df_list <- lapply(df_list, protect_cas)

    wb <- openxlsx::createWorkbook()

    if (is.null(sheet_names)) {
      sheet_names <- names(df_list)
      if (is.null(sheet_names)) {
        sheet_names <- paste0("Sheet", seq_along(df_list))
      }
    }

    for (i in seq_along(df_list)) {
      sheet_name <- sheet_names[i]
      sheet_name <- gsub("[\\[\\]\\*\\?/:]", "_", sheet_name)
      if (nchar(sheet_name) > 31) sheet_name <- substr(sheet_name, 1, 31)

      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, df_list[[i]], rowNames = FALSE, ...)
    }

    file_path <- file.path(output_dir, paste0(base_name, ".xlsx"))
    openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
    return(file_path)
  }

  saved_files <- character()

  # ---- Case 1: x is a data.frame ----
  if (is.data.frame(x)) {
    if (is.null(prefix)) prefix <- "abundance"
    if (format == "csv") {
      file_path <- save_csv(x, prefix, ...)
    } else {
      file_path <- save_xlsx_multi(list(x), prefix, sheet_names = "Sheet1", ...)
    }
    saved_files <- c(saved_files, file_path)
    message("Saved: ", file_path)
    return(invisible(saved_files))
  }

  # ---- Case 2: x is a list ----
  if (!is.list(x)) {
    stop("x must be a data.frame or a list.")
  }

  # ---- Detect list type ----
  is_annotate <- "annotation" %in% names(x) && "abundance_updated" %in% names(x)
  is_filter <- all(c("summary", "abundance", "round1", "round2", "summary_table") %in% names(x))

  if (is_annotate && is_filter) {
    is_annotate <- TRUE
    is_filter <- FALSE
  }

  if (is_annotate) {
    if (is.null(prefix)) prefix <- "annotated"
    available <- c("annotation", "abundance_updated")
    if (is.null(what)) what <- available
    what <- intersect(what, available)
    if (length(what) == 0) {
      warning("No valid components to save for annotate_compounds result.")
      return(invisible(character(0)))
    }

    df_list <- x[what]
    sheet_names <- what

    if (format == "csv") {
      for (comp in what) {
        base_name <- paste0(prefix, "_", comp)
        file_path <- save_csv(x[[comp]], base_name, ...)
        saved_files <- c(saved_files, file_path)
        message("Saved: ", file_path)
      }
    } else {
      base_name <- paste0(prefix, "_all")
      file_path <- save_xlsx_multi(df_list, base_name, sheet_names = sheet_names, ...)
      saved_files <- c(saved_files, file_path)
      message("Saved: ", file_path)
    }
    return(invisible(saved_files))
  }

  if (is_filter) {
    if (is.null(prefix)) prefix <- "filtered"
    available <- c("summary", "abundance", "round1", "round2", "summary_table")
    if (is.null(what)) {
      what <- setdiff(available, "summary")
    } else if (length(what) == 1 && what == "all") {
      what <- available
    }
    what <- intersect(what, available)
    if (length(what) == 0) {
      warning("No valid components to save for filter_by_frequency result.")
      return(invisible(character(0)))
    }

    df_list <- list()
    sheet_names <- character()
    for (comp in what) {
      if (comp %in% names(x)) {
        df <- x[[comp]]
        if (comp == "summary") {
          df <- transpose_summary(df)
        }
        df_list[[comp]] <- df
        sheet_names <- c(sheet_names, comp)
      }
    }

    if (format == "csv") {
      for (i in seq_along(what)) {
        comp <- what[i]
        base_name <- paste0(prefix, "_", comp)
        file_path <- save_csv(df_list[[comp]], base_name, ...)
        saved_files <- c(saved_files, file_path)
        message("Saved: ", file_path)
      }
    } else {
      base_name <- paste0(prefix, "_all")
      file_path <- save_xlsx_multi(df_list, base_name, sheet_names = sheet_names, ...)
      saved_files <- c(saved_files, file_path)
      message("Saved: ", file_path)
    }
    return(invisible(saved_files))
  }

  # ---- Case 3: General list (user-defined) ----
  if (is.null(prefix)) prefix <- "result"

  df_list <- x[vapply(x, is.data.frame, logical(1))]
  if (length(df_list) == 0) {
    warning("No data.frames found in the list to save.")
    return(invisible(character(0)))
  }

  sheet_names <- names(df_list)
  if (is.null(sheet_names)) {
    sheet_names <- paste0("Sheet", seq_along(df_list))
  }

  if (format == "csv") {
    for (i in seq_along(df_list)) {
      base_name <- paste0(prefix, "_", sheet_names[i])
      file_path <- save_csv(df_list[[i]], base_name, ...)
      saved_files <- c(saved_files, file_path)
      message("Saved: ", file_path)
    }
  } else {
    base_name <- paste0(prefix, "_all")
    file_path <- save_xlsx_multi(df_list, base_name, sheet_names = sheet_names, ...)
    saved_files <- c(saved_files, file_path)
    message("Saved: ", file_path)
  }

  invisible(saved_files)
}
