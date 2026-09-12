#' Build a CAS × Sample abundance matrix from extracted peak areas
#'
#' Takes a list of extracted CAS–area vectors (as returned by
#' \code{extract_peak_areas}) and constructs a matrix where rows are CAS
#' compounds and columns are samples. This is the standard abundance matrix
#' used for subsequent annotation and screening.
#'
#' @param area_list A named list of numeric vectors, each vector representing
#'   one sample with CAS numbers as names and total peak areas as values.
#'   Typically the output of \code{extract_peak_areas()}.
#'
#' @return A data frame with the first column \code{"CAS"} and subsequent
#'   columns for each sample (in the order they appear in \code{area_list}).
#'   Missing values are filled with 0.
#'
#' @keywords postvocs
#'
#' @details
#' The function collects all unique CAS numbers from all samples, sorts them,
#' and creates a matrix with CAS as rows and samples as columns. If a sample
#' does not contain a particular CAS, its area is set to 0.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Get paths to example data
#' txt_dir <- system.file("extdata/txt", package = "postvocs")
#' sample_file <- system.file("extdata/SampleID.xlsx", package = "postvocs")
#'
#' # Batch process example files and extract peak areas
#' batch <- batch_process_gcms(txt_dir, sample_file)
#' areas <- extract_peak_areas(batch)
#'
#' # Build abundance matrix (CAS as rows, samples as columns)
#' abund <- build_cas_abundance(areas)
#'
#' # View first few rows and columns
#' head(abund[, 1:5])
#' }
build_cas_abundance <- function(area_list) {

  # ---- Input validation ----
  if (!is.list(area_list) || length(area_list) == 0) {
    stop("area_list must be a non-empty list of named numeric vectors.")
  }

  # Check that all elements are named numeric vectors
  for (i in seq_along(area_list)) {
    item <- area_list[[i]]
    if (!is.numeric(item) || is.null(names(item)) || any(names(item) == "")) {
      stop("Element ", i, " of area_list is not a named numeric vector.")
    }
  }

  # ---- Collect all CAS numbers ----
  all_cas <- unique(unlist(lapply(area_list, names)))
  all_cas <- sort(all_cas)

  # ---- Build matrix (samples as rows, CAS as columns) ----
  sample_names <- names(area_list)
  if (is.null(sample_names) || any(sample_names == "")) {
    # If list has no names, use generic names
    sample_names <- paste0("Sample", seq_along(area_list))
    names(area_list) <- sample_names
  }

  # Create matrix with samples as rows and CAS as columns
  mat <- matrix(0, nrow = length(area_list), ncol = length(all_cas),
                dimnames = list(sample_names, all_cas))

  for (sname in sample_names) {
    cas_vec <- area_list[[sname]]
    mat[sname, names(cas_vec)] <- cas_vec
  }

  # ---- Transpose to have CAS as rows, samples as columns ----
  mat_t <- t(mat)
  df <- as.data.frame(mat_t)
  df <- cbind(CAS = rownames(df), df)
  rownames(df) <- NULL

  # ---- Return ----
  return(df)
}
