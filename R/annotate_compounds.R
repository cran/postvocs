#' Annotate CAS numbers with chemical information
#'
#' This function takes a set of CAS numbers (from an abundance matrix) and
#' retrieves chemical names and properties from user-provided libraries or
#' from online databases (via `webchem`). It supports caching and automatic
#' rate-limiting to comply with PubChem API constraints.
#'
#' @param abundance_data Either:
#'   \itemize{
#'     \item A data.frame with a first column named "CAS"
#'     \item A character string path to a CSV file with the same structure
#'     \item A character string path to an Excel file (.xlsx or .xls) with the same structure
#'   }
#' @param lib_source Character. Source of chemical information:
#'   \itemize{
#'     \item `"user"` - only use `user_lib`; CAS not found become `NA`.
#'     \item `"webchem"` - only query online via `webchem`.
#'     \item `"auto"` - first try `user_lib`, then webchem for missing ones.
#'   }
#'   Default is `"auto"`.
#' @param user_lib Optional data.frame with at least two columns: `CAS` and
#'   `Name`. Used when `lib_source = "user"` or `"auto"`. Must be provided
#'   if `lib_source = "user"`.
#' @param webchem_config A list controlling online queries:
#'   \itemize{
#'     \item `rate_limit` - maximum requests per second (default 5).
#'     \item `chunk_size` - number of CAS per batch (default 100).
#'     \item `cache_file` - path to RDS cache file. If not provided and
#'           `use_cache = TRUE`, defaults to the user cache directory
#'           (e.g., `tools::R_user_dir("postvocs", "cache")/chem_cache.rds`).
#'   }
#' @param use_cache Logical. Whether to use caching. Default is `FALSE`.
#'   If `FALSE`, no cache file is read or written, and `force_retrieve`
#'   is ignored. If `TRUE`, caching is enabled; the cache file path is
#'   printed to the console.
#' @param force_retrieve Logical. If `TRUE` and `use_cache = TRUE`,
#'   ignore existing cache and query all CAS again. Default `FALSE`.
#'   Has no effect when `use_cache = FALSE`.
#'
#' @return A list with two components:
#'   \item{annotation}{A data.frame containing all CAS numbers and their
#'         retrieved information (ID, CAS, Name, MF, MW, IUPAC_Name, SMILES,
#'         InChIKey, InChI, QueryDate, Status, Source).}
#'   \item{abundance_updated}{The original abundance data.frame with two
#'         additional columns: `Compound_Name` and `Annotation_Source`,
#'         inserted after the `CAS` column.}
#'
#' @keywords postvocs
#'
#' @details
#'
#' The function automatically determines the batch strategy based on the number
#' of CAS to query:
#' \itemize{
#'   \item 1-50: direct individual queries (no chunking).
#'   \item 51-1000: chunked queries using `chunk_size`.
#'   \item 1001-5000: chunked queries with smaller chunk size (200) and
#'         extra delay.
#'   \item > 5000: a warning is issued recommending the PubChem PUG Download
#'         service; the function will still attempt chunked queries but may be
#'         slow or rate-limited.
#' }
#'
#' A persistent cache is stored at the location given by `cache_file` (or
#' the default user cache directory) only when `use_cache = TRUE`. The cache
#' is a named list where each key is a CAS number and the value is a list
#' of chemical properties, including a `status` field (`"success"` or
#' `"not_found"`). Cache is updated after each query when enabled.
#'
#' The `Annotation_Source` column indicates the origin of the compound name:
#' \itemize{
#'   \item `"user_lib"` - from the user‑provided library.
#'   \item `"web"` - newly queried from webchem in this session.
#'   \item `"web(cache)"` - retrieved from the cache (previous webchem query).
#' }
#'
#' CAS numbers are automatically cleaned before processing: leading/trailing
#' single quotes are removed and all spaces are stripped.
#'
#' @importFrom dplyr %>% select mutate left_join
#' @importFrom webchem get_cid pc_prop
#' @importFrom utils read.csv txtProgressBar setTxtProgressBar
#' @importFrom tools file_ext R_user_dir
#' @importFrom readxl read_excel
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' \donttest{
#' # Create a small abundance matrix
#' ab <- data.frame(
#'   CAS = c("64-17-5", "67-64-1", "75-05-8"),
#'   Sample1 = c(100, 200, 300),
#'   Sample2 = c(150, 250, 350)
#' )
#'
#' # Create a user-provided compound library
#' user_lib <- data.frame(
#'   CAS = c("64-17-5", "67-64-1"),
#'   Name = c("Ethanol", "Acetone"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Annotate using the user library (no network required)
#' res <- annotate_compounds(ab, lib_source = "user", user_lib = user_lib)
#'
#' # View annotation table
#' head(res$annotation)
#'
#' # Updated abundance matrix with Compound_Name and Annotation_Source
#' head(res$abundance_updated)
#' }
annotate_compounds <- function(
    abundance_data,
    lib_source = c("auto", "user", "webchem"),
    user_lib = NULL,
    webchem_config = list(
      rate_limit = 5,
      chunk_size = 100,
      cache_file = NULL
    ),
    use_cache = FALSE,
    force_retrieve = FALSE
) {

  # ---------- 1. Input validation and setup ----------
  lib_source <- match.arg(lib_source)

  # ---- Read input (supports data.frame, CSV, Excel) ----
  if (is.character(abundance_data) && length(abundance_data) == 1 && file.exists(abundance_data)) {
    ext <- tolower(tools::file_ext(abundance_data))
    if (ext == "csv") {
      ab <- utils::read.csv(abundance_data, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required to read Excel files. Please install it.")
      }
      ab <- readxl::read_excel(abundance_data)
      ab <- as.data.frame(ab, stringsAsFactors = FALSE)
    } else {
      stop("Unsupported file format: ", ext, ". Please use .csv, .xlsx, or .xls")
    }
  } else if (is.data.frame(abundance_data)) {
    ab <- abundance_data
  } else {
    stop("abundance_data must be a data.frame or a path to an existing CSV or Excel file.")
  }

  if (ncol(ab) < 2) stop("abundance_data must have at least two columns (CAS and one sample).")
  if (names(ab)[1] != "CAS") stop("The first column of abundance_data must be named 'CAS'.")

  # ---- Clean CAS numbers: remove quotes and spaces ----
  ab[[1]] <- as.character(ab[[1]])
  ab[[1]] <- gsub("^'|'$", "", ab[[1]])  # Remove leading/trailing single quotes
  ab[[1]] <- gsub(" ", "", ab[[1]])      # Remove all spaces

  # Extract CAS vector for processing
  cas_vec <- ab[[1]]
  cas_vec <- trimws(cas_vec)

  # ---- Remove NA or empty CAS ----
  if (any(cas_vec == "" | is.na(cas_vec))) {
    warning("Empty or NA CAS values found. They will be ignored in annotation.")
    valid_idx <- which(cas_vec != "" & !is.na(cas_vec))
    cas_vec <- cas_vec[valid_idx]
    ab <- ab[valid_idx, ]
  }

  # ---- Deduplicate CAS numbers ----
  n_total_before <- length(cas_vec)
  # Keep first occurrence
  keep_idx <- !duplicated(cas_vec)
  n_duplicates <- sum(!keep_idx)
  if (n_duplicates > 0) {
    cas_vec <- cas_vec[keep_idx]
    ab <- ab[keep_idx, ]
    # Reindex rows
    rownames(ab) <- NULL
  }
  n_total <- length(cas_vec)
  if (n_total == 0) stop("No valid CAS numbers to annotate.")

  # ---- Clean user_lib if provided ----
  if (!is.null(user_lib)) {
    if (!is.data.frame(user_lib)) stop("user_lib must be a data.frame.")
    if (!all(c("CAS", "Name") %in% names(user_lib))) {
      stop("user_lib must contain columns 'CAS' and 'Name' (case-sensitive).")
    }
    user_lib$CAS <- trimws(as.character(user_lib$CAS))
    user_lib$CAS <- gsub("^'|'$", "", user_lib$CAS)  # Remove quotes
    user_lib$CAS <- gsub(" ", "", user_lib$CAS)      # Remove spaces
    user_lib$Name <- as.character(user_lib$Name)
  }

  # ---- Prepare webchem config and cache settings ----
  rate_limit <- webchem_config$rate_limit %||% 5
  chunk_size <- webchem_config$chunk_size %||% 100

  # Determine cache file path (only if use_cache is TRUE)
  if (use_cache) {
    if (is.null(webchem_config$cache_file)) {
      cache_dir <- tools::R_user_dir("postvocs", "cache")
      cache_file <- file.path(cache_dir, "chem_cache.rds")
      if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    } else {
      cache_file <- webchem_config$cache_file
      cache_dir <- dirname(cache_file)
      if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    }
    message("Cache file: ", cache_file)
  } else {
    cache_file <- NULL
    message("Cache disabled.")
  }

  # ---------- 2. Handle user library (if lib_source = "user" or "auto") ----------
  user_lookup <- list()
  if (!is.null(user_lib)) {
    user_lookup <- setNames(user_lib$Name, user_lib$CAS)
  }

  # For "user" mode: only user_lib, no webchem
  if (lib_source == "user") {
    if (is.null(user_lib)) stop("lib_source = 'user' but user_lib is NULL.")
    # Safely match names using %in% to avoid subscript out of bounds
    names_found <- sapply(cas_vec, function(x) {
      if (x %in% names(user_lookup)) user_lookup[[x]] else NA_character_
    })
    sources <- ifelse(cas_vec %in% names(user_lookup), "user_lib", NA_character_)
    annotation <- data.frame(
      ID = seq_along(cas_vec),
      CAS = cas_vec,
      Name = names_found,
      MF = NA_character_,
      MW = NA_real_,
      IUPAC_Name = NA_character_,
      SMILES = NA_character_,
      InChIKey = NA_character_,
      InChI = NA_character_,
      QueryDate = Sys.Date(),
      Status = ifelse(!is.na(names_found), "Found", "Not Found"),
      Source = sources,
      stringsAsFactors = FALSE
    )
    ab_updated <- add_compound_name(ab, annotation)
    message("\n========== Annotation Summary ==========")
    message("Total CAS numbers (after dedup): ", n_total)
    if (n_duplicates > 0) message("Duplicates removed: ", n_duplicates)
    message("Found in user library: ", sum(!is.na(names_found)))
    message("Not found: ", sum(is.na(names_found)))
    return(list(annotation = annotation, abundance_updated = ab_updated))
  }

  # ---------- 3. Load cache (only if enabled) ----------
  cache <- list()
  if (!is.null(cache_file) && !force_retrieve && file.exists(cache_file)) {
    cache <- readRDS(cache_file)
    if (!is.list(cache)) cache <- list()
  }

  # Determine which CAS need querying
  if (lib_source == "auto" && !is.null(user_lib)) {
    found_in_user <- cas_vec %in% names(user_lookup)
    cas_to_query <- cas_vec[!found_in_user & !cas_vec %in% names(cache)]
  } else {
    cas_to_query <- cas_vec[!cas_vec %in% names(cache)]
  }

  cache_hits_before <- cas_vec[cas_vec %in% names(cache)]
  n_cache_hits <- length(cache_hits_before)
  user_hits <- if (!is.null(user_lib)) sum(cas_vec %in% names(user_lookup)) else 0

  if (length(cas_to_query) == 0) {
    message("All CAS numbers already have information. Skipping online queries.")
    annotation <- build_annotation_from_cache_and_user(cas_vec, cache, user_lookup, newly_queried = character(0))
    ab_updated <- add_compound_name(ab, annotation)
    message("\n========== Annotation Summary ==========")
    message("Total CAS numbers (after dedup): ", n_total)
    if (n_duplicates > 0) message("Duplicates removed: ", n_duplicates)
    message("Cache hits: ", n_cache_hits)
    if (!is.null(user_lib)) message("User library hits: ", user_hits)
    message("Not found (NA): ", sum(is.na(annotation$Name)))
    return(list(annotation = annotation, abundance_updated = ab_updated))
  }

  # ---------- 4. Query webchem ----------
  n_query <- length(cas_to_query)
  if (n_query > 5000) {
    warning(
      "More than 5000 CAS numbers to query. This may take a long time and ",
      "could be rate-limited. Consider using PubChem PUG Download service ",
      "or splitting the query into smaller batches."
    )
  }

  if (n_query <= 50) {
    batch_size <- 1
  } else if (n_query <= 1000) {
    batch_size <- chunk_size
  } else if (n_query <= 5000) {
    batch_size <- 200
  } else {
    batch_size <- 200
  }

  batches <- split(cas_to_query, ceiling(seq_along(cas_to_query) / batch_size))
  new_results <- list()

  message("\nStarting chemical annotation via webchem...")
  pb <- utils::txtProgressBar(min = 0, max = n_query, style = 3)

  for (i in seq_along(batches)) {
    batch <- batches[[i]]
    prop_list <- list()

    # Get CID for the batch
    cid_df <- tryCatch(
      webchem::get_cid(batch, from = "cas", match = "first"),
      error = function(e) {
        warning("get_cid failed for batch ", i, ": ", e$message)
        data.frame(query = batch, cid = NA_integer_, stringsAsFactors = FALSE)
      }
    )
    if (nrow(cid_df) == 0) {
      cid_df <- data.frame(query = batch, cid = NA_integer_, stringsAsFactors = FALSE)
    }

    # Separate found and not found
    found_idx <- which(!is.na(cid_df$cid))
    not_found_cas <- batch[is.na(cid_df$cid)]

    # For found, get properties
    if (length(found_idx) > 0) {
      cids <- cid_df$cid[found_idx]
      prop_df <- tryCatch(
        webchem::pc_prop(
          cids,
          properties = c("MolecularWeight", "MolecularFormula", "IUPACName",
                         "CanonicalSMILES", "InChI", "InChIKey")
        ),
        error = function(e) {
          warning("pc_prop failed for batch ", i, ": ", e$message)
          data.frame(cid = cids, stringsAsFactors = FALSE)
        }
      )
      if (!is.data.frame(prop_df) || nrow(prop_df) == 0) {
        prop_df <- data.frame(cid = cids, stringsAsFactors = FALSE)
      } else if (nrow(prop_df) != length(cids)) {
        full_prop <- data.frame(cid = cids, stringsAsFactors = FALSE)
        prop_df <- merge(full_prop, prop_df, by = "cid", all.x = TRUE)
      }

      col_names <- names(prop_df)
      map <- list(
        MW = grep("MolecularWeight|MW", col_names, ignore.case = TRUE, value = TRUE)[1],
        MF = grep("MolecularFormula|MF", col_names, ignore.case = TRUE, value = TRUE)[1],
        IUPAC = grep("IUPACName|IUPAC", col_names, ignore.case = TRUE, value = TRUE)[1],
        SMILES = grep("SMILES|CanonicalSMILES", col_names, ignore.case = TRUE, value = TRUE)[1],
        InChIKey = grep("InChIKey", col_names, ignore.case = TRUE, value = TRUE)[1],
        InChI = grep("InChI$", col_names, ignore.case = TRUE, value = TRUE)[1]
      )
      for (key in names(map)) if (is.na(map[[key]])) map[[key]] <- NA

      for (j in seq_len(nrow(prop_df))) {
        cas <- batch[found_idx[j]]
        if (is.na(cas) || nchar(cas) == 0) next
        prop <- prop_df[j, ]
        prop_list[[cas]] <- list(
          status = "success",
          Name = if (!is.na(map[["IUPAC"]])) prop[[map[["IUPAC"]]]] else NA_character_,
          MF = if (!is.na(map[["MF"]])) prop[[map[["MF"]]]] else NA_character_,
          MW = if (!is.na(map[["MW"]])) as.numeric(prop[[map[["MW"]]]]) else NA_real_,
          IUPAC_Name = if (!is.na(map[["IUPAC"]])) prop[[map[["IUPAC"]]]] else NA_character_,
          SMILES = if (!is.na(map[["SMILES"]])) prop[[map[["SMILES"]]]] else NA_character_,
          InChIKey = if (!is.na(map[["InChIKey"]])) prop[[map[["InChIKey"]]]] else NA_character_,
          InChI = if (!is.na(map[["InChI"]])) prop[[map[["InChI"]]]] else NA_character_,
          QueryDate = Sys.Date()
        )
      }
    }

    # For not found, fill with NA and status = "not_found"
    for (cas in not_found_cas) {
      if (is.na(cas) || nchar(cas) == 0) next
      prop_list[[cas]] <- list(
        status = "not_found",
        Name = NA_character_,
        MF = NA_character_,
        MW = NA_real_,
        IUPAC_Name = NA_character_,
        SMILES = NA_character_,
        InChIKey = NA_character_,
        InChI = NA_character_,
        QueryDate = Sys.Date()
      )
    }

    new_results <- c(new_results, prop_list)
    utils::setTxtProgressBar(pb, length(new_results))

    # Rate limiting
    if (i < length(batches)) {
      delay <- 1 / rate_limit * length(batch)
      Sys.sleep(delay)
    }
  }
  close(pb)

  # ---------- 5. Merge new results into cache (only if enabled) ----------
  if (!is.null(cache_file)) {
    # Ensure cache keys are unique (overwrite existing)
    for (key in names(new_results)) {
      cache[[key]] <- new_results[[key]]
    }
    saveRDS(cache, cache_file)
    message("\nCache saved to: ", cache_file)
  } else {
    message("\nCache not saved (disabled).")
  }

  # ---------- 6. Build annotation data.frame with source info ----------
  newly_queried <- cas_to_query
  annotation <- build_annotation_from_cache_and_user(cas_vec, cache, user_lookup, newly_queried)

  # ---------- 7. Add Compound_Name and Annotation_Source to abundance data ----------
  ab_updated <- add_compound_name(ab, annotation)

  # ---------- 8. Print summary ----------
  n_annotated <- sum(!is.na(annotation$Name))
  n_not_found <- sum(is.na(annotation$Name))
  n_new_queries <- length(cas_to_query)

  message("\n========== Annotation Summary ==========")
  message("Total CAS numbers (after dedup): ", n_total)
  if (n_duplicates > 0) message("Duplicates removed: ", n_duplicates)
  message("Cache hits (already known): ", n_cache_hits)
  if (!is.null(user_lib) && length(user_lookup) > 0) {
    message("User library hits: ", user_hits)
  }
  message("New webchem queries: ", n_new_queries)
  message("Successfully annotated (total): ", n_annotated)
  message("Not found (NA): ", n_not_found)
  if (!is.null(cache_file)) {
    message("Cache file saved to: ", cache_file)
  }

  # ---------- 9. Return ----------
  return(list(
    annotation = annotation,
    abundance_updated = ab_updated
  ))
}

# ---------- Helper functions ----------

#' Build annotation data.frame from cache, user library, and newly queried
#' @keywords internal
build_annotation_from_cache_and_user <- function(cas_vec, cache, user_lookup = list(), newly_queried = character(0)) {
  df <- data.frame(
    ID = seq_along(cas_vec),
    CAS = cas_vec,
    Name = NA_character_,
    MF = NA_character_,
    MW = NA_real_,
    IUPAC_Name = NA_character_,
    SMILES = NA_character_,
    InChIKey = NA_character_,
    InChI = NA_character_,
    QueryDate = as.Date(NA),
    Status = "Not Found",
    Source = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(cas_vec)) {
    cas <- cas_vec[i]
    # 1. Check user library
    if (length(user_lookup) > 0 && cas %in% names(user_lookup) && !is.na(user_lookup[[cas]])) {
      df[i, "Name"] <- user_lookup[[cas]]
      df[i, "Status"] <- "Found (User)"
      df[i, "Source"] <- "user_lib"
      next
    }
    # 2. Check cache
    if (cas %in% names(cache)) {
      info <- cache[[cas]]
      if (is.list(info)) {
        if (!is.null(info$status) && info$status == "not_found") {
          df[i, "Name"] <- NA_character_
          df[i, "Status"] <- "Not Found (Cache)"
          df[i, "Source"] <- "web(cache)"
        } else {
          df[i, "Name"] <- info$Name %||% NA_character_
          df[i, "MF"] <- info$MF %||% NA_character_
          df[i, "MW"] <- info$MW %||% NA_real_
          df[i, "IUPAC_Name"] <- info$IUPAC_Name %||% NA_character_
          df[i, "SMILES"] <- info$SMILES %||% NA_character_
          df[i, "InChIKey"] <- info$InChIKey %||% NA_character_
          df[i, "InChI"] <- info$InChI %||% NA_character_
          df[i, "QueryDate"] <- info$QueryDate %||% as.Date(NA)
          if (!is.na(df[i, "Name"])) {
            df[i, "Status"] <- "Found (Cache)"
          } else {
            df[i, "Status"] <- "Not Found (Cache)"
          }
          if (cas %in% newly_queried) {
            df[i, "Source"] <- "web"
          } else {
            df[i, "Source"] <- "web(cache)"
          }
        }
      }
    }
  }
  return(df)
}

#' Add Compound_Name and Annotation_Source columns to abundance data
#' @keywords internal
add_compound_name <- function(ab, annotation) {
  ab$CAS <- as.character(ab$CAS)
  annotation$CAS <- as.character(annotation$CAS)
  # Deduplicate annotation by CAS (keep first)
  annotation <- annotation[!duplicated(annotation$CAS), ]
  # Merge name and source
  ab_merged <- ab %>%
    dplyr::left_join(annotation[, c("CAS", "Name", "Source")], by = "CAS")
  # Reorder: place Compound_Name and Annotation_Source right after CAS
  cas_col <- which(names(ab_merged) == "CAS")
  name_col <- which(names(ab_merged) == "Name")
  source_col <- which(names(ab_merged) == "Source")
  if (length(name_col) == 1 && length(source_col) == 1) {
    ab_merged <- ab_merged[, c(1, name_col, source_col, setdiff(seq_along(ab_merged), c(1, name_col, source_col)))]
    names(ab_merged)[2] <- "Compound_Name"
    names(ab_merged)[3] <- "Annotation_Source"
  } else {
    if (!"Compound_Name" %in% names(ab_merged)) {
      ab_merged$Compound_Name <- NA_character_
    }
    if (!"Annotation_Source" %in% names(ab_merged)) {
      ab_merged$Annotation_Source <- NA_character_
    }
    cols <- c("CAS", "Compound_Name", "Annotation_Source", setdiff(names(ab_merged), c("CAS", "Compound_Name", "Annotation_Source")))
    ab_merged <- ab_merged[, cols]
  }
  return(ab_merged)
}

# Utility `%||%`
`%||%` <- function(x, y) if (is.null(x)) y else x
