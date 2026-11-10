#' Filter compounds by occurrence frequency across samples and treatment groups
#'
#' This function performs a two‑step filtering of compounds from an abundance
#' matrix. It optionally removes compounds detected in blank samples, then
#' retains compounds based on occurrence frequency across all samples and
#' within treatment groups. It returns a detailed summary table with flags for
#' each compound.
#'
#' @param abundance_data Either a data.frame with a first column named "CAS" and
#'   a second column named "Compound_Name", a character path to a CSV or Excel
#'   file with the same structure, or a list returned by \code{annotate_compounds()}
#'   (containing \code{abundance_updated}).
#' @param sample_group_file Either a data.frame or a character path to the sample
#'   grouping file (CSV or Excel). The first column must be "FileID" (ignored)
#'   and the second column must be "SampleName" (used to match abundance columns).
#' @param sheet Integer or character. Sheet number or name to read from Excel
#'   `sample_group_file`. Default is `1` (first sheet). Ignored if
#'   `sample_group_file` is a data.frame or CSV file.
#' @param group_col Character. Name of the column in `sample_group_file` that
#'   defines treatment groups. Default is `"Combined_Treatment"`.
#' @param threshold_total Numeric. Minimum fraction of total samples in which a
#'   compound must appear to pass the first filter. Default is `0.25` (25 percent).
#' @param threshold_treatment Numeric. Minimum fraction of samples within a
#'   treatment group for the second filter. Default is `0.25` (25 percent).
#' @param blank_indicators Character vector. Columns in `sample_group_file`
#'   that must all equal "0" to identify blank samples. Default is
#'   `c("Species", "Treatment", "Combined_Treatment")`. If some columns are
#'   not present, they are skipped with a warning.
#' @param remove_blank_compounds Logical. If `TRUE`, compounds detected
#'   (area > 0) in any blank sample are removed from all samples. If `FALSE`
#'   (default), blank samples are still excluded from frequency calculations,
#'   but compounds detected in blanks are retained in the abundance matrix.
#'
#' @return A list with five components:
#'   \describe{
#'     \item{summary}{A data frame with filtering statistics.}
#'     \item{abundance}{The full abundance matrix with added frequency columns
#'           and a flag column indicating whether each compound was detected in
#'           blank samples.}
#'     \item{round1}{Data frame of compounds kept in round 1.}
#'     \item{round2}{Data frame of compounds kept in round 2.}
#'     \item{summary_table}{A comprehensive table for all original compounds
#'           with flags.}
#'   }
#'
#' @keywords postvocs
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr filter select mutate left_join group_by summarise
#' @importFrom dplyr bind_rows distinct pull if_all any_of all_of
#' @importFrom dplyr rowwise ungroup c_across
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom utils read.csv head
#' @importFrom readxl read_excel
#' @importFrom rlang sym
#' @importFrom tools file_ext
#' @export
filter_by_frequency <- function(
    abundance_data,
    sample_group_file,
    sheet = 1,
    group_col = "Combined_Treatment",
    threshold_total = 0.25,
    threshold_treatment = 0.25,
    blank_indicators = c("Species", "Treatment", "Combined_Treatment"),
    remove_blank_compounds = FALSE
) {

  # ---------- 1. Input validation and data loading ----------

  if (is.list(abundance_data) && "abundance_updated" %in% names(abundance_data)) {
    abun <- abundance_data$abundance_updated
    if (!is.data.frame(abun)) stop("abundance_updated is not a data.frame.")
  } else if (is.character(abundance_data) && length(abundance_data) == 1 && file.exists(abundance_data)) {
    ext <- tolower(tools::file_ext(abundance_data))
    if (ext == "csv") {
      abun <- utils::read.csv(abundance_data, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required.")
      abun <- as.data.frame(readxl::read_excel(abundance_data), stringsAsFactors = FALSE)
    } else {
      stop("Unsupported file format: ", ext)
    }
  } else if (is.data.frame(abundance_data)) {
    abun <- abundance_data
  } else {
    stop("abundance_data must be a list, data.frame, or file path.")
  }

  if (ncol(abun) < 2) stop("abundance_data must have at least two columns.")
  if (names(abun)[1] != "CAS") stop("First column must be 'CAS'.")
  if (names(abun)[2] != "Compound_Name") {
    warning("Second column is not 'Compound_Name'; adding NA column.")
    abun <- abun %>% dplyr::mutate(Compound_Name = NA_character_, .after = 1)
  }

  # Remove pre-existing frequency columns
  freq_cols_to_remove <- grep("^(N_total|Freq_total|Freq_group_)", names(abun), value = TRUE)
  if (length(freq_cols_to_remove) > 0) {
    abun <- abun %>% dplyr::select(-dplyr::all_of(freq_cols_to_remove))
  }

  has_source <- "Annotation_Source" %in% names(abun)
  exclude_cols <- c("CAS", "Compound_Name")
  if (has_source) exclude_cols <- c(exclude_cols, "Annotation_Source")
  sample_cols <- names(abun)[!names(abun) %in% exclude_cols]
  if (length(sample_cols) == 0) stop("No sample columns found.")

  original_abun <- abun
  original_cas <- abun$CAS
  original_compound <- abun$Compound_Name
  original_source <- if (has_source) abun$Annotation_Source else NA_character_

  # ---------- 2. Read sample_group_file ----------

  if (is.data.frame(sample_group_file)) {
    group_info <- sample_group_file
  } else if (is.character(sample_group_file) && length(sample_group_file) == 1) {
    if (!file.exists(sample_group_file)) stop("sample_group_file not found.")
    ext <- tolower(tools::file_ext(sample_group_file))
    if (ext == "csv") {
      group_info <- utils::read.csv(sample_group_file, stringsAsFactors = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required.")
      group_info <- as.data.frame(readxl::read_excel(sample_group_file, sheet = sheet), stringsAsFactors = FALSE)
    } else {
      stop("Unsupported file format: ", ext)
    }
  } else {
    stop("sample_group_file must be a data.frame or file path.")
  }

  if (ncol(group_info) < 2) stop("group file must have at least two columns.")
  if (names(group_info)[1] != "FileID") warning("First column should be 'FileID' (ignored).")
  if (names(group_info)[2] != "SampleName") warning("Second column should be 'SampleName'.")
  group_info$SampleName <- trimws(as.character(group_info[[names(group_info)[2]]]))

  if (!group_col %in% names(group_info)) {
    stop("Column '", group_col, "' not found in group file.")
  }

  present_indicators <- intersect(blank_indicators, names(group_info))
  if (length(present_indicators) == 0) {
    warning("None of the blank_indicators found. Skipping blank detection.")
    blank_indicators <- NULL
  } else if (length(present_indicators) < length(blank_indicators)) {
    warning("Some blank_indicators missing. Using: ", paste(present_indicators, collapse = ", "))
    blank_indicators <- present_indicators
  }

  # ---------- 3. Trim whitespace ----------

  names(abun)[names(abun) %in% sample_cols] <- trimws(names(abun)[names(abun) %in% sample_cols])
  group_info$SampleName <- trimws(group_info$SampleName)
  orig_cols <- sample_cols

  # ---------- 4. Detect blank samples and flag compounds ----------

  blank_removed_cas <- character(0)
  blank_samples_actual <- character(0)
  in_blank_cas <- character(0)

  if (!is.null(blank_indicators) && length(blank_indicators) > 0) {
    blank_candidates <- group_info %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(blank_indicators), ~ . == "0")) %>%
      dplyr::pull("SampleName")

    if (length(blank_candidates) > 0) {
      blank_samples_actual <- intersect(blank_candidates, orig_cols)

      if (length(blank_samples_actual) > 0) {
        # CAS detected in any blank sample
        in_blank_cas <- abun %>%
          dplyr::filter(dplyr::if_any(dplyr::all_of(blank_samples_actual), ~ . > 0)) %>%
          dplyr::pull("CAS")

        if (remove_blank_compounds && length(in_blank_cas) > 0) {
          n_before <- nrow(abun)
          abun <- abun %>% dplyr::filter(!.data$CAS %in% in_blank_cas)
          blank_removed_cas <- in_blank_cas
          n_removed <- n_before - nrow(abun)
          message("Removed ", n_removed, " compounds that appeared in blank samples.")
        } else if (length(in_blank_cas) > 0) {
          message("Detected ", length(in_blank_cas),
                  " compounds in blank samples (retained; remove_blank_compounds = FALSE).")
        } else {
          message("No compounds found in blank samples.")
        }
      } else {
        message("Blank samples not found in abundance matrix.")
      }
    } else {
      message("No blank samples identified (all indicator columns equal '0').")
    }
  }

  # ---------- 5. Build flag vector (based on original CAS) ----------

  in_blank_flag <- ifelse(original_cas %in% in_blank_cas, "*", "")

  # ---------- 6. Remove blank samples from frequency calculation ----------

  # Keep blank sample columns in "abun" for abundance output
  abun_with_blank <- abun  # for abundance output (with blank columns)

  # For frequency calculation, remove blank columns
  abun_freq_base <- abun
  if (length(blank_samples_actual) > 0) {
    abun_freq_base <- abun_freq_base[, !(names(abun_freq_base) %in% blank_samples_actual), drop = FALSE]
    group_info <- group_info %>% dplyr::filter(!.data$SampleName %in% blank_samples_actual)
  }

  sample_cols_freq <- names(abun_freq_base)[!names(abun_freq_base) %in% exclude_cols]
  if (length(sample_cols_freq) == 0) stop("No non-blank sample columns remain.")

  group_info <- group_info %>% dplyr::filter(.data$SampleName %in% sample_cols_freq)
  if (nrow(group_info) == 0) {
    stop("No non-blank samples match the group file.\n",
         "Abundance columns: ", paste(head(sample_cols_freq, 10), collapse = ", "), "\n",
         "Group SampleNames: ", paste(head(group_info$SampleName, 10), collapse = ", "))
  }

  # ---------- 7. Calculate frequencies (based on non-blank samples) ----------

  abun_long <- abun_freq_base %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(sample_cols_freq),
      names_to = "Sample", values_to = "Area"
    ) %>%
    dplyr::left_join(group_info, by = c("Sample" = "SampleName"))

  freq_total <- abun_long %>%
    dplyr::group_by(.data$CAS, .data$Compound_Name) %>%
    dplyr::summarise(N_total = sum(.data$Area > 0, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Total_Samples = length(sample_cols_freq),
                  Freq_total = .data$N_total / .data$Total_Samples) %>%
    dplyr::select("CAS", "N_total", "Freq_total")

  freq_group <- abun_long %>%
    dplyr::group_by(.data$CAS, .data$Compound_Name, !!rlang::sym(group_col)) %>%
    dplyr::summarise(N_group = sum(.data$Area > 0, na.rm = TRUE),
                     Total_group = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(Freq_group = .data$N_group / .data$Total_group) %>%
    dplyr::select("CAS", "Freq_group", Treatment = !!rlang::sym(group_col))

  group_wide <- freq_group %>%
    tidyr::pivot_wider(names_from = "Treatment", values_from = "Freq_group",
                       values_fill = 0, names_prefix = "Freq_group_")

  # Merge frequencies into abun_freq_base
  abun_freq <- abun_freq_base %>%
    dplyr::left_join(freq_total, by = "CAS") %>%
    dplyr::left_join(group_wide, by = "CAS", suffix = c("", ".group")) %>%
    dplyr::select(-dplyr::ends_with(".group"))

  # Also attach frequencies to abun_with_blank (for abundance output)
  abun_with_blank <- abun_with_blank %>%
    dplyr::left_join(freq_total, by = "CAS") %>%
    dplyr::left_join(group_wide, by = "CAS", suffix = c("", ".group")) %>%
    dplyr::select(-dplyr::ends_with(".group"))

  # Add flag column
  abun_with_blank <- abun_with_blank %>%
    dplyr::mutate(In_Blank = in_blank_flag[match(.data$CAS, original_cas)], .after = if (has_source) 3 else 2)

  abun_freq <- abun_freq %>%
    dplyr::mutate(In_Blank = in_blank_flag[match(.data$CAS, original_cas)], .after = if (has_source) 3 else 2)

  # ---------- 8. Reorder columns ----------

  freq_cols <- c("N_total", "Freq_total", grep("^Freq_group_", names(abun_freq), value = TRUE))

  # abundance: CAS | Compound_Name | Annotation_Source | In_Blank | blank samples | freq | real samples
  base_cols <- c("CAS", "Compound_Name")
  if (has_source) base_cols <- c(base_cols, "Annotation_Source")
  base_cols <- c(base_cols, "In_Blank")

  real_sample_cols <- sample_cols_freq
  blank_sample_cols <- blank_samples_actual[blank_samples_actual %in% names(abun_with_blank)]

  new_order_abun <- c(base_cols, blank_sample_cols, freq_cols, real_sample_cols)
  new_order_abun <- intersect(new_order_abun, names(abun_with_blank))
  abun_with_blank <- abun_with_blank[, new_order_abun]

  # round1/round2: CAS | Compound_Name | Annotation_Source | In_Blank | freq | real samples
  new_order_freq <- c(base_cols, freq_cols, real_sample_cols)
  new_order_freq <- intersect(new_order_freq, names(abun_freq))

  # ---------- 9. Screening ----------

  round1_cas <- abun_freq %>%
    dplyr::filter(.data$Freq_total >= threshold_total) %>%
    dplyr::pull("CAS")
  round1 <- abun_freq %>% dplyr::filter(.data$CAS %in% round1_cas)

  group_freq_cols <- grep("^Freq_group_", names(abun_freq), value = TRUE)

  if (length(group_freq_cols) > 0) {
    round2_candidates <- abun_freq %>%
      dplyr::filter(.data$Freq_total < threshold_total) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(has_group_ok = any(dplyr::c_across(dplyr::all_of(group_freq_cols)) >= threshold_treatment)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(.data$has_group_ok) %>%
      dplyr::select(-"has_group_ok")
  } else {
    round2_candidates <- abun_freq[0, ]
  }

  round2_cas <- round2_candidates %>% dplyr::pull("CAS")
  round2 <- dplyr::bind_rows(round1, round2_candidates) %>%
    dplyr::distinct(.data$CAS, .keep_all = TRUE)

  # Reorder round1/round2
  round1 <- round1[, intersect(new_order_freq, names(round1))]
  round2 <- round2[, intersect(new_order_freq, names(round2))]

  # ---------- 10. Build summary_table ----------

  summary_df <- data.frame(
    CAS = original_cas,
    Compound_Name = original_compound,
    stringsAsFactors = FALSE
  )

  if (has_source) {
    summary_df$Annotation_Source <- original_source
  }

  # Flag column for FALSE mode
  if (!remove_blank_compounds) {
    summary_df$In_Blank <- in_blank_flag
  }

  # Flags
  summary_df$Removed_in_Blank <- ifelse(summary_df$CAS %in% blank_removed_cas, "-", "")
  summary_df$Removed_in_Round1 <- ifelse(summary_df$CAS %in% round1_cas, "", "-")
  summary_df$Added_in_Round2 <- ifelse(summary_df$CAS %in% round2_cas & !summary_df$CAS %in% round1_cas, "+", "")

  # Frequencies: for removed CAS, compute from original data
  if (remove_blank_compounds && length(blank_removed_cas) > 0) {
    # Compute frequencies for removed CAS from original data (excluding blank samples)
    original_freq_base <- original_abun
    if (length(blank_samples_actual) > 0) {
      original_freq_base <- original_freq_base[, !(names(original_freq_base) %in% blank_samples_actual), drop = FALSE]
    }
    original_sample_cols <- names(original_freq_base)[!names(original_freq_base) %in% exclude_cols]

    removed_long <- original_freq_base %>%
      dplyr::filter(.data$CAS %in% blank_removed_cas) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(original_sample_cols),
                          names_to = "Sample", values_to = "Area") %>%
      dplyr::left_join(group_info, by = c("Sample" = "SampleName"))

    removed_freq_total <- removed_long %>%
      dplyr::group_by(.data$CAS) %>%
      dplyr::summarise(N_total = sum(.data$Area > 0, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(Total_Samples = length(original_sample_cols),
                    Freq_total = .data$N_total / .data$Total_Samples)

    removed_freq_group <- removed_long %>%
      dplyr::group_by(.data$CAS, !!rlang::sym(group_col)) %>%
      dplyr::summarise(N_group = sum(.data$Area > 0, na.rm = TRUE),
                       Total_group = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(Freq_group = .data$N_group / .data$Total_group) %>%
      dplyr::select("CAS", "Freq_group", Treatment = !!rlang::sym(group_col))

    removed_group_wide <- removed_freq_group %>%
      tidyr::pivot_wider(names_from = "Treatment", values_from = "Freq_group",
                         values_fill = 0, names_prefix = "Freq_group_")

    removed_freq <- removed_freq_total %>%
      dplyr::left_join(removed_group_wide, by = "CAS")
  } else {
    removed_freq <- data.frame(CAS = character(0), N_total = numeric(0), Freq_total = numeric(0))
  }

  # Merge frequencies: use computed for removed CAS, and abun_freq for the rest
  freq_all <- dplyr::bind_rows(
    abun_freq %>% dplyr::select("CAS", "N_total", "Freq_total", dplyr::any_of(grep("^Freq_group_", names(abun_freq), value = TRUE))),
    removed_freq
  ) %>% dplyr::distinct(.data$CAS, .keep_all = TRUE)

  summary_df <- summary_df %>%
    dplyr::left_join(freq_all, by = "CAS")

  # Reorder summary_df
  base_cols_summary <- c("CAS", "Compound_Name")
  if (has_source) base_cols_summary <- c(base_cols_summary, "Annotation_Source")
  if (!remove_blank_compounds) base_cols_summary <- c(base_cols_summary, "In_Blank")
  flag_cols <- c("Removed_in_Blank", "Removed_in_Round1", "Added_in_Round2")
  order_cols <- c(base_cols_summary, freq_cols, flag_cols)
  order_cols <- intersect(order_cols, names(summary_df))
  summary_df <- summary_df[, order_cols]

  # ---------- 11. Summary statistics ----------

  summary_stats <- data.frame(
    Total_compounds_original = nrow(original_abun),
    Removed_in_blank = length(blank_removed_cas),
    Round1_kept = nrow(round1),
    Round2_additional = nrow(round2) - nrow(round1),
    Final_kept = nrow(round2)
  )

  message("\n========== Filtering Summary ==========")
  message("Total compounds (original): ", nrow(original_abun))
  message("Blank compound removal mode: ", ifelse(remove_blank_compounds, "removed", "retained"))
  message("Removed in blank: ", length(blank_removed_cas))
  message("Compounds after blank removal: ", nrow(abun))
  message("Samples used for frequency calculation: ", length(sample_cols_freq))
  message("Round 1 kept: ", nrow(round1))
  message("Round 2 additional: ", nrow(round2) - nrow(round1))
  message("Final kept: ", nrow(round2))

  invisible(list(
    summary = summary_stats,
    abundance = abun_with_blank,
    round1 = round1,
    round2 = round2,
    summary_table = summary_df
  ))
}
