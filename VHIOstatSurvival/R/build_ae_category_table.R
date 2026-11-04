#' Count adverse events and patients per category and dose cohort
#'
#' For each variable in `vars` (e.g. system organ class, preferred term) the
#' number of AEs (`*_AE`) and distinct patients (`*_n`) is counted overall
#' (`"All"`) and within every level of the cohort column.
#'
#' @param ae Data frame with one row per adverse event.
#' @param vars Named character vector. Names are the labels shown in the
#'   table, values are column names in `ae`, e.g.
#'   `c("System organ class" = "soc")`.
#' @param cohort_col Character. Cohort column; **must be a factor** (its
#'   levels define the strata order). Default `"dose_cohort"`.
#' @param id_col Character. Patient identifier column. Default `"record_id"`.
#' @param include_variable_total Logical. Add a `"Total"` row per variable
#'   block? Default `FALSE`.
#'
#' @return A tibble with columns `variable`, `is_total`, `level` and, for each
#'   stratum `s` in `c("All", levels(cohort))`, `s_AE` and `s_n`. Input for
#'   [render_ae_table()].
#'
#' @seealso [render_ae_table()]
#'
#' @examples
#' set.seed(1)
#' ae <- data.frame(
#'   record_id = sample(1:30, 80, replace = TRUE),
#'   ae_term   = sample(c("Nausea", "Fatigue", "Anemia", "Rash"), 80, replace = TRUE)
#' )
#' ae$soc <- c(Nausea = "Gastrointestinal", Fatigue = "General",
#'             Anemia = "Blood", Rash = "Skin")[ae$ae_term]
#' ae$dose_cohort <- factor(
#'   c("DL1", "DL2", "DL3")[ae$record_id %% 3 + 1],
#'   levels = c("DL1", "DL2", "DL3")
#' )
#' build_ae_category_table(
#'   ae,
#'   vars = c("System organ class" = "soc", "Preferred term" = "ae_term"),
#'   include_variable_total = TRUE
#' )
#'
#' @export
build_ae_category_table <- function(ae,
                                    vars,
                                    cohort_col = "dose_cohort",
                                    id_col = "record_id",
                                    include_variable_total = FALSE) {

  if (!is.factor(ae[[cohort_col]])) {
    stop("`", cohort_col, "` must be a factor.", call. = FALSE)
  }

  cohort_levels <- levels(ae[[cohort_col]])
  strata <- c("All", cohort_levels)

  get_stratum <- function(s) {
    if (s == "All") ae else ae[ae[[cohort_col]] == s & !is.na(ae[[cohort_col]]), ]
  }

  # Summary of the levels (e.g. AE terms) within a variable
  summarise_level <- function(df, var_col) {
    df |>
      dplyr::filter(!is.na(.data[[var_col]])) |>
      dplyr::group_by(level = .data[[var_col]]) |>
      dplyr::summarise(
        n_AE = dplyr::n(),
        n_pt = dplyr::n_distinct(.data[[id_col]]),
        .groups = "drop"
      )
  }

  # Summary of the whole variable (e.g. whole SOC)
  summarise_variable_total <- function(df, var_col) {
    df |>
      dplyr::filter(!is.na(.data[[var_col]])) |>
      dplyr::summarise(
        n_AE = dplyr::n(),
        n_pt = dplyr::n_distinct(.data[[id_col]])
      )
  }

  # One block per variable
  build_var_block <- function(var_name, var_col) {

    all_levels <- if (is.factor(ae[[var_col]])) {
      levels(ae[[var_col]])
    } else {
      unique(stats::na.omit(ae[[var_col]]))
    }

    # Remove levels that do not occur
    all_levels <- all_levels[all_levels %in% ae[[var_col]]]

    out <- tibble::tibble(level = as.character(all_levels))

    for (s in strata) {
      piece <- summarise_level(get_stratum(s), var_col) |>
        dplyr::mutate(level = as.character(.data$level)) |>
        dplyr::rename(
          !!paste0(s, "_AE") := "n_AE",
          !!paste0(s, "_n")  := "n_pt"
        )
      out <- dplyr::left_join(out, piece, by = "level")
    }

    out <- out |>
      dplyr::mutate(
        dplyr::across(-dplyr::all_of("level"), ~ tidyr::replace_na(.x, 0)),
        variable = var_name,
        is_total = FALSE,
        .before  = "level"
      )

    if (include_variable_total) {
      total_row <- tibble::tibble(
        variable = var_name,
        level    = "Total",
        is_total = TRUE
      )
      for (s in strata) {
        total <- summarise_variable_total(get_stratum(s), var_col)
        total_row[[paste0(s, "_AE")]] <- total$n_AE
        total_row[[paste0(s, "_n")]]  <- total$n_pt
      }
      out <- dplyr::bind_rows(total_row, out)
    }

    out
  }

  purrr::imap(vars, ~ build_var_block(.y, .x)) |>
    dplyr::bind_rows()
}
