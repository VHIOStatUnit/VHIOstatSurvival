#' Summarise adverse events by term and grade band
#'
#' Counts the number of distinct patients per adverse event term for all
#' grades, grade 1-2 and grade 3-5, and adds an "Any" row (patients with at
#' least one AE, counted at their highest grade).
#'
#' @param ae Data frame with one row per adverse event.
#' @param N_total Integer. Total number of patients (denominator for %).
#' @param id_col Character. Patient identifier column. Default `"record_id"`.
#' @param term_col Character. AE term column. Default `"ae_term"`.
#' @param grade_col Character. Numeric grade column (1-5). Default `"grade"`.
#'
#' @return A tibble with the term column, formatted columns `All_grade`,
#'   `Grade_1_2`, `Grade_3_4_5` (`"n (x.x%)"`) and numeric percentages
#'   `Perc_all`, `Perc_1_2`, `Perc_3_4_5` (e.g. for plotting). Terms without
#'   any event are dropped; rows are sorted by decreasing frequency.
#'
#' @examples
#' set.seed(1)
#' ae <- data.frame(
#'   record_id = sample(1:30, 80, replace = TRUE),
#'   ae_term   = sample(c("Nausea", "Fatigue", "Anemia", "Rash"), 80, replace = TRUE),
#'   grade     = sample(1:4, 80, replace = TRUE, prob = c(.4, .3, .2, .1))
#' )
#' summarize_ae(ae, N_total = 30)
#'
#' @export
summarize_ae <- function(ae,
                         N_total,
                         id_col    = "record_id",
                         term_col  = "ae_term",
                         grade_col = "grade") {

  # Helper: count distinct patients by grade band
  count_grades <- function(df) {
    df |>
      dplyr::summarise(
        All_grade_n   = dplyr::n_distinct(.data[[id_col]][.data[[grade_col]] %in% 1:5]),
        Grade_1_2_n   = dplyr::n_distinct(.data[[id_col]][.data[[grade_col]] %in% 1:2]),
        Grade_3_4_5_n = dplyr::n_distinct(.data[[id_col]][.data[[grade_col]] %in% 3:5])
      )
  }

  # Per-term rows
  by_term <- ae |>
    dplyr::group_by(dplyr::across(dplyr::all_of(term_col))) |>
    count_grades()

  # "Any" row: each patient at their highest grade
  ae_overall <- ae |>
    dplyr::group_by(dplyr::across(dplyr::all_of(id_col))) |>
    dplyr::slice_max(
      order_by  = .data[[grade_col]],
      n         = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup()

  any_row <- ae_overall |>
    count_grades() |>
    dplyr::mutate(!!term_col := "Any")

  dplyr::bind_rows(by_term, any_row) |>
    dplyr::mutate(
      All_grade   = sprintf("%d (%.1f%%)", .data$All_grade_n,   100 * .data$All_grade_n   / N_total),
      Grade_1_2   = sprintf("%d (%.1f%%)", .data$Grade_1_2_n,   100 * .data$Grade_1_2_n   / N_total),
      Grade_3_4_5 = sprintf("%d (%.1f%%)", .data$Grade_3_4_5_n, 100 * .data$Grade_3_4_5_n / N_total),
      Perc_all    = 100 * .data$All_grade_n   / N_total,
      Perc_1_2    = 100 * .data$Grade_1_2_n   / N_total,
      Perc_3_4_5  = 100 * .data$Grade_3_4_5_n / N_total
    ) |>
    dplyr::filter(.data$All_grade_n > 0) |>
    dplyr::arrange(dplyr::desc(.data$All_grade_n)) |>
    dplyr::select(
      dplyr::all_of(term_col),
      "All_grade", "Grade_1_2", "Grade_3_4_5",
      "Perc_all", "Perc_1_2", "Perc_3_4_5"
    )
}
