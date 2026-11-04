#' Overall AE and patient totals per stratum
#'
#' @param ae Data frame with one row per adverse event.
#' @param strata Character vector of strata (`"All"` plus cohort levels).
#' @param cohort_col Cohort column name.
#' @param id_col Patient identifier column name.
#'
#' @return A tibble with `stratum`, `n_AE`, `n_pt`.
#' @keywords internal
#' @noRd
compute_stratum_totals <- function(ae, strata, cohort_col, id_col) {

  purrr::map_dfr(strata, function(s) {

    df_s <- if (s == "All") {
      ae
    } else {
      ae[ae[[cohort_col]] == s & !is.na(ae[[cohort_col]]), ]
    }

    tibble::tibble(
      stratum = s,
      n_AE    = nrow(df_s),
      n_pt    = dplyr::n_distinct(df_s[[id_col]])
    )
  })
}


#' Render an AE category table as a striped kableExtra (LaTeX) table
#'
#' Takes the output of [build_ae_category_table()] and renders a table with an
#' overall row, one column pair (AE, n) for "Overall" and for each cohort,
#' striped variable blocks and variable labels via [kableExtra::pack_rows()].
#'
#' @param tbl Output of [build_ae_category_table()].
#' @param ae The AE data frame used to build `tbl` (for the totals row).
#' @param cohort_levels Character vector of cohort levels, usually
#'   `levels(ae$dose_cohort)`.
#' @param cohort_col Character. Cohort column. Default `"dose_cohort"`.
#' @param id_col Character. Patient identifier column. Default `"record_id"`.
#' @param stripe_colors Character vector of background colours, recycled over
#'   variable blocks. Default `c("gray!8", "white")`.
#' @param first.col.name Character. Header of the first column. Default
#'   `"Event"`.
#' @param caption Character or `NULL`. Table caption.
#'
#' @return A `kableExtra` table (`knitr_kable` object).
#'
#' @details Intended for PDF output via LaTeX (`booktabs`, `scale_down`,
#'   `repeat_header`).
#'
#' @seealso [build_ae_category_table()]
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
#' tbl <- build_ae_category_table(
#'   ae,
#'   vars = c("System organ class" = "soc", "Preferred term" = "ae_term"),
#'   include_variable_total = TRUE
#' )
#' old <- options(knitr.table.format = "latex")
#' render_ae_table(
#'   tbl, ae,
#'   cohort_levels = levels(ae$dose_cohort),
#'   caption = "Adverse events by dose cohort"
#' )
#' options(old)
#'
#' @export
render_ae_table <- function(tbl,
                            ae,
                            cohort_levels,
                            cohort_col = "dose_cohort",
                            id_col = "record_id",
                            stripe_colors = c("gray!8", "white"),
                            first.col.name = "Event",
                            caption = NULL) {

  strata <- c("All", cohort_levels)

  # Column order
  col_order <- c(
    "level",
    unlist(lapply(strata, function(s) paste0(s, c("_AE", "_n"))))
  )

  display_tbl <- tbl |>
    dplyr::select(dplyr::all_of(col_order))

  # Overall totals row (must be first)
  totals <- compute_stratum_totals(ae, strata, cohort_col, id_col)

  totals_row <- tibble::tibble(level = "Total, N (n)")
  for (s in strata) {
    totals_row[[paste0(s, "_AE")]] <- totals$n_AE[totals$stratum == s]
    totals_row[[paste0(s, "_n")]]  <- totals$n_pt[totals$stratum == s]
  }
  totals_row <- totals_row |>
    dplyr::select(dplyr::all_of(col_order))

  display_tbl <- dplyr::bind_rows(totals_row, display_tbl)

  # Variable blocks (+1 because of the overall row)
  blocks <- tbl |>
    dplyr::mutate(row = dplyr::row_number() + 1) |>
    dplyr::group_by(.data$variable) |>
    dplyr::summarise(
      start = min(.data$row),
      end   = max(.data$row),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$start)

  # Total rows
  total_rows <- which(tbl$is_total) + 1

  # LaTeX table
  k <- kableExtra::kbl(
    display_tbl,
    col.names = c(first.col.name, rep(c("AE", "n"), length(strata))),
    booktabs  = TRUE,
    align     = c("l", rep("c", length(col_order) - 1)),
    escape    = FALSE,
    caption   = caption
  ) |>
    kableExtra::kable_styling(latex_options = c("scale_down", "repeat_header"))

  # Header with Overall + cohorts
  k <- k |>
    kableExtra::add_header_above(
      c(
        " " = 1,
        "Overall" = 2,
        stats::setNames(rep(2, length(cohort_levels)), cohort_levels)
      )
    )

  # Bold column header (row 0) and overall row
  k <- k |> kableExtra::row_spec(0, bold = TRUE)
  k <- k |> kableExtra::row_spec(1, bold = TRUE)

  # Stripe variable blocks BEFORE pack_rows(), which inserts rows
  for (i in seq_len(nrow(blocks))) {
    bg <- stripe_colors[(i - 1) %% length(stripe_colors) + 1]
    k <- k |>
      kableExtra::row_spec(blocks$start[i]:blocks$end[i], background = bg)
  }

  # Bold total rows
  for (r in total_rows) {
    k <- k |> kableExtra::row_spec(r, bold = TRUE)
  }

  # Variable labels LAST (pack_rows() inserts rows)
  for (i in seq_len(nrow(blocks))) {
    k <- k |>
      kableExtra::pack_rows(
        group_label = blocks$variable[i],
        start_row   = blocks$start[i],
        end_row     = blocks$end[i],
        bold        = TRUE,
        italic      = FALSE,
        indent      = TRUE
      )
  }

  k
}
