#' Convert a gtsummary table into a striped kableExtra (LaTeX) table
#'
#' Takes a \pkg{gtsummary} table, converts it with
#' [gtsummary::as_kable_extra()] and shades alternating variable blocks
#' (all rows belonging to one variable share the same background).
#'
#' @param gts_obj A `gtsummary` object, e.g. from [gtsummary::tbl_summary()].
#' @param stripe_color Character. LaTeX/xcolor colour used for shaded blocks.
#'   Default `"gray!10"`.
#' @param caption.table Character or `NULL`. Table caption.
#'
#' @return A `kableExtra` table (`knitr_kable` object) with booktabs style,
#'   font size 8, `scale_down` and `repeat_header` options.
#'
#' @details Intended for PDF output via LaTeX. The `booktabs`, `scale_down`
#'   and `repeat_header` options are LaTeX-specific.
#'
#' @examples
#' old <- options(knitr.table.format = "latex")
#' tbl <- gtsummary::tbl_summary(
#'   gtsummary::trial,
#'   by = trt,
#'   include = c(age, grade, stage)
#' )
#' summary_table_formatter(tbl, caption.table = "Baseline characteristics")
#' options(old)
#'
#' @export
summary_table_formatter <- function(gts_obj,
                                    stripe_color = "gray!10",
                                    caption.table = NULL) {

  if (!inherits(gts_obj, "gtsummary")) {
    stop("`gts_obj` must be a gtsummary object.", call. = FALSE)
  }

  # 1. Row range of each variable
  row_info <- gts_obj$table_body |>
    dplyr::mutate(row_number = dplyr::row_number()) |>
    dplyr::group_by(.data$variable) |>
    dplyr::summarise(
      start = min(.data$row_number),
      end   = max(.data$row_number),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$start) |>
    dplyr::mutate(stripe = dplyr::row_number() %% 2 == 1)

  # 2. Build LaTeX table
  tbl <- gts_obj |>
    gtsummary::as_kable_extra(
      booktabs = TRUE,
      caption  = caption.table
    ) |>
    kableExtra::kable_styling(
      latex_options = c("scale_down", "repeat_header"),
      font_size = 8
    )

  # 3. Shade every second variable block
  for (i in seq_len(nrow(row_info))) {
    if (row_info$stripe[i]) {
      tbl <- tbl |>
        kableExtra::row_spec(
          row_info$start[i]:row_info$end[i],
          background = stripe_color
        )
    }
  }

  tbl
}
