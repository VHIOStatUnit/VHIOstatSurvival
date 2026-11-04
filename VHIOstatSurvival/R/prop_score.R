#' Propensity score weighting with balance diagnostics
#'
#' Estimates propensity score weights for a binary or multi-category treatment
#' using any method and estimand supported by \pkg{WeightIt}, and optionally
#' produces a love plot and a baseline table comparing the weighted and
#' unweighted populations. The table can show either p-values or standardized
#' mean differences (SMD).
#'
#' @details
#' Weights are estimated with [WeightIt::weightit()]. By default
#' (`method = "ps"`, `estimand = "ATE"`, `stabilize = FALSE`) this gives
#' classical inverse probability of treatment weights (IPTW). Other methods
#' (e.g. `"ebal"`, `"cbps"`, `"energy"`) and estimands (e.g. `"ATT"`, `"ATC"`,
#' `"ATO"`) can be selected via `method` and `estimand`; which combinations are
#' valid depends on \pkg{WeightIt} and on the number of treatment levels.
#' Additional arguments (e.g. `link`, `missing`) are passed on through `...`.
#'
#' The number of treatment levels is detected automatically:
#'
#' * **Two levels:** if the chosen method estimates a propensity score (e.g.
#'   `"ps"`/`"glm"`), it is returned in column `PS`. With `smd = FALSE`, group
#'   comparisons in the baseline tables use the Wilcoxon rank-sum test
#'   (continuous) and the chi-squared test (categorical).
#' * **More than two levels:** no `PS` column is returned. Group comparisons
#'   use the Kruskal-Wallis test (continuous) and the chi-squared test
#'   (categorical). `smd = TRUE` is not available.
#'
#' With `smd = TRUE`, the tables contain an overall column and standardized
#' mean differences (via [gtsummary::add_difference()]) instead of p-values;
#' this requires the \pkg{smd} package. Continuous variables are summarised
#' on multiple lines (`"continuous2"`) and dichotomous variables are shown
#' with both levels.
#'
#' The weighted table is built from a [survey::svydesign()] object using the
#' estimated weights, with the corresponding `svy.*` tests from \pkg{gtsummary}.
#' With `smd = FALSE`, rows for factor levels that are empty in all groups are
#' dropped from the tables.
#'
#' @param data A data frame containing the treatment and baseline variables.
#' @param baselinevars Character vector of covariate names used in the
#'   weighting model and shown in the baseline tables.
#' @param treatment Character string, the name of the treatment variable in
#'   `data` (binary or with more than two categories).
#' @param table_var_header Character string used as header of the variable
#'   column in the baseline table. Default `"Characteristic"`.
#' @param table_cols_header Character string used as spanning header over the
#'   treatment columns. Defaults to the name of the treatment variable.
#' @param table_caption Optional character string used as table caption.
#' @param method Character string, the weighting method passed to
#'   [WeightIt::weightit()]. Default `"ps"` (propensity score weighting).
#' @param estimand Character string, the estimand passed to
#'   [WeightIt::weightit()], e.g. `"ATE"` (default), `"ATT"`, `"ATC"`, `"ATO"`.
#' @param focal Optional character string, the focal treatment level. Required
#'   by \pkg{WeightIt} for `"ATT"`/`"ATC"` with more than two treatment levels.
#'   Default `NULL`.
#' @param stabilize Logical. If `TRUE`, stabilized weights are computed
#'   (passed to [WeightIt::weightit()]). Default `FALSE`.
#' @param smd Logical. If `TRUE`, the baseline tables show standardized mean
#'   differences (and an overall column) instead of p-values. Only available
#'   for binary treatments; requires the \pkg{smd} package. Default `FALSE`.
#' @param love.plot.show Logical. If `TRUE`, print a love plot of absolute
#'   standardized mean differences before and after weighting
#'   (threshold 0.1). Default `FALSE`.
#' @param table.show Logical. If `TRUE`, print the merged baseline table
#'   (weighted vs. unweighted). Default `FALSE`.
#' @param table.show.html Logical. If `TRUE`, write the table as raw HTML via
#'   [cat()] (useful for `results = "asis"` chunks in R Markdown).
#'   Default `FALSE`.
#' @param ... Further arguments passed to [WeightIt::weightit()].
#'
#' @return `data` with an additional column `weights_<estimand>` (e.g.
#'   `weights_ATE`, `weights_ATT`) holding the estimated weights and, for
#'   binary treatments with a method that estimates a propensity score, a
#'   column `PS` taken from the fitted weighting model.
#'   Called for its return value; plots and tables are side effects.
#'
#' @seealso [WeightIt::weightit()], [cobalt::love.plot()],
#'   [gtsummary::tbl_summary()], [gtsummary::tbl_svysummary()],
#'   [gtsummary::add_difference()]
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 300
#' d <- data.frame(
#'   age = rnorm(n, 60, 10),
#'   sex = factor(sample(c("F", "M"), n, replace = TRUE))
#' )
#' d$arm  <- factor(rbinom(n, 1, plogis(-0.5 + 0.03 * (d$age - 60))),
#'                  labels = c("Control", "Treated"))
#' d$arm3 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
#'
#' # Binary treatment, default IPTW (ATE), p-values in the table
#' d2 <- prop_score(d, baselinevars = c("age", "sex"), treatment = "arm",
#'                  table_var_header = "Characteristic",
#'                  table_cols_header = "Treatment",
#'                  table_caption = "Baseline characteristics",
#'                  love.plot.show = TRUE, table.show = TRUE)
#' head(d2[, c("PS", "weights_ATE")])
#'
#' # Stabilized weights, standardized mean differences in the table
#' d_smd <- prop_score(d, c("age", "sex"), "arm",
#'                     stabilize = TRUE, smd = TRUE, table.show = TRUE)
#'
#' # Different estimand (ATT) and method (entropy balancing; no PS column)
#' d_eb <- prop_score(d, c("age", "sex"), "arm", method = "ebal",
#'                    estimand = "ATT", table.show = TRUE)
#'
#' # Multi-category treatment
#' d3 <- prop_score(d, baselinevars = c("age", "sex"), treatment = "arm3",
#'                  table.show = TRUE)
#' }
#'
#' @export
prop_score <- function(data,
                       baselinevars,
                       treatment,
                       table_var_header  = "Characteristic",
                       table_cols_header = treatment,
                       table_caption     = NULL,
                       method            = "ps",
                       estimand          = "ATE",
                       focal             = NULL,
                       stabilize         = FALSE,
                       smd               = FALSE,
                       love.plot.show    = FALSE,
                       table.show        = FALSE,
                       table.show.html   = FALSE,
                       ...) {

  # ---- input checks ---------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  missing_vars <- setdiff(c(treatment, baselinevars), names(data))
  if (length(missing_vars) > 0) {
    stop("Variable(s) not found in `data`: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  n_lev <- length(unique(stats::na.omit(data[[treatment]])))
  if (n_lev < 2) {
    stop("`treatment` must have at least two levels.", call. = FALSE)
  }
  if (!is.character(method) || length(method) != 1 ||
      !is.character(estimand) || length(estimand) != 1) {
    stop("`method` and `estimand` must each be a single character string.",
         call. = FALSE)
  }
  if (isTRUE(smd)) {
    if (n_lev != 2) {
      stop("`smd = TRUE` is only available for binary treatments.",
           call. = FALSE)
    }
    if (!requireNamespace("smd", quietly = TRUE)) {
      stop("Package 'smd' is required for `smd = TRUE`.", call. = FALSE)
    }
  }

  wcol         <- paste0("weights_", estimand)
  weight_label <- if (method %in% c("ps", "glm")) "IPTW" else paste0(method, " (", estimand, ")")

  # ---- weights --------------------------------------------------------------
  ps.formula <- stats::reformulate(baselinevars, response = treatment)

  W.out <- WeightIt::weightit(ps.formula, data = data, method = method,
                              estimand = estimand, focal = focal,
                              stabilize = stabilize, ...)
  data[[wcol]] <- W.out$weights

  # propensity score only for binary treatments and methods that estimate one
  if (n_lev == 2 && !is.null(W.out$ps)) {
    data$PS <- W.out$ps
  }

  # ---- love plot ------------------------------------------------------------
  if (isTRUE(love.plot.show)) {
    print(cobalt::love.plot(W.out,
                            thresholds = c(m = .1),
                            abs        = TRUE,
                            var.order  = "unadjusted",
                            line       = TRUE))
  }

  # ---- baseline tables ------------------------------------------------------
  if (isTRUE(table.show) || isTRUE(table.show.html)) {

    if (isTRUE(smd)) {
      tests <- tests_svy <- NULL
    } else if (n_lev == 2) {
      tests     <- list(gtsummary::all_continuous()  ~ "wilcox.test",
                        gtsummary::all_categorical() ~ "chisq.test")
      tests_svy <- list(gtsummary::all_continuous()  ~ "svy.wilcox.test",
                        gtsummary::all_categorical() ~ "svy.chisq.test")
    } else {
      tests     <- list(gtsummary::all_continuous()  ~ "kruskal.test",
                        gtsummary::all_categorical() ~ "chisq.test")
      tests_svy <- list(gtsummary::all_continuous()  ~ "svy.kruskal.test",
                        gtsummary::all_categorical() ~ "svy.chisq.test")
    }

    tbl_type <- list(gtsummary::all_continuous()  ~ "continuous2",
                     gtsummary::all_dichotomous() ~ "categorical")

    tab <- data |>
      gtsummary::tbl_summary(include = dplyr::all_of(baselinevars),
                             by      = dplyr::all_of(treatment),
                             type    = tbl_type,
                             missing = "ifany",
                             missing_text = "Missing") |>
      .ps_style_table(table_var_header, table_cols_header, tests, smd = smd)

    design_PS <- survey::svydesign(ids = ~1, data = data,
                                   weights = stats::reformulate(wcol))

    tab_PS <- design_PS |>
      gtsummary::tbl_svysummary(include = dplyr::all_of(baselinevars),
                                by      = dplyr::all_of(treatment),
                                type    = tbl_type,
                                missing = "ifany",
                                missing_text = "Missing") |>
      .ps_style_table(table_var_header, table_cols_header, tests_svy, smd = smd)

    tab_final <- gtsummary::tbl_merge(
      tbls        = list(tab_PS, tab),
      tab_spanner = c(paste0("**Balanced by ", weight_label, "**"),
                      "**Unbalanced**")
    ) |>
      gtsummary::modify_spanning_header(
        c(dplyr::starts_with("stat_"),
          dplyr::starts_with("p.value"),
          dplyr::starts_with("estimate")) ~
          paste0("**", table_cols_header, "**"),
        level = 2
      )

    if (!is.null(table_caption)) {
      tab_final <- gtsummary::modify_caption(tab_final,
                                             paste0("**", table_caption, "**"))
    }

    if (isTRUE(table.show.html)) {
      cat(gt::as_raw_html(gtsummary::as_gt(tab_final)))
    }
    if (isTRUE(table.show)) {
      print(tab_final)
    }
  }

  data
}


#' Apply the common formatting to a baseline table
#'
#' Internal helper for `prop_score()`: sets headers and bolds labels. With
#' `smd = FALSE` it additionally drops empty factor-level rows and adds N and
#' p-values; with `smd = TRUE` it adds an overall column, N and standardized
#' mean differences, and hides the confidence interval columns.
#'
#' @param tbl A `tbl_summary` or `tbl_svysummary` object.
#' @param var_header,cols_header Character strings for the label column header
#'   and the spanning header of the statistic columns.
#' @param tests List of gtsummary test specifications passed to
#'   [gtsummary::add_p()] (ignored if `smd = TRUE`).
#' @param smd Logical, see [prop_score()].
#'
#' @return A gtsummary table.
#' @keywords internal
#' @noRd
.ps_style_table <- function(tbl, var_header, cols_header, tests = NULL,
                            smd = FALSE) {

  if (isTRUE(smd)) {
    tbl <- tbl |>
      gtsummary::add_overall() |>
      gtsummary::add_n() |>
      gtsummary::add_difference(test = dplyr::everything() ~ "smd")
  }

  tbl <- tbl |>
    gtsummary::modify_header(
      label = paste0("**", var_header, "**"),
      gtsummary::all_stat_cols() ~
        "**{level}**, N = {gtsummary::style_number(n, digits = 1)} ({gtsummary::style_percent(p)}%)"
    ) |>
    gtsummary::modify_spanning_header(
      gtsummary::all_stat_cols() ~ paste0("**", cols_header, "**")
    )

  if (!isTRUE(smd)) {
    tbl <- tbl |>
      gtsummary::modify_table_body(
        \(x) dplyr::filter(
          x,
          # keep all non-level rows so variable names stay;
          # keep level rows only if any stat column is not "0 (0%)"
          .data$row_type != "level" |
            dplyr::if_any(dplyr::starts_with("stat_"),
                          \(v) v != "0 (0%)" & !is.na(v))
        )
      ) |>
      gtsummary::add_n() |>
      gtsummary::add_p(test = tests)
  }

  tbl <- gtsummary::bold_labels(tbl)

  if (isTRUE(smd)) {
    tbl <- gtsummary::modify_column_hide(tbl, columns = dplyr::contains("conf"))
  }

  tbl
}
