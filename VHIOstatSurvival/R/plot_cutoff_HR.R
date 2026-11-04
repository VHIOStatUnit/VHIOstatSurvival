#' Explore Optimal Cutoffs for a Continuous Variable Using Cox Hazard Ratios
#'
#' @description
#' Scans a grid of candidate cutoffs for a continuous variable and, for each
#' cutoff, dichotomises patients into "Low" (< cutoff) and "High" (>= cutoff)
#' groups. For every cutoff a univariable Cox proportional hazards model and a
#' log-rank test are fitted. The function then selects a "best" cutoff
#' according to a user-defined criterion and visualises the distribution of the
#' variable (histogram) together with the hazard ratio and its 95% confidence
#' interval across all cutoffs.
#'
#' @details
#' Candidate cutoffs default to the 5th to 95th percentile of the variable in
#' steps of one percentile. Cutoffs resulting in a group with fewer than
#' `min_n` patients are skipped. Confidence intervals are Wald intervals
#' (`exp(beta +/- 1.96 * se)`).
#'
#' **Caution:** Selecting a cutoff by maximising the HR or minimising the
#' p-value is a data-driven ("minimum p-value") approach. It inflates the
#' type I error and overestimates the effect size. The reported p-values are
#' **not** adjusted for multiple testing and results should be interpreted
#' as exploratory, ideally validated in an independent cohort or with
#' bootstrap/cross-validation.
#'
#' @param data A data frame containing the variable, time and event columns.
#' @param variable Character string. Name of the continuous variable in
#'   `data` for which cutoffs are evaluated.
#' @param time Character string. Name of the column containing follow-up time.
#' @param event Character string. Name of the column containing the event
#'   indicator (1/TRUE = event, 0/FALSE = censored).
#' @param cutoffs Optional numeric vector of candidate cutoffs. If `NULL`
#'   (default), the unique values of the 5th to 95th percentile (step 1%) of
#'   `variable` are used.
#' @param min_n Integer. Minimum number of patients required in each of the
#'   "Low" and "High" groups for a cutoff to be evaluated. Default is `10`.
#' @param direction Character string. Currently only used as a label
#'   (default `"High vs Low"`); the comparison is always High (>= cutoff)
#'   versus Low (< cutoff).
#' @param best_cutoff Character string defining the criterion used to select
#'   the "best" cutoff. One of `"max_HR"` (largest hazard ratio),
#'   `"min_p"` (smallest Cox Wald p-value), `"max_logHR"` (largest log hazard
#'   ratio) or `"max_logrank"` (largest log-rank chi-square statistic).
#'   Default is `"max_HR"`.
#' @param xlim Numeric vector of length 2. x-axis limits for both plots
#'   (on the scale of `variable`). Default is `c(0, 2)`.
#' @param ylim Numeric vector of length 2. y-axis limits of the hazard ratio
#'   plot. Default is `c(0, 5)`.
#'
#' @return A list with four elements:
#' \describe{
#'   \item{plot}{A `patchwork` object: histogram of `variable` (top) and
#'     hazard ratio with 95% CI by cutoff (bottom). The selected cutoff is
#'     marked by a dashed vertical line.}
#'   \item{results}{A tibble with one row per evaluated cutoff and columns
#'     `cutoff`, `n_low`, `n_high`, `HR`, `CI_low`, `CI_high`, `p_value`,
#'     `logHR`, `logrank_chisq` and `logrank_p`.}
#'   \item{best_cutoff}{Numeric. The selected cutoff value.}
#'   \item{best_result}{A one-row tibble with the results for the selected
#'     cutoff.}
#' }
#'
#' @examples
#' \dontrun{
#' library(survival)
#'
#' # Example with the veteran lung cancer data
#' res <- plot_cutoff_HR(
#'   data        = survival::veteran,
#'   variable    = "karno",
#'   time        = "time",
#'   event       = "status",
#'   min_n       = 15,
#'   best_cutoff = "max_logrank",
#'   xlim        = c(20, 100),
#'   ylim        = c(0, 3)
#' )
#'
#' res$plot
#' res$best_cutoff
#' res$best_result
#' head(res$results)
#'
#' # User-defined cutoffs
#' res2 <- plot_cutoff_HR(
#'   data     = survival::veteran,
#'   variable = "karno",
#'   time     = "time",
#'   event    = "status",
#'   cutoffs  = c(40, 50, 60, 70, 80),
#'   xlim     = c(20, 100)
#' )
#' }
#'
#' @importFrom rlang .data
#' @export
plot_cutoff_HR <- function(data,
                           variable,
                           time,
                           event,
                           cutoffs = NULL,
                           min_n = 10,
                           direction = "High vs Low",
                           best_cutoff = c(
                             "max_HR",
                             "min_p",
                             "max_logHR",
                             "max_logrank"
                           ),
                           xlim = c(0, 2),
                           ylim = c(0, 5)) {

  best_cutoff <- base::match.arg(best_cutoff)

  # ---------------------------------------------------------------------------
  # Prepare data
  # ---------------------------------------------------------------------------

  df <- data |>
    dplyr::select(dplyr::all_of(c(variable, time, event))) |>
    dplyr::filter(
      !is.na(.data[[variable]]),
      !is.na(.data[[time]]),
      !is.na(.data[[event]])
    )

  x <- df[[variable]]

  # Candidate cutoffs
  if (base::is.null(cutoffs)) {
    cutoffs <- base::sort(base::unique(
      stats::quantile(
        x,
        probs = base::seq(0.05, 0.95, by = 0.01),
        na.rm = TRUE,
        names = FALSE
      )
    ))
  }

  # ---------------------------------------------------------------------------
  # Run Cox model + log-rank test for every cutoff
  # ---------------------------------------------------------------------------

  results <- base::lapply(cutoffs, function(cut) {

    df_tmp <- df |>
      dplyr::mutate(
        group = base::ifelse(
          .data[[variable]] >= cut,
          "High",
          "Low"
        ),
        group = base::factor(group, levels = c("Low", "High"))
      )

    # Check minimum number of patients in each group
    n_low  <- base::sum(df_tmp$group == "Low")
    n_high <- base::sum(df_tmp$group == "High")

    if (n_low < min_n || n_high < min_n) {
      return(NULL)
    }

    # -------------------------------------------------------------------------
    # Cox model
    # -------------------------------------------------------------------------

    fit <- base::tryCatch(
      survival::coxph(
        stats::as.formula(
          base::paste0(
            "survival::Surv(`", time, "`, `", event, "`) ~ group"
          )
        ),
        data = df_tmp
      ),
      error = function(e) NULL
    )

    if (base::is.null(fit)) {
      return(NULL)
    }

    s <- base::summary(fit)

    # Extract HR and CI
    beta <- stats::coef(fit)[1]
    se   <- base::sqrt(stats::vcov(fit)[1, 1])

    # -------------------------------------------------------------------------
    # Log-rank test
    # -------------------------------------------------------------------------

    logrank <- base::tryCatch(
      survival::survdiff(
        stats::as.formula(
          base::paste0(
            "survival::Surv(`", time, "`, `", event, "`) ~ group"
          )
        ),
        data = df_tmp
      ),
      error = function(e) NULL
    )

    if (base::is.null(logrank)) {
      return(NULL)
    }

    # Log-rank chi-square statistic
    logrank_chisq <- base::as.numeric(logrank$chisq)

    # Corresponding p-value
    logrank_p <- stats::pchisq(
      logrank_chisq,
      df = base::length(logrank$n) - 1,
      lower.tail = FALSE
    )

    tibble::tibble(
      cutoff = cut,
      n_low = n_low,
      n_high = n_high,
      HR = base::exp(beta),
      CI_low = base::exp(beta - 1.96 * se),
      CI_high = base::exp(beta + 1.96 * se),
      p_value = s$coefficients[1, "Pr(>|z|)"],
      logHR = beta,
      logrank_chisq = logrank_chisq,
      logrank_p = logrank_p
    )
  })

  results <- dplyr::bind_rows(results)

  if (base::nrow(results) == 0) {
    base::stop("No valid cutoffs remained. Try reducing min_n.")
  }

  # ---------------------------------------------------------------------------
  # Determine "best" cutoff
  # ---------------------------------------------------------------------------

  if (best_cutoff == "max_HR") {

    best <- results |>
      dplyr::slice_max(HR, n = 1, with_ties = FALSE)

  } else if (best_cutoff == "min_p") {

    best <- results |>
      dplyr::slice_min(p_value, n = 1, with_ties = FALSE)

  } else if (best_cutoff == "max_logHR") {

    best <- results |>
      dplyr::slice_max(logHR, n = 1, with_ties = FALSE)

  } else if (best_cutoff == "max_logrank") {

    best <- results |>
      dplyr::slice_max(logrank_chisq, n = 1, with_ties = FALSE)
  }

  best_value <- best$cutoff

  # ---------------------------------------------------------------------------
  # Histogram
  # ---------------------------------------------------------------------------

  p_hist <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[variable]])) +
    ggplot2::geom_histogram(
      bins = 30,
      fill = "grey80",
      colour = "grey30"
    ) +
    ggplot2::geom_vline(
      xintercept = best_value,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Number of patients",
      title = variable
    ) +
    ggplot2::coord_cartesian(
      xlim = xlim
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      )
    )

  # ---------------------------------------------------------------------------
  # HR plot
  # ---------------------------------------------------------------------------

  p_hr <- ggplot2::ggplot(results, ggplot2::aes(x = cutoff, y = HR)) +

    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey50"
    ) +

    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = CI_low,
        ymax = CI_high
      ),
      width = 0,
      linewidth = 0.5
    ) +

    ggplot2::geom_point(
      size = 2
    ) +

    ggplot2::geom_vline(
      xintercept = best_value,
      linetype = "dashed",
      linewidth = 0.7
    ) +

    ggplot2::labs(
      x = base::paste(variable, "cutoff"),
      y = "Hazard ratio (High vs Low)"
    ) +

    ggplot2::coord_cartesian(
      xlim = xlim,
      ylim = ylim
    ) +

    ggplot2::theme_classic()

  # ---------------------------------------------------------------------------
  # Combine
  # ---------------------------------------------------------------------------

  final_plot <- patchwork::wrap_plots(
    p_hist,
    p_hr,
    ncol = 1,
    heights = c(1, 1.3)
  )

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------

  base::list(
    plot = final_plot,
    results = results,
    best_cutoff = best_value,
    best_result = best
  )
}
