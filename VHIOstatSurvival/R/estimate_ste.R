# R/estimate_ste.R
# Surrogate threshold effect (STE): the minimum treatment effect on a surrogate
# endpoint that reliably predicts a treatment effect on the clinical outcome.


#' Plot a surrogate threshold effect (STE) analysis
#'
#' Draws the trial-level scatter plot of outcome HR versus surrogate HR on log
#' axes, together with the weighted regression line, the 95% prediction
#' interval and the STE (vertical dashed line). Internal helper, called by
#' [estimate_ste()].
#'
#' @param ma A list returned by [estimate_ste()].
#' @param xlabel,ylabel Axis labels (x = surrogate, y = outcome).
#' @param group_var Optional name of a column in `ma$data` used to colour the
#'   points (e.g. tumour type). Ignored if `NULL` or not a column of the data.
#' @param palette Character vector of colours, one per level of `group_var`.
#'   Must be at least as long as the number of levels. Only used when
#'   `group_var` is supplied.
#' @param hr_breaks Numeric vector of hazard ratios used as axis breaks and
#'   labels on both (log10) axes.
#' @param weights Name of the column in `ma$data` mapped to point size
#'   (usually the trial sample size). Defaults to the weight column used in
#'   the fit.
#' @param xlim.min,xlim.max Limits of the x axis (hazard-ratio scale).
#' @param ylim.min,ylim.max Limits of the y axis (hazard-ratio scale).
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @seealso [estimate_ste()]
#' @importFrom rlang .data
#' @noRd
plot_ste <- function(ma,
                     xlabel    = "Hazard ratio (surrogate)",
                     ylabel    = "Hazard ratio (outcome)",
                     group_var = NULL,
                     palette   = c("#F3AB97", "#B58494", "grey80", "#D6DD7F"),
                     hr_breaks = c(0.2, 0.3, 0.5, 0.7, 1.0, 1.5),
                     weights   = ma$weights,
                     xlim.min  = 0.15,
                     xlim.max  = 2.5,
                     ylim.min  = 0.15,
                     ylim.max  = 2.5) {

  plot_data <- ma$data
  plot_data$hr_surrogate <- exp(plot_data[[ma$surrogate]])
  plot_data$hr_outcome   <- exp(plot_data[[ma$outcome]])

  # colour by group_var if supplied and present in the data
  has_group <- !is.null(group_var) && group_var %in% names(plot_data)

  base_aes <- if (has_group) {
    ggplot2::aes(x = .data$hr_surrogate, y = .data$hr_outcome, color = .data[[group_var]])
  } else {
    ggplot2::aes(x = .data$hr_surrogate, y = .data$hr_outcome)
  }

  p <- ggplot2::ggplot(plot_data, base_aes) +

    # prediction interval
    ggplot2::geom_ribbon(
      data = ma$pred,
      ggplot2::aes(x = .data$hr_surrogate, ymin = .data$hr_outcome_lwr, ymax = .data$hr_outcome_upr),
      inherit.aes = FALSE, alpha = 0.15
    ) +

    # regression line
    ggplot2::geom_line(
      data = ma$pred,
      ggplot2::aes(x = .data$hr_surrogate, y = .data$hr_outcome_fit),
      inherit.aes = FALSE, color = "#2C3E50", linewidth = 0.9
    )

  # STE: where the upper prediction limit crosses outcome HR = 1
  if (!is.na(ma$fit_cross)) {
    p <- p + ggplot2::geom_vline(xintercept = ma$fit_cross, color = "tomato",
                                 linewidth = 0.6, linetype = "dashed")
  }

  # points, size mapped to the weights (no legend)
  p <- p + if (has_group) {
    ggplot2::geom_point(ggplot2::aes(size = .data[[weights]]), shape = 16, alpha = 0.85)
  } else {
    ggplot2::geom_point(ggplot2::aes(size = .data[[weights]]), shape = 1, stroke = 0.9, color = "#2C3E50")
  }

  p <- p +
    ggplot2::scale_size_continuous(range = c(2, 8), guide = "none") +
    ggplot2::scale_x_log10(limits = c(xlim.min, xlim.max), breaks = hr_breaks,
                           labels = as.character(hr_breaks)) +
    ggplot2::scale_y_log10(limits = c(ylim.min, ylim.max), breaks = hr_breaks,
                           labels = as.character(hr_breaks))

  if (has_group) p <- p + ggplot2::scale_color_manual(values = palette)

  p +
    ggplot2::labs(
      title    = " ",
      subtitle = sprintf("R\u00b2 (adj.) = %.3f  |  STE = %.2f", ma$r2, ma$ste),
      x        = xlabel,
      y        = ylabel,
      color    = if (has_group) group_var else NULL,
      caption  = "Prediction interval 95%"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = 11),
      plot.caption     = ggplot2::element_text(color = "grey50", size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title       = ggplot2::element_text(size = 18, face = "bold"),
      axis.text        = ggplot2::element_text(size = 15),
      legend.position  = "bottom",
      legend.text      = ggplot2::element_text(size = 14),
      legend.title     = ggplot2::element_text(size = 16, face = "bold"),
      legend.key.size  = ggplot2::unit(1.2, "cm")
    )
}


#' Estimate the surrogate threshold effect (STE)
#'
#' Estimates the surrogate threshold effect from trial-level data using a
#' weighted linear regression of the treatment effect on the clinical outcome
#' (e.g. log HR for OS) against the treatment effect on the surrogate (e.g. log
#' HR for PFS). The STE is the smallest surrogate treatment effect for which
#' the upper limit of the 95% prediction interval for the outcome effect
#' still lies below the null (outcome HR = 1), i.e. the surrogate effect that
#' is reliably predictive of a non-null effect on the outcome.
#'
#' @details
#' The model is `lm(outcome ~ surrogate, weights = <weights>)`, with both
#' effects on the log-HR scale. The 95% prediction interval is evaluated on a
#' grid of 200 surrogate values between `log(min.pred)` and `log(max.pred)`,
#' using the median trial weight as the weight of a new trial.
#'
#' Two closely related STE values are returned: `ste` is the midpoint (on the
#' log scale) of the two grid points enclosing the crossing of the upper
#' prediction limit with 0, `fit_cross` is the linearly interpolated crossing
#' point. The dashed line in the plot is `fit_cross`; the subtitle reports
#' `ste`. If the upper limit does not cross 0 inside the grid, both are `NA`
#' and a warning is issued; widen `min.pred`/`max.pred` in that case.
#'
#' The reported \eqn{R^2} is the **adjusted** \eqn{R^2} of the weighted
#' regression.
#'
#' @param data A data frame with one row per trial/comparison.
#' @param surrogate Name of the column with the log hazard ratio of the
#'   surrogate endpoint.
#' @param outcome Name of the column with the log hazard ratio of the clinical
#'   outcome.
#' @param weights Name of the column with the regression weights, usually the
#'   trial sample size.
#' @param min.pred,max.pred Range of the surrogate hazard ratio (not log
#'   scale) over which the prediction interval is computed.
#' @param plot Logical; if `TRUE` (default) the plot is printed.
#' @param xlabel,ylabel Axis labels (x = surrogate, y = outcome).
#' @param hr_breaks Hazard ratios used as axis breaks and labels.
#' @param group_var Optional name of a column in `data` used to colour the
#'   points (e.g. tumour type). Ignored if `NULL` or not a column of `data`.
#' @param xlim.min,xlim.max,ylim.min,ylim.max Axis limits of the plot on the
#'   hazard-ratio scale.
#' @param palette Colours for the levels of `group_var`; needs at least as many
#'   colours as there are levels. Only used when `group_var` is supplied.
#'
#' @return A list with elements
#' \describe{
#'   \item{surrogate, outcome, weights}{Column names used.}
#'   \item{model}{The fitted [stats::lm()] object.}
#'   \item{coef}{Data frame of coefficients with 95% CI (`fit`, `lwr`, `upr`).}
#'   \item{pred}{Data frame with the prediction grid: `log_surrogate`,
#'     `hr_surrogate`, `fit`/`lwr`/`upr` (log scale) and
#'     `hr_outcome_fit`/`hr_outcome_lwr`/`hr_outcome_upr` (HR scale).}
#'   \item{r2}{Adjusted R-squared of the regression.}
#'   \item{ste, ste_log}{STE as surrogate HR and log HR (grid midpoint).}
#'   \item{fit_cross, fit_cross_log}{STE as surrogate HR and log HR
#'     (interpolated).}
#'   \item{data}{The input data.}
#'   \item{plot}{The [ggplot2::ggplot()] object; only present if
#'     `plot = TRUE`. Use it to restyle or save the figure, e.g. with
#'     [ggplot2::ggsave()].}
#' }
#'
#' @references
#' Burzykowski T, Buyse M (2006). Surrogate threshold effect: an alternative
#' measure for meta-analytic surrogate endpoint validation. *Pharmaceutical
#' Statistics* 5(3):173-186.
#'
#' Buyse M, Molenberghs G, Burzykowski T, Renard D, Geys H (2000). The
#' validation of surrogate endpoints in meta-analyses of randomized
#' experiments. *Biostatistics* 1(1):49-67.
#'
#' @examples
#' set.seed(1)
#' n <- 25
#' trials <- data.frame(N = sample(100:900, n, replace = TRUE))
#' trials$log_PFS <- rnorm(n, log(0.8), 0.3)
#' trials$log_OS  <- 0.05 + 0.7 * trials$log_PFS + rnorm(n, 0, 0.08)
#'
#' res <- estimate_ste(trials, xlabel = "HR PFS", ylabel = "HR OS")
#' res$ste
#' res$r2
#'
#' # save the figure
#' # ggplot2::ggsave("ste.png", res$plot, width = 8, height = 6)
#'
#' @export
estimate_ste <- function(data,
                         surrogate = "log_PFS",
                         outcome   = "log_OS",
                         weights   = "N",
                         min.pred  = 0.15,
                         max.pred  = 3,
                         plot      = TRUE,
                         xlabel    = "Hazard ratio PFS",
                         ylabel    = "Hazard ratio OS",
                         hr_breaks = c(0.2, 0.3, 0.5, 0.7, 1.0, 1.5),
                         group_var = NULL,
                         xlim.min  = 0.15,
                         xlim.max  = 2.5,
                         ylim.min  = 0.15,
                         ylim.max  = 2.5,
                         palette   = c("#F3AB97", "#B58494", "grey80", "#D6DD7F")) {

  # input checks
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  missing_cols <- setdiff(c(surrogate, outcome, weights), names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  # weighted simple linear regression
  .ste_wts <- data[[weights]]
  fml <- stats::as.formula(paste0("`", outcome, "` ~ `", surrogate, "`"))
  reg <- stats::lm(fml, data = data, weights = .ste_wts)

  # coefficients with 95% CI
  coef_tbl <- as.data.frame(cbind(stats::coef(reg), stats::confint(reg)))
  colnames(coef_tbl) <- c("fit", "lwr", "upr")

  # 95% prediction interval on a grid of surrogate values
  x_log  <- seq(from = log(min.pred), to = log(max.pred), length.out = 200)
  x_pred <- stats::setNames(data.frame(x_log), surrogate)
  w_pred <- rep(stats::median(.ste_wts), nrow(x_pred))

  pi_mat <- stats::predict(reg, newdata = x_pred, interval = "prediction", weights = w_pred)

  pred_tbl <- data.frame(
    log_surrogate  = x_log,
    hr_surrogate   = exp(x_log),
    fit            = pi_mat[, "fit"],
    lwr            = pi_mat[, "lwr"],
    upr            = pi_mat[, "upr"],
    hr_outcome_fit = exp(pi_mat[, "fit"]),
    hr_outcome_lwr = exp(pi_mat[, "lwr"]),
    hr_outcome_upr = exp(pi_mat[, "upr"])
  )

  # STE: intersection of the upper prediction limit with the line outcome log HR = 0
  fit_cross_log <- tryCatch(
    stats::approx(x = pred_tbl$upr, y = pred_tbl$log_surrogate, xout = 0)$y,
    error = function(e) NA_real_
  )
  fit_cross_hr <- exp(fit_cross_log)

  # grid-midpoint version of the STE
  above <- pred_tbl[pred_tbl$upr > 0, , drop = FALSE]
  below <- pred_tbl[pred_tbl$upr < 0, , drop = FALSE]

  if (nrow(above) > 0 && nrow(below) > 0) {
    before  <- above[which.min(above$upr), ]
    after   <- below[which.max(below$upr), ]
    ste_log <- (before$log_surrogate + after$log_surrogate) / 2
    ste_hr  <- exp(ste_log)
  } else {
    ste_log <- NA_real_
    ste_hr  <- NA_real_
    warning(sprintf(
      "The prediction interval does not cross outcome HR = 1 for surrogate '%s' within HR [%g, %g]; STE not estimable.",
      surrogate, min.pred, max.pred
    ), call. = FALSE)
  }

  ma <- list(
    surrogate     = surrogate,
    outcome       = outcome,
    weights       = weights,
    model         = reg,
    coef          = coef_tbl,
    pred          = pred_tbl,
    r2            = summary(reg)$adj.r.squared,
    ste           = ste_hr,
    ste_log       = ste_log,
    fit_cross     = fit_cross_hr,
    fit_cross_log = fit_cross_log,
    data          = data
  )

  if (plot) {
    ma$plot <- plot_ste(ma, xlabel = xlabel, ylabel = ylabel, group_var = group_var,
                        palette = palette, hr_breaks = hr_breaks, weights = weights,
                        xlim.min = xlim.min, xlim.max = xlim.max,
                        ylim.min = ylim.min, ylim.max = ylim.max)
    print(ma$plot)
  }

  ma
}
