# =============================================================================
# Customised ggsurvplot
# =============================================================================

#' Customised Kaplan-Meier survival plot with embedded HR/median/survival-rate table
#'
#' @description
#' Wraps \code{survminer::ggsurvplot()} to build a publication-ready
#' Kaplan-Meier plot with a risk table underneath and a compact summary
#' table (hazard ratio, p-value, median survival, and/or survival rate at
#' user-defined time points) annotated directly onto the KM panel.
#'
#' Internally this fits (or accepts pre-computed) survival estimates via
#' \code{survival::survfit()}, fits a Cox proportional hazards model
#' (\code{survival::coxph()}) or an accelerated failure time model
#' (\code{survival::survreg()}), builds the KM plot and risk table via
#' \code{survminer::ggsurvplot()}, restyles both, and stitches everything
#' together with \code{\link{.ggsurvhrplot}}.
#'
#' @param surv.formula A survival formula, e.g. \code{Surv(time, event) ~ arm},
#'   passed to \code{survfit()} and to the Cox/AFT model.
#' @param weights Optional. Case weights passed to \code{survfit()} and to the
#'   Cox/AFT model.
#' @param dat A data.frame containing the variables in \code{surv.formula}.
#' @param table.df Optional. A pre-computed data.frame of HR estimates with
#'   columns \code{HR}, \code{conf.low}, \code{conf.high}, and (if
#'   \code{p.value_df = TRUE}) \code{p.val}. When supplied, no Cox/AFT model
#'   is fit internally -- \code{.ggsurvhrplot()} uses this table directly.
#' @param labels Optional character vector of group labels, used as legend labels and
#'   as the row labels of the summary table. Must be ordered to match the
#'   strata of \code{surv.formula} / \code{table.df}.
#' @param p.value_df Logical. If \code{table.df} is supplied, whether it
#'   contains a \code{p.val} column to show in the table. Default \code{FALSE}.
#' @param labels.table.unadj Logical. If \code{TRUE}, the risk table shown
#'   underneath the KM plot is refit without \code{weights} (i.e. an
#'   unadjusted/unweighted risk table), while the KM curves themselves remain
#'   as fit in step 1. Default \code{FALSE}.
#' @param aft.model Logical. If \code{TRUE}, fits a Weibull AFT model
#'   (\code{survival::survreg()}) instead of a Cox model. Default \code{FALSE}.
#' @param palette Optional colour vector. When NULL, uses ggplot2's default
#'   discrete hue palette for both KM curves and HR-table line segments.
#' @param risk.table.title.font Font spec (size, face, colour) for the risk
#'   table title, passed to \code{\link{.customize_labels}}. Default
#'   \code{c(10, "bold")}.
#' @param risk.table.ytext.font Font spec for the risk table y-axis tick
#'   labels, passed to \code{\link{.customize_labels}}. Default \code{c(10)}.
#' @param survfit.args Named list of additional arguments passed to
#'   \code{survival::survfit()}.
#' @param model.args Named list of additional arguments passed to
#'   \code{survival::coxph()} (or \code{survival::survreg()} if
#'   \code{aft.model = TRUE}).
#' @param tableplot.args Named list of additional arguments passed through to
#'   \code{\link{.ggsurvhrplot}} (e.g. \code{position}, \code{table_width},
#'   \code{median.show}, \code{surv.times}, font-size and line-width
#'   controls).
#'   \code{position}: One of \code{"topright"} (default), \code{"topleft"},
#'   \code{"bottomright"}, \code{"bottomleft"} — the corner of the KM panel
#'   the table is anchored to.
#'   \code{x_margin} Numeric. Horizontal offset of the table, as a fraction of
#'   the x-axis range. Positive shifts the table right, negative shifts it
#'   left. Default \code{0.02}.
#'   \code{y_margin} Numeric. Vertical offset of the table, as a fraction of
#'   the y-axis range. Positive shifts the table up, negative shifts it
#'   down. Default \code{0.02}.
#'   \code{table_width} Numeric. Width of the table, as a fraction of the
#'   x-axis range. Default \code{0.38}.
#'   \code{table_height} Numeric. Height of the table, as a fraction of the
#'   y-axis range. Default \code{0.16}.
#'   \code{hr.table.core.cex} Numeric. Font size (cex) of the table's data
#'   cells (HR/CI/p-value/median/survival-rate). Default \code{0.7}.
#'   \code{hr.table.colhead.cex} Numeric. Font size (cex) of the table's
#'   column headers. Default \code{0.8}.
#'   \code{hr.table.rowhead.cex} Numeric. Font size (cex) of the table's row
#'   labels (group names). Default \code{0.75}.
#'   \code{hr.table.line.lwd} Numeric. Line width of the colour-coded group
#'   segments in the table's leftmost column. Default \code{2.5}.
#'   \code{hr.table.label.width} Optional. A \code{grid::unit} specifying the
#'   width of the row-label (group name) column, e.g. \code{unit(3, "cm")}.
#'   Default \code{NULL} (auto-sized).
#'   \code{surv.times} Optional numeric vector of time points (in the same
#'   time units as the survival formula used to build \code{survfit}). If
#'   supplied, adds one \code{S(t)} column per time point to the
#'   table, giving the Kaplan-Meier survival probability (as a percentage),
#'   extracted from \code{survfit} via
#'   \code{\link{.extract_surv_at_times}}. Requires \code{survfit} to be
#'   supplied. Default \code{NULL} (no survival-rate columns).
#' @param ggsurvplot.args Named list of additional arguments passed through to
#'   \code{survminer::ggsurvplot()}.
#'
#' @return A \code{patchwork} object: the annotated KM plot stacked above the
#'   risk table (or just the annotated KM plot if no risk table is available).
#'
#' @examples
#' \dontrun{
#' surv.plot(
#'   surv.formula = survival::Surv(time, status) ~ arm,
#'   dat          = trial_data,
#'   labels       = c("Control", "Treatment"),
#'   palette      = c("dodgerblue3", "firebrick"),
#'   tableplot.args = list(median.show = TRUE, surv.times = c(6, 12, 24))
#' )
#' }
#'
#' @export
surv_plot <- function(
    surv.formula = NULL,
    weights = NULL,
    dat = NULL,
    table.df = NULL,
    labels = NULL,
    p.value_df = FALSE,
    labels.table.unadj = FALSE,
    aft.model = FALSE,
    palette = NULL,
    auto.HR = TRUE,
    cox.interaction = TRUE,
    risk.table.title.font = c(10, "bold"),
    risk.table.ytext.font = c(10),
    survfit.args = list(),
    model.args = list(),
    tableplot.args = list(),
    ggsurvplot.args = list()
) {

  requireNamespace("survminer", quietly = TRUE)

  if (is.null(surv.formula) || !inherits(surv.formula, "formula")) {
    stop("Please supply 'surv.formula' as a survival formula.")
  }
  if (is.null(dat) || !is.data.frame(dat)) {
    stop("Please supply 'dat' as a data.frame.")
  }

  predictors <- all.vars(surv.formula[[3]])

  # Fit Kaplan-Meier curves
  survfit.base.args <- list(formula = surv.formula, data = dat)
  if (!is.null(weights)) survfit.base.args$weights <- weights

  svf <- do.call(
    survival::survfit,
    c(survfit.base.args, survfit.args)
  )

  n.groups <- if (!is.null(svf$strata)) length(svf$strata) else 1L

  # Use the same colors for the curves and table segments
  if (is.null(palette)) {
    palette <- scales::hue_pal()(max(1L, n.groups))
  }

  if (length(palette) < n.groups) {
    stop("'palette' must contain at least one color for each KM stratum.")
  }

  # Keep the KM legend labels separate from the HR table labels
  km.labels <- labels
  hr.labels <- labels

  # Default labels when the caller does not supply any
  if (is.null(hr.labels) && !is.null(table.df)) {
    if ("Subgroups" %in% names(table.df)) {
      hr.labels <- c("Ref.", as.character(table.df$Subgroups))
    } else if (!is.null(svf$strata)) {
      hr.labels <- names(svf$strata)
    } else {
      hr.labels <- c("Ref.", paste0("Group ", seq_len(nrow(table.df))))
    }
  } else if (is.null(hr.labels) && is.null(table.df)) {
    hr.labels <- if (!is.null(svf$strata)) {
      names(svf$strata)
    } else {
      "Overall"
    }
  }

  # Automatically generate an HR table for formulas with 2+ predictors
  if (is.null(table.df) &&
      isTRUE(auto.HR) &&
      length(predictors) >= 2L &&
      !isTRUE(aft.model)) {

    auto.hr.results <- .make_emmeans_hr_table(
      surv.formula = surv.formula,
      dat = dat,
      weights = weights,
      model.args = model.args,
      interaction = cox.interaction
    )

    table.df <- auto.hr.results$table.df

    # Custom labels are used in the HR table when supplied
    if (!is.null(labels)) {
      hr.labels <- as.character(labels)
    } else {
      hr.labels <- c("Ref.", as.character(table.df$Subgroups))
    }
  }

  # Fit the model if no precomputed HR table was provided/generated
  cox.model <- NULL
  aft.mod <- NULL

  if (is.null(table.df)) {
    model.base.args <- list(formula = surv.formula, data = dat)
    if (!is.null(weights)) model.base.args$weights <- weights

    if (!isTRUE(aft.model)) {
      cox.model <- do.call(
        survival::coxph,
        c(model.base.args, model.args)
      )
    } else {
      aft.mod <- do.call(
        survival::survreg,
        c(model.base.args, list(dist = "weibull"), model.args)
      )
    }
  }

  # Build KM plot. Custom labels go to the legend here.
  ggsurvplot.base.args <- list(
    fit = svf,
    data = dat,
    legend.labs = km.labels,
    palette = palette
  )

  psurv <- do.call(
    survminer::ggsurvplot,
    c(ggsurvplot.base.args, ggsurvplot.args)
  )

  # Style the KM panel
  psurv$plot <- psurv$plot +
    ggplot2::annotate(
      "segment",
      x = 0, xend = 45,
      y = 50, yend = 50,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.2
    ) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(face = "bold"),
      axis.line = ggplot2::element_line(linewidth = 0.2),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_text(
        face = "bold",
        margin = ggplot2::margin(t = 10)
      ),
      axis.title.y = ggplot2::element_text(
        face = "bold",
        margin = ggplot2::margin(r = 10)
      ),
      plot.margin = ggplot2::margin(6, 6, 6, 27)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(linewidth = 0.8)
      )
    )

  # Style risk table
  if (!is.null(psurv$table)) {
    psurv$table <- .customize_labels(
      psurv$table,
      font.title = risk.table.title.font,
      font.ytickslab = risk.table.ytext.font
    )
    psurv$table$theme$plot.title$margin <-
      ggplot2::margin(0, 0, 6, -40)
    psurv$table$theme$text$hjust <- 0
  }

  # Optional unweighted risk table
  if (isTRUE(labels.table.unadj)) {
    svf.unadj <- do.call(survival::survfit, list(formula = surv.formula, data = dat))

    if (!"risk.table" %in% names(ggsurvplot.args)) {
      ggsurvplot.args$risk.table <- TRUE
    }
    if (!"risk.table.title" %in% names(ggsurvplot.args)) {
      ggsurvplot.args$risk.table.title <- "Unadj. No. at Risk"
    }

    unadj.args <- list(
      fit = svf.unadj,
      data = dat,
      legend.labs = km.labels,
      palette = palette
    )

    psurv.unadj <- do.call(
      survminer::ggsurvplot,
      c(unadj.args, ggsurvplot.args)
    )

    if (!is.null(psurv.unadj$table)) {
      psurv.unadj$table <- .customize_labels(
        psurv.unadj$table,
        font.title = risk.table.title.font,
        font.ytickslab = risk.table.ytext.font
      )
      psurv.unadj$table$theme$plot.title$margin <-
        ggplot2::margin(0, 0, 6, -40)
      psurv.unadj$table$theme$text$hjust <- 0
    }

    active.table <- psurv.unadj$table
  } else {
    active.table <- psurv$table
  }

  # HR table labels are passed here
  base.hrplot.args <- list(
    ggsurvplot = psurv,
    ggsurvtable = active.table,
    labels = hr.labels,
    col = palette,
    survfit = svf
  )

  if (!is.null(table.df)) {
    do.call(
      .ggsurvhrplot,
      c(
        base.hrplot.args,
        list(data.table = table.df, p.value = p.value_df),
        tableplot.args
      )
    )
  } else {
    do.call(
      .ggsurvhrplot,
      c(
        base.hrplot.args,
        list(cox_model = cox.model, aft_model = aft.mod),
        tableplot.args
      )
    )
  }
}



### =============================================================================
### Function to add the ggsurvplot KM, table and Cox model together
### =============================================================================

# -----------------------------------------------------------
# Helping function to extract the statistic of the Cox model
# -----------------------------------------------------------

# Automatic HR table using emmeans
.make_emmeans_hr_table <- function(
    surv.formula,
    dat,
    weights = NULL,
    model.args = list(),
    interaction = TRUE
) {

  predictors <- all.vars(surv.formula[[3]])

  if (length(predictors) < 2L) {
    stop("At least two predictors are required for automatic emmeans HRs.")
  }

  response <- paste(deparse(surv.formula[[2]]), collapse = "")

  rhs <- if (isTRUE(interaction)) {
    paste(predictors, collapse = " * ")
  } else {
    paste(predictors, collapse = " + ")
  }

  cox.formula <- stats::as.formula(
    paste(response, "~", rhs),
    env = environment(surv.formula)
  )

  cox.base.args <- list(formula = cox.formula, data = dat)
  if (!is.null(weights)) cox.base.args$weights <- weights

  cox.model <- do.call(
    survival::coxph,
    c(cox.base.args, model.args)
  )

  emm.formula <- stats::as.formula(
    paste("~", paste(predictors, collapse = " * ")),
    env = environment(surv.formula)
  )

  emm <- as.data.frame(
    emmeans::emmeans(cox.model, specs = emm.formula)
  )

  lower.col <- intersect(c("asymp.LCL", "lower.CL"), names(emm))[1]
  upper.col <- intersect(c("asymp.UCL", "upper.CL"), names(emm))[1]

  if (is.na(lower.col) || is.na(upper.col)) {
    stop("Could not identify confidence-interval columns in emmeans output.")
  }
  if (!"emmean" %in% names(emm)) {
    stop("The emmeans output does not contain the expected 'emmean' column.")
  }

  ord.args <- lapply(emm[predictors], function(x) {
    if (is.factor(x)) as.integer(x) else x
  })

  ord <- do.call(order, c(ord.args, list(na.last = TRUE)))
  emm <- emm[ord, , drop = FALSE]

  emm$HR <- exp(emm$emmean)
  emm$conf.low <- exp(emm[[lower.col]])
  emm$conf.high <- exp(emm[[upper.col]])
  emm$Subgroups <- do.call(paste, c(emm[predictors], sep = " + "))

  # The table function adds the reference row itself, so omit the first
  # combination from this table.
  table.df <- emm |>
    dplyr::mutate(
      `HR (95% CI)` = sprintf(
        "%.2f (%.2f, %.2f)",
        HR, conf.low, conf.high
      )
    ) |>
    dplyr::select(
      Subgroups, HR, SE, conf.low, conf.high, `HR (95% CI)`
    ) |>
    dplyr::slice(-1)

  df.plot <- table.df |>
    dplyr::select(Subgroups, HR, conf.low, conf.high)

  list(
    cox.model = cox.model,
    table.df = table.df,
    df.plot = df.plot
  )
}

#' Format a numeric statistic for table display
#'
#' @description
#' Rounds a number to 2 decimal places, or switches to scientific notation
#' (1 significant digit) when the integer part has more than 4 digits.
#'
#' @param x A numeric value.
#'
#' @return A formatted number (rounded numeric, or a character string in
#'   scientific notation for very large values).
#'
#' @keywords internal
#' @noRd
.format_number <- function(x) {
  if (is.na(x)) return("NR")
  digits_before <- nchar(as.character(floor(abs(x))))

  if (digits_before > 4) {
    format(x, scientific = TRUE, digits = 1)
  } else {
    round(x, 2)
  }
}


# -----------------------------------------------------------
# Helping function to extract survival probabilities at
# user-specified time points from a survfit object
# -----------------------------------------------------------

#' Extract survival probabilities at specified time points
#'
#' @description
#' Internal helper that queries a \code{survival::survfit} object at one or
#' more user-supplied time points and returns the survival probability
#' (with 95% CI) for each stratum, formatted as percentages ready for table
#' display.
#'
#' @param survfit A \code{survfit} object (single- or multi-group).
#' @param times Numeric vector of time points at which to evaluate the
#'   survival function.
#'
#' @return A character matrix with one row per stratum (group) and one
#'   column per requested time point. Each cell is a formatted string of the
#'   form \code{"XX.X% (XX.X% - XX.X%)"}, or \code{"NR"} where the
#'   estimate/bound is not reported (e.g. beyond the last observed event).
#'
#' @keywords internal
#' @noRd
.extract_surv_at_times <- function(survfit, times) {

  s <- summary(survfit, times = times, extend = TRUE)

  fmt_pct <- function(x) ifelse(is.na(x), "NR", base::paste0(round(x * 100, 1), "%"))

  has_strata <- !is.null(s$strata)

  if (has_strata) {
    groups <- names(survfit$strata)
    out <- matrix(
      NA_character_,
      nrow = length(groups), ncol = length(times),
      dimnames = list(groups, as.character(times))
    )
    for (g in groups) {
      idx <- which(s$strata == g)
      idx <- idx[order(s$time[idx])]
      # out[g, ] <- base::paste0(
      #   fmt_pct(s$surv[idx]),  " (",
      #   fmt_pct(s$lower[idx]), " - ",
      #   fmt_pct(s$upper[idx]), ")"
      # )

      out[g, ] <- fmt_pct(s$surv[idx])
    }
  } else {
    idx <- order(s$time)
    # out <- matrix(
    #   base::paste0(
    #     fmt_pct(s$surv[idx]),  " (",
    #     fmt_pct(s$lower[idx]), " - ",
    #     fmt_pct(s$upper[idx]), ")"
    #   ),
    #   nrow = 1, dimnames = list("1", as.character(times))
    # )

    out <- matrix(fmt_pct(s$surv[idx]), nrow = 1, dimnames = list("1", as.character(times)))
  }

  out
}


#' Assemble the KM plot with an annotated HR / median / survival-rate table
#'
#' @description
#' Internal workhorse called by \code{\link{surv.plot}}. Takes a fitted
#' \code{survminer::ggsurvplot} object plus either a Cox model, an AFT model,
#' or a pre-computed HR data.frame, builds a compact summary
#' \code{gridExtra::tableGrob} (hazard ratio, p-value, median survival time,
#' and/or survival probability at user-supplied time points), colour-codes
#' it by group, and annotates it directly onto the KM panel at the requested
#' corner. The risk table (if any) is stacked underneath via
#' \code{patchwork}.
#'
#' @param cox_model A fitted \code{coxph} object. Provide exactly one of
#'   \code{cox_model}, \code{aft_model}, or \code{data.table}.
#' @param aft_model A fitted \code{survreg} object (Weibull AFT).
#' @param data.table A data.frame of pre-computed HR estimates with columns
#'   \code{HR}, \code{conf.low}, \code{conf.high}, and (if \code{p.value =
#'   TRUE}) \code{p.val}. The first row is treated as the reference group.
#' @param p.value Logical. When \code{data.table} is supplied, whether it
#'   includes a \code{p.val} column to display. Ignored (always \code{TRUE})
#'   when \code{cox_model} or \code{aft_model} is used. Default \code{TRUE}.
#' @param survfit A \code{survival::survfit} object, required when
#'   \code{median.show = TRUE} and/or \code{surv.times} is supplied, used to
#'   extract median survival times and/or survival probabilities at
#'   specified time points.
#' @param median.show Logical. If \code{TRUE}, adds a \code{Median (95CI)}
#'   column to the table, extracted from \code{survfit}. Default \code{FALSE}.
#' @param HR.show Logical. If \code{TRUE} (default), adds the
#'   \code{HR (95CI)} column (and \code{p-value} column, if applicable) to
#'   the table. If \code{FALSE}, the table shows only \code{Median} and/or
#'   \code{surv.times} columns.
#' @param ggsurvplot A \code{ggsurvplot} object (as returned by
#'   \code{survminer::ggsurvplot()}) whose \code{$plot} panel the summary
#'   table is annotated onto.
#' @param ggsurvtable Optional. A risk table \code{ggplot} to stack beneath
#'   the annotated KM plot. Falls back to \code{ggsurvplot$table} if
#'   \code{NULL}.
#' @param labels Optional character vector of group labels used as table row names.
#'   Must be ordered to match the strata of \code{cox_model} / \code{aft_model}
#'   / \code{data.table} / \code{survfit}.
#' @param col Optional character vector of colours, one per group, used for the
#'   colour-coded line segments in the table's leftmost column.
#' @param position One of \code{"topright"} (default), \code{"topleft"},
#'   \code{"bottomright"}, \code{"bottomleft"} — the corner of the KM panel
#'   the table is anchored to.
#' @param x_margin Numeric. Horizontal offset of the table, as a fraction of
#'   the x-axis range. Positive shifts the table right, negative shifts it
#'   left. Default \code{0.02}.
#' @param y_margin Numeric. Vertical offset of the table, as a fraction of
#'   the y-axis range. Positive shifts the table up, negative shifts it
#'   down. Default \code{0.02}.
#' @param table_width Numeric. Width of the table, as a fraction of the
#'   x-axis range. Default \code{0.38}.
#' @param table_height Numeric. Height of the table, as a fraction of the
#'   y-axis range. Default \code{0.16}.
#' @param hr.table.core.cex Numeric. Font size (cex) of the table's data
#'   cells (HR/CI/p-value/median/survival-rate). Default \code{0.7}.
#' @param hr.table.colhead.cex Numeric. Font size (cex) of the table's
#'   column headers. Default \code{0.8}.
#' @param hr.table.rowhead.cex Numeric. Font size (cex) of the table's row
#'   labels (group names). Default \code{0.75}.
#' @param hr.table.line.lwd Numeric. Line width of the colour-coded group
#'   segments in the table's leftmost column. Default \code{2.5}.
#' @param hr.table.label.width Optional. A \code{grid::unit} specifying the
#'   width of the row-label (group name) column, e.g. \code{unit(3, "cm")}.
#'   Default \code{NULL} (auto-sized).
#' @param surv.times Optional numeric vector of time points (in the same
#'   time units as the survival formula used to build \code{survfit}). If
#'   supplied, adds one \code{S(t)} column per time point to the
#'   table, giving the Kaplan-Meier survival probability (as a percentage),
#'   extracted from \code{survfit} via
#'   \code{\link{.extract_surv_at_times}}. Requires \code{survfit} to be
#'   supplied. Default \code{NULL} (no survival-rate columns).
#'
#' @return A \code{patchwork} object: the annotated KM plot stacked above the
#'   risk table (or just the annotated KM plot if \code{ggsurvtable} and
#'   \code{ggsurvplot$table} are both \code{NULL}).
#'
#' @keywords internal
#' @noRd
.ggsurvhrplot <- function(cox_model    = NULL,
                          aft_model    = NULL,
                          data.table   = NULL,      # alternative to cox/aft model
                          p.value      = TRUE,      # only used when data.table is provided
                          survfit      = NULL,      # survfit object for median/surv.times extraction
                          median.show  = FALSE,     # add a Median (95% CI) column to the table
                          HR.show      = TRUE,      # add the HR and p-value in the table
                          ggsurvplot   = NULL,
                          ggsurvtable  = NULL,
                          labels       = c("Trt A", "Trt B"),
                          col          = c("dodgerblue3", "firebrick"),
                          position     = "topright",
                          x_margin     = 0.02,
                          y_margin     = 0.02,
                          table_width  = 0.38,
                          table_height = 0.16,
                          # Font sizes of the table
                          hr.table.core.cex    = 0.7,   # font size of HR/CI/p-value/median/surv-rate cells
                          hr.table.colhead.cex = 0.8,   # font size of column headers
                          hr.table.rowhead.cex = 0.75,  # font size of row labels (treatment names)
                          # Thickness of the colored group-line segment
                          hr.table.line.lwd    = 2.5,
                          hr.table.label.width = NULL,  # width of the subgroup label column, e.g. unit(3, "cm")
                          # Survival rate at user-defined time points
                          surv.times   = NULL   # e.g. c(6, 12, 24) -> adds one S(t) column per time point
) {

  position <- match.arg(position, c("topright", "topleft", "bottomright", "bottomleft"))

  # =====================================================
  # Input validation
  # =====================================================
  n_sources <- sum(!is.null(cox_model), !is.null(aft_model), !is.null(data.table))
  if (n_sources == 0 && HR.show) stop("Provide one of: cox_model, aft_model, or data.table when HR.show = TRUE.")
  if (n_sources > 1) stop("Provide only one of: cox_model, aft_model, or data.table.")

  if (median.show && is.null(survfit)) {
    stop("A survfit object must be supplied via the 'survfit' argument when median.show = TRUE.")
  }

  if (!is.null(surv.times) && is.null(survfit)) {
    stop("A survfit object must be supplied via the 'survfit' argument when surv.times is specified.")
  }

  # =====================================================
  # Extract statistics — Cox model
  # =====================================================
  if (!is.null(cox_model) && HR.show) {
    cox_summary <- summary(cox_model)
    hr      <- sapply(cox_summary$coefficients[, "exp(coef)"], .format_number)
    ci_low  <- sapply(cox_summary$conf.int[, "lower .95"],     .format_number)
    ci_high <- sapply(cox_summary$conf.int[, "upper .95"],     .format_number)
    pval     <- cox_summary$coefficients[, "Pr(>|z|)"]
    pval_txt <- ifelse(pval < 0.01, "<0.01", round(pval, 2))
    p.value  <- TRUE   # always show p-value from model
  }

  # =====================================================
  # Extract statistics — AFT model
  # =====================================================
  if (!is.null(aft_model) && HR.show) {
    aft_summary <- summary(aft_model)
    names_var   <- names(coef(aft_model))[-1]
    scale       <- aft_model$scale
    coef_val <- vapply(names_var, function(v) coef(aft_model)[v],          numeric(1))
    se       <- vapply(names_var, function(v) sqrt(vcov(aft_model)[v, v]), numeric(1))
    hr      <- sapply(exp(-coef_val / scale),                .format_number)
    ci_low  <- sapply(exp(-(coef_val + 1.96 * se) / scale), .format_number)
    ci_high <- sapply(exp(-(coef_val - 1.96 * se) / scale), .format_number)
    pval     <- vapply(names_var, function(v) aft_summary$table[v, "p"], numeric(1))
    pval_txt <- ifelse(pval < 0.01, "<0.01", round(pval, 2))
    p.value  <- TRUE
  }

  # =====================================================
  # Extract statistics — data.frame input
  # Expected columns: HR, conf.low, conf.high, (p.val optional)
  # First row is treated as reference and prepended automatically
  # =====================================================
  if (!is.null(data.table) && HR.show) {
    required <- c("HR", "conf.low", "conf.high")
    if (!all(required %in% names(data.table))) {
      stop("data.table must contain HR, conf.low and conf.high.")
    }
    if (nrow(data.table) != length(labels) - 1L) {
      stop("data.table must have one row per non-reference group (in strata order).")
    }
    hr      <- sapply(data.table$HR,        .format_number)
    ci_low  <- sapply(data.table$conf.low,  .format_number)
    ci_high <- sapply(data.table$conf.high, .format_number)

    if (p.value) {
      if (!"p.val" %in% names(data.table))
        stop("data.table must contain a 'p.val' column when p.value = TRUE.")
      pval     <- data.table$p.val
      pval_txt <- ifelse(pval < 0.01, "<0.01", round(pval, 2))
    }
  }

  # =====================================================
  # Extract median survival times (if requested)
  # =====================================================
  if (median.show) {
    sv_sum <- summary(survfit)$table

    # summary()$table returns a matrix (multi-group) or named vector (single group)
    if (is.matrix(sv_sum)) {
      med_vals  <- sv_sum[, "median"]
      med_lo    <- sv_sum[, "0.95LCL"]
      med_hi    <- sv_sum[, "0.95UCL"]
    } else {
      # single-group survfit — wrap in a 1-row structure
      med_vals  <- sv_sum["median"]
      med_lo    <- sv_sum["0.95LCL"]
      med_hi    <- sv_sum["0.95UCL"]
    }

    # Helper: format a single median value / CI bound (NA → "NR")
    fmt_med <- function(x) ifelse(is.na(x), "NR", .format_number(x))

    med_txt <- mapply(
      function(m, lo, hi)
        base::paste0(fmt_med(m), " (", fmt_med(lo), " - ", fmt_med(hi), ")"),
      med_vals, med_lo, med_hi
    )
  }

  # =====================================================
  # Extract survival probabilities at user-specified time
  # points (if requested)
  #
  # NB: row order of surv_time_mat follows the strata levels
  # of `survfit`, mirroring the existing assumption used for
  # median.show above — this must match the order of `labels`.
  # =====================================================
  if (!is.null(surv.times)) {
    surv_time_mat <- .extract_surv_at_times(survfit, surv.times)
  }

  # =====================================================
  # Build table data frame
  # Columns are assembled in a fixed left-to-right order:
  #   Median  ->  S(t) at each user time point  ->  HR  ->  p-value
  # =====================================================
  hr_col <- if (HR.show) c("Ref.", base::paste0(hr, " (", ci_low, " - ", ci_high, ")")) else NULL

  cols <- list()

  if (median.show) {
    cols[["Median (95%CI)"]] <- med_txt
  }

  if (!is.null(surv.times)) {
    for (i in seq_along(surv.times)) {
      col_name         <- base::paste0("S(", surv.times[i], ")")
      cols[[col_name]] <- surv_time_mat[, i]
    }
  }

  if (HR.show) {
    cols[["HR (95%CI)"]] <- hr_col
    if (p.value) {
      cols[["p-value"]] <- c("-", pval_txt)
    }
  } else if (length(cols) == 0) {
    # A summary table needs at least one enabled statistic.
    stop("Set median.show = TRUE or supply surv.times when HR.show = FALSE.")
  }

  table_df <- as.data.frame(cols, check.names = FALSE)
  n_groups <- nrow(table_df)
  if (is.null(labels)) labels <- paste0("Group ", seq_len(n_groups))
  if (length(labels) != n_groups) stop("Labels and table rows do not match.")
  if (is.null(col)) col <- scales::hue_pal()(n_groups)
  if (length(col) != n_groups) stop("Colours and table rows do not match.")



  tbl <- gridExtra::tableGrob(
    table_df,
    rows  = labels,
    theme = gridExtra::ttheme_minimal(
      core    = list(fg_params = list(cex = hr.table.core.cex)),
      colhead = list(fg_params = list(fontface = "bold",  cex = hr.table.colhead.cex)),
      rowhead = list(fg_params = list(fontface = "plain", cex = hr.table.rowhead.cex, hjust = 0, x = 0.1))
    )
  )

  if (!is.null(hr.table.label.width)) {
    tbl$widths[1] <- hr.table.label.width
  }
  # =====================================================
  # Inject color segments as extra leftmost gtable column
  # =====================================================
  n_groups <- length(labels)

  line_grobs <- c(
    list(grid::nullGrob()),
    lapply(seq_len(n_groups), function(i) {
      grid::segmentsGrob(
        x0 = grid::unit(0.1, "npc"), x1 = grid::unit(0.9, "npc"),
        y0 = grid::unit(0.5, "npc"), y1 = grid::unit(0.5, "npc"),
        gp = grid::gpar(col = col[i], lwd = hr.table.line.lwd)
      )
    })
  )

  tbl <- gtable::gtable_add_cols(tbl, widths = grid::unit(1.5, "char"), pos = 0)
  for (row_i in seq_along(line_grobs)) {
    tbl <- gtable::gtable_add_grob(
      tbl,
      grobs = line_grobs[[row_i]],
      t = row_i, b = row_i,
      l = 1,     r = 1
    )
  }

  # =====================================================
  # Compute annotation_custom coordinates from axis ranges
  #
  # x_margin / y_margin sign convention (position-independent):
  #   positive x_margin → shift table RIGHT
  #   negative x_margin → shift table LEFT
  #   positive y_margin → shift table UP
  #   negative y_margin → shift table DOWN
  # =====================================================
  built   <- ggplot2::ggplot_build(ggsurvplot$plot)
  x_range <- built$layout$panel_params[[1]]$x.range
  y_range <- built$layout$panel_params[[1]]$y.range

  x_span <- diff(x_range)
  y_span <- diff(y_range)

  xm <- x_span * x_margin   # signed offset in data units
  ym <- y_span * y_margin   # signed offset in data units
  tw <- x_span * table_width
  th <- y_span * table_height

  if (position == "topright") {
    # anchor: top-right corner; positive xm pushes further right, positive ym pushes further up
    xmax <- x_range[2] - (x_span * 0.02) + xm;  xmin <- xmax - tw
    ymax <- y_range[2] - (y_span * 0.02) + ym;  ymin <- ymax - th

  } else if (position == "topleft") {
    # anchor: top-left corner; positive xm pushes further right, positive ym pushes further up
    xmin <- x_range[1] + (x_span * 0.02) + xm;  xmax <- xmin + tw
    ymax <- y_range[2] - (y_span * 0.02) + ym;  ymin <- ymax - th

  } else if (position == "bottomright") {
    # anchor: bottom-right corner; positive xm pushes further right, positive ym pushes further up
    xmax <- x_range[2] - (x_span * 0.02) + xm;  xmin <- xmax - tw
    ymin <- y_range[1] + (y_span * 0.02) + ym;  ymax <- ymin + th

  } else if (position == "bottomleft") {
    # anchor: bottom-left corner; positive xm pushes further right, positive ym pushes further up
    xmin <- x_range[1] + (x_span * 0.02) + xm;  xmax <- xmin + tw
    ymin <- y_range[1] + (y_span * 0.02) + ym;  ymax <- ymin + th
  }

  # =====================================================
  # Inject table and combine with risk table
  # =====================================================
  final_plot <- ggsurvplot$plot +
    ggplot2::annotation_custom(
      grob = tbl,
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax
    )

  bottom <- if (!is.null(ggsurvtable)) ggsurvtable else ggsurvplot$table
  #if (!is.null(bottom)) final_plot / bottom + patchwork::plot_layout(heights = c(5, 1)) else final_plot
  if (!is.null(bottom)) patchwork::wrap_plots(final_plot, bottom, ncol = 1, heights = c(5, 1)) else final_plot
}


# =============================================================================
# A function to customize the labels of the ggsurvplot table
# =============================================================================

#' Apply custom fonts to ggsurvplot risk-table elements
#'
#' @description
#' Internal helper (adapted from \code{survminer:::.ggsurvplot_customize})
#' that applies markdown-aware font specs (via \code{ggtext}) to the title,
#' subtitle, caption, axis titles, and axis tick labels of one or more
#' \code{ggplot} objects — typically the risk table panel of a
#' \code{ggsurvplot} object.
#'
#' @param p A single \code{ggplot} object, or a list of \code{ggplot}
#'   objects (e.g. \code{ggsurvplot$table} when multiple tables are present).
#' @param font.title Font spec (size, face, colour) for the plot title, e.g.
#'   \code{c(10, "bold")}. \code{NULL} leaves it unchanged.
#' @param font.subtitle Font spec for the plot subtitle. \code{NULL} leaves
#'   it unchanged.
#' @param font.caption Font spec for the plot caption. \code{NULL} leaves it
#'   unchanged.
#' @param font.x Font spec for the x-axis title. \code{NULL} leaves it
#'   unchanged.
#' @param font.y Font spec for the y-axis title. \code{NULL} leaves it
#'   unchanged.
#' @param font.xtickslab Font spec for the x-axis tick labels. \code{NULL}
#'   leaves it unchanged.
#' @param font.ytickslab Font spec for the y-axis tick labels. \code{NULL}
#'   leaves it unchanged.
#'
#' @return The restyled \code{ggplot} object, or list of objects, matching
#'   the input type of \code{p}.
#'
#' @keywords internal
#' @noRd
.customize_labels <- function (p, font.title = NULL,
                               font.subtitle = NULL, font.caption = NULL,
                               font.x = NULL, font.y = NULL, font.xtickslab = NULL,
                               font.ytickslab = NULL){
  original.p <- p
  if(ggplot2::is_ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    if (is.null(font)) return(NULL)
    font <- list(size = suppressWarnings(as.numeric(font[1])),
                 face = if (length(font) >= 2L) font[2] else "plain",
                 color = if (length(font) >= 3L) font[3] else "black")
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in seq_along(list.plots)){
    p <- list.plots[[i]]
    if(ggplot2::is_ggplot(p)){
      if (!is.null(font.title)) p <- p + ggplot2::theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + ggplot2::theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + ggplot2::theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + ggplot2::theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + ggplot2::theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + ggplot2::theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + ggplot2::theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(ggplot2::is_ggplot(original.p)) list.plots[[1]]
  else list.plots
}
