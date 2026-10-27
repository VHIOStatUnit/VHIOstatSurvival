# ==============================================================================
# R/cox_subgroup_forest.R
#
# Publication-quality subgroup forest plots from univariable Cox models,
# in the style commonly used for RCT subgroup analyses in NEJM / Lancet / JCO.
#
# Author: Julian Silan
#
# Package requirements (DESCRIPTION):
#   Depends: R (>= 4.1.0)      # native pipe `|>`
#   Imports: broom, dplyr, forestploter, grid, gtable, purrr, stats,
#            survival, tibble, utils
# ==============================================================================


# ------------------------------------------------------------------------------
# Global variable bindings (columns used via dplyr non-standard evaluation).
# Silences "no visible binding for global variable" notes in R CMD check.
# ------------------------------------------------------------------------------
utils::globalVariables(c(
  "indent", "subgroup", "is_first_of_block", "subgroup_text", "block_id",
  "hr", "conf_low", "conf_high", "p_value", "median_a", "median_b",
  "n", "n_event", "covariate_label", "allpatients_n", "covlabel_n",
  "is_first_of_block_covlabel", "is_first_of_block_subgroup",
  "is_first_of_block_N", "n_fmt", "endpoint_label", "covariate_label_key"
))


# ------------------------------------------------------------------------------
# Design note (read this before using the function)
# ------------------------------------------------------------------------------
# This implements the classic *subgroup treatment-effect* forest plot seen in
# trial publications: for every covariate, the plot shows the treatment-arm
# hazard ratio (Arm A vs Arm B) computed separately within each subgroup
# defined by that covariate (e.g. Male / Female), plus an "All patients" row
# for the overall treatment effect. This matches the Cox model
# Surv(time, event) ~ treatment, fit on data filtered to each subgroup, and is
# what makes "Favours Treatment A / Favours Treatment B" axis labels
# meaningful.
#
# Continuous covariates cannot define discrete "subgroups" for that model, so
# they are handled two ways:
#   1. For the forest-plot rows, the covariate is split at its median into
#      two subgroups ("< median" / ">= median"), and the same treatment-effect
#      model is fit within each half - this keeps every row on the plot
#      using the same HR definition (treatment effect), so the log axis and
#      "Favours A / B" labels stay valid throughout.
#   2. Separately, a stratified Cox model
#      Surv(time, event) ~ covariate + strata(treatment)
#      is fit to obtain the HR *per one-unit increase* in the covariate
#      (adjusted for treatment via stratification). This is returned as a
#      supplementary table (`continuous_effects`) rather than as a forest-plot
#      row, because it estimates a different quantity (the covariate's own
#      prognostic effect) than the treatment-effect rows do.
#
# Excluding subgroup rows from the plot (e.g. "Unknown" categories) is fully
# under caller control via `exclude_levels`, keyed by the raw `covariates`
# names (not the display `labels`). Nothing is excluded automatically.
# ------------------------------------------------------------------------------


# ==============================================================================
# Helper functions
# ==============================================================================

#' Validate the core inputs of cox_subgroup_forest()
#' @noRd
.validate_forest_inputs <- function(data, covariates, labels, treatment, time, event) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame (or tibble).", call. = FALSE)
  }
  if (length(covariates) == 0) {
    stop("`covariates` must contain at least one variable name.", call. = FALSE)
  }
  if (length(time) != length(event)) {
    stop(sprintf(
      "`time` (length %d) and `event` (length %d) must have the same length: one time variable per event variable.",
      length(time), length(event)
    ), call. = FALSE)
  }
  if (length(labels) != length(covariates)) {
    stop(sprintf(
      "`labels` (length %d) and `covariates` (length %d) must have the same length: one display label per covariate.",
      length(labels), length(covariates)
    ), call. = FALSE)
  }
  if (!is.character(treatment) || length(treatment) != 1) {
    stop("`treatment` must be a single character string naming the treatment column.", call. = FALSE)
  }

  all_vars <- unique(c(covariates, treatment, time, event))
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "The following variable(s) were not found in `data`: %s",
      paste(missing_vars, collapse = ", ")
    ), call. = FALSE)
  }

  invisible(TRUE)
}

#' Decide whether a covariate should be treated as categorical or continuous
#' @noRd
.detect_covariate_type <- function(data, covariate, continuous_vars = NULL) {
  if (!is.null(continuous_vars) && covariate %in% continuous_vars) return("continuous")
  x <- data[[covariate]]
  if (is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) > 5) return("continuous")
  "categorical"
}

#' Format an HR (95% CI) text string, e.g. "1.34 (1.05-1.72)"
#' @noRd
.format_hr_ci_subgroup <- function(est, lower, upper, digits = 2) {
  # ASCII hyphen (not an en dash) is used deliberately: the base R "pdf"
  # device's default fonts use WinAnsi/Latin-1 encoding, which cannot render
  # typographic dashes and silently corrupts them ("conversion failure").
  # Use a cairo_pdf / Unicode-aware device if you want a true en dash here.

  fmt_num <- function(x, digits) {
    if (is.na(x)) return(NA_character_)
    if (abs(x) > 1000) sprintf(paste0("%.", digits, "e"), x)
    else               sprintf(paste0("%.", digits, "f"), x)
  }

  mapply(function(e, l, u) {
    if (is.na(e) || is.na(l) || is.na(u)) return(NA_character_)
    sprintf("%s (%s-%s)", fmt_num(e, digits), fmt_num(l, digits), fmt_num(u, digits))
  }, est, lower, upper, USE.NAMES = FALSE)
}

#' Format a p-value to `digits` decimals, showing "<0.001" style thresholds
#' @noRd
.format_pvalue_subgroup <- function(p, digits = 3) {
  threshold <- 10^(-digits)
  dplyr::case_when(
    is.na(p)      ~ "NA",
    p < threshold ~ paste0("<", format(threshold, scientific = FALSE)),
    TRUE          ~ formatC(p, digits = digits, format = "f")
  )
}

#' Exclude specific subgroup rows on a per-covariate basis before plotting
#'
#' @param df A results data frame containing `covariate` and `subgroup`
#'   columns (as produced by `.build_covariate_results()`).
#' @param exclude_levels Optional named list keyed by raw covariate name; each
#'   value is a character vector of subgroup values to drop. `NULL` drops
#'   nothing.
#' @noRd
.filter_excluded_subgroups <- function(df, exclude_levels) {
  if (is.null(exclude_levels) || length(exclude_levels) == 0) return(df)

  drop_flag <- purrr::map2_lgl(df$covariate, df$subgroup, function(cov, sg) {
    excl <- exclude_levels[[cov]]
    !is.null(excl) && sg %in% excl
  })

  df[!drop_flag, , drop = FALSE]
}

#' Fit the treatment-effect Cox model Surv(time, event) ~ treatment within one
#' data subset (e.g. one subgroup). Arm order is fixed via `arm_order` so the
#' reported HR is always Arm A vs Arm B (`arm_order[2]` is set as the reference
#' level, so the model coefficient represents `arm_order[1]`).
#' @noRd
.fit_treatment_cox <- function(data, time_var, event_var, treatment, arm_order, weights) {

  needed <- c(time_var, event_var, treatment)
  d <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  d[[treatment]] <- factor(as.character(d[[treatment]]), levels = arm_order)
  d <- d[!is.na(d[[treatment]]), , drop = FALSE]

  n       <- nrow(d)
  n_event <- suppressWarnings(sum(d[[event_var]], na.rm = TRUE))

  if (n == 0 || dplyr::n_distinct(d[[treatment]]) < 2) {
    return(list(model = NULL, n = n, n_event = n_event))
  }

  # Reference = arm_order[2] ("Arm B") so the estimated HR is Arm A vs Arm B.
  d[[treatment]] <- stats::relevel(d[[treatment]], ref = arm_order[2])

  fmla <- stats::as.formula(
    paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", treatment)
  )
  if (is.null(weights)) {
    fit <- tryCatch(survival::coxph(fmla, data = d), error = function(e) NULL)
  } else {
    fit <- tryCatch(survival::coxph(fmla, data = d, weights = get(weights)), error = function(e) NULL)
  }

  list(model = fit, n = n, n_event = n_event)
}

#' Kaplan-Meier median survival time per treatment arm, within one data subset
#' @noRd
.median_survival_by_arm <- function(data, time_var, event_var, treatment, arm_order, weights) {

  out <- tibble::tibble(arm = arm_order, median = NA_real_)

  needed <- c(time_var, event_var, treatment)
  d <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  d[[treatment]] <- factor(as.character(d[[treatment]]), levels = arm_order)
  d <- d[!is.na(d[[treatment]]), , drop = FALSE]
  if (nrow(d) == 0) return(out)

  fmla <- stats::as.formula(
    paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", treatment)
  )
  if (is.null(weights)) {
    fit <- tryCatch(survival::survfit(fmla, data = d), error = function(e) NULL)
  } else {
    fit <- tryCatch(survival::survfit(fmla, data = d, weights = get(weights)), error = function(e) NULL)
  }

  if (is.null(fit)) return(out)

  s <- summary(fit)$table
  if (is.null(nrow(s))) {
    # Single stratum (only one arm present in this subset)
    rn <- if (!is.null(fit$strata)) names(fit$strata) else paste0(treatment, "=", arm_order[1])
    s <- matrix(s, nrow = 1, dimnames = list(rn, names(s)))
  }
  rn  <- gsub(paste0("^", treatment, "="), "", rownames(s))
  med <- unname(s[, "median"])
  out$median[match(rn, out$arm)] <- med
  out
}

#' Assemble one results row (one subgroup x one endpoint) from a fitted model
#' and its matching median-survival table
#' @noRd
.row_from_fit <- function(covariate, covariate_label, subgroup, indent,
                          endpoint, endpoint_label, fit_info, median_tbl,
                          arm_order) {

  fit <- fit_info$model
  if (is.null(fit)) {
    hr <- conf_low <- conf_high <- p_value <- NA_real_
  } else {
    tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
    hr        <- tt$estimate[1]
    conf_low  <- tt$conf.low[1]
    conf_high <- tt$conf.high[1]
    p_value   <- tt$p.value[1]
  }

  med_a <- median_tbl$median[median_tbl$arm == arm_order[1]]
  med_b <- median_tbl$median[median_tbl$arm == arm_order[2]]
  med_a <- if (length(med_a) == 0) NA_real_ else med_a
  med_b <- if (length(med_b) == 0) NA_real_ else med_b

  tibble::tibble(
    covariate       = covariate,
    covariate_label = covariate_label,
    subgroup        = subgroup,
    indent          = indent,
    endpoint        = endpoint,
    endpoint_label  = endpoint_label,
    n               = fit_info$n,
    n_event         = fit_info$n_event,
    median_a        = med_a,
    median_b        = med_b,
    hr              = hr,
    conf_low        = conf_low,
    conf_high       = conf_high,
    p_value         = p_value
  )
}

#' Build the full results table + fitted models for a single covariate,
#' across all requested endpoints.
#' @noRd
.build_covariate_results <- function(data, covariate, label, treatment,
                                     time, event, endpoint_labels,
                                     arm_order, continuous_vars, weights) {

  type <- .detect_covariate_type(data, covariate, continuous_vars)
  x <- data[[covariate]]

  if (type == "categorical") {
    if (is.factor(x)) {
      levels_vec <- levels(droplevels(x))       # preserve existing factor order
    } else {
      levels_vec <- as.character(sort(unique(stats::na.omit(x))))
    }
    subgroup_group <- as.character(x)
  } else {
    med <- stats::median(x, na.rm = TRUE)
    lab_lo <- sprintf("< median (%.1f)", med)
    lab_hi <- sprintf(">= median (%.1f)", med)
    subgroup_group <- ifelse(is.na(x), NA_character_, ifelse(x >= med, lab_hi, lab_lo))
    levels_vec <- c(lab_lo, lab_hi)
  }

  data_work <- data
  data_work$.subgroup_group <- subgroup_group

  models <- list()
  rows   <- list()

  # Subgroup-major ordering: "All patients" first, then each subgroup level,
  # with every endpoint's row for that subgroup immediately adjacent. This is
  # what lets the plotting step blank out repeated subgroup labels so each
  # subgroup name appears to "span" its block of endpoint rows.
  subgroup_levels <- c("All patients", levels_vec)

  for (lev in subgroup_levels) {
    is_overall <- identical(lev, "All patients")
    sub_data <- if (is_overall) {
      data_work
    } else {
      data_work[!is.na(data_work$.subgroup_group) & data_work$.subgroup_group == lev, , drop = FALSE]
    }

    for (ep in seq_along(time)) {
      time_var  <- time[ep]
      event_var <- event[ep]
      ep_label  <- endpoint_labels[ep]

      fit_info <- .fit_treatment_cox(sub_data, time_var, event_var, treatment, arm_order, weights = weights)
      med_tbl  <- .median_survival_by_arm(sub_data, time_var, event_var, treatment, arm_order, weights = weights)

      key <- paste(label, lev, ep_label, sep = " | ")
      models[[key]] <- fit_info$model

      rows[[length(rows) + 1]] <- .row_from_fit(
        covariate = covariate, covariate_label = label, subgroup = lev,
        indent = !is_overall, endpoint = time_var, endpoint_label = ep_label,
        fit_info = fit_info, median_tbl = med_tbl, arm_order = arm_order
      )
    }
  }

  results <- dplyr::bind_rows(rows)

  # --- supplementary per-unit HR for continuous covariates (not plotted) ---
  continuous_summary <- NULL
  if (type == "continuous") {
    continuous_summary <- purrr::map_dfr(seq_along(time), function(ep) {
      time_var  <- time[ep]
      event_var <- event[ep]
      ep_label  <- endpoint_labels[ep]

      needed <- c(time_var, event_var, treatment, covariate)
      d <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
      if (nrow(d) == 0) return(NULL)

      fmla <- stats::as.formula(paste0(
        "survival::Surv(", time_var, ", ", event_var, ") ~ ", covariate,
        " + survival::strata(", treatment, ")"
      ))

      if (is.null(weights)) {
        fit <- tryCatch(survival::coxph(fmla, data = d), error = function(e) NULL)
      } else {
        fit <- tryCatch(survival::coxph(fmla, data = d, weights = get(weights)), error = function(e) NULL)
      }
      if (is.null(fit)) return(NULL)

      tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
      tibble::tibble(
        covariate = covariate, covariate_label = label,
        endpoint = time_var, endpoint_label = ep_label,
        hr = tt$estimate[1], conf_low = tt$conf.low[1], conf_high = tt$conf.high[1],
        p_value = tt$p.value[1], n = nrow(d),
        interpretation = "HR per 1-unit increase, adjusted for treatment (strata)"
      )
    })
  }

  list(type = type, results = results, models = models, continuous_summary = continuous_summary)
}


#' Apply a row-height multiplier to the body rows of a forestploter object
#'
#' @param p A forestploter/gtable object.
#' @param row_height Numeric multiplier. 1 = default height, 1.3 = 30% taller.
#' @noRd
.apply_row_height <- function(p, row_height = 1) {

  if (length(row_height) != 1 || !is.numeric(row_height) ||
      is.na(row_height) || row_height <= 0) {
    stop("`row_height` must be a single positive numeric value.", call. = FALSE)
  }

  # Locate the body rows. forestploter stores the table as a gtable.
  core_fg <- which(p$layout$name == "core-fg")

  if (length(core_fg) == 0) {
    return(p)
  }

  body_rows <- unique(p$layout$t[core_fg])

  # Only enlarge body rows, not the header / title / legend.
  p$heights[body_rows] <- p$heights[body_rows] * row_height

  p
}

#' Build the endpoint colour legend and append it below a forestploter object
#' @noRd
.add_endpoint_legend <- function(p, endpoints, legend_cols, ci_pch) {

  n_legend <- length(endpoints)
  legend_grobs <- list()

  for (i in seq_along(endpoints)) {

    # x position of this legend entry
    x_pos <- (i - 0.5) / n_legend

    # Point
    legend_grobs[[length(legend_grobs) + 1]] <- grid::pointsGrob(
      x = grid::unit(x_pos - 0.035, "npc"),
      y = grid::unit(0.5, "npc"),
      pch = ci_pch,
      size = grid::unit(2.2, "mm"),
      gp = grid::gpar(col = legend_cols[i], fill = legend_cols[i])
    )

    # Label
    legend_grobs[[length(legend_grobs) + 1]] <- grid::textGrob(
      label = endpoints[i],
      x = grid::unit(x_pos - 0.005, "npc"),
      y = grid::unit(0.5, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 9, col = "grey20")
    )
  }

  legend_grob <- grid::gTree(children = do.call(grid::gList, legend_grobs))

  # Additional row below the forest plot
  p <- gtable::gtable_add_rows(p, heights = grid::unit(0.3, "in"), pos = nrow(p))

  # Legend spanning the full width
  gtable::gtable_add_grob(
    p, legend_grob,
    t = nrow(p), l = 1, r = ncol(p),
    name = "manual_endpoint_legend"
  )
}

#' Build a forestploter forest plot object for one covariate's results table
#' @noRd
.build_forest_plot <- function(results, covariate_label, arm_labels, digits, time_unit, endpoint_palette,
                               xlim, ticks_at, length_plot,
                               ci_pch, ci_lty, ci_lwd, ci_Theight, sizes,
                               exclude_levels = NULL,
                               row_height = 1, label_indent = 4) {

  if (length(row_height) != 1 || !is.numeric(row_height) || is.na(row_height) || row_height <= 0) {
    stop("`row_height` must be a single positive numeric value.", call. = FALSE)
  }

  if (is.null(endpoint_palette)) {
    endpoint_palette <- c("#E67E22", "#16A085", "#8E44AD", "#C0392B", "#2980B9")
  }

  # Rows are already ordered subgroup-major (see .build_covariate_results),
  # i.e. every endpoint's row for a given subgroup is adjacent. That lets us
  # blank out the repeated subgroup label so the subgroup name reads as if it
  # "spans" its block of endpoint rows.
  df <- results |>
    .filter_excluded_subgroups(exclude_levels) |>
    dplyr::mutate(
      subgroup_text     = ifelse(indent, paste0(strrep(" ", label_indent), subgroup), subgroup),
      is_first_of_block = subgroup != dplyr::lag(subgroup, default = "\u0001"),
      subgroup_display  = ifelse(is_first_of_block, subgroup_text, ""),
      block_id          = cumsum(is_first_of_block),
      hr_ci        = .format_hr_ci_subgroup(hr, conf_low, conf_high, digits),
      p_fmt        = .format_pvalue_subgroup(p_value),
      median_a_fmt = ifelse(is.na(median_a), "NR", sprintf("%.1f", median_a)),
      median_b_fmt = ifelse(is.na(median_b), "NR", sprintf("%.1f", median_b)),
      n_fmt        = ifelse(is.na(n), "", as.character(n)),
      event_fmt    = ifelse(is.na(n_event), "", as.character(n_event)),
      blank_col    = strrep(" ", length_plot)
    )

  endpoints      <- unique(df$endpoint_label)
  multi_endpoint <- length(endpoints) > 1
  endpoint_col   <- stats::setNames(
    rep(endpoint_palette, length.out = length(endpoints)),
    endpoints
  )

  # No separate "Endpoint" text column: which endpoint a row belongs to is
  # conveyed by the colour of its point/CI (see recolouring + legend below).
  plot_cols <- c("subgroup_display", "n_fmt", "event_fmt", "median_a_fmt",
                 "median_b_fmt", "blank_col", "hr_ci", "p_fmt")

  tbl <- as.data.frame(df[, plot_cols])
  names(tbl) <- c(
    "Subgroup", "N", "Events",
    paste0("Median ", arm_labels[1], "\n(", time_unit, ")"),
    paste0("Median ", arm_labels[2], "\n(", time_unit, ")"),
    " ",
    "HR (95% CI)",
    "p-value"
  )

  ci_col <- which(names(tbl) == " ")

  if (is.null(xlim)) {
    finite_bounds <- c(df$conf_low, df$conf_high, 1)
    finite_bounds <- finite_bounds[is.finite(finite_bounds) & finite_bounds > 0]
    x_range <- range(finite_bounds, na.rm = TRUE)
    xlim <- c(max(0.05, x_range[1] * 0.8), x_range[2] * 1.25)
  }

  # Choose "nice" tick locations, then clip to xlim so forestploter never
  # warns about (or silently drops) ticks that fall outside the axis range.
  if (is.null(ticks_at)) {
    ticks_at <- pretty(xlim, n = 5)
    ticks_at <- ticks_at[ticks_at > 0 & ticks_at >= xlim[1] & ticks_at <= xlim[2]]
    if (length(ticks_at) < 2) ticks_at <- signif(exp(seq(log(xlim[1]), log(xlim[2]), length.out = 4)), 2)
  }

  # Row shading is applied per subgroup *block* (all of that subgroup's
  # endpoint rows share one shade), not per individual row.
  shaded_rows <- which(df$block_id %% 2 == 1)

  ft_theme <- forestploter::forest_theme(
    base_size      = 11,
    ci_pch         = ci_pch,
    ci_col         = endpoint_col[[1]],
    ci_lty         = ci_lty,
    ci_lwd         = ci_lwd,
    ci_Theight     = ci_Theight,
    axis_extend    = c(0.2, 0.05),
    vertline_lwd   = 1,
    vertline_lty   = "dashed",
    vertline_col   = "grey20",
    refline_gp     = grid::gpar(lwd = 1, lty = "dashed", col = "grey30"),
    footnote_gp    = grid::gpar(col = "grey20", fontface = "italic", cex = 0.8),
    core = list(bg_params = list(fill = "white"))   # banding applied manually below
  )

  p <- forestploter::forest(
    tbl,
    est       = df$hr,
    lower     = df$conf_low,
    upper     = df$conf_high,
    ci_column = ci_col,
    sizes     = sizes,
    ref_line  = 1,
    x_trans   = "log",
    xlim      = xlim,
    ticks_at  = ticks_at,
    arrow_lab = c(paste("Favours", arm_labels[1]), paste("Favours", arm_labels[2])),
    title     = covariate_label,
    theme     = ft_theme
  )

  p <- .apply_row_height(p, row_height = row_height)

  # Shade alternating subgroup blocks
  if (length(shaded_rows) > 0) {
    p <- forestploter::edit_plot(p, row = shaded_rows, which = "background",
                                 gp = grid::gpar(fill = "grey95", col = "grey95"))
  }

  # Bold the "All patients" summary row(s)
  bold_rows <- which(df$subgroup == "All patients")
  if (length(bold_rows) > 0) {
    p <- forestploter::edit_plot(p, row = bold_rows, gp = grid::gpar(fontface = "bold"))
  }

  # Vertical separator lines between table columns
  for (col_i in c(0, ci_col - 1, ci_col, length(names(tbl)))) {
    p <- forestploter::add_border(
      p, row = 1:nrow(tbl), col = col_i, part = "body", where = "right",
      gp = grid::gpar(col = "grey30", lwd = 1)
    )
  }

  # Vertical separator lines in header
  for (col_i in c(0, ci_col - 1, ci_col, length(names(tbl)))) {
    p <- forestploter::add_border(
      p, row = 1, col = col_i, part = "header", where = "right",
      gp = grid::gpar(col = "grey30", lwd = 1)
    )
  }

  # Horizontal lines
  for (row_i in c(0, nrow(tbl))) {
    p <- forestploter::add_border(
      p, row = row_i, col = 1:ncol(tbl), part = "body", where = "bottom",
      gp = grid::gpar(col = "grey30", lwd = 1)
    )
  }

  p <- forestploter::add_border(
    p, row = 0, col = 1:ncol(tbl), part = "header", where = "bottom",
    gp = grid::gpar(col = "grey30", lwd = 1)
  )

  # Colour the point/CI for each endpoint, and add a small legend,
  # when there is more than one endpoint.
  if (multi_endpoint) {
    for (ep_label in endpoints) {
      rows_i <- which(df$endpoint_label == ep_label)
      col_i  <- endpoint_col[[ep_label]]
      p <- forestploter::edit_plot(p, row = rows_i, col = ci_col, which = "ci",
                                   gp = grid::gpar(col = col_i, fill = col_i))
    }
    p <- .add_endpoint_legend(p, endpoints, unname(endpoint_col[endpoints]), ci_pch)
  }

  p
}


#' Build a combined forest plot object for all covariates
#' @noRd
.build_combined_forest_plot <- function(results, labels, arm_labels, digits, time_unit,
                                        endpoint_palette, xlim, ticks_at, length_plot,
                                        ci_pch, ci_lty, ci_lwd, ci_Theight, sizes,
                                        exclude_levels = NULL, row_height = 1,
                                        covariate_label_row = FALSE, col_spacing = c(0, 0),
                                        label_indent = 4) {

  if (length(row_height) != 1 || !is.numeric(row_height) || is.na(row_height) || row_height <= 0) {
    stop("`row_height` must be a single positive numeric value.", call. = FALSE)
  }

  if (is.null(endpoint_palette)) {
    endpoint_palette <- c("#E67E22", "#16A085", "#8E44AD", "#C0392B", "#2980B9")
  }

  df_results <- dplyr::bind_rows(results)
  endpoints <- unique(df_results$endpoint_label)

  df <- df_results |>
    .filter_excluded_subgroups(exclude_levels) |>
    dplyr::mutate(allpatients_n = cumsum(subgroup == "All patients")) |>
    dplyr::filter(subgroup != "All patients" | allpatients_n <= length(endpoints)) |>
    dplyr::mutate(covariate_label = dplyr::if_else(subgroup == "All patients", "All patients", covariate_label),
                  subgroup = dplyr::if_else(subgroup == "All patients", "", subgroup)) |>
    dplyr::select(-allpatients_n) |>
    dplyr::group_by(covariate_label) |>
    dplyr::mutate(covlabel_n = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      is_first_of_block_covlabel = covariate_label != dplyr::lag(covariate_label, default = "\u0001"),
      covariate_label_display  = ifelse(is_first_of_block_covlabel,
                                        if (covariate_label_row) "" else paste(covariate_label, "        "), ""),
      is_first_of_block_subgroup = subgroup != dplyr::lag(subgroup, default = "\u0001") |
        covariate_label != dplyr::lag(covariate_label, default = "\u0001"),
      subgroup_display  = ifelse(is_first_of_block_subgroup, subgroup, ""),
      block_id          = cumsum(is_first_of_block_subgroup),
      hr_ci        = .format_hr_ci_subgroup(hr, conf_low, conf_high, 2),
      p_fmt        = .format_pvalue_subgroup(p_value),
      median_a_fmt = ifelse(is.na(median_a), "NR", sprintf("%.1f", median_a)),
      median_b_fmt = ifelse(is.na(median_b), "NR", sprintf("%.1f", median_b)),
      n_fmt        = ifelse(is.na(n), "", as.character(n)),
      is_first_of_block_N = n_fmt != dplyr::lag(n_fmt, default = "\u0001"),
      N_display  = ifelse(is_first_of_block_N, n_fmt, ""),
      event_fmt    = ifelse(is.na(n_event), "", as.character(n_event)),
      blank_col    = strrep(" ", length_plot)
    )

  if (covariate_label_row) {
    # Insert a header row immediately before the first row of each covariate.
    df$is_covariate_header <- FALSE

    # split() on a factor with levels in order of appearance keeps the original
    # covariate order (a plain character split would sort alphabetically).
    cov_factor <- factor(df$covariate_label, levels = unique(df$covariate_label))

    df <- purrr::map_dfr(split(df, cov_factor), function(block) {
      header <- block[1, , drop = FALSE]
      header$covariate_label_display <- unique(block$covariate_label)[1]
      header$subgroup_display <- ""
      header$N_display <- ""
      header$event_fmt <- ""
      header$median_a_fmt <- ""
      header$median_b_fmt <- ""
      header$blank_col <- strrep(" ", length_plot)
      header$hr_ci <- ""
      header$p_fmt <- ""
      header$hr <- NA_real_
      header$conf_low <- NA_real_
      header$conf_high <- NA_real_
      header$is_covariate_header <- TRUE

      block$is_covariate_header <- FALSE
      # The header row now marks the start of the covariate; body rows must not
      # be treated as covariate starts (they would otherwise be bolded).
      block$is_first_of_block_covlabel <- FALSE

      dplyr::bind_rows(header, block)
    })
  } else {
    df$is_covariate_header <- FALSE
  }

  endpoints      <- unique(df$endpoint_label)
  multi_endpoint <- length(endpoints) > 1
  endpoint_col   <- stats::setNames(
    rep(endpoint_palette, length.out = length(endpoints)),
    endpoints
  )

  # No separate "Endpoint" text column: which endpoint a row belongs to is
  # conveyed by the colour of its point/CI (see recolouring + legend below).
  if (covariate_label_row) {
    # Single "Variable" column: covariate label on its own bold header row,
    # subgroup labels indented underneath it in the same column.
    # (col_spacing[1] has no meaning here: there is no Variable/Subgroup gap.)
    df$variable_display <- ifelse(
      df$is_covariate_header,
      df$covariate_label_display,
      ifelse(df$subgroup_display != "",
             paste0(strrep(" ", label_indent), df$subgroup_display), "")
    )
    plot_cols <- "variable_display"
    col_names <- "Variable"
  } else {
    plot_cols <- "covariate_label_display"
    col_names <- "Variable"

    if (col_spacing[1] > 0) {
      df$spacer_var_subgroup <- strrep(" ", col_spacing[1])
      plot_cols <- c(plot_cols, "spacer_var_subgroup")
      col_names <- c(col_names, "")
    }

    plot_cols <- c(plot_cols, "subgroup_display")
    col_names <- c(col_names, "Subgroup")
  }

  if (col_spacing[2] > 0) {
    df$spacer_subgroup_n <- strrep(" ", col_spacing[2])
    plot_cols <- c(plot_cols, "spacer_subgroup_n")
    col_names <- c(col_names, "  ")
  }

  plot_cols <- c(plot_cols, "N_display", "event_fmt", "median_a_fmt",
                 "median_b_fmt", "blank_col", "hr_ci", "p_fmt")
  col_names <- c(col_names, "N", "Events",
                 paste0("Median ", arm_labels[1], "\n(", time_unit, ")"),
                 paste0("Median ", arm_labels[2], "\n(", time_unit, ")"),
                 " ", "HR (95% CI)", "p-value")

  tbl <- as.data.frame(df[, plot_cols])
  names(tbl) <- col_names

  ci_col <- which(names(tbl) == " ")

  if (is.null(xlim)) {
    finite_bounds <- c(df$conf_low, df$conf_high, 1)
    finite_bounds <- finite_bounds[is.finite(finite_bounds) & finite_bounds > 0]
    x_range <- range(finite_bounds, na.rm = TRUE)
    xlim <- c(max(0.05, x_range[1] * 0.8), x_range[2] * 1.25)
  }

  # Choose "nice" tick locations, then clip to xlim so forestploter never
  # warns about (or silently drops) ticks that fall outside the axis range.
  if (is.null(ticks_at)) {
    ticks_at <- pretty(xlim, n = 5)
    ticks_at <- ticks_at[ticks_at > 0 & ticks_at >= xlim[1] & ticks_at <= xlim[2]]
    if (length(ticks_at) < 2) ticks_at <- signif(exp(seq(log(xlim[1]), log(xlim[2]), length.out = 4)), 2)
  }

  # Row shading is applied per *covariate* block (all rows belonging to
  # every second variable share one shade), which visually separates
  # consecutive variables when many covariates are stacked.
  if (covariate_label_row) {
    covariate_block_id <- cumsum(df$is_covariate_header)
  } else {
    covariate_block_id <- cumsum(df$is_first_of_block_covlabel)
  }

  shaded_rows <- which(covariate_block_id %% 2 == 1)

  ft_theme <- forestploter::forest_theme(
    base_size      = 11,
    ci_pch         = ci_pch,
    ci_col         = endpoint_col[[1]],
    ci_lty         = ci_lty,
    ci_lwd         = ci_lwd,
    ci_Theight     = ci_Theight,
    axis_extend    = c(0.2, 0.05),
    vertline_lwd   = 1,
    vertline_lty   = "dashed",
    vertline_col   = "grey20",
    refline_gp     = grid::gpar(lwd = 1, lty = "dashed", col = "grey30"),
    footnote_gp    = grid::gpar(col = "grey20", fontface = "italic", cex = 0.8),
    core = list(bg_params = list(fill = "white"))   # banding applied manually below
  )

  p <- forestploter::forest(
    tbl,
    est       = df$hr,
    lower     = df$conf_low,
    upper     = df$conf_high,
    sizes     = sizes,
    ci_column = ci_col,
    ref_line  = 1,
    x_trans   = "log",
    xlim      = xlim,
    ticks_at  = ticks_at,
    arrow_lab = c(paste("Favours", arm_labels[1]), paste("Favours", arm_labels[2])),
    theme     = ft_theme
  )

  p <- .apply_row_height(p, row_height = row_height)

  if (covariate_label_row) {
    covariate_header_rows <- which(df$is_covariate_header)

    if (length(covariate_header_rows) > 0) {
      p <- forestploter::edit_plot(
        p, row = covariate_header_rows, col = 1, which = "text",
        gp = grid::gpar(fontface = "bold")
      )
    }
  }

  # Shade alternating covariate blocks
  if (length(shaded_rows) > 0) {
    p <- forestploter::edit_plot(p, row = shaded_rows, which = "background",
                                 gp = grid::gpar(fill = "grey95", col = "grey95"))
  }

  # Bold the "All patients" summary row(s)
  allpatients_row <- which(df$covariate_label == "All patients")
  if (length(allpatients_row) > 0) {
    p <- forestploter::edit_plot(p, row = allpatients_row, gp = grid::gpar(fontface = "bold"))
  }

  # Bold the Variable column for the other covariates
  covlabel_rows <- which(df$is_first_of_block_covlabel & df$covariate_label != "All patients")
  if (length(covlabel_rows) > 0) {
    p <- forestploter::edit_plot(p, row = covlabel_rows, col = 1,
                                 gp = grid::gpar(fontface = "bold"))
  }

  # Vertical separator lines between table columns
  for (col_i in c(0, ci_col - 1, ci_col, length(names(tbl)))) {
    p <- forestploter::add_border(
      p, row = 1:nrow(tbl), col = col_i, part = "body", where = "right",
      gp = grid::gpar(col = "grey30", lwd = 1)
    )
  }

  # Vertical separator lines in header
  for (col_i in c(0, ci_col - 1, ci_col, length(names(tbl)))) {
    p <- forestploter::add_border(
      p, row = 1, col = col_i, part = "header", where = "right",
      gp = grid::gpar(col = "grey30", lwd = 1)
    )
  }

  # Horizontal lines
  for (row_i in c(0, nrow(tbl))) {
    p <- forestploter::add_border(
      p, row = row_i, col = 1:ncol(tbl), part = "body", where = "bottom",
      gp = grid::gpar(col = "grey30", lwd = 1)
    )
  }

  p <- forestploter::add_border(
    p, row = 0, col = 1:ncol(tbl), part = "header", where = "bottom",
    gp = grid::gpar(col = "grey30", lwd = 1)
  )

  # Colour the point/CI for each endpoint, and add a small legend,
  # when there is more than one endpoint.
  if (multi_endpoint) {
    for (ep_label in endpoints) {
      rows_i <- which(df$endpoint_label == ep_label)
      col_i  <- endpoint_col[[ep_label]]
      p <- forestploter::edit_plot(p, row = rows_i, col = ci_col, which = "ci",
                                   gp = grid::gpar(col = col_i, fill = col_i))
    }
    p <- .add_endpoint_legend(p, endpoints, unname(endpoint_col[endpoints]), ci_pch)
  }

  p
}


# ==============================================================================
# Main function
# ==============================================================================

#' Publication-quality subgroup forest plots from univariable Cox models
#'
#' Fits, for each requested covariate, a treatment-effect Cox proportional
#' hazards model (\code{Surv(time, event) ~ treatment}) separately within
#' every subgroup defined by that covariate (plus an overall "All patients"
#' row), and renders the results as a publication-style forest plot in the
#' convention used for RCT subgroup analyses.
#' Kaplan-Meier median survival time is also reported per subgroup and
#' treatment arm.
#'
#' @param data A data.frame or tibble containing all variables referenced
#'   below.
#' @param covariates Character vector of covariate names in \code{data} to
#'   analyse (one forest plot is produced per covariate).
#' @param labels Character vector of display labels, one per element of
#'   \code{covariates} (same length and order).
#' @param treatment Character string naming the treatment/arm variable in
#'   \code{data}. Must have exactly two levels.
#' @param time Character vector of one or more survival time variable names.
#' @param event Character vector of event indicator variable names, the same
#'   length as \code{time} (paired by position: \code{time[i]}/\code{event[i]}).
#' @param endpoint_labels Optional character vector of display labels for
#'   each (time, event) pair. Defaults to the \code{time} variable names.
#' @param continuous_vars Optional character vector identifying which
#'   \code{covariates} should be treated as continuous. If \code{NULL},
#'   numeric covariates with more than 5 distinct values are auto-detected
#'   as continuous.
#' @param arm_order Optional length-2 character vector giving the two
#'   treatment-arm levels in the order \code{c("Arm A", "Arm B")}. The
#'   reported hazard ratio is always Arm A vs Arm B (HR < 1 favours Arm A).
#'   Defaults to the levels of \code{treatment} in their existing order
#'   (factor level order, or sorted order otherwise).
#' @param treatment_labels Optional length-2 character vector of display
#'   labels for the two arms (used in column headers and "Favours ..."
#'   annotations). Defaults to the raw \code{arm_order} values when they are
#'   informative (i.e. \code{treatment} is not purely numeric), and to
#'   generic "Treatment A" / "Treatment B" otherwise.
#' @param digits Number of decimal places for HR/CI display. Default 2.
#' @param time_unit Character string describing the time unit for median
#'   survival columns (e.g. "months", "days"). Default "months".
#' @param palette Optional character vector of colours, one per endpoint
#'   (same length as \code{endpoint_labels}). If \code{NULL} (default), a
#'   built-in palette of five colours is used (recycled if necessary).
#' @param xlim Optional length-2 numeric vector with the limits of the
#'   (log-scaled) HR axis. If \code{NULL}, derived from the confidence limits.
#' @param ticks_at Optional numeric vector of HR axis tick positions. If
#'   \code{NULL}, "nice" ticks within \code{xlim} are chosen automatically.
#' @param length_plot Integer. Width of the blank column reserved for the CI
#'   graphic, expressed as a number of spaces. Default 40.
#' @param weights Optional character string naming a weight column in
#'   \code{data}, passed to \code{survival::coxph()} and
#'   \code{survival::survfit()}. Default \code{NULL} (unweighted).
#' @param ci_pch Point symbol for the HR estimate. Default 15 (filled square).
#' @param ci_lty Line type of the confidence interval. Default 1.
#' @param ci_lwd Line width of the confidence interval. Default 3.
#' @param ci_Theight Height of the confidence interval end caps. Default 0.4.
#' @param sizes Size of the HR point estimate. Default 0.6.
#' @param exclude_levels Optional named list for dropping specific
#'   subgroup rows from the forest plots on a per-covariate basis. Names must
#'   match entries of \code{covariates} (the raw variable names - not
#'   \code{labels}); each value is a character vector of \code{subgroup}
#'   values to exclude for that covariate, e.g.
#'   \code{list(sex = "Unknown", ecog_ps = "Unknown")}. Covariates not listed
#'   are unaffected. Nothing is excluded by default (\code{NULL}). This only
#'   affects the plots; \code{results_table} and \code{median_survival_table}
#'   still contain all rows.
#' @param row_height Numeric multiplier controlling the height of body rows.
#'   \code{1} uses the default forestploter row height; values above 1
#'   increase the row height proportionally, e.g. \code{1.3} makes rows
#'   approximately 30% taller.
#' @param covariate_label_row Logical. If \code{TRUE}, each covariate label is
#'   displayed in a separate bold row immediately before its subgroup rows,
#'   and the subgroup labels are shown indented in the same "Variable" column
#'   (there is no separate "Subgroup" column). If \code{FALSE} (default), the
#'   covariate label is displayed in the first column alongside the first
#'   subgroup, with subgroups in their own column.
#' @param col_spacing Numeric vector of length 2 giving the width (in spaces)
#'   of an extra spacer column between the Variable and Subgroup columns, and
#'   between the Subgroup and N columns, respectively. \code{0} (default)
#'   adds no spacer. If \code{covariate_label_row = TRUE} there is no
#'   Subgroup column, so only \code{col_spacing[2]} is used (spacer between
#'   the Variable and N columns) and \code{col_spacing[1]} is ignored.
#' @param label_indent Number of spaces used to indent the subgroup labels
#'   (below the covariate label when \code{covariate_label_row = TRUE}, and
#'   below the "All patients" row in the per-covariate plots). Default 4.
#'
#' @return A list with elements:
#' \describe{
#'   \item{models}{Named list of lists of fitted \code{coxph} model objects
#'     (one nested list per covariate, keyed by "label | subgroup | endpoint").}
#'   \item{all_result}{Named list of the per-covariate results tables.}
#'   \item{results_table}{Tidy data frame with one row per
#'     covariate x subgroup x endpoint, including HR, 95 percent CI, p-value, N,
#'     events, and median survival per arm.}
#'   \item{median_survival_table}{Subset of \code{results_table} with just
#'     the median survival columns, for easy tabulation.}
#'   \item{plots}{Named list of forestploter plot objects, one per covariate
#'     (name = covariate label). Draw with \code{grid::grid.draw()}.}
#'   \item{combined_plots}{A single forestploter plot object containing all
#'     covariates.}
#'   \item{continuous_effects}{Data frame of supplementary per-unit hazard
#'     ratios for continuous covariates (see Details), or \code{NULL} if
#'     none were analysed.}
#'   \item{arm_order}{The two treatment arm levels, in the order used for
#'     the HR calculation (Arm A vs Arm B).}
#'   \item{arm_labels}{The display labels used for the two arms.}
#' }
#'
#' @details
#' For every covariate the treatment-effect HR (Arm A vs Arm B) is estimated
#' separately within each subgroup, so the "Favours ..." axis labels are
#' valid for every row. Continuous covariates are split at their median for
#' the plot; the HR per one-unit increase (stratified by treatment) is
#' returned separately in \code{continuous_effects}. Subgroups are fitted
#' independently: no interaction test is computed.
#'
#' @examples
#' \dontrun{
#' result <- cox_subgroup_forest(
#'   data = my_trial_data,
#'   covariates = c("sex", "age"),
#'   labels = c("Sex", "Age"),
#'   treatment = "arm",
#'   time = "os_time",
#'   event = "os_event",
#'   continuous_vars = "age",
#'   exclude_levels = list(sex = "Unknown")
#' )
#' grid::grid.draw(result$plots[["Sex"]])
#' grid::grid.draw(result$combined_plots)
#' }
#' @export
cox_subgroup_forest <- function(data, covariates, labels, treatment, time, event,
                                endpoint_labels = NULL, continuous_vars = NULL,
                                arm_order = NULL, treatment_labels = NULL,
                                digits = 2, time_unit = "months", palette = NULL,
                                xlim = NULL, ticks_at = NULL, length_plot = 40,
                                weights = NULL, ci_pch = 15, ci_lty = 1,
                                ci_lwd = 3, ci_Theight = 0.4, sizes = 0.6,
                                exclude_levels = NULL,
                                row_height = 1,
                                covariate_label_row = FALSE,
                                col_spacing = c(0, 0),
                                label_indent = 4) {

  .validate_forest_inputs(data, covariates, labels, treatment, time, event)

  if (!is.numeric(label_indent) || length(label_indent) != 1 ||
      is.na(label_indent) || label_indent < 0) {
    stop("`label_indent` must be a single non-negative number.", call. = FALSE)
  }

  if (is.null(endpoint_labels)) endpoint_labels <- time
  if (length(endpoint_labels) != length(time)) {
    stop("`endpoint_labels` must be the same length as `time`.", call. = FALSE)
  }
  if (!is.null(palette) && length(endpoint_labels) != length(palette)) {
    stop("`palette` must be the same length as `endpoint_labels`.", call. = FALSE)
  }

  # --- resolve the two treatment arms and their display labels ---
  trt_raw      <- data[[treatment]]
  was_numeric  <- is.numeric(trt_raw)
  trt_levels   <- if (is.factor(trt_raw)) {
    levels(droplevels(trt_raw))
  } else {
    as.character(sort(unique(stats::na.omit(trt_raw))))
  }
  if (length(trt_levels) != 2) {
    stop(sprintf(
      "`treatment` variable '%s' must have exactly 2 levels for a two-arm subgroup forest plot; found %d (%s).",
      treatment, length(trt_levels), paste(trt_levels, collapse = ", ")
    ), call. = FALSE)
  }

  if (is.null(arm_order)) arm_order <- trt_levels
  if (length(arm_order) != 2 || !all(arm_order %in% trt_levels)) {
    stop("`arm_order` must be a length-2 character vector matching the levels of `treatment`.", call. = FALSE)
  }

  if (!is.null(treatment_labels)) {
    if (length(treatment_labels) != 2) {
      stop("`treatment_labels` must be a length-2 character vector.", call. = FALSE)
    }
    arm_labels <- treatment_labels
  } else if (was_numeric) {
    arm_labels <- c("Treatment A", "Treatment B")
  } else {
    arm_labels <- arm_order
  }

  # --- loop over covariates ---
  all_results <- list()
  all_models  <- list()
  all_plots   <- list()

  continuous_summaries <- list()

  for (i in seq_along(covariates)) {
    cov <- covariates[i]
    lab <- labels[i]

    res <- .build_covariate_results(
      data = data, covariate = cov, label = lab, treatment = treatment,
      time = time, event = event, endpoint_labels = endpoint_labels,
      arm_order = arm_order, continuous_vars = continuous_vars, weights = weights
    )

    all_results[[lab]] <- res$results
    all_models[[lab]]  <- res$models
    if (!is.null(res$continuous_summary)) continuous_summaries[[lab]] <- res$continuous_summary

    all_plots[[lab]] <- .build_forest_plot(
      results = res$results, covariate_label = lab, arm_labels = arm_labels,
      digits = digits, time_unit = time_unit,
      endpoint_palette = palette,
      xlim = xlim, ticks_at = ticks_at, length_plot = length_plot,
      ci_pch = ci_pch, ci_lty = ci_lty, ci_lwd = ci_lwd,
      ci_Theight = ci_Theight, sizes = sizes,
      exclude_levels = exclude_levels,
      row_height = row_height,
      label_indent = label_indent
    )
  }

  combined_plot <- .build_combined_forest_plot(
    results = all_results, labels = labels, arm_labels = arm_labels,
    digits = digits, time_unit = time_unit,
    endpoint_palette = palette,
    xlim = xlim, ticks_at = ticks_at, length_plot = length_plot,
    ci_pch = ci_pch, ci_lty = ci_lty, ci_lwd = ci_lwd,
    ci_Theight = ci_Theight, sizes = sizes,
    exclude_levels = exclude_levels,
    row_height = row_height,
    covariate_label_row = covariate_label_row,
    col_spacing = col_spacing,
    label_indent = label_indent
  )

  results_table <- dplyr::bind_rows(all_results, .id = "covariate_label_key")

  median_survival_table <- results_table |>
    dplyr::select(covariate_label, subgroup, endpoint_label, n, n_event,
                  median_a, median_b)
  names(median_survival_table)[match(c("median_a", "median_b"), names(median_survival_table))] <-
    paste0("median_", arm_labels)

  list(
    models                = all_models,
    all_result            = all_results,
    results_table         = dplyr::select(results_table, -dplyr::all_of("covariate_label_key")),
    median_survival_table = median_survival_table,
    plots                 = all_plots,
    combined_plots        = combined_plot,
    continuous_effects    = if (length(continuous_summaries) > 0) dplyr::bind_rows(continuous_summaries) else NULL,
    arm_order             = arm_order,
    arm_labels            = arm_labels
  )
}
