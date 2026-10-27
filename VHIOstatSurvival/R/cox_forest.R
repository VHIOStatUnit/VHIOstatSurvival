# ==============================================================================
# cox_forest.R
#
# Publication-quality forest plots for univariable and/or multivariable Cox
# proportional hazards models across one or more endpoints.
# ==============================================================================


# ==============================================================================
# Internal helpers
# ==============================================================================

#' @noRd
.validate_cox_inputs <- function(data, covariates, labels, time, event, weights = NULL) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame (or tibble).", call. = FALSE)
  }
  if (length(covariates) == 0) {
    stop("`covariates` must contain at least one variable name.", call. = FALSE)
  }
  if (length(labels) != length(covariates)) {
    stop(sprintf(
      "`labels` (length %d) and `covariates` (length %d) must have the same length.",
      length(labels), length(covariates)
    ), call. = FALSE)
  }
  if (length(time) != length(event)) {
    stop(sprintf(
      "`time` (length %d) and `event` (length %d) must have the same length.",
      length(time), length(event)
    ), call. = FALSE)
  }
  if (!is.null(weights) && (!is.character(weights) || length(weights) != 1)) {
    stop("`weights` must be NULL or a single column name (character).", call. = FALSE)
  }

  all_vars <- unique(c(covariates, time, event, weights))
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "The following variable(s) were not found in `data`: %s",
      paste(missing_vars, collapse = ", ")
    ), call. = FALSE)
  }

  invisible(TRUE)
}

#' Decide whether a covariate is categorical or continuous
#' @noRd
.detect_covariate_type <- function(data, covariate, continuous_vars = NULL) {
  if (!is.null(continuous_vars) && covariate %in% continuous_vars) return("continuous")
  x <- data[[covariate]]
  if (is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) > 5) return("continuous")
  "categorical"
}

#' Coerce to factor, moving `ref_level` to the front if supplied
#' @noRd
.factor_with_ref <- function(x, ref_level = NULL) {
  f <- if (is.factor(x)) droplevels(x) else factor(as.character(x))
  if (!is.null(ref_level)) {
    if (!ref_level %in% levels(f)) {
      stop(sprintf("Reference level '%s' not found among observed levels (%s).",
                   ref_level, paste(levels(f), collapse = ", ")), call. = FALSE)
    }
    f <- stats::relevel(f, ref = ref_level)
  }
  f
}

#' Format "HR (lower-upper)". "Reference" is printed ONLY for reference rows;
#' non-reference rows without an estimate (e.g. too few patients/events) get
#' `na_label`. ASCII hyphen on purpose: the base pdf() device corrupts true
#' en dashes (WinAnsi encoding).
#' @noRd
.format_hr_ci <- function(est, lower, upper, is_reference, digits = 2,
                          na_label = "NA") {
  fmt <- paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)")
  dplyr::case_when(
    is_reference                             ~ "Reference",
    is.na(est) | is.na(lower) | is.na(upper) ~ na_label,
    TRUE                                     ~ sprintf(fmt, est, lower, upper)
  )
}

#' @noRd
.format_pvalue <- function(p, digits = 3) {
  threshold <- 10^(-digits)
  dplyr::case_when(
    is.na(p)      ~ "",
    p < threshold ~ paste0("<", format(threshold, scientific = FALSE)),
    TRUE          ~ formatC(p, digits = digits, format = "f")
  )
}

#' Set excluded levels to NA so complete-case filtering removes them
#' @noRd
.apply_exclude_levels <- function(data, exclude_levels) {
  if (is.null(exclude_levels) || length(exclude_levels) == 0) return(data)
  for (cov in names(exclude_levels)) {
    if (!cov %in% names(data)) next
    lvls_to_drop <- exclude_levels[[cov]]
    data[[cov]][as.character(data[[cov]]) %in% lvls_to_drop] <- NA
  }
  data
}

#' @noRd
.level_counts <- function(data, covariate, event_var, level = NULL) {
  d <- if (is.null(level)) data else data[as.character(data[[covariate]]) == level, , drop = FALSE]
  list(n = nrow(d), n_event = suppressWarnings(sum(d[[event_var]], na.rm = TRUE)))
}

#' @noRd
.rows_from_categorical_term <- function(tidy_fit, model_data, covariate, covariate_label,
                                        event_var, fac_levels, ref_level) {

  non_ref_levels <- setdiff(fac_levels, ref_level)
  term_names <- paste0(covariate, non_ref_levels)

  level_rows <- purrr::map_dfr(seq_along(non_ref_levels), function(i) {
    lev <- non_ref_levels[i]
    trm <- term_names[i]
    hit <- tidy_fit[tidy_fit$term == trm, , drop = FALSE]
    cnt <- .level_counts(model_data, covariate, event_var, lev)
    tibble::tibble(
      covariate = covariate, covariate_label = covariate_label,
      level = lev, is_reference = FALSE,
      n = cnt$n, n_event = cnt$n_event,
      hr = if (nrow(hit) == 1) hit$estimate[1] else NA_real_,
      conf_low = if (nrow(hit) == 1) hit$conf.low[1] else NA_real_,
      conf_high = if (nrow(hit) == 1) hit$conf.high[1] else NA_real_,
      p_value = if (nrow(hit) == 1) hit$p.value[1] else NA_real_
    )
  })

  ref_cnt <- .level_counts(model_data, covariate, event_var, ref_level)
  ref_row <- tibble::tibble(
    covariate = covariate, covariate_label = covariate_label,
    level = ref_level, is_reference = TRUE,
    n = ref_cnt$n, n_event = ref_cnt$n_event,
    hr = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_
  )

  dplyr::bind_rows(ref_row, level_rows)
}

#' @noRd
.row_from_continuous_term <- function(tidy_fit, model_data, covariate, covariate_label,
                                      event_var, increase = 1) {
  hit <- tidy_fit[tidy_fit$term == covariate, , drop = FALSE]
  cnt <- .level_counts(model_data, covariate, event_var, level = NULL)

  tibble::tibble(
    covariate = covariate, covariate_label = covariate_label,
    level = paste0("per ", increase, "-unit increase"),
    is_reference = FALSE,
    n = cnt$n, n_event = cnt$n_event,
    hr = if (nrow(hit) == 1) hit$estimate[1]^increase else NA_real_,
    conf_low = if (nrow(hit) == 1) hit$conf.low[1]^increase else NA_real_,
    conf_high = if (nrow(hit) == 1) hit$conf.high[1]^increase else NA_real_,
    p_value = if (nrow(hit) == 1) hit$p.value[1] else NA_real_
  )
}

#' Fit a Cox model, optionally weighted; returns NULL on failure
#' @noRd
.safe_coxph <- function(fmla, d, weights) {
  tryCatch({
    if (is.null(weights)) {
      survival::coxph(fmla, data = d)
    } else {
      d$.cox_weights <- d[[weights]]
      survival::coxph(fmla, data = d, weights = .cox_weights)
    }
  }, error = function(e) NULL)
}

#' @noRd
.fit_univariable <- function(data, time_var, event_var, covariate, cov_type,
                             ref_level, weights, continuous_increase = 1) {

  needed <- c(time_var, event_var, covariate, weights)
  d <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]

  if (nrow(d) == 0) return(tibble::tibble())

  if (cov_type == "categorical") {
    d[[covariate]] <- .factor_with_ref(d[[covariate]], ref_level)
    fac_levels <- levels(d[[covariate]])
    ref_level  <- fac_levels[1]
  }

  fmla <- stats::as.formula(
    paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", covariate)
  )
  fit <- .safe_coxph(fmla, d, weights)
  if (is.null(fit)) return(tibble::tibble())

  tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

  if (cov_type == "categorical") {
    .rows_from_categorical_term(tt, d, covariate, covariate, event_var, fac_levels, ref_level)
  } else {
    .row_from_continuous_term(tt, d, covariate, covariate, event_var,
                              increase = continuous_increase)
  }
}

#' @noRd
.fit_multivariable <- function(data, time_var, event_var, covariates, cov_types,
                               ref_levels, weights, continuous_increases = NULL) {

  needed <- c(time_var, event_var, covariates, weights)
  d <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]

  if (nrow(d) == 0) return(tibble::tibble())

  fac_level_map <- list()
  for (cov in covariates) {
    if (cov_types[[cov]] == "categorical") {
      d[[cov]] <- .factor_with_ref(d[[cov]], ref_levels[[cov]])
      fac_level_map[[cov]] <- levels(d[[cov]])
    }
  }

  fmla <- stats::as.formula(paste0(
    "survival::Surv(", time_var, ", ", event_var, ") ~ ", paste(covariates, collapse = " + ")
  ))
  fit <- .safe_coxph(fmla, d, weights)
  if (is.null(fit)) return(tibble::tibble())

  tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

  purrr::map_dfr(covariates, function(cov) {
    if (cov_types[[cov]] == "categorical") {
      .rows_from_categorical_term(tt, d, cov, cov, event_var,
                                  fac_level_map[[cov]], fac_level_map[[cov]][1])
    } else {
      increase <- if (!is.null(continuous_increases[[cov]])) continuous_increases[[cov]] else 1
      .row_from_continuous_term(tt, d, cov, cov, event_var, increase = increase)
    }
  })
}

#' Multiply body-row heights of a forestploter object
#' @noRd
.apply_row_height <- function(p, row_height = 1) {

  if (length(row_height) != 1 || !is.numeric(row_height) ||
      is.na(row_height) || row_height <= 0) {
    stop("`row_height` must be a single positive numeric value.", call. = FALSE)
  }

  core_fg <- which(p$layout$name == "core-fg")
  if (length(core_fg) == 0) return(p)

  body_rows <- unique(p$layout$t[core_fg])
  p$heights[body_rows] <- p$heights[body_rows] * row_height
  p
}

#' Long-format results for one model type across all endpoints
#' @noRd
.build_cox_results <- function(data, covariates, labels, time, event, endpoint_labels,
                               cov_types, ref_levels, weights, model_type,
                               continuous_increases = NULL) {

  label_lookup <- stats::setNames(labels, covariates)

  purrr::map_dfr(seq_along(time), function(ep) {
    time_var  <- time[ep]
    event_var <- event[ep]
    ep_label  <- endpoint_labels[ep]

    if (model_type == "univariable") {
      rows <- purrr::map_dfr(covariates, function(cov) {
        increase <- if (!is.null(continuous_increases[[cov]])) continuous_increases[[cov]] else 1
        .fit_univariable(data, time_var, event_var, cov, cov_types[[cov]],
                         ref_levels[[cov]], weights, continuous_increase = increase)
      })
    } else {
      rows <- .fit_multivariable(data, time_var, event_var, covariates, cov_types,
                                 ref_levels, weights, continuous_increases = continuous_increases)
    }

    if (nrow(rows) == 0) return(rows)

    rows |>
      dplyr::mutate(
        covariate_label = unname(label_lookup[covariate]),
        endpoint = time_var, endpoint_label = ep_label,
        model_type = model_type
      )
  })
}

#' Build the forestploter object for one model's results table
#' @noRd
.build_cox_forest_plot <- function(results, digits, endpoint_palette, xlim, ticks_at,
                                   length_plot, ci_pch, ci_lty, ci_lwd, ci_Theight, sizes,
                                   risk_labels, title, row_height, covariate_label_row = FALSE,
                                   col_spacing = c(0, 0), label_indent = 4) {

  if (length(row_height) != 1 || !is.numeric(row_height) || is.na(row_height) || row_height <= 0) {
    stop("`row_height` must be a single positive numeric value.", call. = FALSE)
  }

  if (is.null(endpoint_palette)) {
    endpoint_palette <- c("#E67E22", "#16A085", "#8E44AD", "#C0392B", "#2980B9")
  }

  # Preserve covariate and endpoint order as supplied by the user
  cov_order      <- unique(results$covariate_label)
  endpoint_order <- unique(results$endpoint_label)

  df <- results |>
    dplyr::mutate(covariate_label = factor(covariate_label, levels = cov_order)) |>
    dplyr::arrange(covariate_label, dplyr::desc(is_reference), level,
                   match(endpoint_label, endpoint_order)) |>
    dplyr::group_by(covariate_label, level, is_reference) |>
    dplyr::mutate(endpoint_order = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      is_first_of_block = dplyr::row_number() == 1 | covariate_label != dplyr::lag(covariate_label),
      covariate_display = ifelse(is_first_of_block, as.character(covariate_label), ""),
      level_display     = ifelse(endpoint_order == 1, level, ""),
      n_fmt             = ifelse(endpoint_order == 1 & !is.na(n), as.character(n), ""),
      block_id          = cumsum(is_first_of_block),
      hr_ci             = .format_hr_ci(hr, conf_low, conf_high, is_reference, digits),
      p_fmt             = ifelse(!is_reference & is.na(p_value), "NA", .format_pvalue(p_value)),
      event_fmt         = ifelse(is.na(n_event), "", as.character(n_event)),
      blank_col         = strrep(" ", length_plot)
    )

  # Optional covariate header rows: the covariate label gets its own (bold)
  # row and the levels are shown indented below it in the SAME column.
  if (covariate_label_row) {
    df$is_covariate_header <- FALSE
    df <- purrr::map_dfr(
      split(df, df$covariate_label),
      function(block) {
        header <- block[1, , drop = FALSE]
        header$covariate_display <- as.character(block$covariate_label[1])
        header$label_display     <- as.character(block$covariate_label[1])
        header$level_display     <- ""
        header$n_fmt             <- ""
        header$event_fmt         <- ""
        header$blank_col         <- strrep(" ", length_plot)
        header$hr_ci             <- ""
        header$p_fmt             <- ""
        header$hr                <- NA_real_
        header$conf_low          <- NA_real_
        header$conf_high         <- NA_real_
        header$is_covariate_header <- TRUE

        block$covariate_display  <- ""
        block$label_display      <- ifelse(
          block$level_display == "", "",
          paste0(strrep(" ", label_indent), block$level_display)
        )
        block$is_covariate_header <- FALSE

        dplyr::bind_rows(header, block)
      }
    )
  } else {
    df$is_covariate_header <- FALSE
  }

  endpoints      <- unique(df$endpoint_label)
  multi_endpoint <- length(endpoints) > 1
  endpoint_col   <- stats::setNames(rep(endpoint_palette, length.out = length(endpoints)), endpoints)

  # Table columns (optional spacer columns)
  if (covariate_label_row) {
    # One combined column: covariate label (bold) + indented levels
    plot_cols <- "label_display"
    col_names <- "Variable"
  } else {
    plot_cols <- "covariate_display"
    col_names <- "Variable"

    if (col_spacing[1] > 0) {
      df$spacer_var_subgroup <- strrep(" ", col_spacing[1])
      plot_cols <- c(plot_cols, "spacer_var_subgroup")
      col_names <- c(col_names, "")
    }

    plot_cols <- c(plot_cols, "level_display")
    col_names <- c(col_names, "Subgroup")
  }

  if (col_spacing[2] > 0) {
    df$spacer_subgroup_n <- strrep(" ", col_spacing[2])
    plot_cols <- c(plot_cols, "spacer_subgroup_n")
    col_names <- c(col_names, "  ")
  }

  plot_cols <- c(plot_cols, "n_fmt", "event_fmt", "blank_col", "hr_ci", "p_fmt")
  col_names <- c(col_names, "N", "Events", " ", "HR (95% CI)", "p-value")

  tbl <- as.data.frame(df[, plot_cols])
  names(tbl) <- col_names

  ci_col <- which(names(tbl) == " ")

  if (is.null(xlim)) {
    finite_bounds <- c(df$conf_low, df$conf_high, 1)
    finite_bounds <- finite_bounds[is.finite(finite_bounds) & finite_bounds > 0]
    x_range <- range(finite_bounds, na.rm = TRUE)
    xlim <- c(max(0.05, x_range[1] * 0.8), x_range[2] * 1.25)
  }
  if (is.null(ticks_at)) {
    ticks_at <- pretty(xlim, n = 5)
    ticks_at <- ticks_at[ticks_at > 0 & ticks_at >= xlim[1] & ticks_at <= xlim[2]]
    if (length(ticks_at) < 2) {
      ticks_at <- signif(exp(seq(log(xlim[1]), log(xlim[2]), length.out = 4)), 2)
    }
  }

  # Alternating shading; with header rows the whole block shares one background
  if (covariate_label_row) {
    covariate_block_id <- cumsum(df$is_covariate_header)
    shaded_rows <- which(covariate_block_id %% 2 == 1)
  } else {
    shaded_rows <- which(df$block_id %% 2 == 1)
  }

  ft_theme <- forestploter::forest_theme(
    base_size  = 11,
    ci_pch     = ci_pch,
    ci_col     = if (multi_endpoint) endpoint_col[[1]] else endpoint_palette[[1]],
    ci_lty     = ci_lty,
    ci_lwd     = ci_lwd,
    ci_Theight = ci_Theight,
    axis_extend = c(0.2, 0.05),
    vertline_lwd = 1, vertline_lty = "dashed", vertline_col = "grey20",
    refline_gp = grid::gpar(lwd = 1, lty = "dashed", col = "grey30"),
    footnote_gp = grid::gpar(col = "grey20", fontface = "italic", cex = 0.8),
    core = list(bg_params = list(fill = "white"))
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
    arrow_lab = risk_labels,
    title     = title,
    theme     = ft_theme
  )

  p <- .apply_row_height(p, row_height = row_height)

  # Bold covariate labels
  bold_rows <- if (covariate_label_row) which(df$is_covariate_header) else which(df$is_first_of_block)
  if (length(bold_rows) > 0) {
    p <- forestploter::edit_plot(p, row = bold_rows, col = 1,
                                 gp = grid::gpar(fontface = "bold"))
  }

  if (length(shaded_rows) > 0) {
    p <- forestploter::edit_plot(p, row = shaded_rows, which = "background",
                                 gp = grid::gpar(fill = "grey95", col = "grey95"))
  }

  if (multi_endpoint) {
    for (ep_label in endpoints) {
      rows_i <- which(df$endpoint_label == ep_label & !df$is_reference)
      col_i  <- endpoint_col[[ep_label]]
      if (length(rows_i) > 0) {
        p <- forestploter::edit_plot(p, row = rows_i, col = ci_col, which = "ci",
                                     gp = grid::gpar(col = col_i, fill = col_i))
      }
    }
  }

  sep_gp   <- grid::gpar(col = "grey30", lwd = 1)
  sep_cols <- c(0, ci_col - 1, ci_col, length(names(tbl)))

  # Vertical separators (body + header)
  for (col_i in sep_cols) {
    p <- forestploter::add_border(p, row = seq_len(nrow(tbl)), col = col_i,
                                  part = "body", where = "right", gp = sep_gp)
    p <- forestploter::add_border(p, row = 1, col = col_i,
                                  part = "header", where = "right", gp = sep_gp)
  }

  # Horizontal lines
  for (row_i in c(0, nrow(tbl))) {
    p <- forestploter::add_border(p, row = row_i, col = seq_len(ncol(tbl)),
                                  part = "body", where = "bottom", gp = sep_gp)
  }
  p <- forestploter::add_border(p, row = 0, col = seq_len(ncol(tbl)),
                                part = "header", where = "bottom", gp = sep_gp)

  # Manual legend for multiple endpoints
  if (multi_endpoint) {
    legend_cols <- unname(endpoint_col[endpoints])
    n_legend    <- length(endpoints)
    legend_grobs <- list()

    for (i in seq_along(endpoints)) {
      x_pos <- (i - 0.5) / n_legend

      legend_grobs[[length(legend_grobs) + 1]] <- grid::pointsGrob(
        x = grid::unit(x_pos - 0.035, "npc"),
        y = grid::unit(0.5, "npc"),
        pch = ci_pch,
        size = grid::unit(2.2, "mm"),
        gp = grid::gpar(col = legend_cols[i], fill = legend_cols[i])
      )

      legend_grobs[[length(legend_grobs) + 1]] <- grid::textGrob(
        label = endpoints[i],
        x = grid::unit(x_pos - 0.005, "npc"),
        y = grid::unit(0.5, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(fontsize = 9, col = "grey20")
      )
    }

    legend_grob <- grid::gTree(children = do.call(grid::gList, legend_grobs))

    p <- gtable::gtable_add_rows(p, heights = grid::unit(0.3, "in"), pos = nrow(p))
    p <- gtable::gtable_add_grob(p, legend_grob, t = nrow(p), l = 1, r = ncol(p),
                                 name = "manual_endpoint_legend")
  }

  p
}


# ==============================================================================
# Exported function
# ==============================================================================

#' Forest plots for univariable and/or multivariable Cox models
#'
#' For every covariate in `covariates` (a treatment/arm variable may be
#' included like any other covariate), estimates the hazard ratio for each
#' non-reference level and each requested endpoint, and renders the results as
#' a forest plot. A "Reference" row is shown for the reference level of every
#' categorical covariate. Non-reference levels for which no estimate could be
#' obtained (e.g. too few patients or events) are displayed as "NA".
#'
#' Two model types are supported (see `model_type`):
#' \describe{
#'   \item{univariable}{One Cox model per covariate:
#'     `Surv(time, event) ~ covariate`.}
#'   \item{multivariable}{One Cox model per endpoint containing all covariates:
#'     `Surv(time, event) ~ cov1 + cov2 + ... + covk`.}
#' }
#'
#' @param data A data.frame or tibble.
#' @param covariates Character vector of covariate names in `data`.
#' @param labels Character vector of display labels, same length/order as
#'   `covariates`.
#' @param time Character vector of one or more survival time variable names.
#' @param event Character vector of event indicator variable names, paired by
#'   position with `time`.
#' @param endpoint_labels Optional display labels for each (time, event) pair.
#'   Defaults to `time`.
#' @param continuous_vars Optional character vector of covariates to treat as
#'   continuous. If `NULL`, numeric variables with more than 5 distinct values
#'   are auto-detected as continuous.
#' @param continuous_increase Optional named list (or named vector) giving the
#'   unit increase for the HR of continuous covariates, e.g.
#'   `list(age = 10)` reports the HR per 10 units. Default is 1 unit.
#' @param ref_levels Optional named list giving the reference level of
#'   categorical covariates, e.g. `list(sex = "Female")`. Defaults to the
#'   first factor level.
#' @param exclude_levels Optional named list of levels to exclude per
#'   covariate, e.g. `list(ecog = "Unknown")`. Excluded values are treated as
#'   missing: only that covariate's univariable model loses these patients,
#'   whereas the multivariable model drops them entirely (complete cases).
#' @param model_type One of `"univariable"`, `"multivariable"` or `"both"`.
#'   Default `"univariable"`.
#' @param digits Number of decimals for HR/CI display. Default 2.
#' @param palette Colours per endpoint (only relevant with more than one
#'   endpoint). Must have the same length as `endpoint_labels`.
#' @param xlim Numeric length-2 vector, limits of the HR axis. Computed
#'   automatically if `NULL`.
#' @param ticks_at Numeric vector of axis tick positions. Computed
#'   automatically if `NULL`.
#' @param length_plot Width of the CI column in blank characters. Default 40.
#' @param weights Optional single column name in `data` with case weights
#'   passed to [survival::coxph()]. Rows with missing weights are dropped.
#' @param ci_pch,ci_lty,ci_lwd,ci_Theight,sizes Appearance of the point
#'   estimates and confidence intervals, passed to
#'   [forestploter::forest_theme()] / [forestploter::forest()].
#' @param risk_labels Optional length-2 character vector for the axis arrow
#'   annotations, e.g. `c("Lower risk", "Higher risk")`. `NULL` omits them.
#' @param title Plot title. Defaults to "Univariable Cox regression" /
#'   "Multivariable Cox regression".
#' @param row_height Multiplier for the body row height (`1` = default,
#'   `1.3` = 30 percent taller).
#' @param covariate_label_row Logical. If `TRUE`, a separate bold row with the
#'   covariate label is inserted before each covariate block, and the levels
#'   are shown indented below it in the same column (no separate "Subgroup"
#'   column). Default `FALSE`.
#' @param col_spacing Numeric length-2 vector: width (in blank characters) of
#'   spacer columns between the first and second and between the second and
#'   third table column. `0` adds no spacer. The first element is ignored when
#'   `covariate_label_row = TRUE`, since only one label column exists.
#' @param label_indent Number of spaces used to indent the levels when
#'   `covariate_label_row = TRUE`. Default 4.
#'
#' @return A list with elements:
#' \describe{
#'   \item{results_table}{Long-format tibble, one row per covariate level x
#'     endpoint x model type, with HR, 95% CI, p-value, N and events.}
#'   \item{plots}{Named list with up to two forestploter plot objects,
#'     `$univariable` and/or `$multivariable`.}
#' }
#'
#' @examples
#' \dontrun{
#' lung <- survival::lung
#'
#' lung2 <- lung %>%
#'  mutate(
#'     event_os = status - 1,
#'     time_os  = time,
#'     sex      = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
#'     ph.ecog  = factor(ph.ecog),
#'     treatment = sample(c("Arm A", "Arm B"), n(), replace = TRUE),
#'     time_pfs  = pmin(time_os, time_os * runif(n(), 0.4, 0.9)),
#'     event_pfs = ifelse(time_pfs < time_os, 1, event_os)
#'   ) %>%
#'   filter(!is.na(ph.ecog), !is.na(sex), !is.na(age))
#'
#' ps.formula <- as.formula(paste("treatment", "~", paste(c("sex", "age", "ph.ecog"), collapse = "+")))
#'
#' W.out <- WeightIt::weightit(ps.formula, data = lung2, method = "ps", estimand = "ATE")
#'
#' lung2$weights_ATE <- W.out$weights
#'
#' # Usage of the function
#' res2 <- cox_forest(
#'   data = lung2,
#'   # Specify here all covariates that you want to include
#'   covariates = c("treatment", "sex", "ph.ecog", "age"),
#'   # If you have a continuous variable please specify it here (didn't implement
#'   # the automatic search (next step))
#'   continuous_vars = "age",
#'   # Optional: If you want a x unit increase for a continuous increase specify it here
#'   continuous_increase = list(age = 10),
#'   # Specify here the names you want to have to be displayed
#'   # for each covariate
#'   labels = c("Treatment", "Sex", "ECOG PS", "Age"),
#'   # Specify the time and event column names of the endpoint of interest
#'   time = c( "time_pfs", "time_os"),
#'   event = c("event_pfs", "event_os"),
#'   # Names of each endpoint
#'   endpoint_labels = c("Progression-free Survival", "Overall Survival"),
#'   # Optional: if you want to change the reference level of a covariate
#'   ref_levels = list(sex = "Male"),
#'   # Optional: can exlude labels of each covariate just for the displaying
#'   # if for example only 4 patient are in one group the Cox model will be
#'   # obtained with those patient and group but not displayed in the forestplot
#'   # because it has no meaning
#'   exclude_levels = list(ph.ecog = "3"),
#'   # Specify which model type you want ("univariable", "multivariable", "both")
#'   model_type = "both",
#'   # Colour palette
#'   palette = c("#0F084B", "#A6CFD5"),
#'   # Limits of the x-axis for the forestplot
#'   xlim = c(0.25, 4),
#'   # Ticks at the x-axis
#'   ticks_at = c(0.5, 1, 2),
#'   # Width of the forestplot
#'   length_plot = 60,
#'   # Optional: only if using PS
#'   #weights = "weights_ATE", # optional if PS used
#'   # Argument for the displaying of the CI (check forestmodel::forest for each meaning)
#'   ci_pch = 15,
#'   ci_lty = 1,
#'   ci_lwd = 3,
#'   ci_Theight = 0.4,
#'   sizes = 1,
#'   # Optional: Specify row height
#'   row_height = 1.4,
#'   # Optional: TRUE creates a seperate row for the variable name
#'   covariate_label_row = FALSE,
#'   # Optional: increases the space between the first 3 columns
#'   col_spacing = c(1, 1)
#' )
#'
#' # Ploting univariable model
#' res2$plots$univariable
#'
#' # Ploting multivariable model
#' res2$plots$multivariable
#'
#' # Table with all uni and multivariable results
#' res2$results_table
#' }
#'
#' @importFrom stats complete.cases as.formula relevel setNames
#' @export
cox_forest <- function(data, covariates, labels, time, event, endpoint_labels = NULL,
                       continuous_vars = NULL, continuous_increase = NULL,
                       ref_levels = NULL, exclude_levels = NULL,
                       model_type = c("univariable", "multivariable", "both"),
                       digits = 2, palette = NULL, xlim = NULL, ticks_at = NULL,
                       length_plot = 40, weights = NULL, ci_pch = 15, ci_lty = 1,
                       ci_lwd = 3, ci_Theight = 0.4, sizes = 0.5, risk_labels = NULL,
                       title = NULL, row_height = 1, covariate_label_row = FALSE,
                       col_spacing = c(0, 0), label_indent = 4) {

  model_type <- match.arg(model_type, several.ok = TRUE)
  if ("both" %in% model_type) model_type <- c("univariable", "multivariable")
  model_type <- unique(model_type)

  .validate_cox_inputs(data, covariates, labels, time, event, weights)

  if (is.null(endpoint_labels)) endpoint_labels <- time
  if (length(endpoint_labels) != length(time)) {
    stop("`endpoint_labels` must be the same length as `time`.", call. = FALSE)
  }
  if (!is.null(palette) && length(endpoint_labels) != length(palette)) {
    stop("`palette` must be the same length as `endpoint_labels`.", call. = FALSE)
  }

  continuous_increase <- as.list(continuous_increase)
  if (is.null(ref_levels)) ref_levels <- list()

  data <- .apply_exclude_levels(data, exclude_levels)

  cov_types <- stats::setNames(
    lapply(covariates, function(cov) .detect_covariate_type(data, cov, continuous_vars)),
    covariates
  )

  results_table <- purrr::map_dfr(model_type, function(mt) {
    .build_cox_results(data, covariates, labels, time, event, endpoint_labels,
                       cov_types, ref_levels, weights, mt,
                       continuous_increases = continuous_increase)
  })

  plots <- list()
  for (mt in model_type) {
    res_mt <- dplyr::filter(results_table, .data$model_type == mt)
    if (nrow(res_mt) == 0) next
    plot_title <- if (is.null(title)) {
      if (mt == "univariable") "Univariable Cox regression" else "Multivariable Cox regression"
    } else title
    plots[[mt]] <- .build_cox_forest_plot(
      results = res_mt, digits = digits, endpoint_palette = palette, xlim = xlim,
      ticks_at = ticks_at, length_plot = length_plot, ci_pch = ci_pch, ci_lty = ci_lty,
      ci_lwd = ci_lwd, ci_Theight = ci_Theight, sizes = sizes, risk_labels = risk_labels,
      title = plot_title, row_height = row_height, covariate_label_row = covariate_label_row,
      col_spacing = col_spacing, label_indent = label_indent
    )
  }

  list(
    results_table = results_table,
    plots         = plots
  )
}
