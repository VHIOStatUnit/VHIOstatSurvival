# ==============================================================================
# R/model_forest.R
#
# Publication-quality forest plots for ANY regression model in R that
# broom::tidy() supports (coxph, glm [logistic, Poisson, ...], lm, clogit,
# ...), driven by a model formula.
#
# One row per covariate level (with a "Reference" row for the reference level
# of every categorical covariate), one row per continuous covariate.
# Optionally several formulas (e.g. one per endpoint) can be shown in one plot.
#
# Author: Julian Silan
# Requires: survival, broom, dplyr, purrr, forestploter, tibble, grid, gtable
# ==============================================================================


# ==============================================================================
# Helper functions
# ==============================================================================

#' Coerce a variable to a factor, moving `ref_level` to the front if supplied
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

#' Format an "estimate (CI)" text string, e.g. "1.34 (1.05-1.72)"
#' @noRd
.format_est_ci <- function(est, lower, upper, is_reference, digits = 2) {
  # ASCII hyphen on purpose (base pdf() device mangles en dashes).
  fmt <- paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)")
  dplyr::case_when(
    is_reference                             ~ "Reference",
    is.na(est) | is.na(lower) | is.na(upper) ~ "NA",
    TRUE                                     ~ sprintf(fmt, est, lower, upper)
  )
}


#' Format a p-value, showing "<0.001" style thresholds
#' @noRd
.format_pvalue <- function(p, digits = 3) {
  threshold <- 10^(-digits)
  dplyr::case_when(
    is.na(p)      ~ "",
    p < threshold ~ paste0("<", format(threshold, scientific = FALSE)),
    TRUE          ~ formatC(p, digits = digits, format = "f")
  )
}

#' Row-height multiplier for the body rows of a forestploter object
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


#' Exclude specific levels of specific covariates from analysis, by setting
#' matching values to NA. This lets each covariate's excluded level(s) drop
#' out via the existing complete-case filtering: for univariable models only
#' that covariate's own model is affected, while for the multivariable model
#' the affected patients are dropped from the whole model (since it requires
#' complete cases across all covariates simultaneously).
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

#' Set reference levels in the data before fitting
#' @noRd
.apply_ref_levels <- function(data, ref_levels) {
  if (is.null(ref_levels) || length(ref_levels) == 0) return(data)
  for (cov in names(ref_levels)) {
    if (!cov %in% names(data)) {
      stop(sprintf("`ref_levels`: variable '%s' not found in `data`.", cov), call. = FALSE)
    }
    data[[cov]] <- .factor_with_ref(data[[cov]], ref_levels[[cov]])
  }
  data
}


#' Model frame without relying on re-evaluating fit$call
#' (coxph/lm/glm store the frame in fit$model when fitted with model = TRUE;
#' otherwise fall back to model.frame() with the data supplied explicitly)
#' @noRd
.get_model_frame <- function(fit, data = NULL) {
  mf <- fit$model
  if (is.data.frame(mf)) return(mf)
  tryCatch(
    stats::model.frame(fit),
    error = function(e) {
      if (is.null(data)) stop(e)
      stats::model.frame(fit, data = data)
    }
  )
}


#' Fit the model for one formula
#' @noRd
.fit_model <- function(formula, data, fit_fun, weights = NULL, ...) {
  if (is.null(fit_fun)) stop("`fit_fun`: a fitting model must be specified.", call. = FALSE)
  if (is.character(fit_fun)) fit_fun <- match.fun(fit_fun)


  if (!is.null(weights) && !weights %in% names(data)) {
    stop(sprintf("`weights`: column '%s' not found in `data`.", weights), call. = FALSE)
  }
  # keep the model frame whenever the fitter supports it (coxph, lm, glm, ...);
  # names(match.call(...)$...) avoids forcing `...` (e.g. id = id, cluster = id)
  keep <- "model" %in% names(formals(fit_fun)) &&
    !"model" %in% names(match.call(expand.dots = FALSE)$...)


  # Build the call with `...` left unevaluated; `weights` goes in as the
  # column *symbol* so the fitter resolves it inside `data`.
  cl <- quote(fit_fun(formula = formula, data = data, ...))
  if (keep) cl$model <- TRUE
  if (!is.null(weights)) cl$weights <- as.name(weights)
  eval(cl, envir = environment())

}

#' Defaults (label, exponentiate, reference value) depending on the model class
#' @noRd
.model_defaults <- function(fit) {
  if (inherits(fit, "clogit")) return(list(label = "OR", exponentiate = TRUE))
  if (inherits(fit, "coxph"))  return(list(label = "HR", exponentiate = TRUE))
  if (inherits(fit, "glm")) {
    link <- fit$family$link
    if (identical(link, "logit")) return(list(label = "OR", exponentiate = TRUE))
    if (identical(link, "log"))   return(list(label = "RR", exponentiate = TRUE))
  }
  list(label = "Estimate", exponentiate = FALSE)
}

#' Per-observation event indicator (Surv status or binary glm outcome);
#' NULL if events are not defined for this model
#' @noRd
.model_events <- function(fit, mf) {
  y <- tryCatch(stats::model.response(mf), error = function(e) NULL)
  if (inherits(y, "Surv")) return(as.numeric(y[, ncol(y)]))
  if (inherits(fit, "glm") && identical(fit$family$family, "binomial")) {
    yy <- fit$y
    if (!is.null(yy) && all(yy %in% c(0, 1))) return(as.numeric(yy))
  }
  NULL
}

#' Extract one row per covariate level, including reference categories.
#' Supports single-coefficient numeric transformations (e.g., log(x), I(x^2))
#' and categorical terms with treatment/dummy coding (e.g., factor(x)).
#' Interactions and multi-column transformations are not plotted.
#' @noRd
.extract_model_rows <- function(fit, exponentiate, labels, ci_level, data = NULL) {
  mf <- .get_model_frame(fit, data)
  tt <- broom::tidy(fit, conf.int = TRUE, conf.level = ci_level,
                    exponentiate = exponentiate)
  events <- .model_events(fit, mf)
  xlv <- if (is.list(fit)) fit$xlevels else NULL

  # Reconstruct a standard model matrix from the fitted model frame.
  # Unlike model.matrix.coxph(), this preserves the integer "assign"
  # attribute required to associate coefficients with formula terms.
  terms_obj <- stats::delete.response(stats::terms(fit))
  mm <- tryCatch(
    stats::model.matrix(
      terms_obj,
      data = mf,
      contrasts.arg = fit$contrasts,
      xlev = fit$xlevels
    ),
    error = function(e) {
      stop(
        "Cannot reconstruct the model matrix: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  assign_idx <- attr(mm, "assign")
  if (is.null(assign_idx)) {
    stop("Reconstructed model matrix has no term assignments.", call. = FALSE)
  }
  term_labels <- attr(terms_obj, "term.labels")
  coef_map <- stats::setNames(
    lapply(seq_along(term_labels), function(i) colnames(mm)[assign_idx == i]),
    term_labels
  )

  # Match model-matrix rows to the model frame, including after NA removal.
  if (!identical(rownames(mm), rownames(mf))) {
    idx <- match(rownames(mf), rownames(mm))
    if (anyNA(idx)) {
      stop("Model matrix rows cannot be aligned with the model frame.", call. = FALSE)
    }
    mm <- mm[idx, , drop = FALSE]
  }

  getv <- function(h, nm) {
    if (nrow(h) == 1L && nm %in% names(h)) h[[nm]][1] else NA_real_
  }
  ev_sum <- function(idx) {
    if (is.null(events)) NA_real_ else sum(events[idx], na.rm = TRUE)
  }

  out <- list()
  for (tl in term_labels) {
    if (grepl(":", tl, fixed = TRUE)) {
      warning(sprintf("Interaction term '%s' is not shown.", tl), call. = FALSE)
      next
    }
    if (grepl("^(strata|cluster|frailty|tt)\\(", tl)) next
    if (!tl %in% names(mf)) {
      warning(sprintf("Term '%s' was not found in the model frame; skipping.", tl),
              call. = FALSE)
      next
    }

    x <- mf[[tl]]
    if (is.matrix(x) || is.ordered(x)) {
      warning(sprintf("Matrix-valued or ordered term '%s' is not shown.", tl),
              call. = FALSE)
      next
    }
    lab <- if (!is.null(labels) && tl %in% names(labels)) {
      unname(labels[[tl]])
    } else tl
    coef_names <- coef_map[[tl]]
    if (!length(coef_names)) next

    if (is.factor(x) || is.character(x) || is.logical(x)) {
      lv <- xlv[[tl]]
      if (is.null(lv)) {
        lv <- if (is.factor(x)) levels(x) else if (is.logical(x)) {
          c("FALSE", "TRUE")
        } else sort(unique(as.character(x)))
      }
      if (length(lv) < 2L) next
      xc <- as.character(x)
      non_ref <- lv[-1L]

      # A reference row is meaningful only for ordinary treatment coding.
      # Determine the correct coefficient for each non-reference level from
      # the actual indicator columns, not from assumed coefficient names.
      if (length(coef_names) != length(non_ref)) {
        warning(sprintf("Term '%s' does not use one coefficient per non-reference level; skipping.", tl),
                call. = FALSE)
        next
      }
      term_mm <- mm[, coef_names, drop = FALSE]
      matched <- character(length(non_ref))
      available <- seq_along(coef_names)
      valid <- TRUE
      for (i in seq_along(non_ref)) {
        expected <- as.numeric(xc == non_ref[i])
        candidates <- available[vapply(available, function(j) {
          isTRUE(all.equal(as.numeric(term_mm[, j]), expected,
                           tolerance = 1e-8, check.attributes = FALSE))
        }, logical(1))]
        if (length(candidates) != 1L) {
          valid <- FALSE
          break
        }
        matched[i] <- coef_names[candidates]
        available <- setdiff(available, candidates)
      }
      if (!valid) {
        warning(sprintf("Term '%s' is not treatment/dummy coded; skipping.", tl),
                call. = FALSE)
        next
      }

      ref <- lv[1L]
      ref_idx <- xc == ref
      ref_row <- tibble::tibble(
        covariate = tl, covariate_label = lab, level = ref, level_idx = 1L,
        is_reference = TRUE, n = sum(ref_idx), n_event = ev_sum(ref_idx),
        estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
        p_value = NA_real_
      )
      lev_rows <- purrr::map_dfr(seq_along(non_ref), function(i) {
        hit <- tt[!is.na(tt$term) & tt$term == matched[i], , drop = FALSE]
        idx <- xc == non_ref[i]
        tibble::tibble(
          covariate = tl, covariate_label = lab,
          level = non_ref[i], level_idx = i + 1L,
          is_reference = FALSE, n = sum(idx), n_event = ev_sum(idx),
          estimate = getv(hit, "estimate"), conf_low = getv(hit, "conf.low"),
          conf_high = getv(hit, "conf.high"), p_value = getv(hit, "p.value")
        )
      })
      out[[length(out) + 1L]] <- dplyr::bind_rows(ref_row, lev_rows)

    } else if (is.numeric(x)) {
      if (length(coef_names) != 1L) {
        warning(sprintf("Term '%s' has multiple coefficients; skipping.", tl),
                call. = FALSE)
        next
      }
      hit <- tt[!is.na(tt$term) & tt$term == coef_names[1L], , drop = FALSE]
      if (nrow(hit) != 1L) {
        warning(sprintf("No unique coefficient found for '%s'; skipping.", tl),
                call. = FALSE)
        next
      }
      all_idx <- rep(TRUE, length(x))
      out[[length(out) + 1L]] <- tibble::tibble(
        covariate = tl, covariate_label = lab,
        level = "per 1-unit increase", level_idx = 1L,
        is_reference = FALSE, n = length(x), n_event = ev_sum(all_idx),
        estimate = getv(hit, "estimate"),
        conf_low = getv(hit, "conf.low"),
        conf_high = getv(hit, "conf.high"),
        p_value = getv(hit, "p.value")
      )
    }
  }
  dplyr::bind_rows(out)
}

# ------------------------------------------------------------------------------
# Plot builder
# ------------------------------------------------------------------------------

#' @noRd
.build_model.forest_plot <- function(results, digits, endpoint_palette, xlim, ticks_at,
                               length_plot, ci_pch, ci_lty, ci_lwd, ci_Theight, sizes,
                               risk_labels, title, row_height, covariate_label_row = FALSE,
                               col_spacing = c(0, 0), exponentiate = TRUE,
                               effect_label = "Estimate",
                               label_indent = 4, ci_level = 0.95) {

  if (length(row_height) != 1 || !is.numeric(row_height) || is.na(row_height) || row_height <= 0) {
    stop("`row_height` must be a single positive numeric value.", call. = FALSE)
  }
  if (is.null(endpoint_palette)) {
    endpoint_palette <- c("#E67E22", "#16A085", "#8E44AD", "#C0392B", "#2980B9")
  }
  ref_value <- if (exponentiate) 1 else 0

  cov_order      <- unique(results$covariate_label)
  endpoint_order <- unique(results$endpoint_label)

  df <- results |>
    dplyr::mutate(covariate_label = factor(covariate_label, levels = cov_order)) |>
    dplyr::arrange(covariate_label, level_idx, match(endpoint_label, endpoint_order)) |>
    dplyr::group_by(covariate_label, level, is_reference) |>
    dplyr::mutate(endpoint_order = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      is_first_of_block = dplyr::row_number() == 1 | covariate_label != dplyr::lag(covariate_label),
      covariate_display = ifelse(is_first_of_block, as.character(covariate_label), ""),
      level_display     = ifelse(endpoint_order == 1, level, ""),
      n_fmt             = ifelse(endpoint_order == 1 & !is.na(n), as.character(n), ""),
      block_id          = cumsum(is_first_of_block),
      est_ci            = .format_est_ci(estimate, conf_low, conf_high, is_reference, digits),
      p_fmt             = ifelse(!is_reference & is.na(p_value), "NA", .format_pvalue(p_value)),
      event_fmt         = ifelse(is.na(n_event), "", as.character(n_event)),
      blank_col         = strrep(" ", length_plot)
    )

  show_events <- any(!is.na(df$n_event))

  # Optional covariate label rows
  if (covariate_label_row) {
    df$is_covariate_header <- FALSE
    df <- purrr::map_dfr(
      split(df, df$covariate_label),
      function(block) {
        header <- block[1, , drop = FALSE]
        header$covariate_display <- as.character(block$covariate_label[1])
        header$label_display     <- as.character(block$covariate_label[1])
        header$level_display <- ""
        header$n_fmt <- ""
        header$event_fmt <- ""
        header$blank_col <- strrep(" ", length_plot)
        header$est_ci <- ""
        header$p_fmt <- ""
        header$estimate <- NA_real_
        header$conf_low <- NA_real_
        header$conf_high <- NA_real_
        header$is_reference <- FALSE
        header$is_covariate_header <- TRUE
        block$covariate_display <- ""
        block$label_display      <- ifelse(block$level_display == "", "", paste0(strrep(" ", label_indent), block$level_display))
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
  plot_cols <- c(plot_cols, "n_fmt")
  col_names <- c(col_names, "N")
  if (show_events) {
    plot_cols <- c(plot_cols, "event_fmt")
    col_names <- c(col_names, "Events")
  }
  plot_cols <- c(plot_cols, "blank_col", "est_ci", "p_fmt")
  col_names <- c(col_names, " ", paste0(effect_label, " (", format(100 * ci_level, trim = TRUE, scientific = FALSE), "% CI)"), "p-value")

  tbl <- as.data.frame(df[, plot_cols])
  names(tbl) <- col_names
  ci_col <- which(names(tbl) == " ")

  # Axis limits / ticks
  if (is.null(xlim)) {
    b <- c(df$conf_low, df$conf_high, ref_value)
    b <- b[is.finite(b)]
    if (exponentiate) {
      b <- b[b > 0]
      r <- range(b)
      xlim <- c(max(0.05, r[1] * 0.8), r[2] * 1.25)
    } else {
      r <- range(b)
      pad <- max(diff(r) * 0.1, 1e-6)
      xlim <- c(r[1] - pad, r[2] + pad)
    }
  }
  if (is.null(ticks_at)) {
    ticks_at <- pretty(xlim, n = 5)
    ticks_at <- ticks_at[ticks_at >= xlim[1] & ticks_at <= xlim[2]]
    if (exponentiate) {
      ticks_at <- ticks_at[ticks_at > 0]
      if (length(ticks_at) < 2) {
        ticks_at <- signif(exp(seq(log(xlim[1]), log(xlim[2]), length.out = 4)), 2)
      }
    }
  }

  # Alternating shading (whole covariate block incl. header row)
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
    est       = df$estimate,
    lower     = df$conf_low,
    upper     = df$conf_high,
    sizes     = sizes,
    ci_column = ci_col,
    ref_line  = ref_value,
    x_trans   = if (exponentiate) "log" else "none",
    xlim      = xlim,
    ticks_at  = ticks_at,
    arrow_lab = risk_labels,
    title     = title,
    theme     = ft_theme
  )

  p <- .apply_row_height(p, row_height = row_height)

  # Bold covariate labels
  if (covariate_label_row) {
    hdr_rows <- which(df$is_covariate_header)
    if (length(hdr_rows) > 0) {
      p <- forestploter::edit_plot(p, row = hdr_rows, col = 1, which = "text",
                                   gp = grid::gpar(fontface = "bold"))
    }
  } else {
    lab_rows <- which(df$is_first_of_block)
    if (length(lab_rows) > 0) {
      p <- forestploter::edit_plot(p, row = lab_rows, col = 1,
                                   gp = grid::gpar(fontface = "bold"))
    }
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

  # Vertical separators (body + header)
  sep_cols <- c(0, ci_col - 1, ci_col, length(names(tbl)))
  for (col_i in sep_cols) {
    p <- forestploter::add_border(p, row = 1:nrow(tbl), col = col_i, part = "body",
                                  where = "right", gp = grid::gpar(col = "grey30", lwd = 1))
    p <- forestploter::add_border(p, row = 1, col = col_i, part = "header",
                                  where = "right", gp = grid::gpar(col = "grey30", lwd = 1))
  }

  # Horizontal lines
  for (row_i in c(0, nrow(tbl))) {
    p <- forestploter::add_border(p, row = row_i, col = 1:ncol(tbl), part = "body",
                                  where = "bottom", gp = grid::gpar(col = "grey30", lwd = 1))
  }
  p <- forestploter::add_border(p, row = 0, col = 1:ncol(tbl), part = "header",
                                where = "bottom", gp = grid::gpar(col = "grey30", lwd = 1))

  # Manual legend for multiple models/endpoints
  if (multi_endpoint) {
    legend_cols <- unname(endpoint_col[endpoints])
    n_legend <- length(endpoints)
    legend_grobs <- list()
    for (i in seq_along(endpoints)) {
      x_pos <- (i - 0.5) / n_legend
      legend_grobs[[length(legend_grobs) + 1]] <- grid::pointsGrob(
        x = grid::unit(x_pos - 0.035, "npc"), y = grid::unit(0.5, "npc"),
        pch = ci_pch, size = grid::unit(2.2, "mm"),
        gp = grid::gpar(col = legend_cols[i], fill = legend_cols[i])
      )
      legend_grobs[[length(legend_grobs) + 1]] <- grid::textGrob(
        label = endpoints[i],
        x = grid::unit(x_pos - 0.005, "npc"), y = grid::unit(0.5, "npc"),
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
# Main function
# ==============================================================================

#' Forest plot for any regression model, specified by a formula
#'
#' Fits the model given by \code{formula} (multivariable, i.e. all terms
#' together) and renders the coefficients as a forest plot. Works with any
#' model that has \code{model.frame()}, \code{terms()} and a
#' \code{broom::tidy()} method (coxph, glm, lm, clogit, ...).
#'
#' @param formula A model formula, e.g. \code{Surv(os_time, os_event) ~ arm + age + sex},
#'   or a (preferably named) list of formulas to show several models/endpoints
#'   in one colour-coded plot; the list names are used as legend labels.
#' @param data A data.frame or tibble.
#' @param fit_fun Fitting function or its name (e.g. \code{survival::coxph},
#'   \code{"glm"}, or \code{stats::lm}). Must be specified.
#' @param labels Optional named character vector mapping term labels (including
#'   expressions such as \code{"log(age)"} and \code{"I(age^2)"}) to
#'   display labels, e.g. \code{c(age = "Age", sex = "Sex")}. Unnamed terms
#'   keep their variable name.
#' @param endpoint_labels Optional legend labels when \code{formula} is a list
#'   (overrides the list names).
#' Numeric transformations such as \code{log(x)} and \code{I(x^2)}
#' are displayed per one-unit increase on the transformed scale.
#' Categorical expressions such as \code{factor(x)} are supported with
#' standard treatment/dummy contrasts. Interactions, splines and other
#' multi-coefficient numeric terms are not plotted.
#'
#' @param ref_levels Optional named list of reference levels, e.g.
#'   \code{list(sex = "Female")} (applied to \code{data} before fitting).
#' @param exclude_levels Optional named list of levels to set to NA before
#'   fitting, e.g. \code{list(ecog = "Unknown")}.
#' @param exponentiate \code{NULL} (auto: TRUE for coxph and log/logit-link
#'   glm), or TRUE/FALSE.
#' @param effect_label \code{NULL} (auto: HR / OR / RR / Estimate) or a string.
#' @param ci_level Confidence level strictly between 0 and 1. Default 0.95.
#'   The confidence interval column heading reflects this value.
#' @param digits,palette,xlim,ticks_at,length_plot,ci_pch,ci_lty,ci_lwd,
#'   ci_Theight,sizes,risk_labels,title,row_height,covariate_label_row,
#'   col_spacing Plot appearance arguments, as in \code{cox_forest()}.
#' @param ... Passed on to \code{fit_fun}, e.g. \code{family = binomial},
#'   \code{weights = w}, \code{ties = "efron"}.
#' @param covariate_label_row Logical. If `TRUE`, a separate bold row with the
#'   covariate label is inserted before each covariate block. Default `FALSE`.
#'   covariate label is inserted before each covariate block, and the levels
#'   are shown indented below it in the same column (no separate "Subgroup"
#'   column). Default `FALSE`.
#' @param col_spacing Numeric length-2 vector: ... (the first element is
#' ignored when `covariate_label_row = TRUE`, since only one label column exists)
#' @param label_indent Number of spaces used to indent the levels when
#' `covariate_label_row = TRUE`. Default 4.
#'
#' @return A list with \code{results_table} (long tibble: covariate level x
#'   model, with estimate, CI, p-value, N, events), \code{models} (the fitted
#'   model objects) and \code{plot} (forestploter object).
#'
#' @examples
#' # Transformed numeric predictors in a Cox regression
#'
#' lung <- survival::lung
#'
#' lung2 <- lung %>%
#' mutate(
#'  event_os = status - 1,
#'  time_os  = time,
#'  sex      = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
#'  ph.ecog  = factor(ph.ecog),
#'  treatment = sample(c("Arm A", "Arm B"), n(), replace = TRUE),
#'  time_pfs  = pmin(time_os, time_os * runif(n(), 0.4, 0.9)),
#'  event_pfs = ifelse(time_pfs < time_os, 1, event_os)
#'  ) %>%
#'  filter(!is.na(ph.ecog), !is.na(sex), !is.na(age))
#'
#' res <- model_forest(
#'   survival::Surv(time_os, event_os) ~ log(age) + I(age^2) + factor(sex),
#'   data = lung2,
#'   fit_fun = survival::coxph,
#'   labels = c("log(age)" = "Log age", "I(age^2)" = "Age squared",
#'              "factor(sex)" = "Sex"),
#'   ci_level = 0.90
#' )
#' plot(res$plot)
#'
#' # Logistic regression with a categorical expression
#' res_logit <- model_forest(
#'   vs ~ mpg + factor(am), data = mtcars,
#'   fit_fun = stats::glm, family = stats::binomial()
#' )
#' plot(res_logit$plot)
#'
#' @export
model_forest <- function(formula, data, fit_fun = NULL, labels = NULL,
                         endpoint_labels = NULL,
                         ref_levels = NULL, exclude_levels = NULL,
                         exponentiate = NULL, effect_label = NULL, ci_level = 0.95,
                         weights = NULL, digits = 2, palette = NULL, xlim = NULL,
                         ticks_at = NULL, length_plot = 40, ci_pch = 15, ci_lty = 1,
                         ci_lwd = 3, ci_Theight = 0.4, sizes = 0.5, risk_labels = NULL,
                         title = NULL, row_height = 1, covariate_label_row = FALSE,
                         col_spacing = c(0, 0), label_indent = 4, ...) {

  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  formulas <- if (inherits(formula, "formula")) list(formula) else formula
  if (!is.list(formulas) || !all(vapply(formulas, inherits, logical(1), what = "formula"))) {
    stop("`formula` must be a formula or a list of formulas.", call. = FALSE)
  }

  ep_names <- names(formulas)
  if (is.null(ep_names)) {
    ep_names <- if (length(formulas) == 1) "Model" else paste("Model", seq_along(formulas))
  }

  if (!is.null(endpoint_labels)) {
    if (length(endpoint_labels) != length(formulas)) {
      stop("`endpoint_labels` must have one entry per formula.", call. = FALSE)
    }
    ep_names <- endpoint_labels
  }
  if (!is.null(palette) && length(palette) != length(formulas)) {
    stop("`palette` must have one colour per formula.", call. = FALSE)
  }

  data <- .apply_exclude_levels(data, exclude_levels)
  data <- .apply_ref_levels(data, ref_levels)

  fits <- lapply(formulas, function(f) .fit_model(f, data, fit_fun, weights, ...))
  names(fits) <- ep_names

  defs <- .model_defaults(fits[[1]])
  if (is.null(exponentiate)) exponentiate <- defs$exponentiate
  if (is.null(effect_label)) effect_label <- defs$label

  results_table <- purrr::map_dfr(seq_along(fits), function(i) {
    rows <- .extract_model_rows(fits[[i]], exponentiate, labels, ci_level, data = data)
    if (nrow(rows) == 0) return(rows)
    dplyr::mutate(rows, endpoint_label = ep_names[i])
  })
  if (nrow(results_table) == 0) {
    stop("No plottable terms found in the model(s).", call. = FALSE)
  }

  plot_title <- if (is.null(title)) {
    if (length(formulas) == 1) paste(deparse(formulas[[1]]), collapse = " ") else NULL
  } else title

  p <- .build_model.forest_plot(
    results = results_table, digits = digits, endpoint_palette = palette, xlim = xlim,
    ticks_at = ticks_at, length_plot = length_plot, ci_pch = ci_pch, ci_lty = ci_lty,
    ci_lwd = ci_lwd, ci_Theight = ci_Theight, sizes = sizes, risk_labels = risk_labels,
    title = plot_title, row_height = row_height, covariate_label_row = covariate_label_row,
    col_spacing = col_spacing, exponentiate = exponentiate, effect_label = effect_label,
    label_indent = label_indent, ci_level = ci_level
  )

  list(results_table = results_table, models = fits, plot = p)
}
