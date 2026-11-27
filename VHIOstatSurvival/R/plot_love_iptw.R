#' Love plot of covariate balance before and after IPTW
#'
#' @description
#' Draws a love plot showing the absolute standardized mean difference
#' (|SMD|) of each covariate before (unweighted) and after (weighted)
#' inverse probability of treatment weighting (IPTW). The plot background is
#' divided into four shaded zones that make it easy to judge balance at a
#' glance:
#'
#' * **Balanced**: |SMD| < 0.1
#' * **Minor** imbalance: 0.1 to < 0.25
#' * **Moderate** imbalance: 0.25 to < 0.5
#' * **Severe** imbalance: >= 0.5
#'
#' Each covariate gets one point per sample (unweighted and weighted), and a
#' grey line joins the two so the change after weighting is easy to see.
#'
#' For treatments with more than two levels, the maximum |SMD| across all
#' pairs of treatment groups is shown for each covariate.
#'
#' @details
#' Balance statistics are computed with [cobalt::bal.tab()] using
#' `un = TRUE`, so both unadjusted and adjusted differences are available. The
#' propensity score row (`prop.score`) is removed from the plot.
#'
#' Covariates are taken from the `weightit` object (`w.out$covs`) and shown
#' in the order in which they appear in the model formula, with the first
#' covariate at the top. Factor covariates are split into one row per level,
#' labelled `"<label>: <level>"` (for example `"Sex: Female"`). When one
#' covariate name is the start of another (for example `age` and
#' `age_group`), each row is matched to the longest fitting name. Terms that
#' cannot be matched to an original covariate (for example interactions or
#' transformations such as `I(age^2)`) are kept, shown with their
#' `bal.tab()` name, and placed at the bottom.
#'
#' The x-axis runs from 0 to the larger of 0.6 and 110% of the largest
#' |SMD|, so the "Severe" zone is always visible.
#'
#' @param w.out An object of class `weightit`, as returned by
#'   [WeightIt::weightit()].
#' @param labels Optional named character vector of display labels. The names
#'   must be the covariate names used in the model formula and the values the
#'   labels to show, e.g. `c(sex = "Sex", age = "Age")`. Covariates without a
#'   label keep their original name. If `NULL` (default), all original names
#'   are used.
#' @param title Character string. Title of the plot.
#' @param binary Character string, `"std"` (default) or `"raw"`. How
#'   differences for binary covariates are computed; passed to
#'   [cobalt::bal.tab()].
#' @param continuous Character string, `"std"` (default) or `"raw"`. How
#'   differences for continuous covariates are computed; passed to
#'   [cobalt::bal.tab()].
#'
#' @return A `ggplot` object. It can be printed, extended with further
#'   `ggplot2` layers or themes, or saved with [ggplot2::ggsave()].
#'
#' @seealso [prop_score()] (argument `plot_love_iptw`), [cobalt::bal.tab()], [cobalt::love.plot()],
#'   [WeightIt::weightit()]
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(123)
#'   lung2 <- survival::lung
#'   lung2$sex       <- factor(lung2$sex, levels = c(1, 2),
#'                             labels = c("Male", "Female"))
#'   lung2$ph.ecog   <- factor(lung2$ph.ecog)
#'   lung2$treatment <- sample(c("Arm A", "Arm B"), nrow(lung2), replace = TRUE)
#'   lung2 <- lung2[stats::complete.cases(lung2[, c("sex", "age", "ph.ecog")]), ]
#'
#'   w.out <- WeightIt::weightit(treatment ~ sex + age + ph.ecog,
#'                               data = lung2, method = "ps",
#'                               estimand = "ATT")
#'
#'   plot_love_iptw(w.out,
#'                  labels = c(sex = "Sex", age = "Age", ph.ecog = "ECOG PS"))
#' }
#'
#' @importFrom rlang .data
#' @export
plot_love_iptw <- function(w.out,
                           labels     = NULL,
                           title      = "Covariate balance before and after IPTW",
                           binary     = "std",
                           continuous = "std") {

  # Input checks --------------------------------------------------------------
  if (!inherits(w.out, "weightit")) {
    stop("`w.out` must be an object of class 'weightit' (see WeightIt::weightit()).",
         call. = FALSE)
  }
  if (!is.null(labels) && (is.null(names(labels)) || any(names(labels) == ""))) {
    stop("`labels` must be a named character vector, e.g. c(sex = \"Sex\").",
         call. = FALSE)
  }
  binary     <- match.arg(binary, c("std", "raw"))
  continuous <- match.arg(continuous, c("std", "raw"))

  # 1. Balance table ----------------------------------------------------------
  bt <- cobalt::bal.tab(w.out, un = TRUE,
                        binary = binary, continuous = continuous)

  if (!is.null(bt$Balance)) {
    # binary treatment
    bal    <- bt$Balance
    col_un <- "Diff.Un"
    col_ad <- "Diff.Adj"
    x_lab  <- "Absolute standardized mean difference (|SMD|)"
  } else if (!is.null(bt$Balance.Across.Pairs)) {
    # multi-category treatment: maximum SMD across all pairs of groups
    bal    <- bt$Balance.Across.Pairs
    col_un <- "Max.Diff.Un"
    col_ad <- "Max.Diff.Adj"
    x_lab  <- "Maximum pairwise |SMD|"
  } else {
    stop("Could not find a balance table in the output of cobalt::bal.tab().",
         call. = FALSE)
  }
  bal <- bal[rownames(bal) != "prop.score", , drop = FALSE]

  # 2. Labels -----------------------------------------------------------------
  covs <- colnames(w.out$covs)
  if (is.null(labels)) labels <- stats::setNames(covs, covs)
  missing_lab <- setdiff(covs, names(labels))
  labels[missing_lab] <- missing_lab

  # Match each bal.tab row to its covariate (longest match wins)
  vars <- rownames(bal)
  prefix <- vapply(vars, function(v) {
    hits <- covs[v == covs | startsWith(v, paste0(covs, "_"))]
    if (length(hits) == 0) NA_character_ else hits[which.max(nchar(hits))]
  }, character(1))

  grp <- match(prefix, covs)
  grp[is.na(grp)] <- length(covs) + 1          # unmatched terms go last

  level <- ifelse(is.na(prefix), "", substring(vars, nchar(prefix) + 2))
  lab   <- unname(ifelse(is.na(prefix), vars,
                         ifelse(level == "", labels[prefix],
                                paste0(labels[prefix], ": ", level))))

  # Order follows the formula, first covariate on top
  ord     <- order(grp, seq_along(grp))
  lab_ord <- unique(lab[ord])

  # 3. Long format ------------------------------------------------------------
  df_lp <- data.frame(
    variable = rep(lab, 2),
    sample   = rep(c("Unweighted", "Weighted"), each = nrow(bal)),
    smd      = abs(c(bal[[col_un]], bal[[col_ad]]))
  )
  df_lp$variable <- factor(df_lp$variable, levels = rev(lab_ord))
  df_lp$sample   <- factor(df_lp$sample, levels = c("Unweighted", "Weighted"))

  n_var <- length(lab_ord)
  xmax  <- max(0.6, max(df_lp$smd, na.rm = TRUE) * 1.1)

  x_breaks <- c(0, 0.1, 0.25, 0.5)
  if (xmax > 0.75) x_breaks <- c(x_breaks, seq(0.75, xmax, by = 0.25))

  # Colours
  col_ok  <- "#2E9E5B"
  col_min <- "#E8C33D"
  col_mod <- "#E8873D"
  col_sev <- "#D64545"
  y_lab   <- n_var + 0.9

  # 4. Plot -------------------------------------------------------------------
  ggplot2::ggplot(df_lp, ggplot2::aes(x = .data$smd, y = .data$variable)) +
    # Shaded zones
    ggplot2::annotate("rect", xmin = 0,    xmax = 0.1,  ymin = -Inf, ymax = Inf, fill = col_ok,  alpha = 0.12) +
    ggplot2::annotate("rect", xmin = 0.1,  xmax = 0.25, ymin = -Inf, ymax = Inf, fill = col_min, alpha = 0.15) +
    ggplot2::annotate("rect", xmin = 0.25, xmax = 0.5,  ymin = -Inf, ymax = Inf, fill = col_mod, alpha = 0.14) +
    ggplot2::annotate("rect", xmin = 0.5,  xmax = xmax, ymin = -Inf, ymax = Inf, fill = col_sev, alpha = 0.12) +
    # Cut-off lines
    ggplot2::geom_vline(xintercept = 0.1,  linetype = "dashed", colour = col_ok,  linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = 0.25, linetype = "dashed", colour = col_mod, linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = 0.5,  linetype = "dashed", colour = col_sev, linewidth = 0.6) +
    # Zone labels
    ggplot2::annotate("text", x = 0.05,             y = y_lab, label = "Balanced", colour = col_ok,    fontface = "bold", size = 3.3) +
    ggplot2::annotate("text", x = 0.175,            y = y_lab, label = "Minor",    colour = "#B8962A", fontface = "bold", size = 3.3) +
    ggplot2::annotate("text", x = 0.375,            y = y_lab, label = "Moderate", colour = col_mod,   fontface = "bold", size = 3.3) +
    ggplot2::annotate("text", x = (0.5 + xmax) / 2, y = y_lab, label = "Severe",   colour = col_sev,   fontface = "bold", size = 3.3) +
    # Segment joining before and after
    ggplot2::geom_line(ggplot2::aes(group = .data$variable), colour = "grey65", linewidth = 0.5) +
    # Points
    ggplot2::geom_point(ggplot2::aes(colour = .data$sample, shape = .data$sample), size = 3.2) +
    ggplot2::scale_colour_manual(values = c("Unweighted" = "#8C8C8C", "Weighted" = "#2C5F8A")) +
    ggplot2::scale_shape_manual(values = c("Unweighted" = 16, "Weighted" = 17)) +
    ggplot2::scale_x_continuous(limits = c(0, xmax), breaks = x_breaks,
                                expand = ggplot2::expansion(mult = c(0.01, 0))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.6, 1.3))) +
    ggplot2::labs(x = x_lab,
                  y = NULL, colour = NULL, shape = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position    = "bottom",
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
      plot.title         = ggplot2::element_text(face = "bold"),
      axis.text.y        = ggplot2::element_text(colour = "grey20")
    )
}
