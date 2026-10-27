rm(list = ls(all = TRUE))

devtools::install("~/Desktop/VHIO/R code/Package/VHIOstatSurvival")

library(VHIOstatSurvival)
library(tidyverse)

# ============================================================================ #
#                               Data Preparation                               #
# ============================================================================ #

lung <- survival::lung

lung2 <- lung |>
  mutate(event_os  = status - 1,
         time_os   = time,
         sex       = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
         ph.ecog   = factor(ph.ecog),
         treatment = sample(c("Arm A", "Arm B"), n(), replace = TRUE),
         time_pfs  = pmin(time_os, time_os * runif(n(), 0.4, 0.9)),
         event_pfs = ifelse(time_pfs < time_os, 1, event_os)) |>
  filter(!is.na(ph.ecog), !is.na(sex), !is.na(age))

ps.formula <- as.formula(paste("treatment", "~", paste(c("sex", "age", "ph.ecog"), collapse = "+")))

W.out <- WeightIt::weightit(ps.formula, data = lung2, method = "ps", estimand = "ATE")

lung2$weights_ATE <- W.out$weights



################################################################################
#                               Survival Plot                                  #
################################################################################

surv_plot(surv.formula = Surv(time_os, event_os) ~ treatment,
          dat          = lung2,
          #weights      = lung2$weights_ATE,
          palette      = c("dodgerblue3", "firebrick"),
          labels       = c("Arm A", "Arm B"), 
          labels.table.unadj = FALSE,
          # Increase the text size of risk table title labels
          #risk.table.title.font = c(14, "bold"),   # bigger, bold risk table title
          #risk.table.ytext.font = c(13),           # bigger risk table group labels
          # List of arguments for HR Table
          tableplot.args = list(median.show = TRUE,
                                position    = "topright",
                                x_margin = -0.04,
                                y_margin = -0.05,
                                # bigger HR table cell text
                                hr.table.core.cex    = 0.75,
                                hr.table.colhead.cex = 0.8,
                                hr.table.rowhead.cex = 0.75,
                                hr.table.line.lwd = 2
          ),
          # List of arguments for ggsurvplot
          ggsurvplot.args = list(ylab                  = "Overall survival (%)",
                                 xlab                  = "Time in Months until Death",
                                 tables.theme          = survminer::theme_cleantable(),
                                 title                 = "Test Survival Curve",
                                 pval                  = F, 
                                 pval.method           = F,   
                                 censor                = T,
                                 censor.shape          = 124,
                                 censor.size           = 1.5,
                                 conf.int              = FALSE,
                                 legend                = "none",
                                 legend.title          = "",
                                 surv.median.line      = "hv",
                                 tables.height         = 0.2,
                                 tables.y.text.col     = T,
                                 tables.y.text         = T,
                                 break.time.by         = 365,
                                 size                  = 1,
                                 linewidth             = 1,
                                 risk.table            = TRUE,
                                 risk.table.title      = "No. at Risk",
                                 font.title            = c(12, "bold"),    # title font
                                 font.x                = c(10, "bold"),        # x-axis label font
                                 font.y                = c(10, "bold"),        # y-axis label font
                                 font.tickslab         = c(10),
                                 risk.table.fontsize   = 3,
                                 risk.table.size       = 5,
                                 axes.offset           = T,
                                 risk.table.y.fontsize = 1))


surv_plot(surv.formula = Surv(time_os, event_os) ~ treatment,
          dat          = lung2,
          weights      = lung2$weights_ATE,
          palette      = c("dodgerblue3", "firebrick"),
          labels       = c("Arm A", "Arm B"), 
          labels.table.unadj = FALSE,
          # Increase the text size of risk table title labels
          #risk.table.title.font = c(14, "bold"),   # bigger, bold risk table title
          #risk.table.ytext.font = c(13),           # bigger risk table group labels
          # List of arguments for HR Table
          tableplot.args = list(median.show = TRUE,
                                position    = "topright",
                                x_margin = -0.04,
                                y_margin = -0.05,
                                surv.times = c(365),
                                # bigger HR table cell text
                                hr.table.core.cex    = 0.75,
                                hr.table.colhead.cex = 0.8,
                                hr.table.rowhead.cex = 0.75,
                                hr.table.line.lwd = 2
          ),
          # List of arguments for ggsurvplot
          ggsurvplot.args = list(ylab                  = "Overall survival (%)",
                                 xlab                  = "Time in Months until Death",
                                 tables.theme          = survminer::theme_cleantable(),
                                 title                 = "Test Survival Curve",
                                 pval                  = F, 
                                 pval.method           = F,   
                                 censor                = T,
                                 censor.shape          = 124,
                                 censor.size           = 1.5,
                                 conf.int              = FALSE,
                                 legend                = "none",
                                 legend.title          = "",
                                 surv.median.line      = "hv",
                                 tables.height         = 0.2,
                                 tables.y.text.col     = T,
                                 tables.y.text         = T,
                                 break.time.by         = 365,
                                 size                  = 1,
                                 linewidth             = 1,
                                 risk.table            = TRUE,
                                 risk.table.title      = "No. at Risk",
                                 font.title            = c(12, "bold"),    # title font
                                 font.x                = c(10, "bold"),        # x-axis label font
                                 font.y                = c(10, "bold"),        # y-axis label font
                                 font.tickslab         = c(10),
                                 risk.table.fontsize   = 3,
                                 risk.table.size       = 5,
                                 axes.offset           = T,
                                 risk.table.y.fontsize = 1))


surv_plot(surv.formula = Surv(time_os, event_os) ~ treatment,
          dat          = lung2,
          weights      = lung2$weights_ATE,
          palette      = c("dodgerblue3", "firebrick"),
          labels       = c("Arm A", "Arm B"), 
          labels.table.unadj = FALSE,
          # Increase the text size of risk table title labels
          #risk.table.title.font = c(14, "bold"),   # bigger, bold risk table title
          #risk.table.ytext.font = c(13),           # bigger risk table group labels
          # List of arguments for HR Table
          tableplot.args = list(median.show = TRUE,
                                HR.show = FALSE,
                                position    = "topright",
                                x_margin = -0.04,
                                y_margin = -0.05,
                                surv.times = c(365),
                                # bigger HR table cell text
                                hr.table.core.cex    = 0.75,
                                hr.table.colhead.cex = 0.8,
                                hr.table.rowhead.cex = 0.75,
                                hr.table.line.lwd = 2
          ),
          # List of arguments for ggsurvplot
          ggsurvplot.args = list(ylab                  = "Overall survival (%)",
                                 xlab                  = "Time in Months until Death",
                                 tables.theme          = survminer::theme_cleantable(),
                                 title                 = "Test Survival Curve",
                                 pval                  = F, 
                                 pval.method           = F,   
                                 censor                = T,
                                 censor.shape          = 124,
                                 censor.size           = 1.5,
                                 conf.int              = FALSE,
                                 legend                = "none",
                                 legend.title          = "",
                                 surv.median.line      = "hv",
                                 tables.height         = 0.2,
                                 tables.y.text.col     = T,
                                 tables.y.text         = T,
                                 break.time.by         = 365,
                                 size                  = 1,
                                 linewidth             = 1,
                                 risk.table            = TRUE,
                                 risk.table.title      = "No. at Risk",
                                 font.title            = c(12, "bold"),    # title font
                                 font.x                = c(10, "bold"),        # x-axis label font
                                 font.y                = c(10, "bold"),        # y-axis label font
                                 font.tickslab         = c(10),
                                 risk.table.fontsize   = 3,
                                 risk.table.size       = 5,
                                 axes.offset           = T,
                                 risk.table.y.fontsize = 1))


surv_plot(surv.formula = Surv(time_os, event_os) ~ treatment + sex,
          dat          = lung2,
          weights      = lung2$weights_ATE,
          #palette      = c("dodgerblue3", "firebrick"),
          #labels       = c("Arm A", "Arm B"), 
          labels.table.unadj = FALSE,
          # Increase the text size of risk table title labels
          #risk.table.title.font = c(14, "bold"),   # bigger, bold risk table title
          #risk.table.ytext.font = c(13),           # bigger risk table group labels
          # List of arguments for HR Table
          tableplot.args = list(median.show = TRUE,
                                position    = "topright",
                                x_margin = -0.04,
                                y_margin = -0.05,
                                surv.times = c(365),
                                # bigger HR table cell text
                                hr.table.core.cex    = 0.75,
                                hr.table.colhead.cex = 0.8,
                                hr.table.rowhead.cex = 0.75,
                                hr.table.line.lwd = 2
          ),
          # List of arguments for ggsurvplot
          ggsurvplot.args = list(ylab                  = "Overall survival (%)",
                                 xlab                  = "Time in Months until Death",
                                 tables.theme          = survminer::theme_cleantable(),
                                 title                 = "Test Survival Curve",
                                 pval                  = F, 
                                 pval.method           = F,   
                                 censor                = T,
                                 censor.shape          = 124,
                                 censor.size           = 1.5,
                                 conf.int              = FALSE,
                                 legend                = "none",
                                 legend.title          = "",
                                 surv.median.line      = "hv",
                                 tables.height         = 0.2,
                                 tables.y.text.col     = T,
                                 tables.y.text         = T,
                                 break.time.by         = 365,
                                 size                  = 1,
                                 linewidth             = 1,
                                 risk.table            = TRUE,
                                 risk.table.title      = "No. at Risk",
                                 font.title            = c(12, "bold"),    # title font
                                 font.x                = c(10, "bold"),        # x-axis label font
                                 font.y                = c(10, "bold"),        # y-axis label font
                                 font.tickslab         = c(10),
                                 risk.table.fontsize   = 3,
                                 risk.table.size       = 5,
                                 axes.offset           = T,
                                 risk.table.y.fontsize = 1))


surv_plot(surv.formula = Surv(time_os, event_os) ~ treatment + sex,
          dat          = lung2,
          weights      = lung2$weights_ATE,
          #palette      = c("dodgerblue3", "firebrick"),
          labels       = c("Arm A + Male", "Arm A + Female", "Arm B + Male", "Arm B + Female"), 
          labels.table.unadj = FALSE,
          # Increase the text size of risk table title labels
          #risk.table.title.font = c(14, "bold"),   # bigger, bold risk table title
          #risk.table.ytext.font = c(13),           # bigger risk table group labels
          # List of arguments for HR Table
          tableplot.args = list(median.show = TRUE,
                                position    = "topright",
                                x_margin = -0.04,
                                y_margin = -0.05,
                                surv.times = c(365),
                                # bigger HR table cell text
                                hr.table.core.cex    = 0.75,
                                hr.table.colhead.cex = 0.8,
                                hr.table.rowhead.cex = 0.75,
                                hr.table.line.lwd = 2
          ),
          # List of arguments for ggsurvplot
          ggsurvplot.args = list(ylab                  = "Overall survival (%)",
                                 xlab                  = "Time in Months until Death",
                                 tables.theme          = survminer::theme_cleantable(),
                                 title                 = "Test Survival Curve",
                                 pval                  = F, 
                                 pval.method           = F,   
                                 censor                = T,
                                 censor.shape          = 124,
                                 censor.size           = 1.5,
                                 conf.int              = FALSE,
                                 legend                = "none",
                                 legend.title          = "",
                                 surv.median.line      = "hv",
                                 tables.height         = 0.2,
                                 tables.y.text.col     = T,
                                 tables.y.text         = T,
                                 break.time.by         = 365,
                                 size                  = 1,
                                 linewidth             = 1,
                                 risk.table            = TRUE,
                                 risk.table.title      = "No. at Risk",
                                 font.title            = c(12, "bold"),    # title font
                                 font.x                = c(10, "bold"),        # x-axis label font
                                 font.y                = c(10, "bold"),        # y-axis label font
                                 font.tickslab         = c(10),
                                 risk.table.fontsize   = 3,
                                 risk.table.size       = 5,
                                 axes.offset           = T,
                                 risk.table.y.fontsize = 1))


################################################################################
#                            Survival Forest Plot                              #
################################################################################

# ============================================================================ #
#                          Uni- & Multivariable Cox                            #
# ============================================================================ #


res1 <- cox_forest(
  data = lung2,
  # Specify here all covariates that you want to include
  covariates = c("treatment", "sex", "ph.ecog", "age"),
  # If you have a continuous variable please specify it here (didn't implement
  # the automatic search (next step))
  continuous_vars = "age",
  # Optional: If you want a x unit increase for a continuous increase specify it here
  continuous_increase = list(age = 10),
  # Specify here the names you want to have to be displayed
  # for each covariate
  labels = c("Treatment", "Sex", "ECOG PS", "Age"),
  # Specify the time and event column names of the endpoint of interest
  time = c( "time_pfs", "time_os"),
  event = c("event_pfs", "event_os"),
  # Names of each endpoint
  endpoint_labels = c("Progression-free Survival", "Overall Survival"),
  # Optional: if you want to change the reference level of a covariate
  ref_levels = list(sex = "Male"),
  # Optional: can exlude labels of each covariate just for the displaying 
  # if for example only 4 patient are in one group the Cox model will be 
  # obtained with those patient and group but not displayed in the forestplot
  # because it has no meaning
  exclude_levels = list(ph.ecog = "3"),
  # Specify which model type you want ("univariable", "multivariable", "both")
  model_type = "both",
  # Colour palette
  palette = c("#0F084B", "#A6CFD5"),
  # Limits of the x-axis for the forestplot
  xlim = c(0.25, 4),
  # Ticks at the x-axis
  ticks_at = c(0.5, 1, 2),
  # Width of the forestplot
  length_plot = 60,
  # Optional: only if using PS
  #weights = "weights_ATE", # optional if PS used
  # Argument for the displaying of the CI (check forestmodel::forest for each meaning)
  ci_pch = 15, 
  ci_lty = 1,
  ci_lwd = 3, 
  ci_Theight = 0.4, 
  sizes = 1,
  # Optional: Specify row height 
  row_height = 1.4,
  # Optional: TRUE creates a seperate row for the variable name
  covariate_label_row = FALSE,
  # Only use of covariate_label_row = TRUE,
  # label_indent = 4,
  # Optional: increases the space between the first 3 columns
  col_spacing = c(0, 0)
)

# Ploting univariable model
res1$plots$univariable

# Ploting multivariable model
res1$plots$multivariable

# Table with all uni and multivariable results
res1$results_table




res2 <- cox_forest(
  data = lung2,
  # Specify here all covariates that you want to include
  covariates = c("treatment", "sex", "ph.ecog", "age"),
  # If you have a continuous variable please specify it here (didn't implement
  # the automatic search (next step))
  continuous_vars = "age",
  # Optional: If you want a x unit increase for a continuous increase specify it here
  continuous_increase = list(age = 10),
  # Specify here the names you want to have to be displayed
  # for each covariate
  labels = c("Treatment", "Sex", "ECOG PS", "Age"),
  # Specify the time and event column names of the endpoint of interest
  time = c( "time_pfs", "time_os"),
  event = c("event_pfs", "event_os"),
  # Names of each endpoint
  endpoint_labels = c("Progression-free Survival", "Overall Survival"),
  # Optional: if you want to change the reference level of a covariate
  ref_levels = list(sex = "Male"),
  # Optional: can exlude labels of each covariate just for the displaying 
  # if for example only 4 patient are in one group the Cox model will be 
  # obtained with those patient and group but not displayed in the forestplot
  # because it has no meaning
  exclude_levels = list(ph.ecog = "3"),
  # Specify which model type you want ("univariable", "multivariable", "both")
  model_type = "both",
  # Colour palette
  palette = c("#0F084B", "#A6CFD5"),
  # Limits of the x-axis for the forestplot
  xlim = c(0.25, 4),
  # Ticks at the x-axis
  ticks_at = c(0.5, 1, 2),
  # Width of the forestplot
  length_plot = 60,
  # Optional: only if using PS
  #weights = "weights_ATE", # optional if PS used
  # Argument for the displaying of the CI (check forestmodel::forest for each meaning)
  ci_pch = 15, 
  ci_lty = 1,
  ci_lwd = 3, 
  ci_Theight = 0.4, 
  sizes = 1,
  # Optional: Specify row height 
  row_height = 1.4,
  # Optional: TRUE creates a seperate row for the variable name
  covariate_label_row = TRUE,
  # Only use of covariate_label_row = TRUE,
  label_indent = 4,
  # Optional: increases the space between the first 3 columns
  col_spacing = c(0, 0)
)

# Ploting univariable model
res2$plots$univariable

# Ploting multivariable model
res2$plots$multivariable



# ============================================================================ #
#                                  Subgroup Cox                                #
# ============================================================================ #


res3 <- cox_subgroup_forest(
  data = lung2,
  # Specify here all covariates that you want to include
  covariates = c("sex", "ph.ecog", "age"),
  # Specify here the names you want to have to be displayed
  # for each covariate
  labels = c("Sex", "ECOG PS", "Age"),
  # Treatment variable
  treatment = "treatment",
  # Optional: order of the treatment labels
  # Default: order of treatment variable
  arm_order = c("Arm A", "Arm B"),
  # Optional: rename the labels of the treatment variable
  #           in the order of the arm order
  # treatment_labels = c("A", "B"),
  # Specify the time and event column names of the endpoint of interest
  time = c( "time_os", "time_pfs"),
  event = c("event_os", "event_pfs"),
  # Names of each endpoint
  endpoint_labels = c("Overall Survival", "Progression-free Survival"),
  # Specify the time unit of the survival time
  time_unit = "days",
  # Colour palette
  palette = c("#0F084B", "#A6CFD5"),
  # Optional: can exlude labels of each covariate just for the displaying 
  # if for example only 4 patient are in one group the Cox model will be 
  # obtained with those patient and group but not displayed in the forestplot
  # because it has no meaning
  # exclude_levels = list(ph.ecog = "3"),
  # Limits of the x-axis for the forestplot
  xlim = c(0.25, 4),
  # Ticks at the x-axis
  ticks_at = c(0.5, 1, 2),
  # Width of the forestplot
  length_plot = 60,
  # Optional: only if using PS
  # weights = "weights_ATE", # optional if PS used
  # Argument for the displaying of the CI (check forestmodel::forest for each meaning)
  ci_pch = 15, 
  ci_lty = 1,
  ci_lwd = 3, 
  ci_Theight = 0.4, 
  sizes = 1,
  # Optional: Specify row height 
  row_height = 1.4,
  # Optional: TRUE creates a seperate row for the variable name
  covariate_label_row = TRUE,
  # Only use of covariate_label_row = TRUE,
  label_indent = 4,
  # Optional: increases the space between the first 3 columns
  col_spacing = c(0, 0)
)

res3$plots$Age
res3$plots$`ECOG PS`
res3$plots$Sex

res3$combined_plots



################################################################################
#                             General Forest Plot                              #
################################################################################


res <- model_forest(
  data = lung2,
  # Formula you want to fit
  formula = list(as.formula("Surv(time_os, event_os) ~ treatment * sex")),
  # Function used to fit the model
  fit_fun = survival::coxph,
  # Endpoint label of interest
  endpoint_labels = "Overall Survival",
  # Optional: if you want to change the reference level of a covariate
  ref_levels = NULL, 
  # Optional: can exlude labels of each covariate just for the displaying 
  # if for example only 4 patient are in one group the Cox model will be 
  # obtained with those patient and group but not displayed in the forestplot
  # because it has no meaning
  exclude_levels = NULL,
  # Exponentiate the estimate of the outcome of interest in the forestplot
  exponentiate = TRUE, 
  # Optional: Label of the effect of interest (auto: HR / OR / RR / Estimate) 
  # or a string.
  effect_label = NULL, 
  # Confidence Level
  ci_level = 0.95,
  # Digits 
  digits = 2, 
  # Colour palette
  palette = "#0F084B", 
  # Limits of the x-axis for the forestplot
  xlim = c(0.25, 4),
  # Ticks at the x-axis
  ticks_at = c(0.5, 1, 2),
  # Width of the forestplot
  length_plot = 60, 
  # Argument for the displaying of the CI (check forestmodel::forest for each meaning)
  ci_pch = 15, 
  ci_lty = 1,
  ci_lwd = 3, 
  ci_Theight = 0.4, 
  sizes = 1,
  # Optional: Specify row height 
  row_height = 1.4,
  # Optional: TRUE creates a seperate row for the variable name
  covariate_label_row = TRUE,
  # Only use of covariate_label_row = TRUE,
  label_indent = 4,
  # Optional: increases the space between the first 3 columns
  col_spacing = c(0, 0))

res$plot
res$models
res$results_table


lung2 <- lung2 |> mutate(treatment_bin = as.numeric(as.factor(lung2$treatment))-1)

res_glm <- model_forest(formula = list(as.formula("treatment_bin ~ sex + log(age) + ph.ecog")),
                        data = lung2,
                        fit_fun = glm,
                        family = "binomial",
                        endpoint_labels = "Overall Survival",
                        ref_levels = NULL, 
                        exclude_levels = NULL,
                        exponentiate = TRUE, 
                        effect_label = NULL, 
                        ci_level = 0.95,
                        digits = 2, 
                        palette = "#0F084B", 
                        xlim = c(0.25, 4), 
                        ticks_at = c(0.5, 1, 2),
                        length_plot = 60, 
                        ci_pch = 15, 
                        ci_lty = 1, 
                        ci_lwd = 3,
                        ci_Theight = 0.4, 
                        sizes = 1, 
                        risk_labels = NULL,
                        title = NULL, 
                        row_height = 1.5, 
                        covariate_label_row = FALSE,
                        col_spacing = c(1, 1))

res_glm$plot


set.seed(1)
d <- data.frame(
  resp = rbinom(300, 1, 0.4),
  arm  = factor(sample(c("Control", "Treatment"), 300, TRUE)),
  sex  = factor(sample(c("F", "M"), 300, TRUE)),
  age  = rnorm(300, 60, 10)
)

res <- model_forest(formula = as.formula(resp ~ arm + sex + I(age/10)), 
                    data = d, #family = binomial,
                    fit_fun = glm,
                    covariate_label_row = TRUE,
                    ci_level = 0.9)
plot(res$plot)

res$results_table














