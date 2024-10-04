#' Check for NaN Values in Data Frame
#'
#' This function checks each element of a data frame to determine if it is NaN
#' (Not a Number). It returns a logical matrix where each element is `TRUE` if
#' the corresponding element in the data frame is NaN, and `FALSE` otherwise.
#'
#' @param x A data frame to be checked for NaN values.
#' @return A logical matrix of the same dimensions as the input data frame,
#' with `TRUE` indicating NaN values and `FALSE` indicating non-NaN values.
#'
#' @examples
#' df <- data.frame(a = c(1, NaN, 3), b = c(NaN, 2, NaN))
#' is.nan(df)
#'
#' @export
is.nan.data.frame <- function(x) {
  do.call(cbind, lapply(x, is.nan))
}

#' Format P-Value with Standardized Precision
#'
#' This function formats a p-value to a standardized number of significant
#' digits. If the p-value is smaller than the threshold (10^(-sig_digits)), it is
#' returned as "<0.001" (or another appropriate threshold based on `sig_digits`),
#' otherwise, it is rounded to the specified number of significant digits.
#'
#' @param value A numeric value representing the p-value to be formatted.
#' @param sig_digits An integer specifying the number of significant digits to retain.
#' Default is 3.
#'
#' @return A character string representing the formatted p-value.
#'
#' @details
#' The function takes a p-value and checks if it is smaller than a certain
#' threshold determined by the `sig_digits` parameter. If the p-value is smaller than
#' this threshold, it returns a string indicating that the value is smaller than
#' the threshold. Otherwise, the p-value is returned rounded to the specified
#' number of significant digits.
#'
#' @examples
#' describe_p_value(0.0005) # "<0.001"
#' describe_p_value(0.12345) # "0.123"
#' describe_p_value(0.0001, 4) # "<0.0001"
#'
#' @export
describe_p_value <- function(value, sig_digits = 3) {
  text <- if (value < 10^(-sig_digits)) {
    paste0("<", 10^(-sig_digits))
  } else {
    round(value, sig_digits)
  }
  text
}

#' Calculate Standardized Beta Coefficients from a Linear Mixed Model
#'
#' This function calculates standardized beta coefficients from a fitted linear
#' mixed-effects model (using the `lmer` function from the `lme4` package).
#' Standardized beta coefficients are obtained by scaling the fixed effects
#' coefficients by the standard deviations of the predictor variables and
#' the outcome variable.
#'
#' @param mod A fitted linear mixed-effects model object of class `lmerMod`,
#' typically created with the `lmer` function from the `lme4` package.
#' @return A vector of standardized beta coefficients corresponding to the
#' fixed effects (excluding the intercept) in the model.
#' @importFrom lme4 fixef getME
#' @importFrom stats sd
#' @examples
#' library(lme4)
#' # Fit a linear mixed-effects model
#' mod <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' # Calculate standardized beta coefficients
#' std_beta_lmer(mod)
#' @export
std_beta_lmer <- function(mod) {
  b <- lme4::fixef(mod)[-1]
  x_mat <- lme4::getME(mod, "X")[, -1]
  sd_x <- if (is.vector(x_mat)) {
    sd(x_mat)
  } else {
    apply(x_mat, 2, sd)
  }
  sd_y <- sd(lme4::getME(mod, "y"))
  b * sd_x / sd_y
}

#' Compute Confidence Intervals for Fixed Effects in a Linear Mixed Model
#'
#' This function calculates confidence intervals for the fixed effects of a
#' fitted linear mixed-effects model (of class `rlmerMod`). The intervals are
#' computed based on the estimated coefficients and their standard errors using
#' a normal distribution.
#'
#' @param mod A fitted linear mixed-effects model object of class
#' `rlmerMod`. This is typically created using the `rlmer` function from
#' the `robustlmm` package.
#' @param parm A character vector specifying which parameters (fixed effects)
#' to compute confidence intervals for. If missing, confidence intervals for
#' all fixed effects are computed.
#' @param level A numeric value between 0 and 1 indicating the confidence
#' level for the intervals. The default is 0.95 for 95% confidence intervals.
#' @return A matrix of confidence intervals for the specified fixed effects.
#' Each row corresponds to a fixed effect, with columns representing the lower
#' and upper bounds of the confidence intervals.
#' @importFrom lme4 fixef
#' @importFrom stats vcov qnorm
#' @importFrom Matrix diag
#'
#' @examples
#' library(robustlmm)
#' # Fit a linear mixed-effects model
#' mod <- rlmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' # Compute 95% confidence intervals for all fixed effects
#' confint(mod)
#' # Compute 99% confidence intervals for specific parameters
#' confint(mod, parm = "Days", level = 0.99)
#'
#' @export
confint.rlmerMod <- function(mod, parm, level = 0.95) {
  beta <- lme4::fixef(mod)
  if (missing(parm)) {
    parm <- names(beta)
  }
  se <- sqrt(Matrix::diag(stats::vcov(mod)))
  z <- stats::qnorm((1 + level) / 2)
  c_tab <- cbind(beta - z * se, beta + z * se)
  colnames(c_tab) <- complmrob:::format.perc(
    c((1 - level) / 2, (1 + level) / 2),
    digits = 3
  )
  return(c_tab[parm, ])
}

#' Diagonal plot for linear mixed effects models
#'
#' This function creates a diagonal plot for a linear mixed effects model using
#' the \code{ggplot2} package. The plot includes a scatterplot of the residuals
#' versus the fitted values, a Q-Q plot of the residuals, and a histogram of the
#' residuals.
#'
#' @param model a linear mixed effects model created using the
#' \code{lme4} package
#'
#' @return a \code{ggplot2} object representing the diagonal plot
#'
#' @importFrom stats fitted resid
#' @importFrom ggplot2 aes geom_point
#' @importFrom ggpubr ggarrange ggqqplot gghistogram
#' @importFrom ggsci scale_color_nejm scale_fill_nejm
#'
#' @examples
#' library(lme4)
#' model <- lmer(Sepal.Length ~ Sepal.Width + (1 | Species), data = iris)
#' diag_plot_lmer(model)
#'
#' @export
diag_plot_lmer <- function(model) {
  # Get fitted values and residuals
  fitted_values <- stats::fitted(model)
  resid_values <- stats::resid(model)

  # Create data frame for plotting
  data_f <- data.frame(fitted_values, resid_values)

  # Create scatter plot of fitted values and residuals
  g1 <-
    ggplot2::ggplot(data_f, ggplot2::aes(x = fitted_values, y = resid_values)) +
    ggplot2::geom_point() +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.position = "none"
    ) +
    ggplot2::xlab("Fitted values") +
    ggplot2::ylab("Residuals")

  # Create QQ-plot of residuals
  g2 <-
    ggpubr::ggqqplot(resid(model)) +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.position = "none"
    ) +
    ggplot2::xlab("Theoretical") +
    ggplot2::ylab("Residuals")

  # Create histogram of residuals
  g3 <-
    ggpubr::gghistogram(resid(model), bins = 10) +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.position = "none"
    ) +
    ggplot2::xlab("Residuals") +
    ggplot2::ylab("Count")

  # Arrange plots into a grid
  figure <-
    ggpubr::ggarrange(
      g1,
      ggpubr::ggarrange(
        g2,
        g3,
        labels = c("B", "C"),
        ncol = 2,
        nrow = 1
      ),
      labels = c("A"),
      ncol = 1,
      nrow = 2
    )

  return(figure)
}

#' Standardize a Contrast Table
#'
#' This function standardizes a contrast table from model summary statistics,
#' including the estimates, standard errors, confidence intervals,
#' and test statistics. It adjusts column names and formats the results for
#' easy interpretation, providing concise representations of the differences
#' and test results.
#'
#' @param raw_contrasts_table A contrast table from a model summary, typically
#'   containing columns for estimates, standard errors, confidence intervals,
#'   and test statistics. This object can either contain asymptotic or
#'   confidence limits, which the function will detect and handle accordingly.
#'
#' @return A standardized data frame containing:
#'   \describe{
#'     \item{Contrast}{The name of the contrast tested.}
#'     \item{Difference±SE \[95% CI\]}{Formatted difference estimates with
#'     standard error and confidence intervals.}
#'     \item{Test}{Formatted test statistic (either `t` or `z` statistic),
#'     degrees of freedom, and p-value.}
#'   }
#'
#' @details
#' The function checks if the input table contains asymptotic or
#' regular confidence limits (`asymp.LCL`/`asymp.UCL` or `lower.CL`/`upper.CL`).
#' It also detects whether to use `z.ratio` or `t.ratio` for the test statistic.
#' The final table provides a concise summary of the contrasts, including
#' estimates, standard errors, confidence intervals, test statistics,
#' and p-values, all in a human-readable format.
#'
#' @examples
#' library(emmeans)
#' model <- lm(Sepal.Length ~ Sepal.Width + Species, data = iris)
#' contrast_table <- emmeans(model, pairwise ~ Species)$contrasts
#' result <- describe_contrasts(contrast_table)
#' print(result)
#'
#' @importFrom dplyr rowwise mutate
#' @export
describe_contrasts <- function(raw_contrasts_table) {
  if ("asymp.LCL" %in% names(data.frame(raw_contrasts_table)) ||
    "asymp.LCL" %in% names(
      data.frame(summary(raw_contrasts_table, infer = c(TRUE, TRUE)))
    )
  ) {
    ci_names <- c("asymp.LCL", "asymp.UCL")
    stat_names <- c("z.ratio", "z")
  } else {
    ci_names <- c("lower.CL", "upper.CL")
    stat_names <- c("t.ratio", "t")
  }

  a_table <- data.frame(summary(raw_contrasts_table, infer = c(TRUE, TRUE))) |>
    dplyr::rowwise() |>
    dplyr::mutate(sig_digits = ceiling(log10(1 / SE))) |>
    dplyr::rowwise() |>
    dplyr::mutate(Difference = paste0(
      round(estimate, sig_digits),
      "±",
      round(SE, sig_digits),
      " [",
      round(base::get(ci_names[1]), sig_digits),
      ", ",
      round(base::get(ci_names[2]), sig_digits), "]"
    ), Test = paste0(
      stat_names[2],
      "[", round(df, 1), "]", "=",
      round(base::get(stat_names[1]), 3), ", ",
      rchiro::describe_p_value(p.value)
    ))
  a_table$estimate <- NULL
  a_table$SE <- NULL
  a_table$df <- NULL
  a_table$lower.CL <- NULL
  a_table$upper.CL <- NULL
  a_table$asymp.LCL <- NULL
  a_table$asymp.UCL <- NULL
  a_table$sig_digits <- NULL
  a_table$t.ratio <- NULL
  a_table$z.ratio <- NULL
  a_table$p.value <- NULL
  a_table <- as.data.frame(a_table)

  a_table_names <- names(a_table)
  a_table_names[length(a_table_names) - 1] <- "Difference±SE [95% CI]"
  a_table_names[length(a_table_names)] <- paste0(stat_names[2], "[df], p-value")
  a_table_names[1] <- "Contrast"
  names(a_table) <- a_table_names

  a_table
}

#' Standardize an Estimated Marginal Means (EMM) Table
#'
#' This function standardizes a table of estimated marginal means (EMMs)
#' from model summary statistics. It reformats the EMM estimates,
#' standard errors, confidence intervals, and test statistics into a
#' human-readable format.
#'
#' @param raw_emmeans_table A table of estimated marginal means from a model summary,
#'   typically containing columns for EMM estimates, standard errors, confidence
#'   intervals, and test statistics. The function will detect if the table
#'   contains asymptotic or regular confidence limits.
#'
#' @return A standardized data frame containing:
#'   \describe{
#'     \item{Estimate±SE \[95% CI\]}{Formatted estimates with standard errors
#'     and 95% confidence intervals.}
#'     \item{Test}{Formatted test statistic (either `t` or `z` statistic),
#'     degrees of freedom, and p-value.}
#'   }
#'
#' @details
#' This function checks whether the input table includes asymptotic or regular
#' confidence limits (`asymp.LCL`/`asymp.UCL` or `lower.CL`/`upper.CL`) and
#' whether it should use `z.ratio` or `t.ratio` for the test statistic.
#' It processes the raw table and returns a user-friendly summary with
#' formatted estimates and test statistics.
#'
#' @examples
#'
#' library(emmeans)
#' model <- lm(Sepal.Length ~ Sepal.Width + Species, data = iris)
#' contrast_table <- emmeans(model, ~Species)
#' result <- describe_emmeans(contrast_table)
#' print(result)
#'
#' @importFrom dplyr rowwise mutate
#' @export
describe_emmeans <- function(raw_emmeans_table) {
  if ("asymp.LCL" %in% names(data.frame(raw_emmeans_table)) ||
    "asymp.LCL" %in% names(
      data.frame(summary(raw_emmeans_table, infer = c(TRUE, TRUE)))
    )
  ) {
    ci_names <- c("asymp.LCL", "asymp.UCL")
    stat_names <- c("z.ratio", "z")
  } else {
    ci_names <- c("lower.CL", "upper.CL")
    stat_names <- c("t.ratio", "t")
  }

  a_table <- data.frame(summary(raw_emmeans_table, infer = c(TRUE, TRUE))) |>
    dplyr::rowwise() |>
    dplyr::mutate(sig_digits = ceiling(log10(1 / SE))) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      Estimate = paste0(
        round(emmean, sig_digits),
        "±",
        round(SE, sig_digits),
        " [",
        round(base::get(ci_names[1]), sig_digits),
        ", ",
        round(base::get(ci_names[2]), sig_digits), "]"
      ),
      Test = paste0(
        stat_names[2], "[",
        round(df, 1), "]", "=",
        round(base::get(stat_names[1]), 3), ", ",
        rchiro::describe_p_value(p.value)
      )
    )

  a_table$emmean <- NULL
  a_table$SE <- NULL
  a_table$df <- NULL
  a_table$lower.CL <- NULL
  a_table$upper.CL <- NULL
  a_table$asymp.LCL <- NULL
  a_table$asymp.UCL <- NULL
  a_table$sig_digits <- NULL
  a_table$t.ratio <- NULL
  a_table$z.ratio <- NULL
  a_table$p.value <- NULL
  a_table <- as.data.frame(a_table)

  a_table_names <- names(a_table)
  a_table_names[length(a_table_names) - 1] <- "Estimate±SE [95% CI]"
  a_table_names[length(a_table_names)] <- paste0(stat_names[2], "[df], p-value")
  names(a_table) <- a_table_names

  a_table
}

#' Get Data Frame Column Names Except Specified
#'
#' This function returns the column names of a data frame excluding those specified.
#'
#' @param df A data frame from which column names will be retrieved.
#' @param ex A character vector of column names to exclude.
#'
#' @return A character vector of column names from the data frame, excluding those in \code{ex}.
#'
#' @examples
#' df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
#' names_of_df_except(df, c("b", "c"))
#'
#' @export
names_of_df_except <- function(df, ex) {
  names_all <- base::names(df)
  ex_which <- which(names_all %in% ex)

  names_all[-ex_which]
}

#' Get Data Frame Columns Except those Specified by Name
#'
#' This function returns the columns of a data frame excluding those specified.
#'
#' @param df A data frame from which column names will be retrieved.
#' @param ex A character vector of column names to exclude.
#'
#' @return A data frame with columns  excluding those in \code{ex}.
#'
#' @examples
#' df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
#' columns_of_df_except(df, c("b", "c"))
#'
#' @export
columns_of_df_except <- function(df, ex) {
  df[, df |> names_of_df_except(ex)]
}


#' Separate Baseline Scores in a Data Frame
#'
#' This function separates baseline data in a data frame by separating baseline and non-baseline records.
#'
#' @param df A data frame containing the data.
#' @param time_var A string specifying the name of the time variable (default is "Time").
#' @param baseline_lvl A string specifying the baseline level (default is "Baseline").
#' @param var_name A string specifying the name of the outcome variable (default is "Outcome").
#' @param id_vars A character vector of identifier variables (default is c("PartId", "Group")).
#'
#' @return A data frame that includes the non-baseline records and adds a column for the baseline outcome.
#'
#' @examples
#' df <- data.frame(
#'   PartId = c(1, 1, 2, 2),
#'   Group = c("A", "A", "B", "B"),
#'   Time = c("Baseline", "Post", "Baseline", "Post"),
#'   Outcome = c(10, 12, 20, 22)
#' )
#'
#' result <- separate_baseline_df(df)
#' print(result)
#'
#' @export

separate_baseline_df <- function(df,
                                 time_var = "Time",
                                 baseline_lvl = "Baseline",
                                 var_name = "Outcome",
                                 id_vars = c("PartId", "Group")) {
  df_long <- df
  df_long_pre <- subset(df_long, df[, time_var] == baseline_lvl)
  df_long_pre[, paste0(var_name, "Pre")] <- df_long_pre[, var_name]
  df_long_pre[, var_name] <- NULL
  df_long_pre[, time_var] <- NULL

  df_long_post <- subset(df_long, df[, time_var] != baseline_lvl)

  df_long <- df_long_post |>
    merge(
      df_long_pre,
      by = id_vars
    )

  df_long_pre <- NULL
  df_long_post <- NULL

  df_long
}
