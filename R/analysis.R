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
#' digits. If the p-value is smaller than the threshold (10^(-sigD)), it is
#' returned as "<0.001" (or another appropriate threshold based on `sigD`),
#' otherwise, it is rounded to the specified number of significant digits.
#'
#' @param value A numeric value representing the p-value to be formatted.
#' @param sigD An integer specifying the number of significant digits to retain.
#' Default is 3.
#'
#' @return A character string representing the formatted p-value.
#'
#' @details
#' The function takes a p-value and checks if it is smaller than a certain
#' threshold determined by the `sigD` parameter. If the p-value is smaller than
#' this threshold, it returns a string indicating that the value is smaller than
#' the threshold. Otherwise, the p-value is returned rounded to the specified
#' number of significant digits.
#'
#' @examples
#' standard.p.value(0.0005) # "<0.001"
#' standard.p.value(0.12345) # "0.123"
#' standard.p.value(0.0001, 4) # "<0.0001"
#'
#' @export
standard.p.value <- function(value, sigD = 3) {
  text <- if (value < 10^(-sigD)) {
    paste0("<", 10^(-sigD))
  } else {
    round(value, sigD)
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
#' lmer.std.beta(mod)
#' @export
lmer.std.beta <- function(mod) {
  b <- lme4::fixef(mod)[-1]
  x.mat <- lme4::getME(mod, "X")[, -1]
  sd.x <- if (is.vector(x.mat)) {
    sd(x.mat)
  } else {
    apply(x.mat, 2, sd)
  }
  sd.y <- sd(lme4::getME(mod, "y"))
  b * sd.x / sd.y
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
  ctab <- cbind(beta - z * se, beta + z * se)
  colnames(ctab) <- complmrob:::format.perc(
    c((1 - level) / 2, (1 + level) /
      2),
    digits = 3
  )
  return(ctab[parm, ])
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
#' diag.plot.lmer(model)
#'
#' @export
diag.plot.lmer <- function(model) {
  # Get fitted values and residuals
  fitted.vals <- stats::fitted(model)
  resid.vals <- stats::resid(model)

  # Create data frame for plotting
  data.f <- data.frame(fitted.vals, resid.vals)

  # Create scatter plot of fitted values and residuals
  g1 <-
    ggplot2::ggplot(data.f, ggplot2::aes(x = fitted.vals, y = resid.vals)) +
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
#' @param raw.ctt.table A contrast table from a model summary, typically
#'   containing columns for estimates, standard errors, confidence intervals,
#'   and test statistics. This object can either contain asymptotic or
#'   confidence limits, which the function will detect and handle accordingly.
#'
#' @return A standardized data frame containing:
#'   \describe{
#'     \item{Contrast}{The name of the contrast tested.}
#'     \item{Difference±SE [95% CI]}{Formatted difference estimates with
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
#' result <- standard.ctt.table(contrast_table)
#' print(result)
#'
#' @importFrom dplyr rowwise mutate
#' @export
standard.ctt.table <- function(raw.ctt.table) {
  if ("asymp.LCL" %in% names(data.frame(raw.ctt.table)) ||
    "asymp.LCL" %in% names(
      data.frame(summary(raw.ctt.table, infer = c(T, T)))
    )
  ) {
    ci.names <- c("asymp.LCL", "asymp.UCL")
    stat.names <- c("z.ratio", "z")
  } else {
    ci.names <- c("lower.CL", "upper.CL")
    stat.names <- c("t.ratio", "t")
  }

  aT <- data.frame(summary(raw.ctt.table, infer = c(T, T))) |>
    dplyr::rowwise() |>
    dplyr::mutate(sigD = ceiling(log10(1 / SE))) |>
    dplyr::rowwise() |>
    dplyr::mutate(Difference = paste0(
      round(estimate, sigD),
      "±",
      round(SE, sigD),
      " [",
      round(base::get(ci.names[1]), sigD),
      ", ",
      round(base::get(ci.names[2]), sigD), "]"
    ), Test = paste0(
      stat.names[2],
      "[", round(df, 1), "]", "=",
      round(base::get(stat.names[1]), 3), ", ",
      rchiro::standard.p.value(p.value)
    ))
  aT$estimate <- NULL
  aT$SE <- NULL
  aT$df <- NULL
  aT$lower.CL <- NULL
  aT$upper.CL <- NULL
  aT$asymp.LCL <- NULL
  aT$asymp.UCL <- NULL
  aT$sigD <- NULL
  aT$t.ratio <- NULL
  aT$z.ratio <- NULL
  aT$p.value <- NULL
  aT <- as.data.frame(aT)

  aT.names <- names(aT)
  aT.names[length(aT.names) - 1] <- "Difference±SE [95% CI]"
  aT.names[length(aT.names)] <- paste0(stat.names[2], "[df], p-value")
  aT.names[1] <- "Contrast"
  names(aT) <- aT.names

  aT
}

#' Standardize an Estimated Marginal Means (EMM) Table
#'
#' This function standardizes a table of estimated marginal means (EMMs)
#' from model summary statistics. It reformats the EMM estimates,
#' standard errors, confidence intervals, and test statistics into a
#' human-readable format.
#'
#' @param raw.em.table A table of estimated marginal means from a model summary,
#'   typically containing columns for EMM estimates, standard errors, confidence
#'   intervals, and test statistics. The function will detect if the table
#'   contains asymptotic or regular confidence limits.
#'
#' @return A standardized data frame containing:
#'   \describe{
#'     \item{Estimate±SE [95% CI]}{Formatted estimates with standard errors
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
#' result <- standard.em.table(contrast_table)
#' print(result)
#'
#' @importFrom dplyr rowwise mutate
#' @export
standard.em.table <- function(raw.em.table) {
  if ("asymp.LCL" %in% names(data.frame(raw.em.table)) ||
    "asymp.LCL" %in% names(
      data.frame(summary(raw.em.table, infer = c(T, T)))
    )
  ) {
    ci.names <- c("asymp.LCL", "asymp.UCL")
    stat.names <- c("z.ratio", "z")
  } else {
    ci.names <- c("lower.CL", "upper.CL")
    stat.names <- c("t.ratio", "t")
  }

  aT <- data.frame(summary(raw.em.table, infer = c(T, T))) |>
    dplyr::rowwise() |>
    dplyr::mutate(sigD = ceiling(log10(1 / SE))) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      Estimate = paste0(
        round(emmean, sigD),
        "±",
        round(SE, sigD),
        " [",
        round(base::get(ci.names[1]), sigD),
        ", ",
        round(base::get(ci.names[2]), sigD), "]"
      ),
      Test = paste0(
        stat.names[2], "[",
        round(df, 1), "]", "=",
        round(base::get(stat.names[1]), 3), ", ",
        rchiro::standard.p.value(p.value)
      )
    )

  aT$emmean <- NULL
  aT$SE <- NULL
  aT$df <- NULL
  aT$lower.CL <- NULL
  aT$upper.CL <- NULL
  aT$asymp.LCL <- NULL
  aT$asymp.UCL <- NULL
  aT$sigD <- NULL
  aT$t.ratio <- NULL
  aT$z.ratio <- NULL
  aT$p.value <- NULL
  aT <- as.data.frame(aT)

  aT.names <- names(aT)
  aT.names[length(aT.names) - 1] <- "Estimate±SE [95% CI]"
  aT.names[length(aT.names)] <- paste0(stat.names[2], "[df], p-value")
  names(aT) <- aT.names

  aT
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
#' xnames.data.frame(df, c("b", "c"))
#'
#' @export 
xnames.data.frame <- function(df, ex) {
  names.all <- base::names(df)
  ex.which <- which(names.all %in% ex)
  
  names.all[-ex.which]
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
#' columns.data.frame(df, c("b", "c"))
#'
#' @export 
columns.data.frame <- function(df, ex) {
  df[, df |> xnames.data.frame(ex)]
}