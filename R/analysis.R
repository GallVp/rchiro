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
