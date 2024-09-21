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
