#' Check for NaN Values in Data Frame
#'
#' This function checks each element of a data frame to determine if it is NaN
#' (Not a Number). It returns a logical matrix where each element is `TRUE` if
#' the corresponding element in the data frame is NaN, and `FALSE` otherwise.
#'
#' @param x A data frame to be checked for NaN values.
#' @return A logical matrix of the same dimensions as the input data frame,
#' with `TRUE` indicating NaN values and `FALSE` indicating non-NaN values.
#' @examples
#' df <- data.frame(a = c(1, NaN, 3), b = c(NaN, 2, NaN))
#' is.nan.data.frame(df)
#' # Returns:
#' #      a     b
#' # FALSE  TRUE
#' # TRUE  FALSE
#' # FALSE  TRUE
#' @export
is.nan.data.frame <- function(x) {
  do.call(cbind, lapply(x, is.nan))
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
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod, "X")[, -1], 2, sd)
  sd.y <- sd(getME(mod, "y"))
  b * sd.x / sd.y
}


#' Compute Confidence Intervals for Fixed Effects in a Linear Mixed Model
#'
#' This function calculates confidence intervals for the fixed effects of a
#' fitted linear mixed-effects model (of class `rlmerMod`). The intervals are
#' computed based on the estimated coefficients and their standard errors using
#' a normal distribution.
#'
#' @param object A fitted linear mixed-effects model object of class
#' `rlmerMod`. This is typically created using the `rlmer` function from
#' the `rstanarm` package.
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
#' @examples
#' library(rstanarm)
#' # Fit a linear mixed-effects model
#' mod <- rlmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' # Compute 95% confidence intervals for all fixed effects
#' confint.rlmerMod(mod)
#' # Compute 99% confidence intervals for specific parameters
#' confint.rlmerMod(mod, parm = "Days", level = 0.99)
#' @export
confint.rlmerMod <- function(object, parm, level = 0.95) {
  beta <- fixef(object)
  if (missing(parm)) {
    parm <- names(beta)
  }
  se <- sqrt(diag(vcov(object)))
  z <- qnorm((1 + level) / 2)
  ctab <- cbind(beta - z * se, beta + z * se)
  colnames(ctab) <- stats:::format.perc(
    c((1 - level) / 2, (1 + level) /
      2),
    digits = 3
  )
  return(ctab[parm, ])
}


#' Convert a `rlmerMod` Model to a `glm` Model
#'
#' This function converts a fitted `rlmerMod` model object (from the `rstanarm`
#' package) to a generalized linear model (`glm`) object. This is done by
#'  extracting the fixed effects and variance-covariance matrix from the
#' `rlmerMod` model and creating a new `glm` model with equivalent fixed
#' effects.
#'
#' @param mod A fitted `rlmerMod` model object from the `rstanarm` package.
#' @return A `glm` model object with the following attributes:
#' \itemize{
#'   \item \code{coefficients}: Fixed effects coefficients from the original
#' model.
#'   \item \code{vcov}: Variance-covariance matrix of the fixed effects
#' coefficients.
#'   \item \code{linear.predictors}: Linear predictors for the `glm` model.
#'   \item \code{fitted.values}: Fitted values based on the `glm` model.
#'   \item \code{weights}: Weights used in the `glm` model.
#'   \item \code{residuals}: Residuals calculated for the `glm` model.
#' }
#' The object is assigned a new class `fakeglm` in addition to the standard
#' `glm` class.
#' @details The function extracts the call, fixed effects, and
#' variance-covariance matrix from the `rlmerMod` model. It then constructs
#' a new `glm` object with these components and adjusts its attributes to
#' match the original model's structure.
#' @importFrom stats gaussian glm.control
#' @importFrom lme4 fixef
#' @importFrom stats vcov
#' @examples
#' library(rstanarm)
#' # Fit a linear mixed-effects model
#' mod <- rlmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' # Convert the model to a glm object
#' mod_glm <- rlmerMod.to.glm(mod)
#' # Examine the glm object
#' summary(mod_glm)
#' @export
rlmerMod.to.glm <- function(mod) {
  family <- gaussian()
  link <- family$link
  family <- family$family
  cl <- mod@call
  cl$control <- glm.control(epsilon = 1)
  .m <- match(
    c(
      "formula",
      "family",
      "data",
      "weights",
      "subset",
      "na.action",
      "offset",
      "model",
      "contrasts"
    ),
    names(cl),
    0L
  )
  cl <- cl[c(1L, .m)]
  cl[[1L]] <- as.name("glm")
  cl$formula <- effects:::fixFormula(as.formula(cl$formula))
  mod2 <- eval(cl)
  mod2$coefficients <- lme4::fixef(mod)
  mod2$vcov <- as.matrix(vcov(mod))
  mod2$linear.predictors <- model.matrix(mod2) %*% mod2$coefficients
  mod2$fitted.values <- mod2$family$linkinv(mod2$linear.predictors)
  mod2$weights <-
    as.vector(with(mod2, prior.weights * (
      family$mu.eta(linear.predictors)^2 / family$variance(fitted.values)
    )))
  mod2$residuals <-
    with(mod2, prior.weights * (y - fitted.values) / weights)
  class(mod2) <- c("fakeglm", class(mod2))
  mod2
}

#' Calculate Over dispersion for a Model
#'
#' This function calculates a measure of overdispersion for a fitted model.
#' Overdispersion is assessed using Pearson's chi-squared statistic, which
#' compares the observed and expected values of residuals. This is particularly
#' useful for evaluating models fitted to count data or other types of data
#' where overdispersion might be a concern.
#'
#' @param model A fitted model object from which to calculate overdispersion.
#' This function is typically used with models such as generalized linear
#' models (`glm`) or other types of models where residuals and degrees of
#' freedom can be computed.
#' @return A named vector containing:
#' \itemize{
#'   \item \code{chisq}: Pearson's chi-squared statistic.
#'   \item \code{ratio}: Ratio of the chi-squared statistic to the residual
#' degrees of freedom.
#'   \item \code{rdf}: Residual degrees of freedom of the model.
#'   \item \code{p}: p-value for the chi-squared statistic.
#' }
#' @details The function computes Pearson's chi-squared statistic from the
#' Pearson residuals of the model and then calculates the ratio of this
#' statistic to the residual degrees of freedom. A large ratio indicates
#' potential overdispersion.
#' @importFrom stats df.residual residuals pchisq
#' @examples
#' # Fit a generalized linear model
#' model <- glm(cbind(successes, failures) ~ predictor,
#'   family = binomial,
#'   data = example_data
#' )
#' # Calculate overdispersion
#' overdisp_fun(model)
#' @export
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(
    chisq = Pearson.chisq,
    ratio = prat,
    rdf = rdf,
    p = pval
  )
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
#' @importFrom ggplot2 aes geom_point
#' @importFrom ggpubr ggarrange ggqqplot gghistogram
#' @importFrom ggsci scale_color_nejm scale_fill_nejm
#' @importFrom gridExtra arrangeGrob
#'
#' @examples
#' model <- lmer(Sepal.Length ~ Sepal.Width + (1 | Species), data = iris)
#' diag.plot.lmer(model)
#'
#' @export
diag.plot.lmer <- function(model) {
  # Get fitted values and residuals
  fitted.vals <- fitted(model)
  resid.vals <- resid(model)

  # Create data frame for plotting
  data.f <- data.frame(fitted.vals, resid.vals)

  # Create scatter plot of fitted values and residuals
  g1 <-
    ggplot(data.f, aes(x = fitted.vals, y = resid.vals)) +
    geom_point() +
    scale_color_nejm() +
    scale_fill_nejm() +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) +
    xlab("Fitted values") +
    ylab("Residuals")

  # Create QQ-plot of residuals
  g2 <-
    ggqqplot(resid(model)) +
    scale_color_nejm() +
    scale_fill_nejm() +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) +
    xlab("Theoretical") +
    ylab("Residuals")

  # Create histogram of residuals
  g3 <-
    gghistogram(resid(model), bins = 10) +
    scale_color_nejm() +
    scale_fill_nejm() +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) +
    xlab("Residuals") +
    ylab("Count")

  # Arrange plots into a grid
  figure <-
    ggarrange(
      g1,
      ggarrange(
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

#' Diagnostic Plot for Linear Mixed-Effects Model
#'
#' This function generates diagnostic plots for a linear mixed-effects
#' model (`lme`).
#' The function creates three plots: a residuals vs fitted values plot,
#' a Q-Q plot of the residuals,
#' and a histogram of the residuals. These plots are arranged in a single
#' figure to visually assess the model's performance.
#'
#' @param model An object of class `lme`, representing a linear
#' mixed-effects model.
#'
#' @return A `ggarrange` object containing three diagnostic plots:
#'   - A: Residuals vs Fitted values plot.
#'   - B: Q-Q plot of the residuals.
#'   - C: Histogram of the residuals.
#'
#' @details
#' The function takes a fitted `lme` model object and uses `ggplot2` and
#' related functions to create three diagnostic plots.
#'
#' @examples
#' \dontrun{
#' library(nlme)
#' library(ggplot2)
#' model <- lme(
#'   fixed = distance ~ age, random = ~ 1 | Subject,
#'   data = Orthodont
#' )
#' diag.plot.lme(model)
#' }
#'
#' @import ggplot2 ggpubr
#' @export
diag.plot.lme <- function(model) {
  fitted.vals <- fitted(model)
  resid.vals <- resid(model, type = "pearson")
  data.f <- data.frame(fitted.vals, resid.vals)
  g1 <-
    ggplot(data.f, aes(x = fitted.vals, y = resid.vals)) +
    geom_point() +
    scale_color_nejm() +
    scale_fill_nejm() +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) +
    xlab("Fitted values") +
    ylab("Residuals")
  g2 <-
    ggqqplot(resid(model)) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab("Theoretical") + ylab("Residuals")
  g3 <-
    gghistogram(resid(model), bins = 10) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab("Residuals") + ylab("Count")
  figure <-
    ggarrange(
      g1,
      ggarrange(
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
}

#' Variance Table for Linear Mixed-Effects Model
#'
#' This function calculates and returns a table of variance components
#' for a linear mixed-effects model (`lmer`).
#' The table includes the fixed, random, and error components, as well as
#' their percentages relative to the total variance
#' and the Intraclass Correlation Coefficient (ICC).
#'
#' @param model An object of class `lmerMod`, representing a linear
#' mixed-effects model, typically fitted using the `lme4` package.
#'
#' @return A data frame with the following columns:
#' \item{Component}{The type of variance component: "Fixed", "Random",
#' or "Error".}
#' \item{Variance}{The amount of variance attributable to each component.}
#' \item{Percentage}{The percentage of total variance explained by each
#' component.}
#' \item{ICC}{The Intraclass Correlation Coefficient (ICC) for the random
#' component, if applicable.}
#'
#' @details
#' The function extracts variance components using the
#' `insight::get_variance()` function and calculates the percentage
#' of variance explained by each component. Additionally, the ICC for the
#' random component is calculated.
#'
#' @examples
#' \dontrun{
#' library(lme4)
#' library(insight)
#' model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' var.table.lmer(model)
#' }
#'
#' @import lme4 insight
#' @export
var.table.lmer <- function(model) {
  var <- as.data.frame(unlist(insight::get_variance(model)))
  var$Variance <- var$`unlist(insight::get_variance(model))`
  var$`unlist(insight::get_variance(model))` <- NULL
  var$Component <- row.names(var)
  row.names(var) <- NULL
  var <- var[c(1, 2, 3), c(2, 1)]
  var$Component <- c("Fixed", "Random", "Error")

  var.total <- var[var$Component == "Fixed", 2] + var[var$Component == "Random", 2] + var[var$Component == "Error", 2]
  var$Percentage <- var$Variance / var.total * 100

  var$Variance <- round(var$Variance, 4)
  var$Percentage <- round(var$Percentage, 1)

  var$ICC <- var[2, "Variance"] / (var[2, "Variance"] + var[3, "Variance"]) * 100
  var[c(2, 3), "ICC"] <- NA

  var$Percentage <- round(var$Percentage, 1)
  var$ICC <- round(var$ICC, 1)

  return(var)
}

#' Create a Factor with Custom Levels and Labels
#'
#' This function transforms a vector into a factor, assigns specified levels
#' and labels, and optionally handles missing values by explicitly assigning
#' them a level.
#'
#' @param Data.vector A vector to be converted into a factor.
#' @param levels A vector of levels to be assigned to the factor.
#' @param labels A vector of labels corresponding to the levels.
#' @param explicit.na Logical, indicating whether missing values (NA) should
#' be explicitly assigned to a "Missing" level. Default is TRUE.
#'
#' @return A factor with specified levels and labels, potentially including
#' an explicit "Missing" level for NA values.
#'
#' @details
#' This function allows for the creation of factors from a vector while
#' specifying both levels and labels.
#' If `explicit.na` is set to `TRUE`, missing values are explicitly
#' assigned to a level named "Missing".
#' The function also drops any unused levels from the factor.
#'
#' @examples
#' \dontrun{
#' vec <- c("low", "medium", "high", NA, "low")
#' levels <- c("low", "medium", "high")
#' labels <- c("Low", "Medium", "High")
#' make.factor(vec, levels, labels)
#' }
#'
#' @import forcats
#' @export
make.factor <- function(Data.vector, levels, labels, explicit.na = TRUE) {
  Result.vector <- as.factor(Data.vector)
  Result.vector <- factor(Result.vector, levels = levels, labels = labels)

  if (explicit.na == TRUE) {
    Result.vector <- fct_na_value_to_level(Result.vector, level = "Missing")
  }
  Result.vector <- droplevels(Result.vector)
}

#' Present Factor Variable as a Bar Plot and Summary Table
#'
#' This function generates a bar plot and a summary table for a factor variable
#' from a given data source. The bar plot displays the counts of each level, and
#' the summary table includes both counts and percentages of each level.
#'
#' @param Datasource A data frame containing the factor variable to be analyzed.
#' @param var.name A string representing the name of the factor variable to be
#' plotted and summarized.
#'
#' @return A list with two elements:
#' \item{graph}{A ggplot object representing the bar plot of the factor
#' variable.}
#' \item{table}{A data frame with two columns: `N` (the count of each factor
#' level) and `Percentage` (the percentage of total observations for each
#' level).}
#'
#' @details
#' This function takes a data frame and a factor variable name, and generates
#' both a bar plot (using `ggplot2`) and a summary table that shows the counts
#' and percentages of each level of the factor.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(dplyr)
#' data(iris)
#' result <- present.factor(iris, "Species")
#' result$graph # Display the bar plot
#' result$table # Display the summary table
#' }
#'
#' @import ggplot2 dplyr ggpubr
#' @export
present.factor <- function(Datasource, var.name) {
  gg <- ggplot(Datasource, aes(x = Datasource[, var.name])) +
    geom_bar() +
    scale_color_nejm() +
    scale_fill_nejm() +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) +
    xlab(var.name) +
    ylab("N")

  total <- length(Datasource[, var.name])
  stats <- as.data.frame(Datasource %>% group_by(.dots = var.name) %>% summarise(N = n(), Percentage = round(n() / total * 100, 1)))

  return(list(graph = gg, table = stats))
}

present.continuous <- function(Datasource, var.name, na.rm = FALSE) {
  gHist <- gghistogram(Datasource[, var.name], bins = 7) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  ) + xlab(var.name) + ylab("N")
  gQQ <- ggqqplot(Datasource[, var.name]) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )

  gBP <- ggboxplot(Datasource[, var.name]) + coord_flip() + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  ) + xlab("") + ylab(var.name)

  gg <- ggarrange(ggarrange(gHist, gQQ, nrow = 1, ncol = 2), gBP, nrow = 2, ncol = 1)

  stats <- data.frame(list(Mean = mean(Datasource[, var.name], na.rm = na.rm), Median = median(Datasource[, var.name], na.rm = na.rm), Min = min(Datasource[, var.name], na.rm = na.rm), Max = max(Datasource[, var.name], na.rm = na.rm), SD = sd(Datasource[, var.name], na.rm = na.rm), IQR = IQR(Datasource[, var.name], na.rm = na.rm)))
  return(list(graph = gg, table = stats))
}

Univariate.summary.table <- function(Datasource, exclude.vars, na.rm = FALSE, round.to = 1) {
  stats.table <- data.frame(list(Variable = NA, Category = NA, N = NA, Percentage = NA, Mean = NA, Median = NA, Min = NA, Max = NA, SD = NA, IQR = NA))
  for (var.name in names(Datasource)) {
    if (sum(exclude.vars == var.name) > 0) {
      next
    }

    if (is.factor(Datasource[, var.name])) {
      total <- length(Datasource[, var.name])
      stats <- as.data.frame(Datasource %>% group_by(base::get(var.name)) %>% summarise(N = n(), Percentage = round(n() / total * 100, 1)))
      stats$Mean <- NA
      stats$Median <- NA
      stats$Min <- NA
      stats$Max <- NA
      stats$SD <- NA
      stats$IQR <- NA
      stats$Variable <- var.name
      stats <- stats[, c(10, seq(1, 9))]
      names(stats)[2] <- "Category"
      stats[seq(2, length(stats[, 1])), 1] <- NA
      stats.table <- rbind(stats.table, stats)
    } else {
      stats <- data.frame(list(N = NA, Percentage = NA, Mean = mean(Datasource[, var.name], na.rm = na.rm), Median = median(Datasource[, var.name], na.rm = na.rm), Min = min(Datasource[, var.name], na.rm = na.rm), Max = max(Datasource[, var.name], na.rm = na.rm), SD = sd(Datasource[, var.name], na.rm = na.rm), IQR = IQR(Datasource[, var.name], na.rm = na.rm)))
      stats$Variable <- var.name
      stats$Category <- NA
      stats <- stats[, c(9, 10, seq(1, 8))]
      stats.table <- rbind(stats.table, stats)
    }
  }

  stats.table[, seq(3, 10)] <- round(stats.table[, seq(3, 10)], round.to)
  stats.table <- stats.table[-c(1), ]
  rownames(stats.table) <- NULL

  return(stats.table)
}

Plot.factor.vs.factor <- function(Datasource, factor.name.x, factor.name.y) {
  gg <- ggplot(Datasource, aes(x = Datasource[, factor.name.x], y = Datasource[, factor.name.y])) +
    geom_jitter(width = 0.2, height = 0.2) +
    scale_color_nejm() +
    scale_fill_nejm() +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) +
    xlab(factor.name.x) +
    ylab(factor.name.y)

  return(gg)
}

#' Format P-Value with Standardized Precision
#'
#' This function formats a p-value to a standardized number of significant digits.
#' If the p-value is smaller than the threshold (10^(-sigD)), it is returned as
#' "<0.001" (or another appropriate threshold based on `sigD`), otherwise, it is
#' rounded to the specified number of significant digits.
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
#' \dontrun{
#' standard.p.value(0.0005) # "<0.001"
#' standard.p.value(0.12345) # "0.123"
#' standard.p.value(0.0001, 4) # "<0.0001"
#' }
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

standard.icc.interpret <- function(r) {
  interpret <- if (r < 0.5) {
    "Poor"
  } else {
    if (r < 0.75) {
      "Moderate"
    } else {
      if (r < 0.9) {
        "Good"
      } else {
        "Excellent"
      }
    }
  }

  return(interpret)
}

standard.convergence.interpret <- function(r) {
  interpret <- if (r < 0.3) {
    "Poor"
  } else {
    if (r < 0.6) {
      "Adequate"
    } else {
      "Excellent"
    }
  }

  return(interpret)
}

standard.cor.test <- function(x, y, standard.interpret = T, interpret.func = standard.icc.interpret) {
  result <- cor.test(x, y)
  interpret.str <- if (standard.interpret == T) {
    standard.icc.interpret(result$conf.int[1])
  } else {
    standard.icc.interpret(result$estimate)
  }
  result.str <- paste0(
    "r=", round(result$estimate, 2),
    " 95% CI [", round(result$conf.int[1], 2), ", ", round(result$conf.int[2], 2), "]",
    " t[", result$parameter, "]=", round(result$statistic[[1]], 2),
    ", p=", standard.p.value(result$p.value, sigD = 3)
  )

  return(paste0(result.str, ", ", interpret.str))
}


hedges.g <- function(n1, n2, cohens.d) {
  Vd <- ((n1 + n2) / (n1 * n2)) + ((cohens.d^2) / (2 * (n1 + n2)))

  df <- n1 + n2 - 2

  J <- 1 - (3 / (4 * df - 1))

  ges <- cohens.d * J

  ses <- sqrt((J^2) * Vd)

  return(list("g" = ges, "g.se" = ses))
}

example.rma <- function() {
  n1 <- c(
    85,
    50,
    29,
    75,
    219,
    18
  )

  n2 <- c(
    79,
    51,
    28,
    78,
    219,
    18
  )

  d <- c(
    -0.032673065,
    0.21942446,
    0.294432707,
    0.194900425,
    0.274426269,
    -0.142434924
  )

  calculatedEffs <- hedges.g(n1, n2, d)

  ges <- calculatedEffs$ges
  ses <- calculatedEffs$ses

  slabs <- c(
    "Bossen-2013",
    "Chapman-2018",
    "Lee-2014",
    "Maddison-2015",
    "Wong-2020",
    "Nasseri-2020"
  )

  mdl <- rma(ges, sei = ses, slab = slabs, method = "FE")

  gg <- forest(mdl, xlab = "Std. Mean Difference", header = "Author-Year", mlab = "Fixed Effects Model for All Studies\nHeterogeneity Test: Q[5]=3.99, p=0.55")

  return(list("mdl" = mdl, "gg" = gg))
}

Append.within.changes.exp <- function(Input.df, Pre.post.correlation = 0.5) {
  Output.df <- Input.df

  # Pre-post means and SDs
  Output.df$Exp.post.minus.pre.mean.calc <- Output.df$Exp.post.mean - Output.df$Exp.pre.mean
  Output.df$Exp.post.minus.pre.sd.calc <- sqrt((Output.df$Exp.post.sd)^2 +
    (Output.df$Exp.pre.sd)^2 -
    2 * Pre.post.correlation * Output.df$Exp.post.sd * Output.df$Exp.pre.sd)

  # Pre-post means and SDs selection
  Output.df$Exp.post.minus.pre.mean.sel <- Output.df$Exp.post.minus.pre.mean
  Output.df$Exp.post.minus.pre.mean.sel[which(is.na(Output.df$Exp.post.minus.pre.mean.sel))] <- Output.df$Exp.post.minus.pre.mean.calc[which(is.na(Output.df$Exp.post.minus.pre.mean.sel))]

  Output.df$Exp.post.minus.pre.sd.sel <- Output.df$Exp.post.minus.pre.sd
  Output.df$Exp.post.minus.pre.sd.sel[which(is.na(Output.df$Exp.post.minus.pre.sd.sel))] <- Output.df$Exp.post.minus.pre.sd.calc[which(is.na(Output.df$Exp.post.minus.pre.sd.sel))]

  Output.df[Output.df$Increase...Improvement == "No", "Exp.post.minus.pre.mean.sel"] <- Output.df[Output.df$Increase...Improvement == "No", "Exp.post.minus.pre.mean.sel"] * -1.0

  return(Output.df)
}

Append.within.changes.ctrl <- function(Input.df, Pre.post.correlation = 0.5) {
  Output.df <- Input.df
  # Pre-post means and SDs
  Output.df$Ctrl.post.minus.pre.mean.calc <- Output.df$Ctrl.post.mean - Output.df$Ctrl.pre.mean
  Output.df$Ctrl.post.minus.pre.sd.calc <- sqrt((Output.df$Ctrl.post.sd)^2 +
    (Output.df$Ctrl.pre.sd)^2 -
    2 * Pre.post.correlation * Output.df$Ctrl.post.sd * Output.df$Ctrl.pre.sd)

  # Pre-post means and SDs selection
  Output.df$Ctrl.post.minus.pre.mean.sel <- Output.df$Ctrl.post.minus.pre.mean
  Output.df$Ctrl.post.minus.pre.mean.sel[which(is.na(Output.df$Ctrl.post.minus.pre.mean.sel))] <- Output.df$Ctrl.post.minus.pre.mean.calc[which(is.na(Output.df$Ctrl.post.minus.pre.mean.sel))]

  Output.df$Ctrl.post.minus.pre.sd.sel <- Output.df$Ctrl.post.minus.pre.sd
  Output.df$Ctrl.post.minus.pre.sd.sel[which(is.na(Output.df$Ctrl.post.minus.pre.sd.sel))] <- Output.df$Ctrl.post.minus.pre.sd.calc[which(is.na(Output.df$Ctrl.post.minus.pre.sd.sel))]

  Output.df[Output.df$Increase...Improvement == "No", "Ctrl.post.minus.pre.mean.sel"] <- Output.df[Output.df$Increase...Improvement == "No", "Ctrl.post.minus.pre.mean.sel"] * -1.0

  return(Output.df)
}

Append.between.diffs <- function(Input.df) {
  Output.df <- Input.df

  n1 <- Output.df$Exp.post.n
  n2 <- Output.df$Ctrl.post.n

  S1 <- Output.df$Exp.post.minus.pre.sd.sel
  S2 <- Output.df$Ctrl.post.minus.pre.sd.sel

  Swithin <- sqrt(((n1 - 1) * S1^2 + (n2 - 1) * S2^2) / (n1 + n2 - 2))
  d <- (Output.df$Exp.post.minus.pre.mean.sel - Output.df$Ctrl.post.minus.pre.mean.sel) / Swithin

  g <- hedges.g(n1, n2, d)

  Output.df$Hedges.g <- g$g
  Output.df$Hedges.g.se <- g$g.se

  return(Output.df)
}

standard.ctt.table <- function(raw.ctt.table) {
  if ("asymp.LCL" %in% names(data.frame(raw.ctt.table)) ||
    "asymp.LCL" %in% names(data.frame(summary(raw.ctt.table, infer = c(T, T))))) {
    ci.names <- c("asymp.LCL", "asymp.UCL")
    stat.names <- c("z.ratio", "z")
  } else {
    ci.names <- c("lower.CL", "upper.CL")
    stat.names <- c("t.ratio", "t")
  }

  aT <- data.frame(summary(raw.ctt.table, infer = c(T, T))) %>%
    rowwise() %>%
    mutate(sigD = ceiling(log10(1 / SE))) %>%
    rowwise() %>%
    mutate(Difference = paste0(
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
      round(base::get(stat.names[1]), 3), ", ", standard.p.value(p.value)
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

standard.em.table <- function(raw.em.table) {
  if ("asymp.LCL" %in% names(data.frame(raw.em.table)) ||
    "asymp.LCL" %in% names(data.frame(summary(raw.em.table, infer = c(T, T))))) {
    ci.names <- c("asymp.LCL", "asymp.UCL")
    stat.names <- c("z.ratio", "z")
  } else {
    ci.names <- c("lower.CL", "upper.CL")
    stat.names <- c("t.ratio", "t")
  }

  aT <- data.frame(raw.em.table) %>%
    rowwise() %>%
    mutate(sigD = ceiling(log10(1 / SE))) %>%
    rowwise() %>%
    mutate(
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
        round(base::get(stat.names[1]), 3), ", ", standard.p.value(p.value)
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
