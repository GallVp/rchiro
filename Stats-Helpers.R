# function for is nan
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
  
lmer.std.beta <- function(mod) {
   b <- fixef(mod)[-1]
   sd.x <- apply(getME(mod,"X")[,-1],2,sd)
   sd.y <- sd(getME(mod,"y"))
   b*sd.x/sd.y
}

confint.rlmerMod <- function(object, parm, level = 0.95) {
  beta <- fixef(object)
  if (missing(parm))
    parm <- names(beta)
  se <- sqrt(diag(vcov(object)))
  z <- qnorm((1 + level) / 2)
  ctab <- cbind(beta - z * se, beta + z * se)
  colnames(ctab) <- stats:::format.perc(c((1 - level) / 2, (1 + level) /
                                            2),
                                        digits = 3)
  return(ctab[parm, ])
}

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
      family$mu.eta(linear.predictors) ^ 2 / family$variance(fitted.values)
    )))
  mod2$residuals <-
    with(mod2, prior.weights * (y - fitted.values) / weights)
  class(mod2) <- c("fakeglm", class(mod2))
  mod2
}


overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp ^ 2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(
    chisq = Pearson.chisq,
    ratio = prat,
    rdf = rdf,
    p = pval
  )
}

diag.plot.lmer <- function(model) {
  fitted.vals <- fitted(model)
  resid.vals <- resid(model)
  data.f <- data.frame(fitted.vals, resid.vals)
  g1 <-
    ggplot(data.f, aes(x = fitted.vals, y = resid.vals)) + geom_point() + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Fitted values') + ylab('Residuals')
  g2 <-
    ggqqplot(resid(model)) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Theoretical') + ylab('Residuals')
  g3 <-
    gghistogram(resid(model), bins = 10) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Residuals') + ylab('Count')
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

diag.plot.lme <- function(model) {
  fitted.vals <- fitted(model)
  resid.vals <- resid(model, type = "pearson")
  data.f <- data.frame(fitted.vals, resid.vals)
  g1 <-
    ggplot(data.f, aes(x = fitted.vals, y = resid.vals)) + geom_point() + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Fitted values') + ylab('Residuals')
  g2 <-
    ggqqplot(resid(model)) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Theoretical') + ylab('Residuals')
  g3 <-
    gghistogram(resid(model), bins = 10) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Residuals') + ylab('Count')
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

var.table.lmer <- function(model) {
  var <- as.data.frame(unlist(insight::get_variance(model)))
  var$Variance <-var$`unlist(insight::get_variance(model))`
  var$`unlist(insight::get_variance(model))` <- NULL
  var$Component <- row.names(var)
  row.names(var) <- NULL
  var <- var[c(1,2,3), c(2, 1)]
  var$Component <- c("Fixed", "Random", "Error")
  
  var.total <- var[var$Component == "Fixed", 2] + var[var$Component == "Random", 2] + var[var$Component == "Error", 2]
  var$Percentage <- var$Variance / var.total * 100
  
  var$Variance <- round(var$Variance, 4)
  var$Percentage <- round(var$Percentage, 1)
  
  var$ICC <- var[2, "Variance"] / (var[2, "Variance"] + var[3, "Variance"])*100
  var[c(2, 3), 'ICC'] <- NA
  
  var$Percentage <- round(var$Percentage, 1)
  var$ICC <- round(var$ICC, 1)
  
  return (var)
}

make.factor <- function(Data.vector, levels, labels) {
  Result.vector <- as.factor(Data.vector)
  Result.vector <- factor(Result.vector, levels = levels, labels = labels)
  Result.vector <- fct_explicit_na(Result.vector, na_level = "Missing")
}

present.factor <- function(Datasource, var.name) {
  gg <- ggplot(Datasource, aes(x = Datasource[, var.name])) + geom_bar() + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none") + xlab(var.name) + ylab('N')
  
  total <- length(Datasource[, var.name])
  stats <- as.data.frame(Datasource %>% group_by(.dots = var.name) %>% summarise(N = n(), Percentage = round(n()/total*100, 1)))
  
  return(list(graph = gg, table = stats))
}

present.continuous <- function(Datasource, var.name, na.rm = FALSE) {
  gHist <- gghistogram(Datasource[, var.name], bins = 7) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none") + xlab(var.name) + ylab('N')
  gQQ <- ggqqplot(Datasource[, var.name]) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none")
  
  gBP <- ggboxplot(Datasource[, var.name]) + coord_flip() + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none") + xlab('') + ylab(var.name)
  
  gg <- ggarrange(ggarrange(gHist, gQQ, nrow = 1, ncol = 2), gBP, nrow = 2, ncol = 1)
  
  stats <- data.frame(list(Mean = mean(Datasource[, var.name], na.rm = na.rm), Median = median(Datasource[, var.name], na.rm = na.rm), Min = min(Datasource[, var.name], na.rm = na.rm), Max = max(Datasource[, var.name], na.rm = na.rm), SD = sd(Datasource[, var.name], na.rm = na.rm), IQR = IQR(Datasource[, var.name], na.rm = na.rm)))
  return(list(graph = gg, table = stats))
}

Univariate.summary.table <- function(Datasource, exclude.vars, na.rm = FALSE, round.to = 1) {
  stats.table <- data.frame(list(Variable = NA, Category = NA, N = NA, Percentage = NA, Mean = NA, Median = NA, Min = NA, Max = NA, SD = NA, IQR = NA))
  for (var.name in names(Datasource)) {
    
    if(sum(exclude.vars == var.name) > 0) {
      next
    }
    
    if(is.factor(Datasource[, var.name])) {
      total <- length(Datasource[, var.name])
      stats <- as.data.frame(Datasource %>% group_by(.dots = var.name) %>% summarise(N = n(), Percentage = round(n()/total*100, 1)))
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
  gg<- ggplot(Datasource, aes(x = Datasource[, factor.name.x], y = Datasource[, factor.name.y])) + geom_jitter(width = 0.2, height = 0.2) +
    scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none") + xlab(factor.name.x) + ylab(factor.name.y)
  
  return(gg)
}

standard.p.value <- function(value, sigD = 3) {
  text <- if(value < 10^(-sigD)) {paste0("<", 10^(-sigD))} else {round(value, sigD)}
  text
}

standard.em.table <- function(raw.em.table) {
  aT <- as.data.frame(raw.em.table) %>%
    rowwise %>%
    mutate(sigD = ceiling(log10(1/SE))) %>%
    rowwise %>%
    mutate(Estimate = paste0(round(emmean, sigD),
                             '±',
                             round(SE, sigD),
                             ' [',
                             round(lower.CL, sigD),
                             ', ',
                             round(upper.CL, sigD), ']'))
  aT$emmean <- NULL
  aT$SE <- NULL
  aT$df <- NULL
  aT$lower.CL <- NULL
  aT$upper.CL <- NULL
  aT$sigD <- NULL
  aT <- as.data.frame(aT)

  aT.names <- names(aT)
  aT.names[length(aT.names)] <- "Estimate±SE [95% CI]"
  names(aT) <- aT.names

  aT
}

standard.ctt.table <- function(raw.ctt.table) {
  aT <- as.data.frame(summary(raw.ctt.table, infer = c(T, T))) %>%
    rowwise %>%
    mutate(sigD = ceiling(log10(1/SE))) %>%
    rowwise %>%
    mutate(Difference = paste0(round(estimate, sigD),
                             '±',
                             round(SE, sigD),
                             ' [',
                             round(lower.CL, sigD),
                             ', ',
                             round(upper.CL, sigD), ']'), Test = paste0('t[', round(df, 1), ']', '=', round(t.ratio, 3), ', ', standard.p.value(p.value)))
  aT$estimate <- NULL
  aT$SE <- NULL
  aT$df <- NULL
  aT$lower.CL <- NULL
  aT$upper.CL <- NULL
  aT$sigD <- NULL
  aT$t.ratio <- NULL
  aT$p.value <- NULL
  aT <- as.data.frame(aT)

  aT.names <- names(aT)
  aT.names[length(aT.names)-1] <- "Difference±SE [95% CI]"
  aT.names[length(aT.names)] <- "t[df], p-value"
  aT.names[1] <- "Contrast"
  names(aT) <- aT.names

  aT
}

standard.icc.interpret <- function(r) {
  interpret <- if(r < 0.5) { "Poor" } else{ if(r < 0.75) { "Moderate" }
    else { if(r < 0.9) { "Good" } else { "Excellent" } } }

  return(interpret)
}

standard.convergence.interpret <- function(r) {
  interpret <- if(r < 0.3) { "Poor" } else{ if(r < 0.6) { "Adequate" }
    else { "Excellent" } }

  return(interpret)
}

standard.cor.test <- function(x, y, standard.interpret = T, interpret.func = standard.icc.interpret) {
  result <- cor.test(x, y)
  interpret.str <- if(standard.interpret == T) {standard.icc.interpret(result$conf.int[1])} else {standard.icc.interpret(result$estimate)}
  result.str <- paste0('r=', round(result$estimate, 2),
  ' 95% CI [', round(result$conf.int[1], 2), ', ', round(result$conf.int[2], 2), ']',
  ' t[', result$parameter, ']=', round(result$statistic[[1]], 2),
  ', p=', standard.p.value(result$p.value, sigD = 3))

  return(paste0(result.str, ', ', interpret.str))
}


hedges.g <- function(n1, n2, cohens.d) {
  Vd <- ((n1+n2)/(n1*n2)) + ((cohens.d^2)/(2*(n1+n2)))
  
  df <- n1 + n2 - 2
  
  J <- 1-(3/(4*df-1))
  
  ges <- cohens.d*J
  
  ses <- sqrt((J^2)*Vd)
  
  return(list("g" = ges, "g.se" = ses))
}

example.rma <- function() {
  n1 <- c(85,
          50,
          29,
          75,
          219,
          18)

  n2 <- c(79,
          51,
          28,
          78,
          219,
          18)

  d <- c(-0.032673065,
        0.21942446,
        0.294432707,
        0.194900425,
        0.274426269,
        -0.142434924)

  calculatedEffs <- hedges.g(n1, n2, d)

  ges <- calculatedEffs$ges
  ses <- calculatedEffs$ses

  slabs <- c('Bossen-2013',
            'Chapman-2018',
            'Lee-2014',
            'Maddison-2015',
            'Wong-2020',
            'Nasseri-2020')

  mdl<-rma(ges, sei = ses, slab = slabs, method = 'FE')

  gg <- forest(mdl, xlab = 'Std. Mean Difference', header='Author-Year', mlab = "Fixed Effects Model for All Studies\nHeterogeneity Test: Q[5]=3.99, p=0.55")

  return(list("mdl" = mdl, "gg" = gg))
}

Append.within.changes.exp <- function(Input.df, Pre.post.correlation=0.5) {
  
  Output.df <- Input.df
  
  # Pre-post means and SDs
  Output.df$Exp.post.minus.pre.mean.calc <- Output.df$Exp.post.mean -  Output.df$Exp.pre.mean
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

Append.within.changes.ctrl <- function(Input.df, Pre.post.correlation=0.5) {
  
  Output.df <- Input.df
  # Pre-post means and SDs
  Output.df$Ctrl.post.minus.pre.mean.calc <- Output.df$Ctrl.post.mean -  Output.df$Ctrl.pre.mean
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
  
  Swithin <- sqrt(((n1 - 1)*S1^2 + (n2 - 1)*S2^2) / (n1 + n2 - 2))
  d <- (Output.df$Exp.post.minus.pre.mean.sel - Output.df$Ctrl.post.minus.pre.mean.sel) / Swithin

  g <- hedges.g(n1, n2, d)
  
  Output.df$Hedges.g <- g$g
  Output.df$Hedges.g.se <- g$g.se
  
  return(Output.df)
}