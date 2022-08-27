Power.cross.over.2x2  <- function(simNum, N, eSize, baselineMean, baselineSD , timeCor = 0.7, periodCor = 0.5) {
  
  Sigma               <- matrix(c(1,          timeCor,    periodCor,  periodCor,
                                  timeCor,    1,          periodCor,  periodCor,
                                  periodCor,  periodCor,  1,          timeCor,
                                  periodCor,  periodCor,  timeCor,    1),
                                ncol = 4, nrow = 4)
  
  simData             <- mvrnorm(n=N, mu = c(0, 0, 0, 0), Sigma = Sigma)
  
  simData[, 1]        <- simData[, 1]*baselineSD + baselineMean
  simData[, 2]        <- simData[, 2]*baselineSD + baselineMean
  simData[, 3]        <- simData[, 3]*baselineSD + baselineMean
  simData[, 4]        <- simData[, 4]*baselineSD + baselineMean + eSize
  
  
  dataSource          <- data.frame(factor(c(seq(1, N), seq(1, N))),
                                    factor(c(rep(1, N), rep(2, N))),
                                    c(simData[,1], simData[,3]),
                                    c(simData[,2], simData[,4]))
  names(dataSource)   <- c("PartID", "Condition", "OutcomePre", "Outcome")
  
  tryCatch({
    lmModel           <- lmer(Outcome ~ OutcomePre + Condition + (1|PartID), dataSource)
    
    tp                <- emmeans(lmModel, pairwise~Condition)
    cp                <- as.data.frame(tp$contrasts)
    
    p                 <- cp$p.value
    stat              <- cp$estimate
  }, error = function(error_condition) {
    print('Error occured')
    stat              <- 0
    p                 <- 1
  })
  
  return(c(t=stat, p=p, sig=(p < .05)))
}