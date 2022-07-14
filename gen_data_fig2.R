# generate the data for figure 2

library(DBpower)
library(dplyr)
library(tidyr)
library(magrittr)
library(GBJ)
library(mvtnorm)
library(bindata)

# parameters assumed in paper
betaVal <- 0.25
alpha <- 0.01
betaVec <- c(betaVal, 0, 0, 0, 0)
MAF <- 0.3
n <- 400
sigSqY <- 1
Wg <- diag(rep(sqrt(2 * MAF * (1 - MAF)), 5))
xMat <- matrix(data=1, nrow=n, ncol=1)

# high correlation setting
# make the correlation matrix
corMat5 <- matrix(data=NA, nrow=5, ncol=5)
corMat5[1:2, 1:2] <- 0.7
corMat5[3:5, 3:5] <- 0.7
corMat5[1:2, 3:5] <- 0.5
corMat5[3:5, 1:2] <- 0.5
diag(corMat5) <- 1

# bounds
bjBounds5 <- set_BJ_bounds(alpha = 0.01, J=5)
gbjBounds5 <- set_GBJ_bounds(alpha = 0.01, J=5, sig_vec = corMat5[lower.tri(corMat5)])

# low corelation
# make the correlation matrix
corMat0 <- matrix(data=NA, nrow=5, ncol=5)
corMat0[1:2, 1:2] <- 0.3
corMat0[3:5, 3:5] <- 0.3
corMat0[1:2, 3:5] <- 0.1
corMat0[3:5, 1:2] <- 0.1
diag(corMat0) <- 1

# bounds
bjBounds0 <- bjBounds5
gbjBounds0 <- set_GBJ_bounds(alpha = 0.01, J=5, sig_vec = corMat0[lower.tri(corMat0)])


# calculate the bounds and power, simulation to check
set.seed(0)
allPowerCalc <- c()
allPowerSim <- c()
for (corType in c("low", "high")) {

  # high and low correlation
  if (corType == "high") {
    corMat <- corMat5
    bjBound <- bjBounds5
    gbjBound <- gbjBounds5
  } else {
    corMat <- corMat0
    bjBound <- bjBounds0
    gbjBound <- gbjBounds0
  }
  # eigendecomposition
  eVals <- eigen(corMat)$values
  eVecs <- eigen(corMat)$vectors

  # switch between SS and MP
  for (multType in c("SS", "MP")) {
    tempRes <- data.frame(Cor = corType, Mult = multType, location = 1:5, iBJlower=NA, iBJlowerUnion=NA, iBJupper = NA,
                          GBJlower = NA, GBJupper = NA, iBJpower=NA, GBJpower=NA)
    tempSim <- data.frame(Cor = corType, Mult = multType, location = 1:5, iBJlower=NA, iBJlowerUnion = NA, iBJupper = NA,
                          GBJlower = NA, GBJupper = NA, iBJpower=NA, GBJpower=NA)
  
    # move through signal locations
    for (signal_it in 1:5) {
      # move signal
      tempBeta <- rep(0, 5)
      tempBeta[signal_it] <- betaVal

      # for simulating test statistics
      gMat <- bindata::rmvbin(n=n, margprob = rep(MAF, 5), bincorr = corMat) +
        bindata::rmvbin(n=n, margprob = rep(MAF, 5), bincorr = corMat)
      corMatEmp <- cor(gMat)
      decompEmp <- eigen(corMatEmp)

      # mean vector
      if (multType == "MP") {
        genMean <- sqrt(2 * MAF * (1 - MAF) * n) * tempBeta
        innMean <- sqrt(2 * MAF * (1 - MAF) * n) * diag(1 / sqrt(eVals)) %*% t(eVecs) %*% tempBeta
        innMeanEmp <- innMean
      } else if (multType == "SS") {
        genMean <- sqrt(n / sigSqY) * corMat %*% Wg %*% tempBeta
        innMean <- sqrt(n / sigSqY) * diag(sqrt(eVals)) %*% t(eVecs) %*% Wg %*% tempBeta
        innMeanEmp <- sqrt(n / sigSqY) * diag(sqrt(decompEmp$values)) %*% t(decompEmp$vectors) %*% Wg %*% tempBeta
      }
      # make it numeric
      innMean <- as.numeric(innMean)
      genMean <- as.numeric(genMean)

      # calculate bounds
      innBounds <- calcb2(lower = TRUE, upper = TRUE, muVec = innMean, sigMat = diag(rep(1, 5)), bounds = bjBound)
      genBounds <- calcb2(lower = TRUE, upper = TRUE, muVec = genMean, sigMat = corMat, bounds = gbjBound)

      # calculate exact power
      innPower <- calc_exact_power(bounds = bjBound, sig_mat = diag(rep(1, 5)), muVec = innMean)$power
      genPower <- calc_exact_power(bounds = gbjBound, sig_mat = corMat, muVec = genMean)$power

      # record calculated power
      tempRes$iBJlower[signal_it] <- innBounds$lowerProb
      tempRes$iBJupper[signal_it] <- innBounds$upperProb
      tempRes$iBJpower[signal_it] <- innPower
      tempRes$GBJlower[signal_it] <- genBounds$lowerProb
      tempRes$GBJupper[signal_it] <- genBounds$upperProb
      tempRes$GBJpower[signal_it]  <- genPower

      # simulate bounds and exact power
      simBoundsGen <- sim_b2(lower=TRUE, upper=TRUE, n = 30000, muVec = genMean, sigMat = corMat, bounds = gbjBound)
      simBoundsInn <- sim_b2(lower=TRUE, upper=TRUE, n = 30000, muVec = innMean, sigMat = diag(rep(1, 5)), bounds = bjBound)
      simPowerGen <- sim_power_mvn(n = 30000, muVec = genMean, sigMat = corMat, bounds=gbjBound, test=NULL, alpha = alpha)
      simPowerInn <- sim_power_mvn(n = 30000, muVec = innMean, sigMat = diag(rep(1, 5)), bounds=bjBound, test=NULL, alpha = alpha)
    
      # record simulation power
      tempSim$iBJlower[signal_it] <- simBoundsInn$lowerBound
      tempSim$iBJupper[signal_it] <- simBoundsInn$upperBound
      tempSim$iBJpower[signal_it] <- simPowerInn$boundsPower
      tempSim$GBJlower[signal_it] <- simBoundsGen$lowerBound
      tempSim$GBJupper[signal_it] <- simBoundsGen$upperBound
      tempSim$GBJpower[signal_it]  <- simPowerGen$boundsPower

      # checkpoint
      cat(signal_it, '\n')
    } # end loop through signal location

    # append results
    allPowerCalc <- rbind(allPowerCalc, tempRes)
    allPowerSim <- rbind(allPowerSim, tempSim)

    # checkpoint
    cat(multType, '\n')
  } # end loop through multType
} # end loop through corType

# save
write.table(allPowerCalc, "allPowerCalc.txt", append=F, quote=F, row.names=F, col.names=T)
write.table(allPowerSim, "allPowerSim.txt", append=F, quote=F, row.names=F, col.names=T)


