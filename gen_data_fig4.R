# Generate the data for Figure 4.
# This script is designed to be run on a computing cluster.
# We need to run 59 + 90 = 149 jobs with two input arguments.
# First, we need to run it with leading argument (aID) 1-59 and second argument (Snum) 1.
# Then we need to run it with leading argument (aID) 1-90 and second argument (Snum) 3.


library(bindata)
library(dplyr)
library(magrittr)
library(mvtnorm)
library(GBJ)
library(CompQuadForm)
library(MSKAT)
library(SKAT)
library(ACAT)

args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])


# Change output directory to reflect file structure
outputDir <- "/rsrch3/home/biostatistics/rsun3/finiteSamplePower/simulation/output/"
outName <- paste0("sim_rho50_mp40_ss25_S", Snum, "_aID", aID, ".txt")

# parameters that stay fixed
D <- 40
n <- 400
# num causal
r1 <- 4
r2 <- D - r1
mafVec <- rep(0.3, D)
B <- 500
gammaMinSS <- 0
gammaMaxSS <- 0.25
gammaMinMP <- 0
gammaMaxMP <- 0.4

# the different correlation structures and effect sizes
if (Snum == 1) {
  rho1 <- 0.5
  rho2 <- 0.5
  rho3 <- seq(from=0, to=0.6, by=0.01)[aID]
}  else if (Snum == 3) {
  rho1 <- seq(from=0, to=0.95, by=0.01)[aID]
  rho2 <- rho1
  rho3 <- rho2
  oppCor <- matrix(data=0.4, nrow=D, ncol=D)
  diag(oppCor) <- 1
}

cat("gammaMinSS ", gammaMinSS, "\n")
cat("gammaMaxSS ", gammaMaxSS, "\n")
cat("gammaMinMP ", gammaMinMP, "\n")
cat("gammaMaxMP ", gammaMaxMP, "\n")
cat("rho1 ", rho1, "\n")
cat("rho2 ", rho2, "\n")
cat("rho3 ", rho3, "\n")

# when doing a MP simulation, we need to keep the genotype 
# correlation matrix the same when varying the outcome one, and
# vice versa for SS simulation. this is oppCor, same for all scenarios.
oppCor <- matrix(data=0.4, nrow=D, ncol=D)
diag(oppCor) <- 1

# non-SNP effect sizes are the same for each outcome
betaMat <- rbind(0, matrix(data=1, nrow=2, ncol=D))
gammaMat <- matrix(data=0, nrow=D, ncol=D)

# build correlation matrix
corMat <- matrix(data=rho3, nrow=D, ncol=D)
corMat[1:r1, 1:r1] <- rho1
corMat[(r1+1):D, (r1+1):D] <- rho2
diag(corMat) <- 1

# eigendecomposition
eVals <- eigen(corMat)$values
eVecs <- eigen(corMat)$vectors
if (min(eVals) <= 10^(-10)) {stop("Correlation matrix not positive definite")}
# this matrix multiplied by a p*1 vector of correlated statistics creates
# a p*1 vector of independent test statistics
invSqrtCor <- diag(1 / sqrt(eVals)) %*% t(eVecs)
# this is the inverse of corMat
corMatInv <- eVecs %*%  diag(1 / eVals) %*% t(eVecs)
  
# for simulating genotypes
cprob <- bincorr2commonprob(margprob = mafVec, bincorr = corMat)
rmvbinSigma <- commonprob2sigma(commonprob = cprob)

# function to get the individual test statistics
# we want the first outcome against all the SNPs, and we want
# the first SNP against all the outcomes
getTestStatistics <- function(genoMat, yMatMP, xMat, yVecSS, gVecMP) {

  n <- nrow(yMatMP)
  
  # first outcome against all the SNPs
  projX <- xMat %*% solve(t(xMat) %*% xMat) %*% t(xMat)
  fittedYVecSS <- projX %*% yVecSS
  residVecSS <- yVecSS - fittedYVecSS 
  sigmaSqYHatSS <- sum(residVecSS^2) / (n - ncol(xMat))
  # score statistics - these need to be divided by (varY * denominator in GBJ paper)
  snpScoreStats <- t(genoMat) %*% residVecSS
  # covariance matrix - to get denominator in GBJ paper
  snpCovMat <- t(genoMat) %*% genoMat - t(genoMat) %*% projX %*% genoMat
  snpCovVec <- diag(snpCovMat)
  # correlation matrix - use this to innovate
  # it matches the cor_mat from GBJ::calc_cor_stats
  snpCorMat <- sweep(snpCovMat, MARGIN=2, STATS = sqrt(snpCovVec), FUN="/") %>%
    sweep(., MARGIN=1, STATS = sqrt(snpCovVec), FUN="/")
  # make it symmetric
  snpCorMat[lower.tri(snpCorMat)] = t(snpCorMat)[lower.tri(snpCorMat)]
  # eigenvalues for SKAT
  eValsSS <- eigen(snpCovMat, symmetric = TRUE)$values
  # standardized score statistics - marginally N(0,1) under the null, also matches calc_cor_stats
  snpStandStats <- snpScoreStats / sqrt(snpCovVec * sigmaSqYHatSS)

  # first SNP against all the outcomes
  fittedYMat <- projX %*% yMatMP
  residMat <- yMatMP - fittedYMat 
  sigmaSqYHat <- apply(residMat^2, 2, sum) / (n - ncol(xMat))
  outcomeScoreStats <- as.numeric(t(gVecMP) %*% residMat)
  # can check against function checkMP()
  snpVarMP <- as.numeric(t(gVecMP) %*% gVecMP - t(gVecMP) %*% projX %*% gVecMP)
  outcomeStandStats <- outcomeScoreStats / (sqrt(snpVarMP * sigmaSqYHat))
  # eigenvalues for MP SKAT
  #eValsMP <- eigen(snpCovMat[1, 1] * corMatInv, symmetric = TRUE)$values
  # covariance matrix of Y, should I use residuals?
  outcomeCovMat <- cov(residMat) * ((n-1) / (n - ncol(xMat))) * snpVarMP
  # this solve is for the MP SKAT, where we don't need the scaling term
  yCovMatSolve <- solve(cov(residMat) * ((n-1) / (n - ncol(xMat))))
  yCovMatSolve[lower.tri(yCovMatSolve)] <- t(yCovMatSolve)[lower.tri(yCovMatSolve)]
  outcomeCovVec <- diag(outcomeCovMat)
  # use this to innovate
  outcomeCorMat <- sweep(outcomeCovMat, MARGIN=2, STATS = sqrt(outcomeCovVec), FUN="/") %>%
    sweep(., MARGIN=1, STATS = sqrt(outcomeCovVec), FUN="/")
  # make it symmetric
  outcomeCorMat[lower.tri(outcomeCorMat)] <- t(outcomeCorMat)[lower.tri(outcomeCorMat)]
  # eigenvalues for MP SKAT
  #eValsMP <- eigen(snpCovMat[1, 1] * outcomeCovMatSolve, symmetric = TRUE)$values
  eValsMP <- eigen(snpVarMP * yCovMatSolve, symmetric = TRUE)$values 
  # eigenvalues for MSKAT Q statistic
  #eValsMSKAT <- rep(snpCovMat[1, 1], ncol(yMatMP))
  eValsMSKAT <- rep(snpVarMP, ncol(yMatMP))
  
  # variance component snp-set
  # this is equal to the output from the SKAT package times 2
  nullMod <- SKAT_Null_Model(yVecSS ~ xMat - 1)
  skatOutputNW <- SKAT(Z = genoMat, obj=nullMod, weights=rep(1, ncol(genoMat)))
  skatOutput <- SKAT(Z = genoMat, obj=nullMod)
  vcSS <- sum(snpScoreStats^2) / sigmaSqYHatSS
  # this p-value is very close to what SKAT spits out
  vcSS_pval <- CompQuadForm::davies(q=vcSS, lambda=eValsSS, delta=rep(0,length(eValsSS)), acc = 1e-09, lim = 1e+06)$Qq
  vcMP <- t(outcomeScoreStats) %*% yCovMatSolve %*% yCovMatSolve %*% outcomeScoreStats 
  vcMP_pval <- CompQuadForm::davies(q=vcMP, lambda=eValsMP, delta=rep(0,length(eValsMP)), acc = 1e-09, lim = 1e+06)$Qq
  # this is the Q of mskat
  mskat <- t(outcomeScoreStats) %*% yCovMatSolve %*% outcomeScoreStats 
  # mskat pvalue
  mskat_pval <- CompQuadForm::davies(q=mskat, lambda=eValsMSKAT, delta=rep(0,length(eValsMSKAT)), acc = 1e-09, lim = 1e+06)$Qq
  # can check mskat with this:
  mskatOutput <-  MSKAT(MSKAT.cnull(yMatMP, xMat[,-1]), matrix(data=gVecMP, ncol=1), W = rep(1, 1))
  
  return(list(outcomeStandStats = outcomeStandStats, snpStandStats = snpStandStats, sigmaSqYHatSS = sigmaSqYHatSS,
              snpCorMat = snpCorMat, outcomeCorMat = outcomeCorMat,
              vcSS = vcSS, vcSS_pval = vcSS_pval, vcMP = vcMP, vcMP_pval= vcMP_pval,
              skatOutput = skatOutput, skatOutputNW = skatOutputNW, mskatOutput = mskatOutput, mskat_pval = mskat_pval))
}

# save the generalized and innovated test statistics for each case
genSS <- matrix(data=NA, nrow=B, ncol=D)
innSS <- matrix(data=NA, nrow=B, ncol=D)
genMP <- matrix(data=NA, nrow=B, ncol=D)
innMP <- matrix(data=NA, nrow=B, ncol=D)

# loop
set.seed(1000 + 1000*Snum + aID)
resultsDF <- data.frame(rho1=rep(rho1, B), rho2=rho2, rho3=rho3, iBJ_SS=NA, iHC_SS=NA, VC_SS=NA,
                        ACAT_SS=NA, GBJ_SS=NA, GHC_SS=NA, iBJ_MP=NA, iHC_MP=NA, VC_MP=NA,
                        ACAT_MP=NA, GBJ_MP=NA, GHC_MP=NA, mskat=NA, skat=NA, skatNW=NA, mymskat=NA)
for (sim_it in 1:B) {
 
  # new effect sizes every time
  gammaMat[1:r1, 1:r1] <- runif(n=r1^2, min=gammaMinSS, max=gammaMaxSS)

  # simulate D genotypes for SS
  genoMat <- rmvbin(n=n, margprob = mafVec, sigma=rmvbinSigma) + 
    rmvbin(n=n, margprob = mafVec, sigma=rmvbinSigma)
  # simulate genotype for MP
  gVecMP <- rbinom(n=n, size=2, prob=mafVec[1])

  # simulate other covariates
  xMat <- cbind(1, rnorm(n), rbinom(n, size=1, prob=0.5))
  
  # truth for SNP-set
  eMat <- rmvnorm(n=n, mean=rep(0, D), sigma=oppCor)
  yOrigSS <- xMat %*% betaMat + genoMat %*% gammaMat + eMat
  # center the y
  yMeansSS <- apply(yOrigSS, 2, mean)
  yMatSS <- sweep(yOrigSS, MARGIN=2, STATS = yMeansSS, FUN="-")
  yVecSS <- yMatSS[, 1] - mean(yMatSS[, 1])

  # truth for MP
  gammaVec <- rep(0, D)
  gammaVec[1:r1] <- runif(n=r1, min=gammaMinMP, max = gammaMaxMP) 
  eMat2 <- rmvnorm(n=n, mean=rep(0, D), sigma=corMat)
  yOrigMP <- xMat %*% betaMat + gVecMP %*% t(gammaVec) + eMat2
  yMeansMP <- apply(yOrigMP, 2, mean)
  yMatMP <- sweep(yOrigMP, MARGIN=2, STATS = yMeansMP, FUN="-")

  # calculate test statistics
  testStats <- getTestStatistics(genoMat = genoMat, yMatMP = yMatMP, xMat = xMat, 
                                 yVecSS = yVecSS, gVecMP = gVecMP) 
  
  # SNP-set tests
  ssGBJ <- GBJ(test_stats = testStats$snpStandStats, cor_mat = testStats$snpCorMat)
  ssGHC <- GHC(test_stats = testStats$snpStandStats, cor_mat = testStats$snpCorMat)
  decompSNP <- eigen(testStats$snpCorMat, symmetric = TRUE)
  SNPinvRoot <- sweep(t(decompSNP$vectors), MARGIN=1, STATS=sqrt(decompSNP$values), FUN="/")
  iSNPstats <- SNPinvRoot %*% testStats$snpStandStats
  ssiBJ <- BJ(test_stats = iSNPstats, cor_mat = diag(rep(1, D)))
  ssiHC <- HC(test_stats = iSNPstats, cor_mat = diag(rep(1, D)))
  ssACAT <- ACAT(1 - pchisq(as.numeric(testStats$snpStandStats)^2, df=1))

  # MP tests
  mpGBJ <- GBJ(test_stats = testStats$outcomeStandStats, cor_mat = testStats$outcomeCorMat)
  mpGHC <- GHC(test_stats = testStats$outcomeStandStats, cor_mat = testStats$outcomeCorMat)
  decompOutcome <- eigen(testStats$outcomeCorMat, symmetric = TRUE)
  outcomeInvRoot <- sweep(t(decompOutcome$vectors), MARGIN=1, STATS=sqrt(decompOutcome$values), FUN="/")
  iOutcomeStats <- outcomeInvRoot %*% testStats$outcomeStandStats
  mpiBJ <- BJ(test_stats = iOutcomeStats, cor_mat = diag(rep(1, D)))
  mpiHC <- HC(test_stats = iOutcomeStats, cor_mat = diag(rep(1, D)))
  mpACAT <- ACAT(1 - pchisq(as.numeric(testStats$outcomeStandStats)^2, df=1))
  
  # record
  resultsDF$iBJ_SS[sim_it] <- ssiBJ$BJ_pvalue
  resultsDF$iHC_SS[sim_it] <- ssiHC$HC_pvalue
  resultsDF$VC_SS[sim_it] <- testStats$vcSS_pval
  resultsDF$ACAT_SS[sim_it] <- ssACAT
  resultsDF$GBJ_SS[sim_it] <- ssGBJ$GBJ_pvalue
  resultsDF$GHC_SS[sim_it] <- ssGHC$GHC_pvalue
  resultsDF$skat[sim_it] <- testStats$skatOutput$p.value 
  resultsDF$skatNW[sim_it] <- testStats$skatOutputNW$p.value 
  resultsDF$iBJ_MP[sim_it] <- mpiBJ$BJ_pvalue
  resultsDF$iHC_MP[sim_it] <- mpiHC$HC_pvalue
  resultsDF$VC_MP[sim_it] <- testStats$vcMP_pval
  resultsDF$ACAT_MP[sim_it] <- mpACAT
  resultsDF$GBJ_MP[sim_it] <- mpGBJ$GBJ_pvalue
  resultsDF$GHC_MP[sim_it] <- mpGHC$GHC_pvalue
  resultsDF$mskat[sim_it] <- testStats$mskatOutput$p.value[1]
  resultsDF$mymskat[sim_it] <- testStats$mskat_pval

  # checkpoint
  cat(sim_it)
}

# write
setwd(outputDir)
write.table(resultsDF, outName, row.names=F, col.names=T, sep='\t', quote=F)





