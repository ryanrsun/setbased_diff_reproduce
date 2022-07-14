# Generate data for Figure 3
# This script is designed to be run on a computing cluster.
# We need to run 59 + 90 = 149 jobs with two input arguments.
# First, we need to run it with leading argument (aID) 1-59 and second argument (Snum) 1.
# Then we need to run it with leading argument (aID) 1-90 and second argument (Snum) 3.

library(combinat)
library(GBJ)
library(dplyr)
library(mvtnorm)
library(magrittr)
library(data.table)
library(DBpower)
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# change the output directories based on your file structure
outDir <- "/rsrch3/home/biostatistics/rsun3/finiteSamplePower/simulation/output_bounds_theory"
outName <- paste0("bounds_theory_MP_Gen_SS25_MP40_S", Snum, "_aID", aID, ".txt")

# Doing MP or SS
SS <- FALSE
# Doing Gen or Inn
Gen <- TRUE

D <- 40
n <- 400
alpha <- 0.0001
# num causal
r1 <- 4
r2 <- D - r1
# MAF of SNPs
mafValue <- 0.3
gammaValSS <- 0.125
gammaVecSS <- rep(0, D)
gammaVecSS[1:r1] <- gammaValSS
gammaValMP <- 0.2
gammaVecMP <- rep(0, D)
gammaVecMP[1:r1] <- gammaValMP

# build correlation matrix
rho1 <- 0.5
rho2 <- 0.5
rho3 <- seq(from=0, to=0.90, by=0.01)[aID]
corMat <- matrix(data=rho3, nrow=D, ncol=D)
if (Snum == 1) {
  corMat[1:r1, 1:r1] <- rho1
  corMat[(r1+1):D, (r1+1):D] <- rho2
  diag(corMat) <- 1
} else if (Snum == 3) {
  diag(corMat) <- 1
}

# for SS statistics
covMatSNP <- diag(rep(sqrt(mafValue * (1 - mafValue) * 2), D)) %*% corMat %*% diag(rep(sqrt(mafValue * (1 - mafValue) * 2), D))

# eigendecomposition
eVals <- eigen(corMat)$values
eVecs <- eigen(corMat)$vectors

# use these for theory
if (SS) {
  #muInn <- sqrt(n) * diag(sqrt(eVals)) %*% t(eVecs) %*% gammaVec
  muGen <- as.numeric(sqrt(n) * covMatSNP %*% gammaVecSS) / rep(sqrt(2 * mafValue * (1 - mafValue)), 40)
  transformMat <- diag(1 / sqrt(eVals)) %*% t(eVecs)
  muInn <- as.numeric(transformMat %*% muGen)
} else {
  muInn <- as.numeric(sqrt(2 * mafValue * (1 - mafValue) * n) * diag(1 / sqrt(eVals)) %*% t(eVecs) %*% gammaVecMP)
  # generalized case signals
  muGen <- as.numeric(sqrt(2 * mafValue * (1 - mafValue) * n) * gammaVecMP)
}
muInn
muGen

# hold results
resultsDF <- data.frame(rho3 = rho3, innLower = NA, innUpper = NA, genLower = NA, genUpper = NA)

# if doing Gen
if (Gen) {
  gbjBoundsSim <- set_GBJ_bounds(alpha = alpha, J=D, sig_vec = corMat[lower.tri(corMat)])
  gbjBoundsSim
  GBJstatSim <- GBJ(test_stats = gbjBoundsSim, cor_mat = corMat)
  GBJstatSim

  # do gen
  resultsDF$genLower[1] <- calcb2(lower=TRUE, upper=FALSE, muVec = muGen, sigMat = corMat, bounds=gbjBoundsSim)$lowerProb
  resultsDF$genUpper[1] <- calcb2(lower=FALSE, upper=TRUE, muVec = muGen, sigMat = corMat, bounds=gbjBoundsSim)$upperProb
  
} else {
  # bounds for BJ
  bjBoundsSim <- set_BJ_bounds(alpha = alpha, J=D)
  bjBoundsSim
  BJstatSim <- BJ(test_stats = bjBoundsSim, cor_mat = diag(rep(1, D)))
  BJstatSim

  # do inn
  resultsDF$innLower[1] <- calcb2(lower=TRUE, upper=FALSE, muVec = muInn, sigMat = diag(rep(1, D)), bounds= bjBoundsSim)$lowerProb 
  resultsDF$innUpper[1] <- calcb2(lower=FALSE, upper=TRUE, muVec = muInn, sigMat = diag(rep(1, D)), bounds= bjBoundsSim)$upperProb 
} 

# write it
setwd(outDir)
write.table(resultsDF, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')

