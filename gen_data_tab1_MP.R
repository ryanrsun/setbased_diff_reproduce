# Generate data for Table 1
# This script is designed to be run on a computing cluster.
# We need to run 18 jobs with one input argument.
# Run it with leading argument (aID) 1-18.

# Note that much of this script is used to clean the raw GTEx data, which we cannot share.
# We have commented out those parts and instead generate scrambled data for use.

library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(magrittr)
library(dplyr)
library(tidyr)
library(data.table)
library(LungCancerAssoc)
library(tibble)
library(GBJ)
library(SKAT)
source("ACAT.R")

# do one expr at a time
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])

# output
individualOutname <- paste0("mp_blood_cis_10kbsnp_2mbexpr_individual", aID, ".txt")
setbasedOutname <- paste0("mp_blood_cis_10kbsnp_2mbexpr_setbased", aID, ".txt")

# location of the GTEx_v7 folder
rootDir <- "/rsrch3/home/biostatistics/rsun3/GTEx_v7/"

# read genotype information from disk and hold it in a smaller format (somehow).
# the gds was constructed from plink files, which were constructed from vcf files by Andy Shi,
# he used vcftools to remove indels, sites with MAF < 0.05, and sex chromosomes.
read_gds <- function(gdsName) {
  gds <- GdsGenotypeReader(gdsName, YchromCode=24L, XYchromCode=25L)
 
  # all information packed into the file by GTEx 
  scanID <- getScanID(gds)
  family_info <- getVariable(gds, "sample.annot/family")
  father_info <- getVariable(gds, "sample.annot/father")
  mother_info <- getVariable(gds, "sample.annot/mother")
  # sex must be recorded as M/F/NA - it appears to be all "" in this dataset
  sex <- getVariable(gds, "sample.annot/sex")
  sex[sex == ""] <- NA
  phenotype <- getVariable(gds, "sample.annot/phenotype")

  # have to put the annotation data in a special structure, that special structure will go inside
  # another special structure, all this to satisfy the file container requirements
  scanAnnot <- data.frame(scanID, father_info, mother_info, sex, phenotype, stringsAsFactors = FALSE) %>%
    ScanAnnotationDataFrame()

  # some information about the SNPs, positions + alleles
  snpID <- getSnpID(gds)
  chromosome <- getChromosome(gds)
  position <- getPosition(gds)
  alleleA <- getAlleleA(gds)
  alleleB <- getAlleleB(gds)
  rsID <- getVariable(gds, "snp.rs.id")
  snpAnnot <- data.frame(snpID, chromosome, position, rsID, alleleA, alleleB,
                          stringsAsFactors = FALSE) %>%
    SnpAnnotationDataFrame(YchromCode=24L, XYchromCode=25L)

  # final object to be returned
  genoData <- GenotypeData(gds, scanAnnot = scanAnnot, snpAnnot = snpAnnot)
  # return
  return(list(genoData = genoData, gds = gds))
}


# clean and subset the genotype and gene expression data
clean_and_subset <- function(exprGeneNames, tissue, genoData, snpsToGet, snpRange=TRUE,
                             covarDF, exprMeta, exprDF, INT=TRUE) {

  # need the snpIDs of the SNPs that we specified
  # if snpRange=TRUE, then snpsToGet should have Chr, Start, End for each range
  # if snpRange=FALSE, then snpsToGet should have Chr, Position for each SNP 
  snpDF <- c()
  uniqueChrs <- unique(snpsToGet$Chr)
  for (chr_it in 1:length(uniqueChrs)) {
    tempChr <- uniqueChrs[chr_it]
    tempToGet <- snpsToGet %>% filter(Chr == tempChr)
    if (snpRange) {
      # specify a range of positions, get all snps in those positions
      for (row_it in 1:nrow(tempToGet)) {
        tempStart <- tempToGet$start[row_it]
        tempEnd <- tempToGet$end[row_it]
        tempDF <- snpPosDF %>% filter(chromosome == tempChr & 
                                     position <= tempEnd & position >= tempStart)  
        snpDF <- rbind(snpDF, tempDF)
      }
    } else {
      # specify SNPs precisely by position 
      tempPos <- tempToGet %>% select(Position) %>% unlist(.)
      tempDF <- snpPosDF %>% filter(chromosome == tempChr & position %in% tempPos)
      snpDF <- rbind(snpDF, tempDF)
    }
  }

  # can't find it
  if (nrow(snpDF) == 0) {
   return(-1)
  }

  # pull the genotypes from the gds 
  gMat <- getGenotypeSelection(genoData, snpID = snpDF$snpID)
  
  # check names
  if (class(gMat) == "matrix") { 
    stopifnot(rownames(gMat) == as.character(snpDF$snpID))
  } else {
    tempNames <- names(gMat)
    gMat <- matrix(data=gMat, nrow=1)
    colnames(gMat) <- tempNames
  }
  # use rsIDs
  rownames(gMat) <- snpDF$rsID

  # extract expression relevant tissue samples
  tissueExprSamps <- exprMeta %>% filter(Note == tissue)
  
  # overlap genotype and expression subjects
  overlapIDs <- intersect(tissueExprSamps$subj_id, colnames(gMat))
  #message(sprintf("Number of subjects: %d", length(overlapIDs)))

  # get the overlap data, pick a single tissue sample if there are multiple
  tissueExprOverlap <- tissueExprSamps %>% filter(subj_id %in% overlapIDs) %>%
    # even if the subject is in the tissue metadata and genotype, the actual gene reads might not
    # be in the exprDF file for whatever reason (although it will be in the transcript file, just not gene file...)
    filter(Sample %in% colnames(exprDF)) %>%
    unique(by = "subj_id") 

  # after this, exprSelected is p*(n+2) where rows are genes and columns are name/description of gene + sample ID
  # gMatSelected is n*q where rows are subjects and columns are snpID
  exprSelected <- exprDF %>% as.data.frame(.) %>%
    select(Name, Description, all_of(tissueExprOverlap$Sample)) %>%
    filter(Description %in% exprGeneNames) %>%
    # sometimes two columns for the same gene
    distinct_at("Description", .keep_all = TRUE) %>% 
    column_to_rownames(var = "Description") %>% 
    select(-Name) %>%
    t(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(var = "Sample") %>%
    merge(., tissueExprOverlap %>% select(Sample, subj_id), by="Sample") %>%
    select(-Sample)

  # GTEx did TMM (we won't do) and also INT for each gene across samples
  if (INT) {
    n <- nrow(exprSelected)
    # there are actually no NAs in the entire gene tpm file, so we don't have to worry about that
    saveID <- exprSelected$subj_id
    exprINT <- exprSelected %>% select(-subj_id)
    for (gene_it in 1:ncol(exprINT)) {
      tempCol <- exprINT[, gene_it]
      # use the blom offset of k=3/8
      tempTrans <- qnorm( (rank(tempCol) - 3/8) / (n + 0.25) )
      exprINT[, gene_it] <- tempTrans
    }
    # add back subject ID
    exprINT <- exprINT %>% mutate(subj_id = saveID)
  } else {
    exprINT <- exprSelected
  }

  # select and transpose genotype matrix
  # works even if one SNP 
  genotypeSubToSelect <- colnames(gMat)[which(colnames(gMat) %in% tissueExprOverlap$subj_id)]
  gMatSelected <- as.data.frame(gMat) %>% select(all_of(genotypeSubToSelect)) %>%
    t(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(var = "subj_id") %>%
    # can't have numbers be columns for lm()
    set_colnames(paste0("rs_", colnames(.)))
  IDcol <- which(colnames(gMatSelected) == "rs_subj_id")
  colnames(gMatSelected)[IDcol] <- "subj_id"

  # covariates
  # may have slightly more subjects, but we will merge in next step
  covarSelected <- covarDF %>% column_to_rownames(var="ID") %>%
    t(.) %>%
    as.data.frame(.) %>% 
    rownames_to_column(var="subj_id")
  # there's a gene named C2
  cIdx <- which(colnames(covarSelected) == "C2")
  colnames(covarSelected)[cIdx] <- "Cov2"

  # merge
  allDat <- merge(gMatSelected, exprINT, by="subj_id") %>%
    merge(covarSelected, by="subj_id")

  return(list(allDat = allDat, snpDF = snpDF, exprINT = exprINT, gMatSelected = gMatSelected))
}

# do the linear regression
MPreg <- function(covarNames, genoNames, exprNames, allDat) {
 
  nGeno <- length(genoNames)

  # string for regression formula
  formulaRoot <- paste0(" ~ ", paste(covarNames, collapse=" + "))
   
  # two loops, outer is expression, inner is SNP
  resultsDF <- c() 
  setbasedDF <- c() 
  for (snp_it in 1:nGeno) {

    tempSNP <- genoNames[snp_it]
    tempGeno <- allDat %>% select(all_of(tempSNP)) %>% unlist(.)
    tempMAF <- mean(tempGeno, na.rm = TRUE) / 2
    
    # check for NA 
    NAidx <- which(is.na(tempGeno))
    if (length(NAidx) > 0) {
      tempGeno[NAidx] <- rbinom(n=length(NAidx), size=2, prob=tempMAF)
    }

    # all outcomes
    rawY <- allDat %>% select(all_of(exprNames)) %>%
      as.matrix(.)
    # filter out ones that don't vary
    lengthUnique <- function(x) {length(unique(x))}
    numUnique <- apply(rawY, 2, lengthUnique)
    yMat <- rawY[, which(numUnique > 10)]

    # design matrix
    xMat <- cbind(1, allDat %>% select(all_of(covarNames))) %>%
      as.matrix(.)
    projX <- xMat %*% solve(t(xMat) %*% xMat) %*% t(xMat)
    # residuals 
    fittedYMat <- projX %*% yMat
    residMat <- yMat - fittedYMat
    # estimated variance 
    sigmaSqYHat <- apply(residMat^2, 2, sum) / (nrow(xMat) - ncol(xMat))
    # variance-type quantity involving projection matrix 
    Ppart <- as.numeric(sum(tempGeno^2) - t(tempGeno) %*% projX %*% tempGeno)

    # one SNP against all the outcomes
    outcomeScoreStats <- as.numeric(t(tempGeno) %*% residMat)
    # can check against function checkMP()
    outcomeStandStats <- outcomeScoreStats / (sqrt(rep(Ppart, length(sigmaSqYHat)) * sigmaSqYHat))

    # covariance matrix of Y
    n <- nrow(residMat) 
    outcomeCovMat <- cov(residMat) * ((n-1) / (n - ncol(xMat))) * Ppart
    # this solve is for the MP SKAT, where we don't need the scaling term
    #outcomeCovMatSolve <- solve(outcomeCovMat / Ppart)
    fullCovVec <- diag(outcomeCovMat)
    # use this to innovate
    outcomeCorMat <- sweep(outcomeCovMat, MARGIN=2, STATS = sqrt(fullCovVec), FUN="/") %>%
      sweep(., MARGIN=1, STATS = sqrt(fullCovVec), FUN="/")
    # make it symmetric
    outcomeCorMat[lower.tri(outcomeCorMat)] <- t(outcomeCorMat)[lower.tri(outcomeCorMat)]
   
    # prune as necessary
    donePrune <- FALSE
    prunedCor <- outcomeCorMat
    prunedCov <- outcomeCovMat
    prunedScoreStats <- outcomeScoreStats
    prunedStandStats <- outcomeStandStats
    iters <- 0
    while(!donePrune) {
      keepIdx <- 1:nrow(prunedCor)
      
      # decorrelate for innovation
      decomp <- eigen(prunedCor, symmetric = TRUE)
      if (decomp$values[length(decomp$values)] > 10^(-5)) {
        donePrune <- TRUE
        corInvRoot <- sweep(t(decomp$vectors), MARGIN = 1, STATS = sqrt(decomp$values), FUN = "/")
        next
      }

      # put in NA
      diag(prunedCor) <- NA
      iters <- iters + 1
      # prune until prunedCor is invertible
      highestCorrs <- which(abs(prunedCor) == max(abs(prunedCor), na.rm = TRUE), arr.ind = TRUE)
      keepIdx <- keepIdx[-highestCorrs[1, 1]]
      prunedCor <- prunedCor[keepIdx, keepIdx]
      prunedCov <- prunedCov[keepIdx, keepIdx] 
      prunedScoreStats <- prunedScoreStats[keepIdx]
      prunedStandStats <- prunedStandStats[keepIdx]
      # remove NA for decomp
      diag(prunedCor) <- 1
    }
     
    # this solve is for MP SKAT, we don't need scaling term
    outcomeCovMatSolve <- solve(prunedCov / Ppart) 
    
    # eigenvalues for MP SKAT
    eValsMP <- eigen(Ppart * outcomeCovMatSolve, symmetric = TRUE)$values
    # eigenvalues for MSKAT Q statistic
    eValsMSKAT <- rep(Ppart, ncol(yMat))
    vcMP <- t(prunedScoreStats) %*% outcomeCovMatSolve %*% outcomeCovMatSolve %*% prunedScoreStats
    vcMP_pval <- CompQuadForm::davies(q=vcMP, lambda=eValsMP, delta=rep(0,length(eValsMP)), acc = 1e-09, lim = 1e+06)$Qq
    # this is the Q of mskat
    mskat <- t(prunedScoreStats) %*% outcomeCovMatSolve %*% prunedScoreStats
    # mskat pvalue
    mskat_pval <- CompQuadForm::davies(q=mskat, lambda=eValsMSKAT, delta=rep(0,length(eValsMSKAT)), acc = 1e-09, lim = 1e+06)$Qq
      
    # record 
    tempResults <- data.frame(Gene=rep(NA, length(prunedStandStats)), SNP=NA, MAF=NA, testStat=NA, pval=NA)
    tempResults$Gene <- names(prunedStandStats)
    tempResults$SNP <- tempSNP
    tempResults$MAF <- tempMAF
    tempResults$testStat <- prunedStandStats
    tempResults$pval <- 1 - pchisq(prunedStandStats^2, df=1)

    # set-based 
    gbjOutput <- GBJ(test_stats = prunedStandStats, cor_mat = prunedCor)
    ghcOutput <- GHC(test_stats = prunedStandStats, cor_mat = prunedCor)
    acatOutput <- ACAT(1 - pchisq(as.numeric(prunedStandStats)^2, df=1)) 
    iStats <- as.numeric(corInvRoot %*% prunedStandStats)
    ibjOutput <- BJ(test_stats = iStats, cor_mat = diag(rep(1, length(iStats))))
    ihcOutput <- HC(test_stats = iStats, cor_mat = diag(rep(1, length(iStats))))

    # about the correlation matrix
    diag(prunedCor) <- NA
    medianRho <- median(prunedCor, na.rm=TRUE)
    medianAbsRho <- median(abs(prunedCor), na.rm=TRUE)

    # record 
    tempSetbased <- data.frame(SNP = tempSNP, GBJ=gbjOutput$GBJ_pvalue, GHC=ghcOutput$GHC_pvalue,
                              ACAT=acatOutput, SKAT=vcMP_pval, iBJ=ibjOutput$BJ_pvalue,
                              iHC=ihcOutput$HC_pvalue, nExprOrig=length(outcomeStandStats),
                              nExprFinal=length(iStats), medianRho=medianRho, medianAbsRho=medianAbsRho)
    
    # append
    resultsDF <- rbind(resultsDF, tempResults)
    setbasedDF <- rbind(setbasedDF, tempSetbased)
  }
  # return 
  return(list(individualDF = resultsDF, setbasedDF = setbasedDF))
}

#----------------------------------------------------_#
# start analysis
gdsName <- paste0(rootDir, "GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_SNPs_maf.gds")

# read from gds on disk
#retGds <- read_gds(gdsName)
#genoData <- retGds$genoData

# make SNP position df
#snpPosDF <- getSnpAnnotation(genoData) %>%
#  pData() 

# read expression metadata
#exprMeta <- fread(paste0(rootDir, "ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_metrics.tsv"))

# read expression dataset - blood
#exprDF <- fread(paste0(rootDir, "ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_WholeBlood_Gene.txt"))

# read covariates
#covarFilename <- paste0(rootDir, "Covariates/Whole_Blood.v7.covariates.txt")
#covarDF <- fread(covarFilename)

# convert GTEx IDs to subject IDs
#subjIDs <- strsplit(exprMeta$Sample, "-") %>%
#  sapply(function(x) {paste0(x[1], "-", x[2])})
#exprMeta$subj_id <- subjIDs

# read gene location information
load(file="ensembl_refgene_hg19_20180109.rda")
geneRanges <-  data.table(ensembl_refgene_hg19_20180109)

# loop through the 18 loci given in Table 2 of McKay et al.
tab2DF <- data.frame(RS = c("rs_1_77967507_T_A_b37", "rs_6_167376466_G_C_b37", " rs_8_27344719_G_A_b37", 
                            "rs_13_32970031_T_C_b37", 
                            "rs_15_47576969_C_T_b37", "rs_15_78857986_C_G_b37", "rs_19_41353107_T_C_b37 ", 
                            "rs_3_189357199_G_T_b37",
                            "rs_5_1285974_C_A_b37", "rs_8_32410110_G_A_b37", "rs_9_21830157_A_G_b37", 
                            "rs_10_105687632_A_C_b3", 
                            "rs_11_118125625_T_C_b37", "rs_15_49376624_T_G_b37", "rs_20_62326579_G_T_b37", 
                            "rs_6_31434111_A_G_b37",
                            "rs_12_998819_G_C_b37", "rs_22_29130012_T_C_b37"),
                     Gene = c("FUBP1", "RNASET2", "CHRNA2", "BRCA2", "SEMA6D", 
                              "CHRNA5", "CYP2A6", "TP63", "TERT", "NRG1", "CDKN2A",
                              "OBFC1", "AMICA1", "SECISBP2L", "RTEL1", "HCP5", "RAD52", "CHEK2"),
                     BP = c(77967507, 167376466, 27344719, 32970031, 47576969, 78857986,
                            41353107, 189357199, 1285974, 32410110, 21830157, 105687632,
                            118125625, 49376624, 62326579, 31434111, 998819, 29130012),
                     Chr = c(1, 6, 8, 13, 15, 15, 19, 3, 5, 8, 9, 10, 11, 15, 20, 6, 12, 22))
setbasedDF <- c()
cisLength <- 1000000
sentinelWindow <- 0
locusToDo <- as.character(tab2DF$Gene[aID])
# outer loop through each SNP locus
for (snp_locus_it in 1:length(locusToDo)) {
 
  tempSNPlocus <- as.character(locusToDo[snp_locus_it])
  tab2Minus <- tab2DF %>% filter(Gene == tempSNPlocus) 
  snpsToGet <- tab2DF %>% filter(Gene == tempSNPlocus) %>%
    mutate(start = BP - sentinelWindow, end = BP + sentinelWindow)
  
  # inner loop through the expr locus
  for (expr_locus_it in 1:nrow(tab2Minus)) {
   
    tempExprLocus <- as.character(tab2Minus$Gene[expr_locus_it]) 
    tempChr <- tab2Minus$Chr[expr_locus_it]
    tempBP <- tab2Minus$BP[expr_locus_it]

    # expression gene names
    exprTab <- geneRanges %>% filter(Chr == tempChr) %>%
      filter(tempBP - txEnd <= cisLength & txStart - tempBP <= cisLength) %>%
      filter(Notes <= 1)

    # expression gene names
    #exprWanted <- as.character(unlist(exprTab$HGNC_name))

    # clean and merge expression, genotypes, covariates
    #bloodCleaned <- clean_and_subset(exprGeneNames = exprWanted, tissue = "Whole Blood", genoData = genoData, 
    #                            snpsToGet = snpsToGet, snpRange=TRUE, covarDF = covarDF, 
    #                            exprMeta = exprMeta, exprDF = exprDF, INT=TRUE) 
    
    

    # found the SNP?
    #if (class(bloodCleaned) == "numeric") {next}
    # only use those SNPs that were p<10^-5 in the association analysis
    #snpInfo <- bloodCleaned$snpDF %>%
    #  merge(., sigSNPs, by=c("chromosome", "position"))

    # names of genes and SNPs
    #exprNames <- colnames(bloodCleaned$exprINT %>% select(-subj_id))
    #genoNames <- colnames(bloodCleaned$gMatSelected %>% select(-subj_id))
    #genoNames <-  paste0("rs_", as.character(snpInfo$rsID))

    # names of covariates
    #covarNames <- unlist(covarDF$ID)
    # there is a gene named C2
    #covarNames[which(covarNames == "C2")] <- "Cov2"

    # checkpoint
    #cat("Starting number ", expr_locus_it, "\n")
    #cat("Number of genes: ", length(exprNames), "\n")
    #cat("Number of SNPs: ", length(genoNames), "\n")

    #------------------------------------------#
    # NEW FOR SCRAMBLED DATA
    covarNames <- c("cov1", "cov2")
    genoNames <- c("rs_1_77967507_T_A_b37")
    exprNames <- paste0("Gene", 1:10)
    allDat <- cbind(rbinom(n=400, size=2, prob=0.5), rnorm(400), rnorm(400),
                    matrix(data=rnorm(400*10), nrow=400))
    colnames(allDat) <- c("rs_1_77967507_T_A_b37", "cov1", "cov2", paste0("Gene", 1:10))
    bloodCleaned <- list(allDat = data.frame(allDat))
    #------------------------------------------#
    
    # run eQTL analysis
    mpOutput <- MPreg(covarNames = covarNames, genoNames = genoNames, 
                      exprNames = exprNames, allDat = bloodCleaned$allDat)
    
    setbasedDF <- rbind(setbasedDF, mpOutput$setbasedDF %>% 
                        mutate(snpLocus = tempSNPlocus, exprLocus = tempExprLocus))
  }
}

write.table(setbasedDF, setbasedOutname, append=F, quote=F, row.names=F, col.names=T, sep='\t')






