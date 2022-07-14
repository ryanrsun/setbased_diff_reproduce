# For each of the 18 loci reported in the McKay paper, test if the region around the
# sentinel SNP is associated with the gene expression at that SNP in blood.
# The McKay paper did it in lung.

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
source("/rsrch3/home/biostatistics/rsun3/R/3.6.0/ACAT.R")

# do one expr at a time
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# output
sentinelBuffer <- 125000
outputDir <- "/rsrch3/home/biostatistics/rsun3/updatedRuns/analysis/output"
snpOutname <- paste0("snpset_blood_individual_cis_sigsnps_250kb_minus4.txt")
setbasedOutname <- paste0("snpset_blood_setbased_cis_sigsnps_250kb_minus4.txt")

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

# prune SNPs in high correlation to make them more manageable size
prune_high_corr <- function(dat, threshold = 0.999) {
  datNew <- dat
  keepNames <- colnames(datNew)
  cormat <- cor(datNew)
  diag(cormat) <- NA
  iters <- 0
  while (max(abs(cormat), na.rm = TRUE) > threshold) {
    iters <- iters + 1
    highest_corrs <- which(abs(cormat) == max(abs(cormat), na.rm = TRUE), arr.ind = TRUE)
    datNew <- datNew[, -highest_corrs[1, 1]]
    keepNames <- keepNames[-highest_corrs[1, 1]]
    cormat <- cor(datNew)
    diag(cormat) <- NA
  }
  #message("Finished pruning after ", iters, " iterations")
  return(list(datNew = datNew, iters=iters, keepNames = keepNames))
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
 
  # pull the genotypes from the gds 
  gMat <- getGenotypeSelection(genoData, snpID = snpDF$snpID)
  # check names
  stopifnot(rownames(gMat) == as.character(snpDF$snpID))
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
  genotypeSubToSelect <- colnames(gMat)[which(colnames(gMat) %in% tissueExprOverlap$subj_id)]
  gMatSelected <- as.data.frame(gMat) %>% select(all_of(genotypeSubToSelect)) %>%
    t(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(var = "subj_id") %>%
    # can't have numbers be columns for lm()
    set_colnames(paste0("rs_", colnames(.)))
  IDcol <- which(colnames(gMatSelected) == "rs_subj_id")
  colnames(gMatSelected)[IDcol] <- "subj_id"

  #if (pruneSNPs) {
  #  gMatClean <- prune_high_corr(gMatClean)
  #}

  # covariates
  # may have slightly more subjects, but we will merge in next step
  covarSelected <- covarDF %>% column_to_rownames(var="ID") %>%
    t(.) %>%
    as.data.frame(.) %>% 
    rownames_to_column(var="subj_id")

  # merge
  allDat <- merge(gMatSelected, exprINT, by="subj_id") %>%
    merge(covarSelected, by="subj_id")

  return(list(allDat = allDat, snpDF = snpDF, exprINT = exprINT, gMatSelected = gMatSelected))
}

# do the linear regression
eqtlReg <- function(covarNames, genoNames, exprNames, allDat, setBased=FALSE) {
 
  nGeno <- length(genoNames)
  nExpr <- length(exprNames)
  nCovar <- length(covarNames)

  # string for regression formula
  formulaRoot <- paste0(" ~ ", paste(covarNames, collapse=" + "))
   
  # two loops, outer is expression, inner is SNP
  setbasedDF <- c()
  individualDF <- c() 
  for (expr_it in 1:length(exprNames)) {

    # make sure expression is not identically 0
    tempGene <- exprNames[expr_it]
    tempExpr <- allDat %>% select(all_of(tempGene)) %>% unlist(.) 
    if (length(unique(tempExpr)) == 0) {next}

    # have to remove those SNPs with 0 MAF, impute mean for NA
    genoMat <- allDat %>% select(all_of(genoNames)) %>% as.matrix(.)
    genoMAFs <- apply(genoMat, 2, mean, na.rm=TRUE) / 2
    cleanGeno <- genoMat[, which(genoMAFs > 0)] 
    cleanMAFs <- apply(cleanGeno, 2, mean, na.rm=TRUE) / 2
    for (col_it in 1:ncol(cleanGeno)) {
      NAidx <- which(is.na(cleanGeno[, col_it]))
      if (length(NAidx) > 0) {
        cleanGeno[NAidx, col_it] <- rbinom(n=length(NAidx), size=2, prob=cleanMAFs[col_it])
      }
    }

    # have to prune the cleanGeno until the cor_mat is invertible for innovation
    donePrune <- FALSE
    pruneThreshold <- 0.96
    while(!donePrune) {
      pruneThreshold <- pruneThreshold - 0.01
      # have to prune the cleanGeno so that the cor_mat is invertible for innovation
      prunedOutput <- prune_high_corr(dat=cleanGeno, threshold = pruneThreshold) 
      prunedGeno <- prunedOutput$datNew
      keepNames <- colnames(prunedOutput$datNew)
      # null model
      nullMod <- glm(paste0(exprNames[1], formulaRoot), data=allDat)
      # find test stats and cor_mat
      scoreOutput <- calc_score_stats(null_model = nullMod, factor_matrix = prunedGeno, link_function="linear")
      filledCorMat <- scoreOutput$cor_mat
      diag(filledCorMat) <- 1 
      decomp <- eigen(filledCorMat)
      # test smallest eigenvalue
      if (decomp$values[length(decomp$values)] >= 10^(-5)) {
       donePrune <- TRUE
      }
    } 

    # run tests
    gbjOutput <- GBJ(test_stats = scoreOutput$test_stats, cor_mat = scoreOutput$cor_mat)
    ghcOutput <- GHC(test_stats = scoreOutput$test_stats, cor_mat = scoreOutput$cor_mat)
    acatOutput <- ACAT(1 - pchisq(as.numeric(scoreOutput$test_stats)^2, df=1)) 
    skatNull <- SKAT_Null_Model(as.formula(paste0(exprNames[1], formulaRoot)), data=allDat, out_type="C")
    skatOutput <- SKAT(Z = prunedGeno, obj = skatNull, weights=rep(1, ncol(prunedGeno)))
    # decorrelate for innovation 
    corInvRoot <- sweep(t(decomp$vectors), MARGIN=1, STATS=sqrt(decomp$values), FUN="/")
    iStats <- corInvRoot %*% scoreOutput$test_stats
    ibjOutput <- BJ(test_stats = iStats, cor_mat = diag(rep(1, length(iStats))))
    ihcOutput <- HC(test_stats = iStats, cor_mat = diag(rep(1, length(iStats))))

    tempIndividual <- data.frame(SNP = keepNames, testStat = scoreOutput$test_stats)
    tempSetbased <- data.frame(Gene = tempGene, GBJ=gbjOutput$GBJ_pvalue, GHC=ghcOutput$GHC_pvalue,
                               ACAT=acatOutput, SKAT=skatOutput$p.value, iBJ=ibjOutput$BJ_pvalue,
                               iHC=ihcOutput$HC_pvalue, nSNPorig=nGeno, nSNPfinal=length(iStats),
                               medRho = median(scoreOutput$cor_mat, na.rm=TRUE), medAbsRho = median(abs(scoreOutput$cor_mat), na.rm=TRUE))
  
    individualDF <- rbind(individualDF, tempIndividual)
    setbasedDF <- rbind(setbasedDF, tempSetbased) 
  } 
  
  # return
  return(list(individualDF=individualDF, setbasedDF=setbasedDF))
}


#----------------------------------------------------_#
# start analysis
gdsName <- paste0(rootDir, "GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_SNPs_maf.gds")

# read from gds on disk
retGds <- read_gds(gdsName)
genoData <- retGds$genoData

# make SNP position df
snpPosDF <- getSnpAnnotation(genoData) %>%
  pData() 

# read expression metadata
exprMeta <- fread(paste0(rootDir, "ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_metrics.tsv"))

# read expression dataset - lung
exprDF <- fread(paste0(rootDir, "ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_Blood_Gene.txt"))

# read covariates
# should be Lung.v7.covariates.txt or Whole_Blood.v7.covariates.txt
covarFilename <- paste0(rootDir, "Covariates/Whole_Blood.v7.covariates.txt")
covarDF <- fread(covarFilename)

# convert GTEx IDs to subject IDs
subjIDs <- strsplit(exprMeta$Sample, "-") %>%
  sapply(function(x) {paste0(x[1], "-", x[2])})
exprMeta$subj_id <- subjIDs

# load the McKay table
setwd('/rsrch3/home/biostatistics/rsun3/updatedRuns/analysis/')
sigSNPs <- fread("sigSNPsMinus2.txt") %>%
  filter(pvalue < 0.0001) %>%
  select(RS, Chr, BP) %>%
  set_colnames(c("rs_number", "chromosome", "position")) %>%
  distinct()

# loop through the 18 loci given in Table 2 of McKay et al.
# note that McKay misspelled CDKN2A, called JAML AMICA1 instead, and used MHC as a gene even though it's a region.
# for the MHC gene, I put it into the HCP5 gene in my lung cancer paper (it's 525 BP past the transcription end site).
tab2DF <- data.frame(RS = c("rs71658797", "rs6920364", "rs11780471", "rs11571833", 
                            "rs66759488", "rs55781567", "rs56113850", "rs13080835",
                            "rs7705526", "rs4236709", "rs885518", "rs11591710", 
                            "rs1056562", "rs77468143", "rs41309931", "rs116822326",
                            "rs7953330", "rs17879961"),
                     Gene = c("FUBP1", "RNASET2", "CHRNA2", "BRCA2", "SEMA6D", 
                              "CHRNA5", "CYP2A6", "TP63", "TERT", "NRG1", "CDKN2A",
                              "OBFC1", "AMICA1", "SECISBP2L", "RTEL1", "HCP5", "RAD52", "CHEK2"),
                     BP = c(77967507, 167376466, 27344719, 32972626, 47577451, 78857986,
                            41353107, 189357199, 1285974, 32410110, 21830157, 105687632,
                            118125625, 49376624, 62326579, 31434111, 998819, 29121087),
                     Chr = c(1, 6, 8, 13, 15, 15, 19, 3, 5, 8, 9, 10, 11, 15, 20, 6, 12, 22))
# just do one expr at a time
#exprToDo <- tab2DF[aID, ]
#setwd("/rsrch3/home/biostatistics/rsun3/github/LungCancerAssoc/Data")
#load(file="ensembl_refgene_hg19_20180109.rda")
#geneRanges <-  data.table(ensembl_refgene_hg19_20180109)
#promoterLen <- 5000
singleSNPdf <- c()
setbasedDF <- c()
# outer loop through each gene
for (gene_it in 1:nrow(tab2DF)) {

  # expression gene names
  tempGene <- as.character(tab2DF$Gene[gene_it])
  #tab2Minus <- tab2DF %>% filter(Gene != tempGene)
  
  # skip the MHC
  # if (tempGene == "HCP5") {next}

  # inner loop through each locus
  for (locus_it in 1:1) {
    
    # which SNPs to test
    tempLocus <- as.character(tab2DF$Gene[gene_it])
    snpsToGet <- tab2DF %>% filter(Gene == tempLocus) %>% mutate(start = BP - sentinelBuffer, end = BP + sentinelBuffer) 
    #snpsToGet <- geneRanges %>% filter(HGNC_name == tempLocus & Notes <= 1) %>%
    #  mutate(start = txStart - promoterLen, end = txEnd + promoterLen)

    # clean and merge expression, genotypes, covariates
    bloodCleaned <- clean_and_subset(exprGeneNames = tempGene, tissue = "Whole Blood", genoData = genoData, 
                                snpsToGet = snpsToGet, snpRange=TRUE, covarDF = covarDF, 
                                exprMeta = exprMeta, exprDF = exprDF, INT=TRUE) 

    # only use those SNPs that were p<10^-5 in the association analysis
    snpInfo <- bloodCleaned$snpDF %>%
      merge(., sigSNPs, by=c("chromosome", "position"))

    # names of genes and SNPs
    exprNames <- colnames(bloodCleaned$exprINT %>% select(-subj_id))
    # genoNames <- colnames(bloodCleaned$gMatSelected %>% select(-subj_id))
    genoNames <- paste0("rs_", as.character(snpInfo$rsID))
    # names of covariates
    covarNames <- unlist(covarDF$ID)

    # run eQTL analysis
    eqtlOutput <- eqtlReg(covarNames = covarNames, genoNames = genoNames, 
                      exprNames = exprNames, allDat = bloodCleaned$allDat)
    
    # merge with SNP info
    outputTot <- merge(eqtlOutput$individualDF, bloodCleaned$snpDF %>% set_colnames(c("gdsID", "chr", "BP", "SNP", "gdsA", "gdsB")) %>%
                  mutate(SNP = paste0("rs_", SNP)), by="SNP") %>%
      mutate(nchars = nchar(as.character(SNP))) %>%
      mutate(Ref = substr(SNP, nchars - 6, nchars - 6)) %>%
      mutate(Alt = substr(SNP, nchars - 4, nchars - 4)) %>%
      mutate(misAllele = ifelse(gdsA == Ref, 1, 0)) %>% 
      mutate(Locus = tempLocus, GeneExpr = tempGene)

      singleSNPdf <- rbind(singleSNPdf, outputTot)
      setbasedDF <- rbind(setbasedDF, eqtlOutput$setbasedDF)
    # checkpoint
    cat(tempLocus) 
  }
 
  # checkpoint 
  cat(gene_it)
}

setwd(outputDir)
write.table(singleSNPdf, snpOutname, append=F, quote=F, row.names=F, col.names=T, sep='\t')
write.table(setbasedDF, setbasedOutname, append=F, quote=F, row.names=F, col.names=T, sep='\t')






