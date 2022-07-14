# Recreate Figure 1
# Will need data from protected sources (we provide scrambled versions in repository)

# Change to current working directory
setwd('/users/rsun3/desktop/reproduce')

library(ggplot2)
library(data.table)
library(cowplot)
library(dplyr)
library(magrittr)

# Load protected data
mpRes <- fread("mp_blood_scrambled.txt") 
ssRes <- fread("snpset_blood_scrambled.txt")
sigSNPs <- fread("sigSNPsMinus2_scrambled.txt") %>%
  dplyr::select(RS, Chr, BP, pvalue, subset) %>%
  mutate(chrBP = paste0(Chr, BP))


# RAD52 data
rad52SS <- ssRes %>% filter(Locus == "RAD52") %>%
  mutate(chrBP = paste0(chr, BP)) %>%
  merge(sigSNPs %>% filter(subset == "Squam"), by="chrBP") %>%
  dplyr::select(chr, BP.x, testStat, pvalue) %>%
  set_colnames(c("Chr", "BP", "testStat", "LCp")) %>%
  mutate(bloodP = 1 - pchisq(testStat^2, df=1)) %>%
  mutate(logLC = -log10(LCp)) %>%
  mutate(colorLC = logLC / 15)

# create RAD52 plot
plotSS <- ggplot(data=rad52SS, aes(x=BP, y=-log10(bloodP), color = colorLC)) + 
  geom_point() + 
  scale_color_gradient(low="grey", high="darkred", name="SNP\nAssociation\nwith Lung\nCancer", 
                       labels=c(expression('p=10'^-4), expression('p=10'^-12)), breaks=c(0.27, 0.8)) + 
  xlab("SNP Location (BP)") + ylab("-log10(Association with RAD52 Expression)") + 
  theme_cowplot()


# data for TERT plot
load("ensembl_refgene_hg19_20180109.rda")
geneRanges <-  data.table(ensembl_refgene_hg19_20180109)
tertMP <- mpRes %>% filter(exprLocus == "TERT") %>%
  # the one used in data analysis
  filter(BP == 1285974) %>%
  dplyr::select(Gene, pval, snpLocus, exprLocus) %>% 
  merge(geneRanges %>% filter(Notes <= 0) %>% dplyr::select(txStart, txEnd, HGNC_name) %>% 
                                set_colnames(c("start", "end", "Gene")), by="Gene") %>%
  mutate(logP = -log10(pval)) %>%
  mutate(mid = (start + end) / 2) %>%
  mutate(ylab = ifelse(Gene == "TRIP13", logP - 0.14, logP - 0.08)) %>%
  mutate(Gene = ifelse(Gene == "CLPTM1L", "", Gene)) %>%
  mutate(Gene = ifelse(Gene == "MRPL36", "", Gene)) %>%
  mutate(Gene = ifelse(Gene == "NKD2", "", Gene)) %>%
  mutate(Gene = ifelse(Gene == "AHRR", "", Gene)) %>%
  mutate(mid = ifelse(Gene == "TERT", mid + 25000, mid)) %>%
  mutate(mid = ifelse(Gene == "PDCD6", mid + 25000, mid)) 
mpPlotDat <- rbind(tertMP %>% dplyr::select(Gene, logP, start) %>% set_colnames(c("Gene", "logP", "BP")) %>% mutate(Pos = "start"),
                   tertMP %>% dplyr::select(Gene, logP, end) %>% set_colnames(c("Gene", "logP", "BP")) %>% mutate(Pos = "end"))

# plot TERT
plotMP <- ggplot(data=mpPlotDat, aes(x=BP, y=logP)) + 
  geom_point(shape = 15) + 
  geom_text(data = tertMP, aes(x = mid, y = ylab, label = Gene), inherit.aes = FALSE) + 
  xlab("Gene Location (BP)") + ylab("-log10(Association with SNP)") + 
  geom_segment(aes(x = start, y = logP, xend = end, yend = logP), data=tertMP) +
  theme_cowplot()

# put together
plot_grid(plotSS, plotMP, ncol=2, labels=c("A", "B"), label_size = 20)






