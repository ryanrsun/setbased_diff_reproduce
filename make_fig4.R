# Recreate Figure 3
# Will need to first generate the data (using provided scripts).
# Once data has been generated, use this file to recreate the plots

# Change to current working directory
setwd('/users/rsun3/desktop/reproduce')

library(dplyr)
library(magrittr)
library(tidyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(stats)

#-----------------------------------------------------------------------#
# simulation 1 - block correlation structure

# read files
fname_root <- 'sim_rho50_mp40_ss25_S1_aID'
resultsBlock <- read.table(paste(fname_root, '1.txt', sep=''), header=T)
for (i in 2:54) {
  temp_fname <- paste(fname_root, i, '.txt', sep='')
  temp_file <- tryCatch(read.table(temp_fname, header=T), warning=function(w) w, error=function(e) e)
  resultsBlock <- rbind(resultsBlock, temp_file)
}

# convert to power
tempAlpha <- 0.0001
powerBlockMP <- data.frame(rho3 = sort(unique(resultsBlock$rho3)), GBJ=NA, GHC=NA, iBJ=NA, iHC=NA, ACAT=NA, VC=NA, VCme=NA)
powerBlockSS <- data.frame(rho3 = sort(unique(resultsBlock$rho3)), GBJ=NA, GHC=NA, iBJ=NA, iHC=NA, ACAT=NA, VC=NA, VCme=NA)
for (row_it in 1:nrow(powerBlockMP)) {
  tempRho3 <- powerBlockMP$rho3[row_it]
  tempDat <- resultsBlock %>% filter(rho3 == tempRho3)
  # MP
  powerBlockMP$GBJ[row_it] <- length(which(tempDat$GBJ_MP < tempAlpha)) / nrow(tempDat)
  powerBlockMP$GHC[row_it] <- length(which(tempDat$GHC_MP < tempAlpha)) / nrow(tempDat)
  powerBlockMP$iBJ[row_it] <- length(which(tempDat$iBJ_MP < tempAlpha)) / nrow(tempDat)
  powerBlockMP$iHC[row_it] <- length(which(tempDat$iHC_MP < tempAlpha)) / nrow(tempDat)
  powerBlockMP$ACAT[row_it] <- length(which(tempDat$ACAT_MP < tempAlpha)) / nrow(tempDat)
  powerBlockMP$VC[row_it] <- length(which(tempDat$mskat < tempAlpha)) / nrow(tempDat)
  powerBlockMP$VCme[row_it] <- length(which(tempDat$VC_MP < tempAlpha)) / nrow(tempDat)
  # SS
  powerBlockSS$GBJ[row_it] <- length(which(tempDat$GBJ_SS < tempAlpha)) / nrow(tempDat)
  powerBlockSS$GHC[row_it] <- length(which(tempDat$GHC_SS < tempAlpha)) / nrow(tempDat)
  powerBlockSS$iBJ[row_it] <- length(which(tempDat$iBJ_SS < tempAlpha)) / nrow(tempDat)
  powerBlockSS$iHC[row_it] <- length(which(tempDat$iHC_SS < tempAlpha)) / nrow(tempDat)
  powerBlockSS$ACAT[row_it] <- length(which(tempDat$ACAT_SS < tempAlpha)) / nrow(tempDat)
  powerBlockSS$VC[row_it] <- length(which(tempDat$skat < tempAlpha)) / nrow(tempDat)
  powerBlockSS$VCme[row_it] <- length(which(tempDat$VC_SS < tempAlpha)) / nrow(tempDat)
}

# need to make it a longer data frame for ggplot
powerBlockSSlong <- powerBlockSS %>%
  select(-VCme) %>%
  pivot_longer(!rho3, names_to="Test", values_to="Power") 
powerBlockMPlong <- powerBlockMP %>%
  select(-VCme) %>%
  pivot_longer(!rho3, names_to="Test", values_to="Power") 

# plot MP
simBlockMP <- ggplot(powerBlockMPlong, aes(x=rho3, y=Power, color=Test, linetype=Test)) +
  geom_smooth(method="loess", se=FALSE, span=0.2) +
  xlab(expression(paste(rho[3]))) + ylab("Power") + 
  ylim(c(0, 1)) +  xlim(c(0, 0.59)) + 
  theme_cowplot() + 
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) 

# plot SS
simBlockSS <- ggplot(powerBlockSSlong, aes(x=rho3, y=Power, color=Test, linetype=Test)) +
  geom_smooth(method="loess", se=FALSE, span=0.2) +
  xlab(expression(paste(rho[3]))) + ylab("Power") + 
  ylim(c(0, 1)) + xlim(c(0, 0.59)) + 
  theme_cowplot() + 
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) 


#-----------------------------------------------------------------------#
# simulation 2
# exchangeable correlation always


fname_root <- 'sim_rho50_mp40_ss25_S3_aID'
resultsExch<- read.table(paste(fname_root, '1.txt', sep=''), header=T)
for (i in 2:90) {
  # New file name, load it
  temp_fname <- paste(fname_root, i, '.txt', sep='')
  temp_file <- tryCatch(read.table(temp_fname, header=T), warning=function(w) w, error=function(e) e)
  
  # Missing?
  if (length(class(temp_file)) > 1) {
    cat("Missing ", i, '\n')
    next
  }
  
  resultsExch <- rbind(resultsExch, temp_file)
}


# convert to power
tempAlpha <- 0.0001
powerExchMP <- data.frame(rho3 = sort(unique(resultsExch$rho3)), GBJ=NA, GHC=NA, iBJ=NA, iHC=NA, ACAT=NA, VC=NA, VCme=NA)
powerExchSS <- data.frame(rho3 = sort(unique(resultsExch$rho3)), GBJ=NA, GHC=NA, iBJ=NA, iHC=NA, ACAT=NA, VC=NA, VCme=NA)
for (row_it in 1:nrow(powerExchMP)) {
  tempRho3 <- powerExchMP$rho3[row_it]
  tempDat <- resultsExch %>% filter(rho3 == tempRho3)
  # MP
  powerExchMP$GBJ[row_it] <- length(which(tempDat$GBJ_MP < tempAlpha)) / nrow(tempDat)
  powerExchMP$GHC[row_it] <- length(which(tempDat$GHC_MP < tempAlpha)) / nrow(tempDat)
  powerExchMP$iBJ[row_it] <- length(which(tempDat$iBJ_MP < tempAlpha)) / nrow(tempDat)
  powerExchMP$iHC[row_it] <- length(which(tempDat$iHC_MP < tempAlpha)) / nrow(tempDat)
  powerExchMP$ACAT[row_it] <- length(which(tempDat$ACAT_MP < tempAlpha)) / nrow(tempDat)
  powerExchMP$VC[row_it] <- length(which(tempDat$mskat < tempAlpha)) / nrow(tempDat)
  powerExchMP$VCme[row_it] <- length(which(tempDat$VC_MP < tempAlpha)) / nrow(tempDat)
  # SS
  powerExchSS$GBJ[row_it] <- length(which(tempDat$GBJ_SS < tempAlpha)) / nrow(tempDat)
  powerExchSS$GHC[row_it] <- length(which(tempDat$GHC_SS < tempAlpha)) / nrow(tempDat)
  powerExchSS$iBJ[row_it] <- length(which(tempDat$iBJ_SS < tempAlpha)) / nrow(tempDat)
  powerExchSS$iHC[row_it] <- length(which(tempDat$iHC_SS < tempAlpha)) / nrow(tempDat)
  powerExchSS$ACAT[row_it] <- length(which(tempDat$ACAT_SS < tempAlpha)) / nrow(tempDat)
  powerExchSS$VC[row_it] <- length(which(tempDat$skat < tempAlpha)) / nrow(tempDat)
  powerExchSS$VCme[row_it] <- length(which(tempDat$VC_SS < tempAlpha)) / nrow(tempDat)
}


# need to make it a longer data frame for ggplot
powerExchSSlong <- powerExchSS %>%
  select(-VCme) %>%
  pivot_longer(!rho3, names_to="Test", values_to="Power") 
powerExchMPlong <- powerExchMP %>%
  select(-VCme) %>%
  pivot_longer(!rho3, names_to="Test", values_to="Power") 

# plot MP
simExchMP <- ggplot(powerExchMPlong, aes(x=rho3, y=Power, color=Test, linetype=Test)) +
  geom_smooth(method="loess", se=FALSE, span=0.2) +
  xlab(expression(paste(rho[3]))) + ylab("Power") + 
  ylim(c(0, 1)) +  xlim(c(0, 1)) + 
  theme_cowplot() + 
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))


# plot SS
simExchSS <- ggplot(powerExchSSlong, aes(x=rho3, y=Power, color=Test, linetype=Test)) +
  geom_smooth(method="loess", se=FALSE, span=0.2) +
  xlab(expression(paste(rho[3]))) + ylab("Power") + 
  ylim(c(0, 1)) + xlim(c(0, 1)) + 
  theme_cowplot() + 
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) 



# put together
simsPlotsBlock <- plot_grid(simBlockSS +  theme(legend.position="none")  +
                              theme(axis.title = element_text(size=14)),
                            simBlockMP +  theme(legend.position="none")  +
                              theme(axis.title = element_text(size=14)),
                            labels=c("A", "B"), label_size=24, ncol=2)
simsLegend <- get_legend(simBlockSS +
                           theme(legend.direction="vertical", legend.justification="center",
                                 legend.box.just="bottom"))
plot_grid(simsPlotsBlock, simsLegend, ncol=2, rel_widths=c(1, 0.15))


simsPlotsExch <- plot_grid(simExchSS +  theme(legend.position="none") +
                             theme(axis.title = element_text(size=14)),
                           simExchMP +  theme(legend.position="none") +
                             theme(axis.title = element_text(size=14)),
                           labels=c("A", "B"), label_size=24, ncol=2)
plot_grid(simsPlotsExch, simsLegend, ncol=2, rel_widths=c(1, 0.15))



