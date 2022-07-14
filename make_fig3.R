# Recreate Figure 3
# Will need to first generate the data (using provided scripts).
# Once data has been generated, use this file to recreate the plots

# Change to current working directory
setwd('/users/rsun3/desktop/reproduce')

library(dplyr)
library(tidyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)


# load block MP files
allFnames <- list.files()
blockMPinnFnames <- allFnames[which(substr(allFnames, 1, 37) == "bounds_theory_MP_Inn_SS25_MP40_S1_aID")]
blockMPinn <- fread(blockMPinnFnames[1])
for (file_it in 2:length(blockMPinnFnames)) {
  tempRes <- fread(blockMPinnFnames[file_it])
  blockMPinn <- rbind(blockMPinn, tempRes)
}

allFnames <- list.files()
blockMPgenFnames <- allFnames[which(substr(allFnames, 1, 37) == "bounds_theory_MP_Gen_SS25_MP40_S1_aID")]
blockMPgen <- fread(blockMPgenFnames[1])
for (file_it in 2:length(blockMPgenFnames)) {
  tempRes <- fread(blockMPgenFnames[file_it])
  blockMPgen <- rbind(blockMPgen, tempRes)
}

# put them together, pivot longer
blockMP <- blockMPgen %>% select(-innLower, -innUpper) %>%
  merge(., blockMPinn %>% select(-genLower, -genUpper), by="rho3") %>%
  arrange(rho3) %>%
  set_colnames(c("rho3", "GBJ lower", "GBJ upper", "iBJ lower", "iBJ upper")) %>%
  pivot_longer(., cols=!(rho3), names_to = "Bound")


# plot
boundsBlockMP <- ggplot(blockMP, aes(x=rho3, y=value, color=Bound, linetype=Bound)) +
  geom_smooth(method="loess", se=FALSE, span=0.2) +
  xlab(expression(paste(rho[3]))) + ylab("Power (Multiple Outcomes)") +
  ylim(c(0, 1)) + xlim(c(0, 0.59)) +
  scale_linetype_manual(values=c("solid", "solid", "dotted", "dotted")) +
  scale_color_manual(values=c("black", "orange", "black", "orange")) +
  theme_cowplot() +
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))


# load block SS files
allFnames <- list.files()
blockSSinnFnames <- allFnames[which(substr(allFnames, 1, 37) == "bounds_theory_SS_Inn_SS25_MP40_S1_aID")]
blockSSinn <- fread(blockSSinnFnames[1])
for (file_it in 2:length(blockSSinnFnames)) {
  tempRes <- fread(blockSSinnFnames[file_it])
  blockSSinn <- rbind(blockSSinn, tempRes)
}

allFnames <- list.files()
blockSSgenFnames <- allFnames[which(substr(allFnames, 1, 37) == "bounds_theory_SS_Gen_SS25_MP40_S1_aID")]
blockSSgen <- fread(blockSSgenFnames[1])
for (file_it in 2:length(blockSSgenFnames)) {
  tempRes <- fread(blockSSgenFnames[file_it])
  blockSSgen <- rbind(blockSSgen, tempRes)
}

# put them together, pivot longer
blockSS <- blockSSgen %>% select(-innLower, -innUpper) %>%
  merge(., blockSSinn %>% select(-genLower, -genUpper), by="rho3") %>%
  arrange(rho3) %>%
  set_colnames(c("rho3", "GBJ lower", "GBJ upper", "iBJ lower", "iBJ upper")) %>%
  pivot_longer(., cols=!(rho3), names_to = "Bound")


# plot
boundsBlockSS <- ggplot(blockSS, aes(x=rho3, y=value, color=Bound, linetype=Bound)) +
  geom_smooth(method="loess", se=FALSE, span=0.2) +
  xlab(expression(paste(rho[3]))) + ylab("Power (Multiple Exp. Factors)") +
  ylim(c(0, 1)) + xlim(c(0, 0.59)) +
  scale_linetype_manual(values=c("solid", "solid", "dotted", "dotted")) +
  scale_color_manual(values=c("black", "orange", "black", "orange")) +
  theme_cowplot() +
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))


# load exch MP files
allFnames <- list.files()
exchMPinnFnames <- allFnames[which(substr(allFnames, 1, 37) == "bounds_theory_MP_Inn_SS25_MP40_S3_aID")]
exchMPinn <- fread(exchMPinnFnames[1])
for (file_it in 2:length(exchMPinnFnames)) {
  tempRes <- fread(exchMPinnFnames[file_it])
  exchMPinn <- rbind(exchMPinn, tempRes)
}

allFnames <- list.files()
exchMPgenFnames <- allFnames[which(substr(allFnames, 1, 37) == "bounds_theory_MP_Gen_SS25_MP40_S3_aID")]
exchMPgen <- fread(exchMPgenFnames[1])
for (file_it in 2:length(exchMPgenFnames)) {
  tempRes <- fread(exchMPgenFnames[file_it])
  exchMPgen <- rbind(exchMPgen, tempRes)
}

# put them together, pivot longer
exchMP <- exchMPgen %>% select(-innLower, -innUpper) %>%
  merge(., exchMPinn %>% select(-genLower, -genUpper), by="rho3") %>%
  arrange(rho3) %>%
  set_colnames(c("rho3", "GBJ lower", "GBJ upper", "iBJ lower", "iBJ upper")) %>%
  pivot_longer(., cols=!(rho3), names_to = "Bound")

# plot
boundsExchMP <- ggplot(exchMP, aes(x=rho3, y=value, color=Bound, linetype=Bound)) +
  geom_smooth(method="loess", se=FALSE, span=0.2) +
  xlab(expression(paste(rho[1], "=", rho[2], "=", rho[3]))) + ylab("Power (Multiple Outcomes)") +
  ylim(c(0, 1)) + xlim(c(0, 0.9)) +
  scale_linetype_manual(values=c("solid", "solid", "dotted", "dotted")) +
  scale_color_manual(values=c("black", "orange", "black", "orange")) +
  theme_cowplot() +
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

# load exch SS files
allFnames <- list.files()
exchSSinnFnames <- allFnames[which(substr(allFnames, 1, 37) == "bounds_theory_SS_Inn_SS25_MP40_S3_aID")]
exchSSinn <- fread(exchSSinnFnames[1])
for (file_it in 2:length(exchSSinnFnames)) {
  tempRes <- fread(exchSSinnFnames[file_it])
  exchSSinn <- rbind(exchSSinn, tempRes)
}

allFnames <- list.files()
exchSSgenFnames <- allFnames[which(substr(allFnames, 1, 37) == "bounds_theory_SS_Gen_SS25_MP40_S3_aID")]
exchSSgen <- fread(exchSSgenFnames[1])
for (file_it in 2:length(exchSSgenFnames)) {
  tempRes <- fread(exchSSgenFnames[file_it])
  exchSSgen <- rbind(exchSSgen, tempRes)
}

# put them together, pivot longer
exchSS <- exchSSgen %>% select(-innLower, -innUpper) %>%
  merge(., exchSSinn %>% select(-genLower, -genUpper), by="rho3") %>%
  arrange(rho3) %>%
  set_colnames(c("rho3", "GBJ lower", "GBJ upper", "iBJ lower", "iBJ upper")) %>%
  pivot_longer(., cols=!(rho3), names_to = "Bound")

# plot
boundsExchSS <- ggplot(exchSS, aes(x=rho3, y=value, color=Bound, linetype=Bound)) +
  geom_smooth(method="loess", se=FALSE, span=0.2) +
  xlab(expression(paste(rho[1], "=", rho[2], "=", rho[3]))) + ylab("Power (Multiple Exp. Factors)") +
  ylim(c(0, 1)) + xlim(c(0, 0.9)) +
  scale_linetype_manual(values=c("solid", "solid", "dotted", "dotted")) +
  scale_color_manual(values=c("black", "orange", "black", "orange")) +
  theme_cowplot() +
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

#---------------------------------------------------------------------------#
# put bounds together

# save
boundsPlotsBlock <- plot_grid(boundsBlockSS +  theme(legend.position="none")  +
                                theme(axis.title = element_text(size=14)),
                           boundsBlockMP +  theme(legend.position="none") +
                             theme(axis.title = element_text(size=14)),
                           labels=c("A", "B"), label_size=24, ncol=2)
boundsLegend <- get_legend(boundsBlockSS +
                               theme(legend.direction="vertical", legend.justification="center",
                                     legend.box.just="bottom"))

plot_grid(boundsPlotsBlock, boundsLegend, ncol=2, rel_widths=c(1, 0.2))


boundsPlotsExch <- plot_grid(boundsExchSS +  theme(legend.position="none")  +
                               theme(axis.title = element_text(size=14)),
                              boundsExchMP +  theme(legend.position="none")  +
                               theme(axis.title = element_text(size=14)),
                             labels=c("C", "D"), label_size=24, ncol=2)
plot_grid(boundsPlotsExch, boundsLegend, ncol=2, rel_widths=c(1, 0.2))


