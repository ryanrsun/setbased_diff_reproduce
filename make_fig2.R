# plot figure 2 using data from gen_data_fig2.R

library(tidyr)
library(magrittr)
library(ggplot2)
library(ggpattern)

setwd('/users/rsun3/desktop')
allPowerCalc <- fread("allPowerCalc.txt")

# high, SS plot
plotDat5SS <- allPowerCalc %>%
  filter(Cor == "high" & Mult == "SS") %>%
  select(-Cor, -Mult, -iBJlowerUnion) %>%
  pivot_longer(., cols=!(location), names_to="Bound", values_to="Power")
plot5SS <- ggplot(plotDat5SS, aes(x=location, y=Power, fill = Bound, pattern = Bound)) +
  geom_bar_pattern(stat = "identity",
                   position = position_dodge(),
                   color = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_fill_manual(values = colorRampPalette(c("#0066CC","#FFFFFF","#FF8C00"))(6)[c(3,2,1,4,5,6)]) +
  scale_pattern_manual(values = c(iBJlower = "stripe", iBJupper = "stripe", iBJpower="stripe",
                                  GBJlower = "none", GBJupper = "none", GBJpower = "none")) +
  ylim(c(0, 1)) +
  xlab("Signal Location - High Correlation, Multiple Explanatory Factors") +
  theme_cowplot() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

#-------------------------------------------------------------------------------------#

# low, SS
plotDat0SS <-  allPowerCalc %>%
  filter(Cor == "low" & Mult == "SS") %>%
  select(-Cor, -Mult, -iBJlowerUnion) %>%
  pivot_longer(., cols=!(location), names_to="Bound", values_to="Power")
plot0SS <- ggplot(plotDat0SS, aes(x=location, y=Power, fill = Bound, pattern = Bound)) +
  geom_bar_pattern(stat = "identity",
                   position = position_dodge(),
                   color = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_fill_manual(values = colorRampPalette(c("#0066CC","#FFFFFF","#FF8C00"))(6)[c(3,2,1,4,5,6)]) +
  scale_pattern_manual(values = c(iBJlower = "stripe", iBJupper = "stripe", iBJpower="stripe",
                                  GBJlower = "none", GBJupper = "none", GBJpower = "none")) +
  ylim(c(0, 1)) +
  xlab("Signal Location - Low Correlation, Multiple Explanatory Factors") +
  theme_cowplot() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

#-------------------------------------------------------------------------------------#
# high, MP

plotDat5MP <-  allPowerCalc %>%
  filter(Cor == "high" & Mult == "MP") %>%
  select(-Cor, -Mult, -iBJlowerUnion) %>%
  pivot_longer(., cols=!(location), names_to="Bound", values_to="Power")
plot5MP <- ggplot(plotDat5MP, aes(x=location, y=Power, fill = Bound, pattern = Bound)) +
  geom_bar_pattern(stat = "identity",
                   position = position_dodge(),
                   color = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_fill_manual(values = colorRampPalette(c("#0066CC","#FFFFFF","#FF8C00"))(6)[c(3,2,1,4,5,6)]) +
  scale_pattern_manual(values = c(iBJlower = "stripe", iBJupper = "stripe", iBJpower="stripe",
                                  GBJlower = "none", GBJupper = "none", GBJpower = "none")) +
  ylim(c(0, 1)) +
  xlab("Signal Location - High Correlation, Multiple Outcomes") +
  theme_cowplot() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

#-------------------------------------------------------------------------------------#
# low, MP

plotDat0MP <- allPowerCalc %>%
  filter(Cor == "low" & Mult == "MP") %>%
  select(-Cor, -Mult, -iBJlowerUnion) %>%
  pivot_longer(., cols=!(location), names_to="Bound", values_to="Power")
plot0MP <- ggplot(plotDat0MP, aes(x=location, y=Power, fill = Bound, pattern = Bound)) +
  geom_bar_pattern(stat = "identity",
                   position = position_dodge(),
                   color = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_fill_manual(values = colorRampPalette(c("#0066CC","#FFFFFF","#FF8C00"))(6)[c(3,2,1,4,5,6)]) +
  scale_pattern_manual(values = c(iBJlower = "stripe", iBJupper = "stripe", iBJpower="stripe",
                                  GBJlower = "none", GBJupper = "none", GBJpower = "none")) +
  ylim(c(0, 1)) +
  xlab("Signal Location - Low Correlation, Multiple Outcomes") +
  theme_cowplot() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

#-------------------------------------------------------------------------------------#
# put it all together

boundsPlot <- plot_grid(plot0SS + theme(legend.position="none"),
                      plot0MP + theme(legend.position="none"),
                      plot5SS + theme(legend.position="none"),
                      plot5MP + theme(legend.position="none"),
                      labels=c("A", "B", "C", "D"), label_size=24, ncol=2)
boundsPlot_legend <- get_legend(plot0MP+
                               theme(legend.direction="horizontal", legend.justification="center",
                                     legend.box.just="bottom"))
plot_grid(boundsPlot, boundsPlot_legend, ncol=1, rel_heights=c(1, 0.2))


