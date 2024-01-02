# This script is created for the plotting for fastqc/multiqc results
# Last modified: 15 March 2023

library(tidyverse)

#### Plot0A Per base quality ####

temp =  read_tsv("./Input/fastqc_per_base_sequence_quality_plot.tsv")

# Renaming the files
colnames(temp) = c("Pos",paste0("ATAC",1:8),paste0("WGBS",1:4),paste0("RNA",1:8))

# Convert to long dataframe
temp = temp %>% gather(key = "variable",value = "value",-Pos)

# Allowing all same tech to be in one colour
temp$colour = gsub('[[:digit:]]+', '',temp$variable)
temp=na.omit(temp)

# Plot
plot0A = ggplot(temp,aes(x = Pos,y=value),) +
  geom_line(aes(group = variable,color = colour),show.legend = FALSE) +
  scale_x_continuous(limits = c(0,150),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,41),expand = c(0,0)) +
  xlab("Position in read(bp)") +
  ylab("Phred Score") +
  theme_bw() +
  scale_color_manual(values = c("#1f88b4","#ff7f0e","#2ca02c"))

#### Plot0B Per sequence quality####

# Per sequence quality
temp = read_tsv("./Input/fastqc_per_sequence_quality_scores_plot.tsv")
# Renaming the files
colnames(temp) = c("Score",paste0("ATAC",1:8),paste0("WGBS",1:4),paste0("RNA",1:8))
# Convert to long dataframe
temp = temp %>% gather(key = "variable",value = "value",-Score)
temp$colour = gsub('[[:digit:]]+', '',temp$variable)

plot0B = ggplot(temp,aes(x = Score,y=value)) +
  geom_line(aes(group = variable,color = colour),show.legend = FALSE) +
  scale_x_continuous(limits = c(0,40),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,150000000),expand = c(0,0)) +
  xlab("Phred Score") +
  ylab("Number of reads") +
  theme_bw() +
  scale_color_manual(values = c("#1f88b4","#ff7f0e","#2ca02c"))



#### Plot for Supplementary ####
ggpubr::ggarrange(plot0A,plot0B)
ggsave(dpi=900,"fastqc.png")