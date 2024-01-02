# This script is created for the plotting for fastqc/multiqc results
# Last modified: 15 March 2023

library(tidyverse)
library(scales)

#### Plot0A Per base quality ####

temp =  read_tsv("./Input/fastqc_per_base_sequence_quality_plot.tsv")

# Renaming the files
colnames(temp) = c("Pos",paste0("ATAC",1:8),paste0("WGBS",1:4),paste0("RNA",1:8))

# Convert to long dataframe
temp = temp %>% gather(key = "variable",value = "value",-Pos)

# Allowing all same tech to be in one colour
temp$colour = gsub('[[:digit:]]+', '',temp$variable)
temp=na.omit(temp)

# New section to remove WGBS
temp = temp[temp$colour != "WGBS",]

# Plot
plot0A = ggplot(temp,aes(x = Pos,y=value),) +
  geom_line(aes(group = variable,color = colour),show.legend = FALSE) +
  scale_x_continuous(limits = c(0,100.5),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,37),expand = c(0,0)) +
  xlab("Position in read(bp)") +
  ylab("Phred Score") +
  theme_bw() +
  scale_color_manual(values = c("#E47D3A","#6668A3"))
  # scale_color_manual(values = c("#ff7f0e","#1f88b4","#2ca02c"))

ggsave(plot = plot0A,filename = "Plot0A.png",width = 6,height = 3)

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
  scale_y_continuous(labels = unit_format(unit = "",scale = 1e-6),limits = c(0,150000000),expand = c(0,0)) +
  xlab("Phred Score") +
  ylab("Number of reads (Million)") +
  theme_bw() +
  scale_color_manual(values = c("purple","darkolivegreen3","#2ca02c"))

ggsave(plot = plot0B,filename = "Plot0B.png",width = 6,height = 3)

#### Plot0C Plot gc content ####

temp = read_tsv("./Input/fastqc_per_sequence_gc_content_plot.tsv")
# Renaming the files
colnames(temp) = c("Score",paste0("ATAC",1:8),paste0("WGBS",1:4),paste0("RNA",1:8))
# Convert to long dataframe
temp = temp %>% gather(key = "variable",value = "value",-Score)
temp$colour = gsub('[[:digit:]]+', '',temp$variable)

temp = temp[temp$colour != "WGBS",]

plot0C = ggplot(temp,aes(x = Score,y=value)) +
  geom_line(aes(group = variable,color = colour)) +
  scale_x_continuous(limits = c(0,100.5),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,6),expand = c(0,0)) +
  xlab("GC content (%)") +
  ylab("Percentage of reads") +
  labs(color = "") +
  theme_bw() +
  scale_color_manual(values = c("#E47D3A","#6668A3")) +
  theme(legend.text = element_text(size=16),legend.title = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(size = 15)))

legend = ggpubr::get_legend(plot0C)

plot0C = plot0C + guides(color = "none")

ggsave(plot = plot0C,filename = "Plot0C.png",width = 6,height = 3)

#### Output for dissertation ####

ggpubr::ggarrange(plot0A,plot0B)
ggsave(dpi=1200,"Figure3-1.png")

#### Figure need to be replotted ####
# z = data.frame(s = c("Untreat Replicate 1","Untreat Replicate 2","Treated Replicate 1","Treated Replicate 2"),v = c(93.36,93.58,92.73,92.76))
# x = ggplot(data=z, aes(x=v, y=s)) +
#   geom_bar(stat="identity", fill="#6668A3") +
#   theme_bw() + 
#   scale_x_continuous(limits = c(0,100),expand = c(0,0)) +
#   xlab("Mapping percentage (%)") +
#   ylab("") 
# 
# ggpubr::ggarrange(plot0A,plot0C,x,legend,labels = c("A","B","C",""))
# ggpubr::ggarrange(ggpubr::ggarrange(plot0A,plot0C,labels = c("A","B")),
#                   ggpubr::ggarrange(x,legend,labels = c("C",""),widths = c(3,2)),
#                   nrow = 2)
# ggsave(dpi = 1200,height = 4.4,width = 6.3,"x.png")


