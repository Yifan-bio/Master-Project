# Script for all plotting involved in section 5 differential accessibility analysis
# Last modified: 17 March 2023

library("tidyverse")
library("gridExtra")

#### Plot5A Scatter plot comparing conditions ####

ScatterPlot_preparation = function(DiffBind_df,L2FC,FDR) {
  DiffBind_df$Condition = "Unchange"
  DiffBind_df$Condition[DiffBind_df$Fold > 2 & DiffBind_df$FDR < 0.05] = "Opening"
  DiffBind_df$Condition[DiffBind_df$Fold < -2 & DiffBind_df$FDR < 0.05] = "closing"
}
x = ATAC_report


ggplot(data = x,aes(x = Conc_Treated,y=Conc_Control,color = Condition)) +
  geom_point() + 
  scale_color_manual(values=c("#4995C6","#B9181A","#C0C0C0")) + 
  theme_bw() +
  theme(axis.title = element_text(size = 18),axis.text = element_text(size = 13),legend.position = "right") +
  labs(x = "Chromatin accessibility after PMA treatment",y = "Chromatin accessibility before PMA treatment") +
  scale_y_continuous(limits = c(0,12.2),expand = c(0,0)) +
  scale_x_continuous(limits = c(0,12.2),expand = c(0,0))

ggsave(dpi = 1200,filename = "Plot5A.png",height = 10,width = 13.1)


#### Plot5B Venn Diagram for DiffBind ####

library("DiffBind")

dba.plotVenn(DBA = ATAC_contrast,mask = ATAC_contrast$masks$PMA)

#### Other plots ####
# Barplot for peak overlaps
x = data.frame(name = c("Unique","Two replicate","Three replicate","All replicate"),number = olap.rate)
ggplot(x,aes(x = name,y=number)) + geom_bar(stat = "identity") +
  labs(x="Number of Overlap Across Peaksets",y="Peak Count") +
  theme_bw() +
  scale_y_continuous(expand = c(0,0),limits = c(0,25000)) +
  theme(axis.text = element_text(size=15),axis.title = element_text(size=18))

ggsave("Overlap_allcon.png",dpi = 900,height = 10,width = 17)
