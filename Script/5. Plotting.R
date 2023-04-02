# Script for all plotting involved in section 5 differential accessibility analysis
# Last modified: 17 March 2023

library("tidyverse")
library("gridExtra")

#### Plot5A Scatter plot comparing conditions ####

ScatterPlot_preparation = function(DiffBind_df,L2FC,FDR) {
  DiffBind_df$Condition = "Unchange"
  DiffBind_df$Condition[DiffBind_df$Fold > 2 & DiffBind_df$FDR < 0.05] = "Opening"
  DiffBind_df$Condition[DiffBind_df$Fold < -2 & DiffBind_df$FDR < 0.05] = "Closing"
  return(DiffBind_df)
}
x = ScatterPlot_preparation(ATAC_report)


Plot5A = ggplot(data = x,aes(x = Conc_Treated,y=Conc_Control,color = Condition)) +
  geom_point() + 
  scale_color_manual(values=c("#4995C6","#B9181A","#C0C0C0")) + 
  theme_bw() +
  theme(axis.title = element_text(size = 18),axis.text = element_text(size = 13),legend.position = "right",legend.text = element_text(size = 18),legend.title = element_text(size = 18)) +
  labs(x = "Chromatin accessibility after PMA treatment",y = "Chromatin accessibility before PMA treatment") +
  scale_y_continuous(limits = c(0,12.2),expand = c(0,0)) +
  scale_x_continuous(limits = c(0,12.2),expand = c(0,0))

ggsave(dpi = 1200,filename = "Plot5A.png",height = 10,width = 13.1)


#### Plot5B Venn Diagram for DiffBind ####

library("DiffBind")

dba.plotVenn(DBA = ATAC_contrast,
             mask = ATAC_contrast$masks$PMA)

#### Plot5C PCA plot ####

dba.plotPCA(DBA = ATAC_analyze,contrast = 1,label = DBA_CONDITION)

#### Plot5D MA plot ####

dba.plotMA(DBA = ATAC_analyze)

#### Plot5E Box plot ####

dba.plotBox(ATAC_analyze)

#### Plot5F Heatmap ####

dba.plotProfile(ATAC_analyze)

#### Plot5G region annotations ####
Plot5G = plot_annotation(
  annotated_regions = annotated,
  annotation_order = c('hg38_genes_promoters','hg38_genes_introns','hg38_genes_exons','hg38_custom_Genhancer') ,
  x_label = 'Genomic regions',
  y_label = 'Gene count') +
  scale_y_continuous(limits = c(0,5000),expand = c(0,0)) +
  theme(axis.title = element_text(size = 18),axis.text = element_text(size = 15),axis.title.x=element_blank())

ggpubr::ggarrange(Plot5A,NULL,Plot5G,ncol = 1,heights = c(2,0.1,1),labels = c("A","","B"))
ggsave("x.png",width = 13,height = 13.5)
