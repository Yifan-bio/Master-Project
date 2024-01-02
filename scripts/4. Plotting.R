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


# Plot5A = ggplot(data = x,aes(x = Conc_Treated,y=Conc_Control,color = Condition)) +
Plot5A = ggplot(data = x,aes(x = Fold,y=-log10(FDR),color = Condition)) +
  geom_point() + 
  scale_color_manual(values=c("#B9181A","#4995C6","#C0C0C0")) + 
  theme_bw() +
  theme(axis.title = element_text(size = 18),axis.text = element_text(size = 13),legend.position = "right",legend.text = element_text(size = 18),legend.title = element_text(size = 18)) +
  labs(x = expression("Change in Chromatin Accessibility (Log"[2]*" Fold Change)"),y = expression("-log"[10]*"(FDR)"),color = "Region status") +
  scale_y_continuous(limits = c(0,75),expand = c(0,0)) +
  scale_x_continuous(limits = c(-7,7),expand = c(0,0))

ggsave(dpi = 1200,filename = "Plot5A.png",height = 8,width = 11.1)


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

#### Plot5G region annotations (Two version) ####
Plot5G = plot_annotation(
  annotated_regions = annotated,
  annotation_order = c('hg38_genes_promoters','hg38_genes_introns','hg38_genes_exons','hg38_custom_Genhancer') ,
  x_label = 'Genomic regions',
  y_label = 'Gene count') +
  scale_y_continuous(limits = c(0,5000),expand = c(0,0)) +
  theme(axis.title = element_text(size = 18),axis.text = element_text(size = 15),axis.title.x=element_blank())

Plot5G_preparation = function(Annotation){
  
  temp = as.data.frame(annotated)
  
  temp = subset_order_tbl(tbl = temp,col = "annot.type",
                          col_order = c('hg38_genes_promoters','hg38_genes_introns','hg38_genes_exons','hg38_custom_Genhancer'))
  temp = dplyr::distinct(dplyr::ungroup(temp),across(c("seqnames", "start", "end", "annot.type")),.keep_all = TRUE)
  temp$annot.type <- factor(temp$annot.type, levels = c("promoters", "Genhancer", "exons", "introns"))
  temp <- temp[order(temp$annot.type), ]
  temp <- temp[!duplicated(temp[, c("seqnames","start", "end")]), ]
  temp <- temp %>% count(annot.type)
  return(temp)
}

Plot5G = Plot5G_preparation(Annotation = annotated)

Plot5G = data.frame(annot.type = c(paste0(Plot5G$annot.type),"intergenic"),n = c(Plot5G$n,(nrow(Diff.bind.sig) - sum(Plot5G$n))))
Plot5G$annot.type[Plot5G$annot.type == "Genhancer"] = "enhancer"

Plot5G = ggplot(Plot5G, aes(y = "",x = n, fill = annot.type)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values=c("pink","#E6CEAC","#A0EEE1","#F4606C","#431A28")) +
  theme_bw() +
  labs(fill = "Annotation Region") +
   theme(axis.line = element_blank(),
         axis.text.x=element_text(size = 12),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="bottom",
         panel.background=element_blank(),
         panel.border=element_blank(),
         legend.text = element_text(size = 18),
         legend.title = element_text(size = 18),
         #panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank()
         )
  
ggsave("Plot5G.png",width = 13,height = 5)

#### Plot 5H ####

plot5H = read.delim("./greatExportAll.tsv")
plot5H = plot5H[1:15,]
plot5H$Desc = str_to_sentence(plot5H$Desc)

plot5H = 
  ggplot(data=plot5H, aes(x=ObsGenes, y=Desc)) +
  geom_bar(stat="identity",fill = "#E47D3A") +
  scale_y_discrete(position = "left") +
  labs(x = "Gene Count") +
  theme(legend.position="none",
        axis.title.y=element_blank(), 
        axis.title.x=element_text(size = 13),
        axis.text.y= element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.line=element_line(color="gray"),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
  )
ggsave("plot5H.png",height = 4,width = 9,dpi = 1200)
