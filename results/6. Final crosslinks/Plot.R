

#### Plotting voocano plot with DEGs with RE####
RNA$image = "discard"
RNA$image[RNA$Row.names %in% ATAC_RNA_RE$gene_id | RNA$gene_name %in% ATAC_RNA_RE$gene_name] = "label"
RNA = RNA %>% arrange(image)

library(gridExtra)
RNA$padj = ifelse(RNA$padj == 0,1e-303,RNA$padj)
ggplot(RNA,aes(x=log2FoldChange,y=-log10(padj))) + 
  geom_point(aes(color = RNA$image)) +
  theme_classic() +                                                                  # make it bw rather then grey background
  theme(text = element_text(size = 20),legend.position = "none",rect = element_rect(fill = "transparent")) +
  scale_colour_manual(values = c("#C0C0C0","blue4")) +                       # Change the colour of dots
  labs(x = expression("Change in Gene Expression (Log"[2]*" Fold Change)"),y = expression("-log"[10]*"(p adjusted)"),color = "Gene status") + #naming
  scale_y_continuous(expand = c(0,0),limits = c(0,310)) + 
  ggrepel::geom_text_repel(data = subset(RNA,RNA$log2FoldChange > 7.5 & RNA$image == "label" & RNA$baseMean > 100),
                           aes(label = gene_name),
                           nudge_x = 15 - subset(RNA,RNA$log2FoldChange > 7.5 & RNA$image == "label" & RNA$baseMean > 100)$log2FoldChange,
                           segment.color = "grey50",
                           direction     = "y") +
  ggrepel::geom_text_repel(data = subset(RNA,RNA$log2FoldChange < -5 & RNA$image == "label" & RNA$baseMean > 100),
                           aes(label = gene_name),
                           nudge_x = -15 + subset(RNA,RNA$log2FoldChange > 7.5 & RNA$image == "label" & RNA$baseMean > 100)$log2FoldChange,
                           segment.color = "grey50",
                           direction     = "y")
ggsave("VPlot.png",height = 6,width = 8,dpi = 1200)


# New plot
ORA_GO = function(sig_gene,term = 'ALL') {
  genelist = sig_gene$log2FoldChange
  names(genelist) = sig_gene$Row.names
  genelist = sort(genelist,decreasing = TRUE)
  
  gene = names(genelist)[abs(genelist) > 2]
  
  res <- enrichGO(gene = gene, 
                  OrgDb = 'org.Hs.eg.db',
                  ont = term,
                  keyType = 'ENSEMBL',
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)
  return(res)
}
temp = ORA_GO(RNA_sigAcc)
s = c( "ERK1 and ERK2 cascade","leukocyte migration",
       "leukocyte proliferation","cell chemotaxis",
       "regulation of inflammatory response",
       "cellular response to lipopolysaccharide","cellular response to molecule of bacterial origin")
temp1 = as.data.frame(temp)
temp1 = temp1[temp1$Description %in% s,]
plot1 = ggplot(temp1,aes(y = Description,x = Count)) + geom_bar(stat = "identity",fill = "#C0C0C0") +
  theme_bw() + scale_x_continuous(limits = c(0,23),expand = c(0,0)) + 
  theme(axis.title.y = element_blank(),axis.text = element_text(size = 14)) + 
  labs(x = "Gene count")


s = c( "ERK1 and ERK2 cascade","leukocyte migration",
           "leukocyte proliferation",
           "regulation of inflammatory response")
plot2 = cnetplot(x = temp,showCategory = s,circular = TRUE, colorEdge = TRUE,node_label = "gene") + 
  scale_size_continuous(guide = "none") + 
  theme(legend.text = element_text(size = 13),legend.position = "bottom",legend.title = element_blank())
  

ggpubr::ggarrange(plot1,plot2,widths = c(1,2),labels = c("A","B"),nrow = 2,heights =c(1,2))
ggsave("test.png",width = 10,height = 14,dpi = 1200,bg = "white")
