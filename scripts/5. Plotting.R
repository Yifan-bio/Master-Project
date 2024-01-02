# Venn diagram

temp = list(Enhancer = unique(enhancer_ATAC$gene_name),Promoter = unique(promoter_ATAC$gene_name))
Plot6A = ggvenn(temp,fill_color = c("#AC92EB","#ED5564")) 


#### Plot5G region annotations (Two version) ####

temp = as.data.frame(ATAC)
temp$con = "Intergenic"
temp$con[temp$name %in% intron_ATAC$annot.name] = "Intron"
temp$con[temp$name %in% exon_ATAC$annot.name] = "Exon"
temp$con[temp$name %in% promoter_ATAC$annot.name] = "Promoter"
temp$con[temp$name %in% enhancer_ATAC$annot.name] = "Enhancer"
temp = temp %>% count(con)

Plot6B = ggplot(temp, aes(y = "",x = n, fill = con)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values=c("#AC92EB","#4FC1E8","#A0D568","#FFCE54","#ED5564")) +
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
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
  )





ggpubr::ggarrange(Plot6B,Plot6A,nrow = 2,labels = c("A","B"))
ggsave("Plot.png",dpi = 900,width = 8,height = 8)






x = list(geneHancer = c(1:74440),H3K27ac = c(73019:131150))

ggvenn(x) + 