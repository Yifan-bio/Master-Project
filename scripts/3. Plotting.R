# Script for all plotting invovled in section 4 differential expression analysis
# Last modified: 17 March 2023

library(tidyverse)
library(gridExtra)

#### Plot4A Volcano Plot ####

# Modifying the file so a pretty volcano plot can be made
Volcano_plot_modification = function(resdata,L2FC,padj_value){
  
  # Adding the conditions for adding of colours
  resdata$con = "non-significant"
  resdata$con[resdata$log2FoldChange > L2FC & resdata$padj < padj_value] = "upregulated"
  resdata$con[resdata$log2FoldChange < -(L2FC) & resdata$padj < padj_value] = "downregulated"
  
  # Deseq2 reports all padj < 1e-303 as 0 so need to be corrected
  resdata$padj = ifelse(resdata$padj == 0,1e-303,resdata$padj)
  
  return(resdata)
}

resdata = Volcano_plot_modification(resdata = resdata,L2FC = 2,padj_value = 0.05)

# Plotting
ggplot(resdata,aes(x=log2FoldChange,y=-log10(padj),color = con)) + geom_point() +
  theme_bw() +                                                                  
  theme(text = element_text(size = 20),legend.position = "none",rect = element_rect(fill = "transparent")) +
  scale_colour_manual(values = c("#B9181A","#C0C0C0","#4995C6")) +                       # Change the colour of dots
  labs(x = expression("Change in Gene Expression (Log"[2]*" Fold Change)"),y = expression("-log"[10]*"(p adjusted)"),color = "Gene status") + #naming
  scale_y_continuous(expand = c(0,0),limits = c(0,310))                                         # Remove the gap between the resdatas and x margin

ggsave("Plot4A.png",width = 9,height = 5,dpi = 1200)

# Cleaning up
remove(Volcano_plot_modification)

#### Plot4B Upset Plot ####

library(ggupset)

UpsetPlot_preparation = function(Enrich_result,GO_terms) {
  
  # Selecting the terms to be used in upset plot
  temp = as.data.frame(Enrich_result)
  temp = temp[temp$Description %in% GO_terms,]
  
  # Splitting the genes into separate columns
  temp_max = max(temp$Count)
  temp[c(paste0("V",1:temp_max))] <- str_split_fixed(temp$geneID,"/",temp_max)
  temp = temp %>% select(Description,V1:paste0("V",temp_max))
  
  # Converting the format into upsetplot required format
  temp = temp %>% pivot_longer(!Description,names_to = "name",values_to = "value")
  temp = temp %>% filter(!value == "") # Remove some noises
  temp = temp %>% select(!name)
  temp = temp %>% group_by(value) %>% summarise(Pathways = list(Description))
  
  return(temp)
}

# Go term for plotting upsetplot
GO_terms = c("inflammatory response","leukocyte migration",
             "nucleosome assembly","nuclear division",
             "regulation of mitotic cell cycle","chemotaxis",
             "DNA replication-dependent chromatin assembly")

x = UpsetPlot_preparation(Enrich_result = ORA,GO_terms = GO_terms)

# PLotting
ggplot(temp,aes(x=Pathways)) + 
  geom_bar() + 
  theme_bw() +
  scale_x_upset() +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())

ggsave(filename = "Plot4B.png",plot = x,dpi = 1200,width = 6.8,height = 4)

#### Plot4C Transcript TPM distribution ####

x = txi.salmon[["abundance"]]
x = data.frame(x)
colnames(x) = c("Untreat_1","Untreat_2","Treat_1","Treat_2")

# Converting to long dataframe and remove NULL values
x$gene = row.names(x)
x = x %>% pivot_longer(!gene,names_to = "replicate",values_to = "abundance")
x = x %>% filter(abundance != 0)

# Prepring for facet
x$condition = (do.call('rbind', strsplit(as.character(x$replicate),'_',fixed=TRUE)))[,1]
x$rep = (do.call('rbind', strsplit(as.character(x$replicate),'_',fixed=TRUE)))[,2]

ggplot(x,aes(x = log10(abundance),y = -0.02)) +
  # Horizontal boxplot
  ggstance::geom_boxploth(aes(fill = condition),width = 0.03) +
  # density plot
  geom_density(aes(x = log10(abundance)),inherit.aes = F) +
  facet_grid(rows = vars(rep),cols = vars(condition)) +
  scale_fill_discrete() +
  geom_hline(yintercept = 0)+
  labs(y = "Proportion of genes",x = "log10(TPM)")

ggsave("Plot4C.png",width = 7,height = 5)

#### Plot4D Biomarkers ####

x = resdata_2 %>% select(gene_name,V1:V4)
colnames(x) = c("Gene","Untreat1","Untreat2","Treated1","Treated2")

# Select gene to plot
x = x[x$Gene %in% c("CD14","CSF1","ICAM1","ITGAM","SPP1","CSF1R"),]
x = x %>% pivot_longer(!Gene,names_to = "Rep",values_to = "abundance")
x$rep = gsub('[[:digit:]]+', '',x$Rep)
x$mark[x$Gene %in% c("ICAM1","ITGAM","SPP1")] = "Migration marker gene"
x$mark[x$Gene %in% c("CD14","CSF1","CSF1R")] = "Macrophage activation gene"
# x = aggregate(x$abundance, list(x$Gene,x$rep), FUN=mean) 
# 
# x$mark[x$Group.1 %in% c("ICAM1","ITGAM","SPP1")] = "Migration marker gene"
# x$mark[x$Group.1 %in% c("CD14","CSF1","CSF1R")] = "Macrophage activation gene"
# 
# ggplot(x,aes(x = Group.2,y = x)) + 
#   geom_bar(stat = "identity") +
#   facet_wrap(~ Group.1,scales = "free") + 
#   theme_bw() +
#   labs(x=NULL,y= "Normalised Gene Count") +
#   theme(axis.title.y = element_text(size = 16),axis.text = element_text(size = 14)) 


ggplot(x,aes(x = reorder(rep,+abundance),y = abundance)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge",show.legend = FALSE) +
 #  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
  # ggsignif::stat_signif(comparisons=list(c("Treated", "Untreat")),
  #                       test="t.test",
  #                       map_signif_level = T,
  #                       annotations = c("*")) +  
  # ggsignif::geom_signif(comparisons=list(c("Treated", "Untreat")),
  #                       test="t.test",
  #                       annotations = c("*")) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", ref.group = "Untreat",
                     hide.ns = T,vjust = 0.5,
                     symnum.args = list(cutpoints = c(0,0.3,Inf),symbols = c("*","ns"))) +
  facet_nested_wrap(mark ~ Gene,
                    scales = "free",
                    strip =  strip_nested(background_x = elem_list_rect(fill = c("white")),
                                          text_x = element_text(size = 12),
                                          text_y = element_text(size = 12))) +

  theme_bw() +
  labs(x=NULL,y= "Normalised Gene Count") +
  theme(axis.title.y = element_text(size = 16),axis.text = element_text(size = 14))

ggsave("Plot4D.png",width = 7.5,height = 6.5,dpi = 1200)

#### Plot4E baseMean distribution ####

ggplot(resdata_2,aes(x = log2(baseMean))) + 
  geom_histogram(bins = 250) + 
  theme_bw() + 
  geom_vline(xintercept = log2(20),show.legend = T) +
  scale_x_continuous(expand = c(0.01,0)) +
  scale_y_continuous(expand = c(0,0))

ggsave("Plot4E.png",width = 8,height = 3.3)

#### Plot4F Barplot ####

plot4F = ORA@result
plot4F = plot4F[1:15,]
plot4F$Description = str_to_sentence(plot4F$Description)

plot4F = 
  ggplot(data= plot4F, aes(x= Count, y= Description)) +
  geom_bar(stat="identity",fill = "#6668A3") +
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

ggsave("plot4F.png",height = 4,width = 9,dpi = 1200)

#### Plot for Main Script ####

# P