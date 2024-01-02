# This script were created for plots for combining accessibility with regulators or DEGs
# Last modified: 17 March 2023

library(tidyverse)
library(gridExtra)
library(ggrepel)

#### Plot7A Voocano plot of DEGs with regulators with accessibility changes ####

Plot7A_preparation = function(RNA,RNA_with_regulators) {
  
  # Labeling the point for colour
  temp = RNA
  temp$image = "Insig"
  temp$image[temp$gene_name %in% RNA_with_regulators$gene_name] = "label"
  temp$image[abs(temp$log2FoldChange) < 2 | RNA$padj > 0.05] = "Insig"
  temp = temp %>% arrange(image)

  # Fixing some bias from DESeq2
  temp$padj = ifelse(temp$padj == 0,1e-303,temp$padj)

  return(temp)
}

temp = Plot7A_preparation(RNA = RNA,RNA_with_regulators = Reg_ATAC)#[Reg_ATAC$name == "enhancer",])

# Plotting
plot7A = ggplot(temp,aes(x=log2FoldChange,y=-log10(padj))) + 
  geom_point(aes(color = image)) +
  theme_bw() + #theme_classic() +                                                                  # make it bw rather then grey background
  theme(text = element_text(size = 20),legend.position = "none",rect = element_rect(fill = "transparent")) +
  scale_colour_manual(values = c("#C0C0C0","#4995C6")) +                       # Change the colour of dots
  labs(x = expression("Change in Gene Expression (Log"[2]*" Fold Change)"),y = expression("-log"[10]*"(p adjusted)"),color = "Gene status") + #naming
  scale_y_continuous(expand = c(0,0),limits = c(0,310))  +
  ggrepel::geom_text_repel(data = subset(temp,temp$log2FoldChange > 7.5 & temp$image == "label" & temp$baseMean > 100),
                           aes(label = gene_name),
                           nudge_x = 15 - subset(temp,temp$log2FoldChange > 7.5 & temp$image == "label" & temp$baseMean > 100)$log2FoldChange,
                           segment.color = "grey50",
                           direction     = "y") +
  ggrepel::geom_text_repel(data = subset(temp,temp$log2FoldChange < -5 & temp$image == "label" & temp$baseMean > 100),
                           aes(label = gene_name),
                           nudge_x = -15 + subset(temp,temp$log2FoldChange > 7.5 & temp$image == "label" & temp$baseMean > 100)$log2FoldChange,
                           segment.color = "grey50",
                           direction     = "y")
ggsave(filename = "Plot7A.png",plot = plot7A,height = 6,width = 10,dpi = 1200)

#### Plot7B Bar plot of ORA ####

GO_term = c("inflammatory response","positive regulation of cell migration",
            "leukocyte cell-cell adhesion","chemotaxis",
            "ERK1 and ERK2 cascade")

plot7B = heatplot(ORA_RNA_sigAcc,showCategory = GO_term) + 
           theme(panel.grid.minor.x = element_line(color = 2,
                                                   size = 0.25,
                                                   linewidth = 1))

Plot7B_preparation = function(GO_term,ORA){
  
  temp = as.data.frame(ORA)
  temp = temp[temp$Description %in% GO_term,]
  return(temp)
}

temp = Plot7B_preparation(GO_term = GO_term,ORA = ORA_RNA_sigAcc)

plot7B = ggplot(temp,aes(y = Description,x = Count)) + geom_bar(stat = "identity",fill = "#C0C0C0") +
  theme_bw() + scale_x_continuous(limits = c(0,40),expand = c(0,0)) + 
  theme(axis.title.y = element_blank(),axis.text = element_text(size = 14)) + 
  labs(x = "Gene count")

ggsave(filename = "Plot7B.png",plot = plot7B,height = 3,width = 10,dpi = 1200)

#### Plot7C cnetplot for ORA####

GO_term = c("inflammatory response","positive regulation of cell migration",
            "leukocyte cell-cell adhesion")

plot7C = cnetplot(x = ORA_RNA_sigAcc,showCategory = GO_term,circular = TRUE, colorEdge = TRUE,node_label = "gene") + 
  scale_size_continuous(guide = "none") + 
  theme(legend.text = element_text(size = 15),legend.position = "bottom",legend.title = element_blank())

ggsave(filename = "Plot7C.png",plot = plot7C,height = 8,width = 9,dpi = 1200)


#### Plot7D heatmap ####

temp = merge(RNA_sigAcc,ATAC_RNA_RE,by.x = "gene_name",by.y = "gene_name")
temp = temp %>% select(gene_name,log2FoldChange,annot.Fold,con.y)
temp = temp %>% distinct(gene_name,log2FoldChange,annot.Fold,con.y)

GO_term = c("ERK1 and ERK2 cascade","leukocyte migration",
            "leukocyte proliferation",
            "regulation of inflammatory response")

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
  
  return(temp)
}
x = UpsetPlot_preparation(Enrich_result = ORA_RNA_sigAcc,GO_terms = GO_term)

temp = temp[temp$gene_name %in% x$value,]
temp$n = paste0(temp$gene_name,temp$con.y)
temp = temp %>% group_by(n) %>% mutate(Fold = mean(annot.Fold)) %>% ungroup()
temp$Fold = ifelse(abs(temp$Fold) < 2, temp$annot.Fold, temp$Fold)
temp = temp %>% distinct(gene_name,log2FoldChange,Fold,con.y)
temp = temp %>% pivot_longer(names_to = "gene",values_to = "L2FC",cols = -gene_name)

write_excel_csv(temp,"test.csv")

temp = read_delim("test.csv",delim = ";")
temp = temp %>% pivot_longer(cols = -gene_name)
temp$value[temp$value < -2] = -2
temp$value[temp$value > 2] = 2
ggplot(temp, aes(name,gene_name, fill = value)) +
  theme_bw() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  geom_tile(color = "black") +
  scale_fill_gradient(high = "#4995C6", low = "#B9181A") +
  labs(x = NULL,y=NULL) +
  theme(legend.position = "none")
# temp = temp %>% dplyr::arrange(desc(con.y)) %>% distinct(gene_name,log2FoldChange,annot.Fold,.keep_all = T)

#### Plot 7E ####

temp = list(Enhancer = unique(ATAC_RNA_RE$gene_name[ATAC_RNA_RE$con == "enhancer"]),promoter = unique(ATAC_RNA_RE$gene_name[ATAC_RNA_RE$con == "promoter"]))
ggvenn(temp)


#### Plot 7F Circular heatmap ####

library(circlize)
library(ComplexHeatmap)

# Start
x = data.frame(gene = gene,promoter = "0",enhancer = "0",expression = "0")

# Test if the known region is accessible in at least one condition
x$promoter[x$gene %in% promoter_peak$gene_name] = "2"
x$enhancer[x$gene %in% enhancer_peak$gene_name] = "2"
# If the region is experiencing a increase in accessibility
x$enhancer[x$gene %in% enhancer_ATAC$gene_name[enhancer_ATAC$annot.Fold > 2]] = "4"
x$promoter[x$gene %in% promoter_ATAC$gene_name[promoter_ATAC$annot.Fold > 2]] = "4"
# If the region is experiencing a decrease in accessibility
x$enhancer[x$gene %in% enhancer_ATAC$gene_name[enhancer_ATAC$annot.Fold < -2]] = "-4"
x$promoter[x$gene %in% promoter_ATAC$gene_name[promoter_ATAC$annot.Fold < -2]] = "-4"
# If the region contains a increase and a decrease in the cis-element
x$enhancer[x$gene %in% enhancer_ATAC$gene_name[enhancer_ATAC$annot.Fold > 2] & x$gene %in% enhancer_ATAC$gene_name[enhancer_ATAC$annot.Fold < -2]] = "-2"
x$promoter[x$gene %in% promoter_ATAC$gene_name[promoter_ATAC$annot.Fold > 2] & x$gene %in% promoter_ATAC$gene_name[promoter_ATAC$annot.Fold < -2]] = "-2"

# If the gene is upregulated
x$expression[x$gene %in% sig$gene_name[sig$log2FoldChange > 2]] = "4"
# If gene down regulated
x$expression[x$gene %in% sig$gene_name[sig$log2FoldChange < -2]] = "-4"
x$expression[x$gene %in% c("CCL3L1")] = "4"

# Labeling the sectors
x$group = "C"
x$group[x$promoter == "0" | x$enhancer == "0"] = "B"
x$group[x$enhancer == "-2"] = "A"

x$group[x$enhancer == "4" & x$expression == "-4"] = "A"
x$group[x$enhancer == "-4" & x$expression == "4"] = "A"
x$group[x$promoter == "-4" & x$expression == "4"] = "A"
x$group[x$promoter == "4" & x$expression == "-4"] = "A"

x$promoter = as.numeric(x$promoter)
x$enhancer = as.numeric(x$enhancer)
x$expression = as.numeric(x$expression)
x$group = as.character(x$group)
x = x %>% column_to_rownames("gene")
split = factor(toupper(x$group), levels = toupper(letters[1:3]))
data = data.matrix(x)

png("./figure7F.png",width = 20,height = 25,units = "cm",res = 2000)

circos.clear()
circos.par(gap.after = c(2,2,38))

circos.heatmap(data[,1:2],
               col = colorRamp2(c(-4,-2,0,2,4), c( "#B9181A", "goldenrod","black","lavender","#4995C6")),
               rownames.side = "outside",
               rownames.cex = 1.3,
               cell.border = "black",
               split = split)

circos.heatmap(data[,3],
               col = colorRamp2(c(-2,2), c("#B9181A", "#4995C6")),
               cell.border = "black",
               split = spilt)

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + convert_y(2, "mm"), 
              paste0("Group ",CELL_META$sector.index),
              facing = "bending.inside", cex = 1.4,
              adj = c(0.5, -11))
})


circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 3) { # the last sector
    circos.text(x = CELL_META$cell.xlim[2] + convert_x(9, "mm"),
                y = 0.5,
                labels = "Expression", 
                cex = 0.8, 
                facing = "bending.inside")
    circos.text(x = CELL_META$cell.xlim[2] + convert_x(9.4, "mm"),
                y = 1.3,
                labels = "Enahncer", 
                cex = 0.8, 
                facing = "bending.inside")
    circos.text(x = CELL_META$cell.xlim[2] + convert_x(9.8, "mm"), 
                y = 1.9,
                labels = "Promoter", 
                cex = 0.8,
                facing = "bending.inside")
    
  }
}, bg.border = NA)

legend(x = "bottom",
       inset = 0,
       title = NULL,
       legend = c("Multiple enhancer with opposite change","No enhancer / Closed promoter","Accessible promoter","Decrease expression / Chromatin closing","Increase expression / Chromatin opening"),
       fill = c( "goldenrod","black","lavender","#B9181A","#4995C6"),
       cex = 1.1,
       box.lty=0,
       xpd = TRUE,
       ncol = 2,bg = NA)

dev.off()

#### Export for Main Script ####

ggpubr::ggarrange(plot1,plot2,widths = c(1,2),labels = c("A","B"),nrow = 2,heights =c(1,2))
ggsave("test.png",width = 10,height = 14,dpi = 1200,bg = "white")

ggpubr::ggarrange(plot1,plot2,widths = c(1,2),labels = c("A","B"),nrow = 2,heights =c(1,2))
ggpubr::ggarrange(plot7A,plot7B,plot7C,ncol = 1,labels = c("A","B","C"),heights = c(1,1,3))
ggsave("test.png",height = 17,width = 12)
