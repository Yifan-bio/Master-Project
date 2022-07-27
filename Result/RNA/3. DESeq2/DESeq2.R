library(tidyverse)
library(scales) # needed for oob parameter
library(viridis)
# Preparing data
files = c('./Input/0hr_rep1/quant.sf',
          './Input/0hr_rep2/quant.sf',
          './Input/2hr_rep1/quant.sf',
          './Input/2hr_rep2/quant.sf',
          './Input/24hr_rep1/quant.sf',
          './Input/24hr_rep2/quant.sf')
txi <- tximport::tximport(files,type = "salmon", txOut = F,tx2gene = tx2gene,ignoreAfterBar = T) 
table <- data.frame(condition = factor(rep(c("untreated","2hr","24hr"),each = 2)))
rownames(table) <- colnames(txi$counts)
dds <- DESeq2::DESeqDataSetFromTximport(txi,table,~ condition)
dds <- dds[(rowSums(DESeq2::counts(dds)) > 10),]
dds <- DESeq2::DESeq(dds,test = 'Wald') 
remove(files,table)

# ------------------------------- Plot3A --------------------------------------
# Plotting the gene count between two condition and marking the uniquely expressed genes

plot3A = function(txi,sigGene,figure_name = NULL){
  x = data.frame(txi$counts)
  x = x[rowSums(x) > 0,] # Remove non-expressed genes
  # Isolating unique genes
  uni = x[rowSums(x > 0) == 1,]
  # Summarising the two replicates
  x$untreat = (x$X1 + x$X2)
  x$treat = (x$X3 + x$X4)
  x = x %>% select(untreat,treat)
  # label the uniquely expressed genes
  x$uni = row.names(x) %in% row.names(uni)
  # Rename the genes for resdata after DESeq2
  x$name = gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1",rownames(x))
  sig = sigGene %>% filter(abs(log2FoldChange) > 1)
  # Label based on uniquely expressed, significantly regulated and other
  x$con = x$name %in% sig$Row.names
  x$con[x$uni == "TRUE"] = "Unique"
  
  ggplot(x,aes(x = log10(untreat),y=log10(treat),color = con)) + geom_point() +
    theme_bw() + labs(color = "Sig Reg(FDR+p)",title = figure_name) + 
    scale_x_continuous(expand = c(.001,.001)) +
    scale_y_continuous(expand = c(.001,.001))
}

plot3A(txi = txi.salmon_24hr,sigGene = sig_24hr,figure_name = "24hr Comparing the gene count at two condition")

# ------------------------------- Plot3B --------------------------------------
# Volcano PLot
Plot3B = function(resdata){
  x = resdata %>% select(log2FoldChange,padj) 
  x$padj = ifelse(x$padj == 0,1e-300,x$padj)
  x$con = "non-significant"
  x$con[x$log2FoldChange > 1 & x$padj < 0.05] = "upregulated"
  x$con[x$log2FoldChange < -1 & x$padj < 0.05] = "downregulated"
  z = ggplot(x,aes(x=log2FoldChange,y=-log10(padj),color = con)) + geom_point() +
    scale_colour_manual(values = c("green","grey","red")) +                       # Change the colour of dots
    theme_bw() +                                                                  # make it bw rather then grey background
    labs(x = "Change in Gene Expression",y = "P-adjustment value (-log10(p-adj))",color = "Gene status") + #naming
    scale_y_continuous(expand = c(0,10)) +                                         # Remove the gap between the results and x margin
    geom_vline(xintercept = 1,linetype="dotted") + geom_hline(yintercept = -log10(0.05),linetype="dotted") + geom_vline(xintercept = -1,linetype="dotted")
  return(z)
}

Plot3B(resdata = resdata_24hr)

# ------------------------------- Plot3C --------------------------------------
# PCA plot (two options)
rld <- rlog(dds, blind=TRUE)
plotPCA(rld) + theme_bw()

#vsdata <- vst(dds, blind=FALSE)
#plotPCA(vsdata)

# ------------------------------- Plot3D --------------------------------------
# Plotting Gene counts and significant regulated genes in each replicate

Plot3D = function(txi,sig) {
  x = as.data.frame(txi$counts)
  colnames(x) = c("Untreat_1","Untreat_2","Treat_1","Treat2")
  
  x$name = gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1",rownames(x))
  # To long dataframe
  x = x %>% pivot_longer(!name,names_to = "con",values_to = "count")
  # Remove non-expressed genes
  x = x %>% filter(count != 0)
  
  # Identify significant genes
  x$sig = x$name %in% sig$Row.names
  
  x = ggplot(x,aes(x=log10(count),fill=sig)) +
        geom_histogram() + 
        facet_wrap(~con)
  return(x)
}
x1 = Plot3D(txi = txi.salmon_24hr,sig = sig_24hr)
x2 = Plot3D(txi = txi.salmon_2hr,sig = sig_2hr)

ggpubr::ggarrange(x1,x2,common.legend = TRUE,labels = c("24 Hour","2 Hour"))

# ----------------------------------- Plot 3E ----------------------------------
library("VennDiagram")
# Venn diagram for significantly regulated genes between the two conditions
x = list(Later = as.character(sig_24hr$Row.names),Early = as.character(sig_2hr$Row.names))
venn.diagram(x = x,filename = "Plot3E.png")

# ---------------------------------- Plot3F ------------------------------------
# MAplot

x = as.data.frame(resIHW_2hr)
x$significant <- ifelse(x$padj < .05, "Significant", NA)
ggplot(x, aes(x = baseMean,y = log2FoldChange, color = significant)) +
  geom_point() +
  scale_y_continuous(limits=c(-3, 3), oob=squish) +
  scale_x_log10() + 
  geom_hline(yintercept = 0, colour="tomato1", size=2) + 
  labs(x="mean of normalized counts", y="log fold change") + 
  scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + 
  theme_bw()
