# This script is for checking the accessibility of regulaotry elements of DEGs

# Loading all required datasets
enhancer = readRDS("../6.Regulators/active_enhancer_DEGs.rds")
promoter = readRDS("../6.Regulators/promoter_DEGs.rds")
ATAC = readRDS("../5.Differential_accessibility/DiffBind/DiffBind.rds")
RNA = readRDS("../4.Differential_expression/resdata.rds")

#### Preparation of datasets for annotation ####

# renamng for easier in later stage
colnames(promoter) = c("chr","start","end","width","strand","gene_id","gene_name","tx_name")
colnames(enhancer) = c("chr","start","end","Genhancer","gene_name","tech","gene_id")

# Rename the rows as annotatr wont go wrong
row.names(promoter) = 1:length(promoter$chr)
row.names(enhancer) = 1:length(enhancer$chr)
row.names(ATAC) = 1:length(ATAC$Chr)


#### Search for intersection ####

library(annotatr)
library(GenomicRanges)
library(tidyverse)

ATAC_enhancer = annotatr::annotate_regions(GRanges(enhancer),GRanges(ATAC))
ATAC_enhancer = as.data.frame(ATAC_enhancer)
ATAC_promoter = annotatr::annotate_regions(GRanges(promoter),GRanges(ATAC))
ATAC_promoter = as.data.frame(ATAC_promoter)
ATAC_RNA_RE = rbind(ATAC_promoter %>% select(seqnames:end,gene_name,gene_id),ATAC_enhancer %>% select(seqnames:end,gene_name,gene_id))

#### Isolate DEG results with regulators experiencing accessibility changes ####

RNA_sigAcc = RNA[RNA$Row.names %in% ATAC_RNA_RE$gene_id | RNA$gene_name %in% ATAC_RNA_RE$gene_name,]

#### Over-representation analysis of DEGs with RE ####

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

ORA_RNA_sigAcc = ORA_GO(RNA_sigAcc)

#### Exporting important results for next step ####

temp = rbind(ATAC_enhancer %>% select(annot.seqnames:annot.end),ATAC_promoter %>% select(annot.seqnames:annot.end))
temp = temp %>% distinct(annot.seqnames,annot.start,annot.end)
saveRDS(temp,"DiffBindPeak_with_DEG.rds")                                      