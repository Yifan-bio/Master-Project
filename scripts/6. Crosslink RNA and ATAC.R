
library(tidyverse)
library(clusterProfiler)
# Load informations
enhancer_ATAC = readRDS("../5.Regulator_accessibility/enhancer_ATAC.rds")
promoter_ATAC = readRDS("../5.Regulator_accessibility/promoter_ATAC.rds")
RNA = readRDS("../3.Differential_expression/resdata.rds")
RNA$gene_name[RNA$gene_name == "CCL3L1"] = "CCL3L3"
sig = RNA[abs(RNA$log2FoldChange) > 2 & RNA$padj < 0.05,]


# Combiniing the information of pormoter and enhnacer into one
Reg_ATAC = rbind(promoter_ATAC %>% select(annot.seqnames:annot.end,gene_name,annot.Fold) %>% mutate(name ="promoter"),enhancer_ATAC %>% select(annot.seqnames:annot.end,gene_name,annot.Fold) %>% mutate(name = "enhancer"))
Reg_ATAC = Reg_ATAC %>% distinct(annot.seqnames,gene_name,annot.Fold,name,.keep_all = T)

#### Isolate DEG results with regulators experiencing accessibility changes ####

Reg_RNA_ATAC = Reg_ATAC[Reg_ATAC$gene_name %in% sig$gene_name,]
Reg_ATAC_RNA = sig[sig$gene_name %in% Reg_ATAC$gene_name,]

#### Over-representation analysis of DEGs with RE ####

ORA_GO = function(sig_gene,term = 'ALL',background = NULL) {
  genelist = sig_gene$log2FoldChange
  names(genelist) = sig_gene$Row.names
  genelist = sort(genelist,decreasing = TRUE)
  
  gene = names(genelist)[abs(genelist) > 2]
  
  res <- enrichGO(gene = gene,
                  universe = background,
                  OrgDb = 'org.Hs.eg.db',
                  ont = term,
                  keyType = 'ENSEMBL',
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)
  return(res)
}

ORA_RNA_sigAcc = ORA_GO(sig_gene = Reg_ATAC_RNA,background = RNA$Row.names)


#### Isolate genes ####

extract_gene = function(GO_term,term) {
  temp = GO_term@result
  temp = temp[temp$Description == term,]
  temp = strsplit(as.character(temp$geneID), "/")[[1]]
  return(temp)
}

gene = extract_gene(GO_term = ORA_RNA_sigAcc,term = "inflammatory response")

annotate_region = function(region,event){
  if(class(region) == "CompressedGRangesList") {
    region = as.data.frame(region)
    region = GRanges(region)
  }
  if(class(event) == "data.frame") {
    row.names(event) = 1:nrow(event)
    event = GRanges(event)
  }
  result = annotatr::annotate_regions(regions = region,annotations = event)
  result = as.data.frame(result)
  return(result)
}

enhancer_peak = annotate_region(region = readRDS("../../../Paper/Paper 1/THP-1/R/annotation/THP-1_enhancer.rds"),event = rtracklayer::import.bed("../../../Paper/Paper 1/THP-1/R/HMMRATAC/merged.bed"))
enhancer_peak = enhancer_peak[enhancer_peak$gene_name %in% gene,]

promoter_peak = annotate_region(region = readRDS("../../../Paper/Paper 1/THP-1/R/annotation/promoter.rds"),event = rtracklayer::import.bed("../../../Paper/Paper 1/THP-1/R/HMMRATAC/merged.bed"))
promoter_peak = promoter_peak[promoter_peak$gene_name %in% gene,]

# ListOfGene = function(GO_term,ORA) {
#   temp = as.data.frame(ORA)
#   temp = temp[temp$Description %in% GO_term,]
#   temp = paste0(temp$geneID[1],"/",temp$geneID[2],"/",temp$geneID[3])
#   temp = unique(str_split(temp,"/"))
#   return(temp)
# }
# 
# GO_term = c("inflammatory response","positive regulation of cell migration",
#             "leukocyte cell-cell adhesion")
# 
# ListOfGene_ATAC = Reg_RNA_ATAC[Reg_RNA_ATAC$gene_name %in% ListOfGene(GO_term = GO_term,ORA = ORA_RNA_sigAcc)[[1]],]
# ListOfGene_RNA = Reg_ATAC_RNA[Reg_ATAC_RNA$gene_name %in% ListOfGene(GO_term = GO_term,ORA = ORA_RNA_sigAcc)[[1]],]
