# This script is to identify regulators for DEGs
library(tidyverse)
library(rtracklayer)

#### Importing DEGs from section 4 ####

resdata = readRDS("../4.Differential_expression/resdata.rds")
sig_rna = resdata[resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 2,]


#### Identify promoters using gtf file ####

promoter_gtf = function(gtf_file) {
  
  # Importing the gtf file
  temp <- rtracklayer::import.gff(gtf_file)
  temp <- as.data.frame(temp)

  # To only include the transcripts from the file and the stable id
  temp <- temp %>% filter(type == "transcript")
  temp$gene_id <- gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1",temp$gene_id)
  
  # Define promoter regions for positive strand
  temp1 <- temp[temp$strand == "+",]
  temp1$promoter_start <- temp1$start - 2000
  temp1$promoter_end <- temp1$start + 1000
  temp1$width = "3000"
  temp1 = temp1 %>% dplyr::select(seqnames,promoter_start,promoter_end,width,strand,gene_id,gene_name,transcript_name)

  # Define promoter regions for negative strand
  temp2 <- temp[temp$strand == "-",]
  temp2$promoter_start <- temp2$end - 1000
  temp2$promoter_end <- temp2$end + 2000
  temp2$width <- "3000"
  temp2 <- temp2 %>% dplyr::select(seqnames,promoter_start,promoter_end,width,strand,gene_id,gene_name,transcript_name)

  # Combine the files into the full promoter
  temp <- rbind(temp1,temp2)
  return(temp)
}

promoter = promoter_gtf("../support_doc/Gencode/gencode.V40.annotation.gtf")


#### Promoters for DEGs ####

# Select transcript records for the differentially expressed genes from section 4
promoter_DEGs <- promoter[promoter$gene_id %in% sig_rna$Row.names,]

saveRDS(promoter_DEGs,file = "promoter_DEGs.rds")
write_delim(promoter_DEGs,"promoter_DEGs.txt")


#### Identification of active enhancers for DEGs ####
# Some genes share multiple gene name so we will need to also use ensembl id as a confirmation

library("AnnotationDbi")
library("org.Hs.eg.db")
library("tidyverse")

enhancer_DEGs = function(active_enhancer_bed) {
  # Import the files
  enhancer = read.delim(active_enhancer_bed,header = FALSE,col.names = c("chr","start","end","Genhancer","gene_name","tech"))
  # Getting the ensembl gene number
  enhancer$ensid = mapIds(x = org.Hs.eg.db,
                          keys = enhancer$gene_name, 
                          column = "ENSEMBL",
                          keytype = "SYMBOL",
                          multiVals = "first")
  # Removing duplicates from bedtools merge
  enhancer = enhancer %>% distinct(chr,start,end,Genhancer,gene_name,tech,ensid)
  return(enhancer)
}
enhancer = enhancer_DEGs(active_enhancer_bed = "./active_enhancer.bed")

#### Locating active enhancer for DEGs ####

# Getting the enhancer for DEGs
active_enhancer_DEGs = enhancer[enhancer$gene_name %in% sig_rna$gene_name | enhancer$ensid %in% sig_rna$Row.names,]

saveRDS(active_enhancer_DEGs,file = "active_enhancer_DEGs.rds")
