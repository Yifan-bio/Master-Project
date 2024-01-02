# This script is to identify regulators for DEGs
library(tidyverse)
library(rtracklayer)
library(plyranges)
library(GenomicFeatures)
library(GenomicRanges)

#### Identify promoters ####

promoter_gtf = function(gtf_file) {
  
  # Importing the gtf file
  temp <- rtracklayer::import.gff(gtf_file)
  
  temp = temp[(elementMetadata(temp)[,"type"] == "transcript")]
  
  # Isolate TSS from the transcript
  start(temp) <- ifelse(strand(temp) == "+",
                        yes = start(temp),
                        no = end(temp))
  end(temp) <- ifelse(strand(temp) == "+",
                      yes = start(temp),
                      no = end(temp))

  # Getting the promoter regions
  temp = plyranges::flank_upstream(x = temp,width = 3000)
  temp = plyranges::shift_downstream(x = temp,shift = 1000)

  return(temp)
}

promoter = promoter_gtf("../support_doc/Gencode/gencode.v40.annotation.gtf")

#### Identify exon and intron ####

gene_gtf = GenomicFeatures::makeTxDbFromGFF("../support_doc/Gencode/gencode.v40.annotation.gtf")

exon <- exonsBy(gene_gtf,by = "tx")
intron <- GenomicFeatures::intronsByTranscript(gene_gtf)

#### Identify enhancers ####
# Some genes share multiple gene name so we will need to also use ensembl id as a confirmation

active_enhancer = function(geneHancer,H3K27ac) {
  
  # Importing Enhancer regions
  enhancer = read.delim(geneHancer,header = FALSE,col.names = c("chr","start","end","Genhancer","gene_name","tech"))
  enhancer = GRanges(enhancer)
  
  # Importing H3K27ac markers
  marker = rtracklayer::import.bed(H3K27ac)
  #marker = read.delim(marker,header = FALSE,col.names = c("chr","start","end"))
  
  # Identify active enhancer
  active_enhancer = IRanges::subsetByOverlaps(x = enhancer,ranges = marker)
  
  return(active_enhancer)
}

active_enhancer = active_enhancer(geneHancer = "./Genhancer/geneHancer_format fix.bed",H3K27ac = "./H3K27ac/H3K27ac.bed")

#### Importing chroamtin accessibilty changes ####

ATAC = readRDS("../5.Differential_accessibility/DiffBind.rds")
rownames(ATAC) = 1:nrow(ATAC)
ATAC$name = paste0(ATAC$Chr,".",ATAC$Start)
ATAC = GRanges(ATAC)

#### Accessibility of regions ####

annotate_region = function(region,event){
  if(class(region) == "CompressedGRangesList") {
    region = as.data.frame(region)
    region = GRanges(region)
  }
  result = annotatr::annotate_regions(regions = region,annotations = ATAC)
  result = as.data.frame(result)
  return(result)
}

enhancer_ATAC = annotate_region(active_enhancer,ATAC)
promoter_ATAC = annotate_region(promoter,ATAC)
intron_ATAC = annotate_region(intron,ATAC)
exon_ATAC = annotate_region(exon,ATAC)

