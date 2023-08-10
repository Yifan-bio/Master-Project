# This script is put together to allow the running of differential expression analysis using RNA-seq
# Last modification: 14 March 2023

library(tidyverse)

#### Importing information from gtf file ####

gtf_file = ""

library(rtracklayer)
library(plyranges)
library(GenomicFeatures)
library(GenomicRanges)

# Creating the tx2gene variable for RNA-seq
# tximport uses the tx2gene object to perform trasncript to gene conversion

tx2gene_creation = function(gtf_file) {
  # Intaking the gtf file as binary
  temp1 = GenomicFeatures::makeTxDbFromGFF(file = gtf_file)
  # Extracting the required information
  temp2 = AnnotationDbi::keys(x = temp1,
                              keytype = 'TXNAME')
  # Selecting the required information
  temp = AnnotationDbi::select(temp1,
                               keys = temp2,
                               columns = 'GENEID',
                               keytype = 'TXNAME')
  return(temp)
}

tx2gene = tx2gene_creation(gtf_file = gtf_file)

# Creating regions with promoter

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

promoter = promoter_gtf(gtf_file = gtf_file)

# Identify exon and intron

gene_gtf = GenomicFeatures::makeTxDbFromGFF(gtf_file)

exon <- exonsBy(gene_gtf,by = "tx")
intron <- GenomicFeatures::intronsByTranscript(gene_gtf)

remove(gene_gtf)

#### Importing enhancer regions ####
# Some genes share multiple gene name so we will need to also use ensembl id as a confirmation

enhancer = ""
marker = ""

active_enhancer = function(geneHancer,histone = NULL) {
  
  # Importing Enhancer regions
  enhancer = read.delim(geneHancer,header = FALSE,col.names = c("chr","start","end","Genhancer","gene_name","tech"))
  enhancer = GRanges(enhancer)
  
  # Get active enhancer using histone marker
  if (is.null(histone) = FALSE) {
    # Importing histone marker locations that want to look at
    marker = rtracklayer::import.bed(histone)
    
    # Identify active enhancer
    active_enhancer = IRanges::subsetByOverlaps(x = enhancer,ranges = marker)
  
    return(active_enhancer)
  }

  return(enhancer)

}

active_enhancer = active_enhancer(geneHancer = enhancer,histone = marker)


#### Function ####

# Annotating genes based on gtf files
gene_annot <- function(result = result,dds = dds,gtf_file) {
  
  # Creating a ID to symbol conversion table
  gtf_file <- rtracklayer::import(gtf_file)
  gtf_file <- as.data.frame(gtf_file)
  gtf_file <- unique(gtf_file[,c("gene_id","gene_name")])
  
  # Putting all important information into a single object
  res <- merge(as.data.frame(result),
               as.data.frame(counts(dds,normalized=TRUE)),
               by = "row.names",
               sort = FALSE)
  # Putting the gene symbol with the gene id
  res <- merge(x = res,
               y = gtf_file,
               by.x = "Row.names",
               by.y = "gene_id")
  
  # Getting gene stable version number
  res$Row.names <- gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1",res$Row.names) 
  
  return(res)
}

#### RNA-seq analysis ####

library("GenomicFeatures")
library("tximport")
library("DESeq2")

# Listing the location of each salmon output count file (quant.sf files)
files = c('./Input/0hr_rep1/quant.sf',
          './Input/0hr_rep2/quant.sf',
          './Input/24hr_rep1/quant.sf',
          './Input/24hr_rep2/quant.sf')

# Importing salmon files into R and combining the transcript counts into gene counts
txi.salmon <- tximport(files = files,
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

# Creating a table to specific the condition of each file.
sampleTable <- data.frame(condition = factor(rep(c("control", "treated"),each = 2)))

# Connecting the conditions with the file imported using tximport.
rownames(sampleTable) <- colnames(txi.salmon$counts) 

# Importing to DESeq2 format
dds <- DESeqDataSetFromTximport(txi = txi.salmon, 
                                colData = sampleTable,
                                design = ~ condition)

# Running the differential expression analysis 
dds <- DESeq(object = dds,
             test = 'Wald')

# Getting the results out of DEseq2
result = lfcShrink(dds = dds,
                   coef = "condition_treated_vs_control",
                   type = "apeglm"
)

resdata = gene_annot(result = result,dds = dds,gtf_file = gtf_file)

resdata = resdata[resdata$baseMean > 20,]

sig = resdata[resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 2,]


#### ATAC-seq ####

library(DiffBind)
library(tidyverse)
library(annotatr)

# Putting in the metadata
sampleTbl = data.frame(sampleID = c('control1','control2','treated1','treated2'),
                       Cell = 'THP1',
                       Factor = 'PMA',
                       Condition= c('Control','Control','Treated','Treated'),
                       Replicate = c("1","2","1","2"),
                       bamReads = c('./Input/bam/untreat_rep1.rmChrM.dedup.filter.bam',
                                    './Input/bam/untreat_rep2.rmChrM.dedup.filter.bam',
                                    './Input/bam/treat_rep1.rmChrM.dedup.filter.bam',
                                    './Input/bam/treat_rep2.rmChrM.dedup.filter.bam'),
                       Peaks = c('./Input/HMM_removeHC/untreat_rep1.gappedPeak',
                                 './Input/HMM_removeHC/untreat_rep2.gappedPeak',
                                 './Input/HMM_removeHC/treat_rep1.gappedPeak',
                                 './Input/HMM_removeHC/treat_rep2.gappedPeak'),
                       ScoreCol = 13,
                       LowerBetter = FALSE
)

# Importing the data
ATAC <- dba(sampleSheet=sampleTbl)

# Perform count on each peak using the BAM file
ATAC_count <- dba.count(ATAC,minOverlap = 2,
                        score = DBA_SCORE_NORMALIZED,
                        bUseSummarizeOverlaps = TRUE
)

# Differential accessibility analysis
ATAC_contrast <- dba.contrast(ATAC_count, categories=DBA_CONDITION,minMembers = 2)

# Perform differential affinitiy analysis
ATAC_analyze <- dba.analyze(ATAC_contrast,method = DBA_ALL_METHODS)

# Extracting the results
ATAC_report <- dba.report(DBA = ATAC_analyze,
                          DataType = DBA_DATA_FRAME,
                          method = DBA_DESEQ2,
                          contrast = 1,
                          th = 1)

#Extracting usable regions
Diff.bind.sig = ATAC_report %>% dplyr::filter(abs(Fold) > 2 & FDR < 0.01)
rownames(Diff.bind.sig) = 1:nrow(Diff.bind.sig)

#### Genomic annotation creation ####

library("rtracklayer")

gtf_file = "../../../"
