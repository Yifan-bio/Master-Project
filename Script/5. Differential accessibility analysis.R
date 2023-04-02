# This script is for the differential accessibily analysis using ATAC-seq
# Last modified: 14 March 2023

library(DiffBind)

#### Creating a dataframe for the input to DESeq2 algorithm ####


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

#### Differential Accessibility Analysis using DiffBind wrapped around DESeq2 ####

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

#### Extracting usable regions ####

library(tidyverse)

Diff.bind.sig = ATAC_report %>% dplyr::filter(abs(Fold) > 2 & FDR < 0.01)
rownames(Diff.bind.sig) = 1:nrow(Diff.bind.sig)

#### Annotating the regions ####

library(annotatr)
# Preping annotations
read_annotations(con="../6.Regulators/active_enhancer.bed",genome = "hg38",name="Genhancer",format = "bed")
annotatr_annotations=build_annotations(genome = 'hg38',annotations = c('hg38_custom_Genhancer','hg38_basicgenes'))

# Annotating the results
annotated = annotate_regions(regions = GRanges(Diff.bind.sig),
                             annotations = annotatr_annotations, 
                             ignore.strand = TRUE,
                             quiet = FALSE) 
annotated = as.data.frame(annotated)
