# This file is for running differential methylation on DiffBind peaks
# Last modified: 17 March 2023

library(methylKit)

#### Intaking peaks regions ####

ATAC = readRDS("../5.Differential_accessibility/DiffBind/DiffBind_reaks.rds")
ATAC = GenomicRanges::GRanges(ATAC)

#### Running differential methylation in peak regions ####

list=list("../3.WGBS_preprocessing/4. MethylDackel/0hr_dedup_CpG.methylKit",
          "../3.WGBS_preprocessing/4. MethylDackel/24hr_dedup_CpG.methylKit")

# read the files to a methylRawList object: myobj (kept the file as it is)
myobj = methRead(list,
                 sample.id=list("ctrl1","treat1"),
                 assembly="hg38",
                 pipeline = list(fraction=FALSE,chr.col=2,start.col=3,end.col=3,coverage.col=5,strand.col=4,freqC.col=6),
                 treatment=c(0,1),
                 context="CpG",
                 mincov = 10
)

# Putting together region informations
myobj <- methylKit::regionCounts(myobj,ATAC)

#Merging samples
meth <- methylKit::unite(myobj)                   # destrand merge reads from both strand allow better coverage 

#Identify DMR/DMC
myDiff=calculateDiffMeth(meth,test = "fast.fisher",adjust = "fdr")

#### Combining the accessibility changes ####

library(tidyverse)

ATAC = readRDS("../5.Differential_accessibility/DiffBind/DiffBind_reaks.rds")
myDiff = data.frame(myDiff)
myDiff = merge(x = myDiff,
               y = ATAC,
               by.x = c("chr","start"),
               by.y = c("Chr","Start"))

myDiff = myDiff %>% dplyr::select(chr:meth.diff,Fold:FDR)
colnames(myDiff) = c("chr","start","end","strand","meth.pvalue","meth.qvalue","meth.diff","Acc.Fold","Acc.p-value","Acc.FDR")

#### Verify if there is a DEGs regulated regulator ####

ATAC_reg = readRDS("../7.RE_accessibility/DiffBindPeak_with_DEG.rds")
myDiff$con = "B"
myDiff$con[myDiff$chr %in% ATAC_reg$annot.seqnames & myDiff$start %in% ATAC_reg$annot.start] = "A"
