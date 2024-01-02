# Script for all plotting involved in section 2 ATAC pre-processing
# Last modified: 19 March 2023

library(tidyverse)
library(scales)
library(ggExtra)
library(annotatr)
library(GenomicRanges)

#### Importing inforamtion ####

# This file was create manually, make a command
Filter = read_tsv("./3. filter/Filter.txt")

#### Plot2A Filtered reads ####

# Change the table into a long dataframe
temp = Filter %>% gather(key = "variable",value = "value",-Name)

Plot2A = ggplot(temp,
                aes(x=value,y=variable,fill= factor(Name,levels = c("Trimmed","Unmapped","MT read","Duplicate read","Filter","Final usable")))) + 
  geom_bar(stat = "identity",position="stack") +
  scale_x_continuous(labels = unit_format(unit = "",scale = 1e-6),expand = c(0,0),limits = c(0,1.5e8)) +
  scale_y_discrete(labels = c("Treat Replicate 1","Treat Replicate 2","Untreat Replicate 1","Untreat Replicate 2")) +
  labs(x = "Read Pairs (Million)",y="",fill = "") +
  theme_bw() +
  scale_fill_manual(values=c("#3a3a3a","#aeaeae","#3066be","#3384d5","#b4c5e4","#19CAAD"))

ggsave(plot = Plot2A,filename = "Plot2A.png",width = 9,height = 4.9,dpi = 1200)

#### Plot2B Peak score and size distribution ####

# Creating plot for the score and size distribution. 
Create_plot1A=function(file_dir){
  peak=Peak_input(file_dir)
  p = ggplot(peak,aes(x=score,y=log10(width))) + geom_point() + geom_smooth(method = lm,se=F) +
    scale_x_continuous(limits = c(0,700),expand = c(0,0)) + scale_y_continuous(limits = c(2,5),expand = c(0,0))
  p = ggExtra::ggMarginal(p,type = "boxplot")
  return(p)
}

T2 = Create_plot1A("./4. HMMRATAC/treat_rep2_peaks_peaks.gappedPeak")
T1 = Create_plot1A("./4. HMMRATAC/treat_rep1_peaks_peaks.gappedPeak")
C1 = Create_plot1A("./4. HMMRATAC/untreat_rep1_peak_peaks.gappedPeak")
C2 = Create_plot1A("./4. HMMRATAC/untreat_rep2_peak_peaks.gappedPeak")

Plot2B = ggpubr::ggarrange(C1,C2,T1,T2,labels = c("Untreat 1","Untreat 2","Treat 1","Treat 2"),label.x = 0.5)

ggsave(plot = Plot2B,filename = "Plot2B.png",bg = "white")

#### Plot2C Number of peaks ####

# Showing number of peaks detected for each sample
URep1 = rtracklayer::import.bed(con = "./4. HMMRATAC/untreat_rep1_peak_peaks.gappedPeak")
URep2 = rtracklayer::import.bed(con = "./4. HMMRATAC/untreat_rep2_peak_peaks.gappedPeak")
TRep1 = rtracklayer::import.bed(con = "./4. HMMRATAC/treat_rep1_peaks_peaks.gappedPeak")
TRep2 = rtracklayer::import.bed(con = "./4. HMMRATAC/treat_rep2_peaks_peaks.gappedPeak")

# Common across replicates in same conditions
UCommon = as.data.frame(annotatr::annotate_regions(URep1,URep2))
TCommon = as.data.frame(annotatr::annotate_regions(TRep1,TRep2))

# Common across all replicates
find_common_all = function(file1,file2,file3,file4) {
  file = annotatr::annotate_regions(file1,file2)
  file = annotatr::annotate_regions(file,file3)
  file = annotatr::annotate_regions(file,file4)
  file = as.data.frame(file)
  return(file)
}

UCommon_all_rep1 = find_common_all(URep1,URep2,TRep1,TRep2)
UCommon_all_rep2 = find_common_all(URep2,URep1,TRep1,TRep2)
TCommon_all_rep1 = find_common_all(TRep1,TRep2,URep1,URep2)
TCommon_all_rep2 = find_common_all(TRep2,TRep1,URep1,URep2)

df = data.frame(common_all = c(nrow(UCommon_all_rep1 %>% count(seqnames,start)),
                               nrow(UCommon_all_rep2 %>% count(seqnames,start)),
                               nrow(TCommon_all_rep1 %>% count(seqnames,start)),
                               nrow(TCommon_all_rep2 %>% count(seqnames,start))
                               ),
                untreat_common = c((nrow(UCommon %>% count(seqnames,start)) - nrow(UCommon_all_rep1 %>% count(seqnames,start))),
                                   (nrow(UCommon %>% count(annot.seqnames,annot.start)) - nrow(UCommon_all_rep2 %>% count(seqnames,start))),
                                   0,
                                   0
                                   ),
                treat_common = c(0,
                                 0,
                                 (nrow(TCommon %>% count(seqnames,start)) - nrow(TCommon_all_rep1 %>% count(seqnames,start))),
                                 (nrow(TCommon %>% count(annot.seqnames,annot.start)) - nrow(TCommon_all_rep2 %>% count(seqnames,start)))
                                 ),
                condition_unique = c(length(URep1) - nrow(UCommon %>% count(seqnames,start)),
                                     length(URep2) - nrow(UCommon %>% count(annot.seqnames,annot.start)),
                                     length(TRep1) - nrow(TCommon %>% count(seqnames,start)),
                                     length(TRep2) - nrow(TCommon %>% count(annot.seqnames,annot.start))
                                     ),
                condition = c("Untreat Replicate 1",
                              "Untreat Replicate 2",
                              "Treat Replicate 1",
                              "Treat Replicate 2")
                )
df = readRDS("Plot3C_dataframe.rds")

df = pivot_longer(df,!condition,names_to = "name",values_to = "value")

plot2C = ggplot(df,aes(y=condition,x=value,fill= factor(name,levels = c("condition_unique","untreat_common","treat_common","common_all")))) + 
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Number of peaks",y="",fill = "") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0),limits = c(0,65000)) +
  scale_fill_manual(values=c("#cddeaa","#68b88e","#41b349","#1a6840"),labels = c("Unique to Replicates","Common in Untreated","Common in Treated","Common across all"))
  
ggsave("Plot2C.png",width = 10,height = 4)

#### Plot2D Venn diagram ####

library(ggvenn)
df = list(Rep1 = 1:41556,Rep2 = 13608:42586)
ggvenn(df,fill_color = c("grey","white"))

#### ####
x = rbind(TCommon_all_rep1,TCommon_all_rep2,UCommon_all_rep1,UCommon_all_rep2)
x = x %>% select(seqnames:end)
rtracklayer::export.bed(x,"Common_peaks.bed")
#THEN bedtools merge


ggpubr::ggarrange(Plot2A,plot2C,nrow = 2,labels = c("A","B"))
ggsave("Plot.png",dpi = 1200,height = 5.2,width = 6.6)
