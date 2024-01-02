# Script for all plotting involved in section 1 RNA pre-processing
# Last modified: 19 March 2023

library(tidyverse)
library(scales)
#### Import information from salmon output ####

# The info file is still manual, need to create a script for it!!!
salmon_meta = read.csv("./1. salmon/salmon.csv",sep = ";",row.names = 1) 

file_list = c('./1. salmon/0hr_rep1/quant.sf',
              './1. salmon/0hr_rep2/quant.sf',
              './1. salmon/24hr_rep1/quant.sf',
              './1. salmon/24hr_rep2/quant.sf')

#### Plot1A Mapping rate and read component ####

# Plotting log info from salmon output
x = salmon_meta %>% select(Total:number.of.fragments.filtered.vm,Orphan)
# Calculating unmapped reads due to not recorded in file and Paired end mapped reads
x$Unrecorded = x$Total - x$Mapped - x$Decoy - x$Dovetail - x$number.of.fragments.filtered.vm
x$PairedMap = x$Mapped - x$Orphan
x = x %>% dplyr::select(PairedMap,Orphan,Unrecorded,Decoy,Dovetail,number.of.fragments.filtered.vm)

colnames(x) = c("Paired Mapped","Orphan Recovered","No Match","Decoy","Dovetail","Quality cutoff")
rownames(x) = c("Untreat Replicate 1","Untreat Replicate 2","Treat Replicate 1","Treat Replicate 2")
# Coverting into long dataframe
x$name = rownames(x)
x = x %>% pivot_longer(!name,names_to = "ReadType",values_to = "read")

Plot1A = ggplot(x,aes(y = read,x=name,fill= factor(ReadType,levels = c("No Match","Decoy","Dovetail","Quality cutoff","Orphan Recovered","Paired Mapped")))) + 
  geom_bar(stat = 'identity',position = "stack") +
  scale_y_continuous(labels = unit_format(unit = "",scale = 1e-6),expand = c(0,0),limits = c(0,1e8)) +
  coord_flip() +
  labs(x = NULL,y="Number of read pairs (Million)",fill = "Read type") +
  theme_bw() +
  scale_fill_manual(values=c("#F4606C","#ECAD9E","#E6CEAC","#D1BA74","#A0EEE1","#19CAAD"))

ggsave(plot = Plot1A,filename = "Plot1A.png",width = 9,height = 4)

#### Plot1B TPM versus length of transcript ####

Plot1B = function(quant_file,Length_cut = FALSE,TPM_cut = FALSE,name_level = FALSE) {
  x = readr::read_delim(quant_file)
  x = x %>% filter(TPM != 0)
  # To rename the Gencode name
  x$Name = (do.call('rbind', strsplit(as.character(x$Name),'|',fixed=TRUE)))[,2]
  x$Name = gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1",x$Name)
  
  # Adding limit to each axis
  if (TPM_cut != FALSE){
    x$TPM[x$TPM > TPM_cut] = TPM_cut
  }
  if (Length_cut != FALSE) {
    x$Length[x$Length > Length_cut] = Length_cut
  }
  
  plot = ggplot(x,aes(x = Length,y=TPM)) + geom_point() + scale_alpha_identity() +
    theme_bw() + labs(x = "Length of Transcript",y = "TPM of each transcript")  
  
  # Adding name to each plot, the name level means the level after split by "/":
  # So for "./name" then it level=2,"././name" then level=3
  if (name_level != FALSE){
    name = strsplit(as.character(quant_file),'/',fixed=TRUE)[[1]][name_level]
    plot = plot + labs(title = name)
  }
  return(plot)
}

plot_list = list(Plot1B(file_list[1],Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3),
                 Plot1B(file_list[2],Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3),
                 Plot1B(file_list[3],Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3),
                 Plot1B(file_list[4],Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3))

Plot1B = ggpubr::ggarrange(plotlist = plot_list,nrow = 2,ncol = 2)
ggsave(plot = Plot1B,filename = "Plot1B.png",width = 9,height = 6)

#### Plot1C ####



