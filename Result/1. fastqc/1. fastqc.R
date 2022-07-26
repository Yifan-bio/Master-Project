# All of my script for fastqc is based on the data that is downloadable through multiQC
# So all you need to prep before hand is a meta data file

library(tidyverse)
library(scales)

info = read.csv("./meta.csv",sep = ";")

# You should have a file with the accession number and condition

# ------------------------- Plot 1A --------------------------------------------
# Plotting Duplication, Read Length, Seq number
# This part ypu first configure the first table on multiQC to select  Duplication, Read Length, Seq number
# Then just copy table and paste it

# Data preperation
x = read_tsv("./Input/firstplot.txt")

plot1A = function(input,meta) {
  input$`% Dups` = as.numeric(gsub("%$","",input$`% Dups`))                               # Remove the percentage sign and convert to numeric
  input$`Read Length` = as.numeric(gsub(" bp","",input$`Read Length`))                    # Remove ` bp` and convert to numeric 
  input$`Sample Name` = (do.call('rbind', strsplit(as.character(input$`Sample Name`),'_',fixed=TRUE)))[,1] # Isolate only the accession number from files names in multiqc
  input = input %>% group_by(`Sample Name`) %>% mutate_each(funs(mean),`% Dups`) %>% distinct(`% Dups`,.keep_all = T) # Getting the average duplication between two reads
  input = merge(input,info,by.x = "Sample Name",by.y = 'Accession')             # binding the dataframe with the info we need
  
  # Plotting individually
  plot1 = ggplot(input,aes(y = `% Dups`,x=Rep,fill = Rep)) +                          # ggplot Duplication
    geom_bar(stat = "identity") +
    facet_wrap(~ Time,ncol = 1) +
    ylim(0,100) +
    theme_bw() +
    labs(x = NULL,y = "Percentage of duplication (%)")

  plot2 = ggplot(input,aes(y = `Read Length`,x=Rep,fill = Rep)) +                     # ggplot Read Length
    geom_bar(stat = "identity") +
    facet_wrap(~ Time,ncol = 1) +
    ylim(0,100) +
    theme_bw() +
    labs(x = NULL,y = "Read Length (bp)")

  plot3 = ggplot(input,aes(y = `M Seqs`,x=Rep,fill = Rep)) +                          # ggplot sequence number in millions
    geom_bar(stat = "identity") +
    facet_wrap(~ Time,ncol = 1) +
    ylim(0,100) +
    theme_bw() +
    labs(x = NULL,y = "Total number of sequencing reads (million)")

  return(plot1)
  output = ggpubr::ggarrange(plot1,plot2,plot3,common.legend = TRUE,legend = "right",nrow = 1)                            # Combining the plots
  return(output)
}

plot1A(input = x,meta = info)

# ------------------------------Plot 1B-----------------------------------------
# Plotting fastqc duplication (number of reads + % duplicated)
# This file is downloaded at <Sequence Counts> section

x = read_tsv("./Input/fastqc_sequence_counts_plot.tsv")
# Renaming using the accession number
x$Category <- (do.call('rbind', strsplit(as.character(x$Category),'_',fixed=TRUE)))[,1]
# Convert to long dataframe
x = x %>% gather("type","number",-Category)
# Using the meta data to allow the grid to be produced to separate time
x = merge(x,info,by.x = "Category",by.y = "Accession")
ggplot(x,aes(y=number,x=Category,fill = forcats::fct_rev(type))) + # forcats reorders the plot using ampped as main
  geom_bar(stat="identity") +
  facet_grid(~ Time,scales = "free_x",space = "free_x",switch = "x") + # this code split each group in to each column (I have not idea for coor_flip)
  labs(x = NULL,y = 'Number b of Fragments (million)',fill = 'Fragment Status') +
  scale_y_continuous(labels = unit_format(unit = "",scale = 1e-6),expand =expansion(mult = c(0.01,0.01))) +
  scale_fill_manual(values = c(colors()[619],colors()[618])) + 
  theme_bw()

# --------------------------------- Plot 1C ------------------------------------
# Plotting fastqc per base quality

x =  read_tsv("./Input/fastqc_per_base_sequence_quality_plot.tsv")
# Following steps rename all column using accession number which is the first phrase separate by "_"
z = colnames(x)
z = as.data.frame(do.call('rbind', strsplit(as.character(z),'_',fixed=TRUE)))   # Isolating the phrases
colnames(x) = paste0(z$V1)                                                      # Changing the column names to accession number
# This code combines the column with same accession number (R1 and R2) and value is the mean
x = as.data.frame(lapply(split(seq_len(ncol(x)),colnames(x)),function(cis) rowMeans(x[cis])))
# Convert to long dataframe
z = x %>% gather(key = "variable",value = "value",-Position..bp.)
# Plot
ggplot(z,aes(x = Position..bp.,y=value)) +
  geom_line(aes(color = variable)) +
  ylim(0,40)+
  xlab("Position in read(bp)") +
  ylab("Phred Score") +
  theme_bw()

# --------------------------------- Plot 1D ------------------------------------
# Plotting fastqc N content

x = read_tsv("./Input/fastqc_per_base_n_content_plot.tsv")

z = colnames(x)
z <- data.frame(do.call('rbind', strsplit(as.character(z),'_',fixed=TRUE)))
colnames(x) = paste0(z$X1,"(",z$X6,")")

x = x %>% gather(key = "variable",value = "value",-`Position in Read (bp)(Position in Read (bp))`)
ggplot(x,aes(x = `Position in Read (bp)(Position in Read (bp))`,y=value)) +
  geom_line(aes(color = variable)) +
  ylim(0,100)+
  xlab("Position in read(bp)") +
  ylab("Percentage of N content") +
  theme_bw() +
  theme(legend.position = "none")

# Remove all datas
remove(x,z,plot1A,info)
