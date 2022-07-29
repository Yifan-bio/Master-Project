
# ------------------------------ Plot1A ---------------------------------------
# Plotting Duplication, Read Length, Seq number, GC content

x = read_tsv("./Input/raw/general_stats.txt")

Plot1A = function(x) {
  
  x$`% Dups` = as.numeric(gsub("%$","",x$`% Dups`))                               # Remove the percentage sign and convert to numeric
  x$`% GC` = as.numeric(gsub("%$","",x$`% GC`))                               # Remove the percentage sign and convert to numeric
  x$`Read Length` = as.numeric(gsub(" bp","",x$`Read Length`))                    # Remove ` bp` and convert to numeric 
  x$Time = (do.call('rbind', strsplit(as.character(x$`Sample Name`),'_',fixed=TRUE)))[,3] # Isolate only the accession number from files names in multiqc
  x$Rep = (do.call('rbind', strsplit(as.character(x$`Sample Name`),'_',fixed=TRUE)))[,4]
  x$`Sample Name` = paste0((do.call('rbind', strsplit(as.character(x$`Sample Name`),'_',fixed=TRUE)))[,3],
                           "_",(do.call('rbind', strsplit(as.character(x$`Sample Name`),'_',fixed=TRUE)))[,4])
  x = x %>% group_by(`Sample Name`) %>% mutate_each(funs(mean),`% Dups`) %>% distinct(`% Dups`,.keep_all = T) # Getting the average duplication between two reads

# Plotting individually
plot1 = ggplot(x,aes(y = `% Dups`,x=Rep,fill = Rep)) +                          # ggplot Duplication
  geom_bar(stat = "identity") +
  facet_wrap(~ Time,ncol = 1) +
  ylim(0,100) +
  theme_bw() +
  labs(x = NULL,y = "Percentage of duplication (%)")

plot2 = ggplot(x,aes(y = `Read Length`,x=Rep,fill = Rep)) +                     # ggplot Read Length
  geom_bar(stat = "identity") +
  facet_wrap(~ Time,ncol = 1) +
  ylim(0,100) +
  theme_bw() +
  labs(x = NULL,y = "Read Length (bp)")

plot3 = ggplot(x,aes(y = `% GC`,x=Rep,fill = Rep)) +                          # ggplot sequence number in millions
  geom_bar(stat = "identity") +
  facet_wrap(~ Time,ncol = 1) +
  ylim(0,100) +
  theme_bw() +
  labs(x = NULL,y = "Percentage of GC (%)")

plot4 = ggplot(x,aes(y = `M Seqs`,x=Rep,fill = Rep)) +                          # ggplot sequence number in millions
  geom_bar(stat = "identity") +
  facet_wrap(~ Time,ncol = 1) +
  ylim(0,150) +
  theme_bw() +
  labs(x = NULL,y = "Total number of sequencing reads (million)")

ggpubr::ggarrange(plot1,plot2,plot3,plot4,common.legend = TRUE,legend = "right",nrow = 1) 
}

Plot1A(x)
ggsave("Plot1A.png",dpi = 900)

# ------------------------------ Plot1B ---------------------------------------
# Plotting N content after + before trim

raw = read_tsv("./Input/raw/fastqc_per_base_n_content_plot.tsv")
raw = raw %>% pivot_longer(cols = !`Position in Read (bp)`,names_to = "sample",values_to = "percent",names_prefix = "ATAC_THP1_")
raw$con = "raw"

trim = read_tsv("./Input/trim/fastqc_per_base_n_content_plot.tsv")
trim = trim %>% pivot_longer(cols = !`Position in Read (bp)`,names_to = "sample",values_to = "percent",names_prefix = "ATAC_THP1_")
trim$con = 'trim'

x = rbind(raw,trim)
ggplot(x,aes(x = `Position in Read (bp)`,y=percent,group = sample)) +
  geom_line(aes(color = con)) +
  ylim(0,2)+
  xlab("Position in read(bp)") +
  ylab("Percentage of N content") +
  theme_bw() +
  theme(legend.position = "none")
ggsave("Plot1B.png",dpi = 909)

# ------------------------------ Plot1C ---------------------------------------
# Plotting Per base quality

raw = read_tsv("./Input/raw/fastqc_per_base_sequence_quality_plot.tsv")
raw = raw %>% pivot_longer(cols = !`Position (bp)`,names_to = "sample",values_to = "percent",names_prefix = "ATAC_THP1_")
raw$con = "raw"

trim = read_tsv("./Input/trim/fastqc_per_base_sequence_quality_plot.tsv")
trim = trim %>% pivot_longer(cols = !`Position (bp)`,names_to = "sample",values_to = "percent",names_prefix = "ATAC_THP1_")
trim$con = 'trim'

x = rbind(raw,trim)
ggplot(x,aes(x = `Position (bp)`,y=percent,group = sample)) +
  geom_line(aes(color = con)) +
  ylim(0,36)+
  xlab("Position in read(bp)") +
  ylab("Percentage of N content") +
  theme_bw() +
  theme(legend.position = "none")
ggsave("Plot1C.png",dpi = 909)

