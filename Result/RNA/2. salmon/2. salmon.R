# Script for the data visualisation for salmon results

# Inputing Datas
salmon_meta = read.csv("./Input/salmon.csv",sep = ";",row.names = 1)                      # Reading file and creating the dataframe
file_list = c('./Input/0hr_rep1/quant.sf',
              './Input/0hr_rep2/quant.sf',
              './Input/2hr_rep1/quant.sf',
              './Input/2hr_rep2/quant.sf',
              './Input/24hr_rep1/quant.sf',
              './Input/24hr_rep2/quant.sf')
txi <- tximport(file_list,type = "salmon", txOut = F,tx2gene = tx2gene,ignoreAfterBar = T)

# ------------------------- Plot 2A ----------------------------
# Plotting Mapping rate
ggplot(salmon_meta,aes(x = row.names(salmon_meta),y = Percentage.mapped)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample Name assigned",y = "Percentage of reads mapped from the library") +
  scale_y_continuous(limits = c(0,100),expand = c(0,0)) + 
  theme_bw() +
  labs(title = "% Reads Mapped")

# -------------------------------- Plot 2B -------------------------------------
# Plotting log info from salmon output
x = salmon_meta %>% select(Total:number.of.fragments.filtered.vm,Orphan)
# Calculating unmapped reads due to not recorded in file and Paired end mapped reads
x$Unrecorded = x$Total - x$Mapped - x$Decoy - x$Dovetail - x$number.of.fragments.filtered.vm
x$PairedMap = x$Mapped - x$Orphan
x = x %>% dplyr::select(PairedMap,Orphan,Unrecorded,Decoy,Dovetail,number.of.fragments.filtered.vm)
# Coverting into long dataframe
x$name = rownames(x)
x = x %>% pivot_longer(!name,names_to = "ReadType",values_to = "read")

ggplot(x,aes(y = read,x=name,fill = ReadType)) + 
  geom_bar(stat = 'identity',position = "fill") + # geom_bar(stat = 'identity) # This makes it a value not percentage
  coord_flip() +
  labs(x = "Sample",y="Number of reads within the replicate",fill = "Read type") +
  theme_bw() 

# -------------------------------- Plot 2C -------------------------------------
# Plotting raw count and compare it to transcript length (if effective length then can be extracted from tximport)

Plot2B = function(quant_file,Length_cut = FALSE,TPM_cut = FALSE,name_level = FALSE) {
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
# Plot2B("./Input/0hr_rep1/quant.sf",Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3)

plot_list = list(Plot2B(file_list[1],Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3),
                 Plot2B(file_list[2],Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3),
                 Plot2B(file_list[3],Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3),
                 Plot2B(file_list[4],Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3),
                 Plot2B(file_list[5],Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3),
                 Plot2B(file_list[6],Length_cut = 2.5e4,TPM_cut = 1e4,name_level = 3))

ggarrange(plotlist = plot_list,nrow = 3,ncol = 2)


# -------------------------------- Plot 2D -------------------------------------
# Plotting the TPM in replicates (density plot + box plot)

x = txi[["abundance"]]
x = data.frame(x)
colnames(x) = c("Untreat_1","Untreat_2","Early_1","Early_2","Later_1","Later_2")

# Converting to long dataframe and remove NULL values
x$gene = row.names(x)
x = x %>% pivot_longer(!gene,names_to = "replicate",values_to = "abundance")
x = x %>% filter(abundance != 0)

# Prepring for facet
x$condition = (do.call('rbind', strsplit(as.character(x$replicate),'_',fixed=TRUE)))[,1]
x$rep = (do.call('rbind', strsplit(as.character(x$replicate),'_',fixed=TRUE)))[,2]

ggplot(x,aes(x = log10(abundance),y = -0.02)) +
  # Horizontal boxplot
  ggstance::geom_boxploth(aes(fill = condition),width = 0.03) +
  # density plot
  geom_density(aes(x = log10(abundance)),inherit.aes = F) +
  facet_grid(rows = vars(rep),cols = vars(condition)) +
  scale_fill_discrete() +
  geom_hline(yintercept = 0)+
  labs(y = "Proportion of genes",x = "log10(TPM)")