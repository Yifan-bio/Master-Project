# Regulatory Element
> This section records the Regulatory element identified for each Differential expressed genes

* [Package version](#package-version)
* [Identification of active enhancers](#identification-of-active-enhancers)
    + [Input and Output](#input-and-output)
    + [Command](#command)
    + [Result](#result)
* [Identification of promoters](#identification-of-promoters)
    + [Input and Output](#input-and-output)
    + [Command](#command)
    + [Result](#result)
* [Reference](#reference)

## Package version
`tidyverse v1.3.1` `bedtools`


## Identification of active enhancers

We want to pull the list of regulatory element from the whole database to THP-1 cell specific. So we combine H3K27ac marker with enhancer to identify THP-1 active enhancers. However, there isn't any H3K27ac ChIP seq for 24 hour 100ng/ml PMA THP-1 cells so we only included 5 replicates of untreated THP-1 cell ChIP-seq results.

### Input and Output

* Input file
    - Genhancer interaction file
        - Collected the Double Elite file from UCSC browser to provide confident
    - H3K27ac marker from ChIP-atlas
        - GSM3514953; GSM3514952; GSM2544236; GSM2544237; GSM4083810
    - List of DEGs from section 4
* Output file
    - List of active enhancers for DEGs in THP-1 cells
    - List of active enhancers in THP-1 cells
    - H3K27ac marker in untreated THP-1 cells

### Command

```sh
# merging all bed files from ChIP-atlas into one file
cat $.bed > H3K27ac_temp1.bed
sort -k1,1 -k2,2n H3K27ac_temp1.bed > H3K27ac_temp2.bed
bedtools merge -i H3K27ac_temp2.bed > H3K27ac.bed
rm H3K27ac_temp1.bed H3K27ac_temp2.bed
```

```R
# Identification of enhancers for DEGs
library(tidyverse)

temp <- read.delim("/path/to/GeneInteraction_doubleElite.bed")
# Select enhancers for DEGs
temp <- temp[temp$geneName %in% sig$hgnc_symbol,]
# Remove metadata (Crucial)
temp <- dplyr::select(.data = temp,geneHancerChrom:geneHancerStrand,geneName,geneAssociationMethods)
# Export for bedtools
write_delim(x = temp,file = "DEG_enhancer.txt",delim = "\t",col_names = F)
```

```sh
# Combine the marker with enhancer
bedtools -a DEG_enhancer.txt -b H3K27ac.bed -wa > DEG_active_enhancer.bed
```

### Result



## Identification of promoters

Promoter are defined as the regions around TSS. As first exon shows a highly similar role as promoter, we will include both directions. The promoter regions will therefor be defined as +2kb -1kb from TSS.

### Input and Output

* Input file
    - List of DEGs from section 4
    - gtf file
        - Use the same file version for file
* Output file
    - List of promoters for DEGs

### Command

```R
library(tidyverse)
library(rtracklayer)
# Importing the gtf file
temp <- rtracklayer::import.gff("./gencode.V40.annotation.gtf")
temp <- as.data.frame(temp)

# Select transcripts for DEGs
## To only include the transcripts from the file
temp <- temp %>% filter(type == "transcript")
## Keep only the stable ID
temp$gene_id <- gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1",temp$gene_id)
## Select transcript records for the differentially expressed genes from section 4
temp <- temp[temp$gene_id %in% sig$Row.names,]

# Define promoter regions for positive strand
temp1 <- temp[temp$strand == "+",]
temp1$promoter_start <- temp1$start - 2000
temp1$promoter_end <- temp1$start + 1000
temp1$width = "3000"
temp1 = temp1 %>% select(seqnames,promoter_start,promoter_end,width,strand,gene_id,gene_name,transcript_name)

# Define promoter regions for negative strand
temp2 <- temp[temp$strand == "-",]
temp2$promoter_start <- temp2$end - 1000
temp2$promoter_end <- temp2$end + 2000
temp2$width <- "3000"
temp2 <- temp2 %>% select(seqnames,promoter_start,promoter_end,width,strand,gene_id,gene_name,transcript_name)

# Combine the files into the full promoter
temp <- rbind(temp1,temp2)
write.delim(temp,"promoter_DEGs.txt")
```

### Result

## Reference


