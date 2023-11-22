# Transcriptome analysis
> Repo containing details of transcriptome analysis

* [Package version](#package-version)
* [Transcriptome quantification](#transcriptome-quantification)
    + [Result](#result)
* [Differential expression](#differential-expression)
    + [Command](#command)
    + [Result](#result)
* [Enrichment analysis](#enrichment-analysis)
    + [Command](#command)
    + [Result](#result)    
* [Reference](#reference)

## Package version
`salmon v1.8.0`
`DESeq2 v1.34.0`
`tximport v1.22.0`
`biomaRt v2.50.3`
`tidyverse v1.3.1`

## Transcriptome quantification

Salmon were selected as the aligner for this study

> Conducted in terminal 

### Command

```sh
# Generating partial decoy for salmon
generateDecoyTranscriptome.sh -a $gtffile -o ${dir} -j 8 -g ${genome.fa} -t ${transcript.fa} 
# Creating salmon index
salmon index -t ${transcripts.fa} -i $transcripts_index --decoys ${decoys.txt} -k 31
```

```sh
# salmon quantification
salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --seqBias --recoverOrphans -o $output
```

### Result

The mapping of results using Salmon ranges between 50% and 60% (Plot 1A). We noticed that many reads were removed due to quality issues or were not recorded within the index. In section 0, we suggested a low likelihood of foreign contamination. Therefore, considering these results together, we believe that the contamination we are observing here is possibly DNA or ncRNA not recorded in the index. This is a possible outcome of rRNA depletion, the method used for RNA extraction, so there isn't strong evidence to worry about quality issues for the RNA sequencing alignment result.

<br />
<p align="center">
  <img width="700" src="./Figure/Plot1A.png">
</p>

_**Plot1A. Alignment rate and component of RNA-seq library.** The alignment rate of each sample varies between 50% to 60%. Paired Mapped reads and Orphan Recovered reads are considered as usable reads (The two green)._

<br />

According to the theory of RNA-seq, long transcripts should have more RNA-seq library reads after read fragmentation. However, there are also studies that suggest that long transcripts are generally expressed at a low level due to the conservation of energy. Therefore, we plotted Plot1C to examine the distribution of TPM versus transcript length. We found that the overall distribution appears highly similar.

<br />
<p align="center">
  <img width="700" src="./Figure/Plot1B.png">
</p>

_**Plot1B. TPM versus transcript length.**_

## Differential expression
> Conducted in R studio

Using DESeq2 to compare the change in gene expression due to PMA treatment.

### Command
Importing all transcriptome quantification from salmon into a gene-level quanitification with a dataframe format as input for DEseq2.

```R
library("tximport")

# Listing the path to each file.
files <- c('../Input/0hr_rep1/quant.sf',
           '../Input/0hr_rep2/quant.sf',
           '../Input/24hr_rep1/quant.sf',
           '../Input/24hr_rep2/quant.sf')
          
# Covert tx into gene          
txi.salmon <- tximport(files = files,
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = TRUE)
                      
# Providing condition name to each replicates
sampleTable <- data.frame(condition = factor(rep(c("control", "treated"),each = 2)))
rownames(sampleTable) <- colnames(txi.salmon$counts)  
```

Running DESeq2 for differential expression detection.

```R
library("DESeq2")
library("IHW")

# Importing the result from tximport to DESeq2
dds <- DESeqDataSetFromTximport(countData = txi.salmon,
                                colData = sampleTable,
                                design = ~ condition) 

# Perform Differential expression analysis using Wald test
dds <- DESeq(object = dds,
             test = 'Wald')

# Extracting results from the analysis, specifiy the significant cut is 0.05
result <- results(object = dds,
                  filterFun = ihw,
                  alpha = 0.05)

# Performing l2FC shrinkage to remove bias.
result <- lfcShrink(dds = dds,
                    coef = "condition_treated_vs_control",
                    type="apeglm",
                    res = result)
```

At this point, we will be done with DESeq2. Now it to make it easier for intepretation

```R
library(GenomicFeatures)
library(biomaRt)

# Only keep the stable gene ID
result$Row.names <- gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1",result$Row.names) 
# Getting the gene number of name
G_list = getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","hgnc_symbol",'entrezgene_id'),
               values = result$Row.names,
               mart= biomart)
# Adding the name and DESeq2 result into one dataframe
result = merge(result,G_list,by.x="Row.names",by.y="ensembl_gene_id") 
# Putting the accession number into transcripts without a official name
result$hgnc_symbol = ifelse(result$hgnc_symbol =='',result$Row.names,result$hgnc_symbol)  
# Remove duplicate recordings
result = result %>% distinct(Row.names,hgnc_symbol,log2FoldChange,.keep_all = T)

# Certain genes were not recprded in the biomart, so we use two databases
result$id = mapIds(org.Hs.eg.db,keys = result$Row.names,column = 'ENTREZID',keytype = 'ENSEMBL',multiVals = 'first')
# Combine the two results into one column
for(i in 1:nrow(result)) {
  if(is.na(result[i,'entrezgene_id'])){
    result[i,'entrezgene_id'] = result[i,'id']
  } 
}
```

Now adding in cutoff values to determine the final results

```R
library(tidyverse)
DEG = result %>% dplyr::filter(padj <0.05 & baseMean > 20 & abs(log2FoldChange) > 2) 

result$con = "non-significant"
result$con[result$log2FoldChange > 2 & result$padj < 0.05] = "upregulated"
result$con[result$log2FoldChange < -2 & result$padj < 0.05] = "downregulated"
```

### Result

<br />
<p align="center">
  <img width="700" src="./Figure/Plot4E.png">
</p>

To shows the number of differential expressed genes and the distribution, we will be plotting a volcanoi lot for DEGs to show there parameters. 

<br />
<p align="center">
  <img width="700" src="./Figure/Plot4C.png">
</p>

<br />
<p align="center">
  <img width="700" src="./Figure/Plot4A.png">
</p>

_**Figure 3.2. Volcano Plot.** L_




## Enrichment analysis

Ads now we have a bunch of DEGs, we will like to obtain the biological meaning of these DEGs through running enrichment analysis. But kegg pathways shows a limited of understanding on the impact so we will only be focusing on Gene ontology of the genes.

### Command

### Result

<br />
<p align="center">
  <img width="700" src="./Figure/Plot4B.png">
</p>


## Reference

