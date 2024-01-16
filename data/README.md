# All data invovled

## Supporting datasets

> This section records all supporting files for the study.

- All reference file were downloaded from the Gencode project using the GRCh38 human reference and transcript version 40
    - [GRCh38 primary reference](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz)
    - [GRCh38 version 40 transcript reference](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.transcripts.fa.gz)
    - [GRCh38 version 40 comprehensive primary annotation(gtf)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.primary_assembly.annotation.gtf.gz)

- Files involved in the ATAC-seq preprocessing
    - [ENCODE blacklist project version 2](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz) can be downloaded from the [Boyle-Lab](https://github.com/Boyle-Lab/Blacklist)
    - [GeneHancer database](https://genome.ucsc.edu/cgi-bin/hgTables) can be downloaded from the [UCSC genome browser](https://genome.ucsc.edu/index.html). 
        - Unfortunately, the full list cannot be directly downloaded from the browser. So we created a [bed file]() including all base pairs in the main chromosomes and use the UCSC browser to download records fall under these regions.
- ChIP atlas file
    - The histone marker file that are involved in this study includes five replicates from three indepdent studies:
        - [GSM3514952](https://chip-atlas.org/view?id=SRX5141106)
        - [GSM3514953](https://chip-atlas.org/view?id=SRX5141108)
        - [GSM2544236](https://chip-atlas.org/view?id=SRX2655451)
        - [GSM2544237](https://chip-atlas.org/view?id=SRX2655452)
        - [GSM4083810](https://chip-atlas.org/view?id=SRX6866214)

## Analysed sequencing samples

> The link to obtain the samples from the ENA database

All samples are obtained from the bioProject [PRJNA533829](https://www.ebi.ac.uk/ena/browser/view/PRJNA533829)
- ATAC sequencing sample
    - SRR8932925 (0 hour replicate 1)
    - SRR8932926 (24 hour replicate 1)   
    - SRR8932927 (0 hour replicate 2)     
    - SRR8932928 (24 hour replicate 2)
- RNA sequencing sample
    - SRR8932929 (24 hour replicate 2)
    - SRR8932930 (0 hour replicate 1)
    - SRR8932934 (24 hour replicate 1)
    - SRR8932935 (0 hour replicate 2)

## Quality Assessed samples

## Link to Google drive storing intermediate results.
