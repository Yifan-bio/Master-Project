# Accessibility of Regulatory Element
> This section records the accessibility status of the regulatory elements of DEGs

* [Package version](#package-version)
* [Accessibility of regulatory elements](#accessibility-of-regulatory-elements)
    + [Input and Output](#input-and-output)
    + [Command](#command)
    + [Result](#result)
* [Reference](#reference)

## Package version
`tidyverse v1.3.1` `bedtools`

## Accessibility of regulatory elements

### Input and Output

* Input file
    - List of promoters for DEGs
    - List of active enhancer for DEGs
    - List of Differential accessibile regions from section 5
* Output file
    - Regulatory element of DEGs with accessible changes

### Command

```R
# Input all required files
enhancer <- read.delim("./DEG_enhancer.txt",header = F)
promoter <- read.delim("./promoter_DEGs.txt")
ATAC_sig <- read.delim



```

### Result

Creating a plot to show the enrichment of DEGs with accessibility change in regulators


<br />
<p align="center">
  <img width="700" src="./Figure/Plot7A.png">
</p>

<br />
<p align="center">
  <img width="700" src="./Figure/Plot7B.png">
</p>

<br />
<p align="center">
  <img width="700" src="./Figure/Plot7C.png">
</p>

## Reference

