# RNA sequencing preprocessing
> Preprocessing of RNA-sequencing including the quantification of transcriptome.

* [Package version](#package-version)
* [Transcriptome quantification](#transcriptome-quantification)
* [](#queried-dataset)

## Package version
`salmon v1.8.0`
`hisat2`

## Transcriptome quantification

### Command
The study used salmon as transcriptome quantifier.

```sh
# salmon quantification
salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --seqBias --recoverOrphans -o $output
```

A second alignment were performed for the reads to identify unmapped reads source using hisat2
```sh
# trimming

# hisat2 alignment
```

### Plotting

## Reference
