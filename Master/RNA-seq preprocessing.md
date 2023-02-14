# RNA sequencing preprocessing
> Preprocessing of RNA-sequencing including the quantification of transcriptome.

* [Package version](#package-version)
* [Transcriptome quantification](#transcriptome-quantification)
    + [Command](#command)
* [Reference](#reference)

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
trim_galore --paired --output_dir $trim_dir $R1 $R2
# hisat2 alignment
hisat2 -p 8 -1 $R1 -2 $R2 -S ./Output.sam
```

### Plotting

## Reference
