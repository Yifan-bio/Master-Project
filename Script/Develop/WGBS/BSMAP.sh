trim_galore --paired

bsmap -a $R1 -b $R2 -d $ref -o $output.sam
samtools view -bS $output.sam -o output.bam
samtools sort output.bam
samtools view -b chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18  chr19 chr20 chr21 chr22 chrX chrY
picard Markduplicates I O M
MethylDackel
