# fixing HMM results
awk -F"\t" -v OFS="\t" '$4 != "HighCoveragePeak_0" {print $0}' ./DiffBind/Input/untreat_rep1_peak_peaks.gappedPeak > ./Peak/Input/untreat_rep1.gappedPeak
awk -F"\t" -v OFS="\t" '$4 == "HighCoveragePeak_0" {print $0}' ./DiffBind/Input/untreat_rep1_peak_peaks.gappedPeak > ./Peak/Input/untreat_rep1_HighCoverage.gappedPeak

# finding overlapping peaks
bedtools intersect -wa -a ./HMM/treat_rep2.gappedPeak -b ./HMM/treat_rep1.gappedPeak > ./treat1.bed
bedtools intersect -wa -a ./HMM/treat_rep1.gappedPeak -b ./HMM/treat_rep2.gappedPeak > ./treat2.bed
cat treat1.bed treat2.bed > treat.bed
sort -k1,1 -k2,2n treat.bed > treat.peak
bedtools merge -i treat.peak > treat.bed

# getting TSS Followed by promoter region
zcat gencode.vM19.annotation.gtf.gz | awk 'OFS="\t" {if ($3=="transcript") {if ($7 == "+") {print $1,$4-1,$4,$12,".",$7} else {print $1,$5-1,$5,$12,".",$7}}}' | tr -d '";' | sort -k1,1V -k2,2n > gencode.vM19.annotation.tss.bed
awk -v OFS="\t" '{print $1,$2-999,$2,$4,$5,$6}' /mnt/e/gencode.v40.annotation.tss.bed > /mnt/e/gencode.v40.annotation.promoter.bed

# Getting overlap peaks overlap promoter regions
bedtools intersect -a /mnt/e/gencode.v40.annotation.promoter.bed -b /mnt/e/ATAC/Peak/Input/treat.bed > /mnt/e/ATAC/Peak/Input/treat_promoter.bed

# Getting tpm for each condition and do some calculation
awk -v OFS="\t" '{print $1,$6,$7,$6+$7,($6+$7)/2}' iso.txt > treat_tpm.txt
# Now on to r script peak_annot.R