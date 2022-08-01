
index=
read1=
read2=
output_dir=

function trim() {
    trim_galore --paired --retain_unpaired --output_dir $output_dir $read1 $read2
}

function Align() {
    bowtie2 --very-sensitive --end-to-end -p 6 --dovetail --no-discordant -X 2000 -k 5 -t -x $index  -1 $read1 -2 $read2 -U $unpaired1,$unpaired2 -S align.sam
}

function MT() {
    samtools idxstats align.sam | cut -f 1 | grep -v MT | xargs samtools view -b align.sam > MT.bam
}

function Dedup() {
    picard MarkDuplicates Input=MT.bam Output=dedup.bam METRICS_FILE=metrics.txt REMOVE_DUPLICATES=true
}

function filter() {
    samtools view -h -@ 6 dedup.bam | grep -v MT | samtools sort -O bam -o dedup_MT.bam
    samtools view -h -b -@ 6 -q 30 -F 1804 -f 2 -o out.bam dedup_MT.bam
}

function HMM_Peak_Caller() {

}

function footprint_analysis() {

}
