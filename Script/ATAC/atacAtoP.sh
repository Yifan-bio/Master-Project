# ! /bin/bash
# 2 August 2022

output_dir=
index=/mnt/f/support_doc/Gencode/salmon_pa_index/trancripts_index
dir=/mnt/f/Dataset/162nM_PMA_RNA_THP1
input=/mnt/f/Dataset/162nM_PMA_RNA_THP1/run.txt
genome_fasta=""

function find_read_pair() {
    # Separating the file into read1 and read2
    # Didnt use the $fastq as parameter as then it will only use one of the reads rather then both
    for file in $fastq;
    do
        if [[ $file == *"_R1"* ]]; then
            R1=$file
            echo "Read 1 file: $R1"
        fi
        if [[ $file == *"_R2"* ]]; then
            R2=$file
            echo "Read 2 file: $R2"
        fi
    done
}

function trim_QC() {
    # Trim the reads
    trim_galore --paired --retain_unpaired --output_dir $output_dir $R1 $R2
    # QC the trimmed reads
    fastqc $output_dir/*.fq.gz -o $output_dir
}

function Align() {
    bowtie2 --very-sensitive --end-to-end -p 4 --dovetail --no-mixed -X 2000 -t -x $index -1 $R1 -2 $R2 > ${line}.bowtie2.log | samtools sort -@ 4 -O bam -o ${line}.bam
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





while IFS= read -r line; do
    # Separating the file into read1 and read2
    fastq=`find $dir -maxdepth 4 -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz" \) -print | grep -i $line*`
    find_read_pair
    o="$out_dir/$line"
    echo "The full command of salmon been executed is as follows:
    salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --seqBias --posBias -o $o"
    # salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --seqBias --posBias -o $o
done < $input
