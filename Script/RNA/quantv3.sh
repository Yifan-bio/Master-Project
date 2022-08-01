
# ! /bin/bash
# 1 August 2022

out_dir=.
index=/mnt/f/support_doc/Gencode/salmon_pa_index/trancripts_index
dir=/mnt/f/Dataset/162nM_PMA_RNA_THP1
input=/mnt/f/Dataset/162nM_PMA_RNA_THP1/run.txt

function find_ena_pair() {
        # Separating the file into read1 and read2
    # Didnt use the $fastq as parameter as then it will only use one of the reads rather then both
    for file in $fastq;
    do
        if [[ $file == *"_1"* ]]; then
            R1=$file
            echo "Read 1 file: $R1"
        fi
        if [[ $file == *"_2"* ]]; then
            R2=$file
            echo "Read 2 file: $R2"
        fi
    done
}

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

while IFS= read -r line; do
    # Separating the file into read1 and read2
    fastq=`find $dir -maxdepth 4 -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz" \) -print | grep -i $line*`
    find_read_pair
    o="$out_dir/$line"
    echo "The full command of salmon been executed is as follows:
    salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --recoverOrphans --gcBias --seqBias --posBias -o $o"
    # salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --recoverOrphans --seqBias --posBias -o $o
done < $input
