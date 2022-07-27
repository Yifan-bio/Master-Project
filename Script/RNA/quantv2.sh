# ! /bin/bash
# 27 July 2022



out_dir=.
index=/media/studentsgh129/project/support_doc/Gencode/salmon_pa_index/trancripts_index
dir=/media/studentsgh129/project/Dataset/162nM_PMA_RNA_THP1
input=/media/studentsgh129/project/Dataset/162nM_PMA_RNA_THP1/run.txt

while IFS= read -r line; do
    # Finding all files that ends with .fastq.gz and contains the accession number as name
    fastq=`find $dir -maxdepth 4 -type f -name "*.fastq.gz" -print | grep -i $line*`
    # Separating the file into read1 and read2
    for file in $fastq
    do
        if [[ $file == *"R1"* ]]; then
            R1=$file
            echo "Read 1 file: $R1"
            continue
        fi
        if [[ $file == *"R2"* ]]; then
            R2=$file
            echo "Read 2 file: $R2"
      fi
    done
    o="$out_dir/$line"
    echo "The full command of salmon been executed is as follows:
    BAM
    BAM
    BAM
    salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --seqBias --posBias -o $o"
    salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --seqBias --posBias -o $o
done < $input
