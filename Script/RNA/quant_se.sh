while IFS= read -r line; do
    # Separating the file into read1 and read2
    fastq=`find $dir -maxdepth 4 -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz" \) -print | grep -i $line*`
    # file does not contain _R1 or _R2
    if [[ $fastq == *"_R1"* ]] || [[ $fastq == *"_R2"* ]]; then
        echo "The file seems to be paired end reads"
    else
        o="$out_dir/$line"
        echo "The full command of salmon been executed is as follows:
        salmon quant -i $index -l A -1 $Read -p 8 --validateMappings --gcBias --seqBias --posBias -o $o"
        # salmon quant -i $index -l A -1 $Read -p 8 --validateMappings --gcBias --seqBias --posBias -o $o
    fi
    fi
done < $input
