
# Find all files under the currect directory and sub directory with the .fastq.gz,.fastq,fq,fq.gz extension
function all_fq() {
    find . -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz"
}

# Runs all the fastqc on the files found in the all_fq function
for file in $(all_fq);
do
    echo $file
    fastqc "$file"
done
