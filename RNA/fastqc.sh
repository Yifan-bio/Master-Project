# search for all file end with .fastq.gz, .fa.gz, .fa or fastq in the current directory and it subdirectories
function fastqc_all() {
    echo "fastqc_all"
    for file in `find . -name "*.fastq.gz" -o -name "*.fa.gz" -o -name "*.fa" -o -name "*.fastq"`
    do
        echo "fastqc $file"
        fastqc $file
    done
}
