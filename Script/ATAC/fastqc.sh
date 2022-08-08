#!/bin/bash
# 8 August

usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-d <string>		Direcotry of the files that want to be run"		
	echo -e "-o <string>		Directory of the output files"
    echo -e "-t <string>		Number of threads to be used (default: 4) (each threads uses 250MB)"
	exit 1
}

# default variables
dir="."
outdir='.'
t=4

while getopts ":d:o:t:" op; do
	case $op in
        d) dir=${OPTARG} ;;
        o) outdir=${OPTARG} ;;
        t) threads=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# check necessary parameters
# if [[ -z $dir ]]; then
# 	echo -e "Guess which ONE is required"
# 	usage
# 	exit -1
# fi
#
#Get absolute file path, so users can use relative/absolute as they like.
[[ ${dir} != "" ]] && Index=`realpath ${dir}`
[[ ${outdir} != "" ]] && Index=`realpath ${outdir}`

# First, output the information of the run
echo "Using $t threads to run fastqc for the following files:"

# Find all files under the currect directory and sub directory with the .fastq.gz,.fastq,fq,fq.gz extension
function all_fq() {
    find $dir -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz"
}

# Creating the variable for putting all of the files to run
list=""

# Runs all the fastqc on the files found in the all_fq function
for file in $(all_fq);
do
    echo "$file"
    list="$list $file"
done

echo "$list"
fastqc -t $t -o $outdir $list

exit 0;