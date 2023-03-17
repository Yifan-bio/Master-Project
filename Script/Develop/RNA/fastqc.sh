#!/bin/bash
# 23 September

##################################################################################
#
#     To run fastqc on files, this script uses the multi thread function to save time
#
#####################################################################################

# default variables
threads=4
depth=1

# Input paramters
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-d <string>		Direcotry of input files"		
	echo -e "-o <string>		Directory for the output files (default: input directory(-d))"
    echo -e "-t <string>		Number of threads to be used (default: 4)
	                Each thread requires 250MB of RAM memory and process 1 file (e.g. 10 threads will use 2500MB and process 10 file simultaneously)"
	echo -e "-h			Print this help message and exit"
	echo -e "-m <string>		Maxdepth of the file to look for (default: 2; the current directory and it direct subfolders)"
	exit 1
}



while getopts ":d:o:t:m:h:" op; do
	case $op in
        d) dir=${OPTARG} ;;
        o) outdir=${OPTARG} ;;
        t) threads=${OPTARG} ;;
		m) depth=${OPTARG} ;;
		h) usage ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# check necessary parameters
if [[ -z $dir ]]; then
	echo -e "Guess which ONE is required"
	usage
	exit -1
fi

# Output file in the same directory as the input file
if [ -z "$outdir" ]; then
	outdir=$dir
fi

#Get absolute file path, so users can use relative/absolute as they like.
[[ ${dir} != "" ]] && dir=`realpath ${dir}`
[[ ${outdir} != "" ]] && outdir=`realpath ${outdir}`

# First, output the information of the run
echo "Using $threads threads to run fastqc for the following files:"

# Find all files under the currect directory and sub directory with the .fastq.gz,.fastq,fq,fq.gz extension
function all_fq() {
    find $dir -maxdepth $depth -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz"
}

# Creating the variable for putting all of the files to run
list=""

# Runs all the fastqc on the files found in the all_fq function
for file in $(all_fq);
do
    echo "$file"
    list="$list $file"
done

echo "############################################################"
echo "############################################################"
echo "############################################################"

echo "full command: "
echo "fastqc -t $threads -o $outdir $list"

echo "############################################################"
echo "############################################################"
echo "############################################################"

echo "Rest is log from fastqc:"

fastqc -t $threads -o $outdir $list

exit 0;

# # Find all files under the currect directory and sub directory with the .fastq.gz,.fastq,fq,fq.gz extension
# function all_fq() {
#     find . -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz"
# }

# # Runs all the fastqc on the files found in the all_fq function
# for file in $(all_fq);
# do
#     echo $file
#     fastqc "$file"
# done
