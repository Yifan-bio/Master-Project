#!/bin/bash
# 17 March 2023

##################################################################################
#
#     To run fastqc on files, this script uses the multi thread function to save time
#
#####################################################################################


##################################################################################
#
# Parameter setups
#
##################################################################################

# default variables
threads=4
depth=1

# help message
usage() {
    echo -e "Please note, this script is devlepoed under fastqc v0.11.9"
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-d <string>		Direcotry of input files"		
	echo -e "-o <string>		Directory for the output files (default: input directory(-d))"
    echo -e "-t <string>		Number of threads to be used (default: 4)"
	echo -e "-h			Print this help message and exit"
	echo -e "-m <string>		Maxdepth of the file to look for (default: 2; the current directory and it direct subfolders)"
	exit 1
}

# Configuring in the parameters from the input in command lines
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

# Check if fastqc packages exists
if command -v fastqc > /dev/null; then
    echo "Found fastqc command"
else
    echo "fastqc were not detected"
	exit -1
fi    

# check necessary parameters
if [[ -z $dir ]]; then
	echo -e "The directory of input files in missing"
	usage
	exit -1
fi

# Output file in the same directory as the input file if not specified
if [ -z "$outdir" ]; then
	outdir=$dir
fi

#Get absolute file path, so users can use relative/absolute as they like.
[[ ${dir} != "" ]] && dir=`realpath ${dir}`
[[ ${outdir} != "" ]] && outdir=`realpath ${outdir}`

##################################################################################
#
# Main
#
##################################################################################

# First, output the information of the run
echo "Using $threads threads, each using $memory, to run fastqc for the following files:"

# Creating a list to store the files that need to be processed through fastqc
list=""

# A function to locate all fastq or fq files under the queried directory 
function all_fq() {
    find $dir -maxdepth $depth -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz"
}

# Adding all the fastq / fq files found in the queired directory to the running list
for file in $(all_fq);
do
    echo "$file"
    list="$list $file"
done

# Running fastqc with multiple files in parallel
fastqc -t $threads -o $outdir $list

# Generate multiqc result for all the fastqc file if multiqc exists
if command -v multiqc > /dev/null; then
    multiqc $outdir
	echo "Multiqc result generated"
else
    echo "End of fastqc analysis. Multiqc were not detected."
fi    


exit 0;
