# ! /bin/bash
# 5 August 2022
##########################################################################
echo "#########################################################################"
echo "RNA sequencing pseudoalignment analysis"
echo "#########################################################################"
##########################################################################
#set -o pipefail
#export LC_ALL=C
#
# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-i <string>		Index used for salmon"		
	echo -e "-l <string>		sample list file that need to run"
    echo -e "-f <string>		Directory of the cDNA files (or can be 3 directry higher)"
	echo -e "-o <string>		Output directory with all result files"
    echo -e "-s <string>		Directory to salmon binary file if it is not in the PATH"
	exit 1
}

#Initiate parameters with NULL
index=""
sample_list=""
dir=""
WDIR="."
salmon="salmon"

while getopts ":i:l:o:f:s:" op; do
	case $op in
		i) index=${OPTARG} ;;
		l) sample_list=${OPTARG} ;;
		o) WDIR=${OPTARG} ;;
        f) dir=${OPTARG} ;;
        s) salmon=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))
#
# check necessary parameters
if [[ -z $index ]] || [[ -z $sample_list ]] || [[ -z $dir ]]; then
	echo -e "Guess which 3 is required"
	usage
	exit -1
fi
#
#Get absolute file path, so users can use relative/absolute as they like.
[[ ${index} != "" ]] && Index=`realpath ${index}`
[[ ${sample_list} != "" ]] && GenomeFasta=`realpath ${sample_list}`
[[ ${dir} != "" ]] && cDNA_file=`realpath ${dir}`
[[ ${WDIR} != "" ]] && WDIR=`realpath ${WDIR}`
#
# make directories if not exist and enter working directory.
[[ ! -d ${WDIR} ]] && mkdir -p ${WDIR}
cd ${WDIR}

function find_ena_pair() {
        # Separating the file into read1 and read2
    # Didnt use the $fastq as parameter as then it will only use one of the reads rather then both
    for file in $fastq;
    do
        if [[ $file == *"_R1."* ]] || [[ $file == *"_2."* ]]; then
            R1=$file
            echo "Read 1 file: $R1"
        fi
        if [[ $file == *"_R2."* ]] || [[ $file == *"_2."* ]]; then
            R2=$file
            echo "Read 2 file: $R2"
        fi
    done
}

function run_salmon() {
    echo "The full command of salmon been executed is as follows:
    $salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --seqBias --posBias -o $o"
    #$salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --seqBias --posBias -o $o
}

while IFS= read -r line; do
    # Separating the file into read1 and read2
    fastq=`find $dir -maxdepth 4 -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz" \) -print | grep -i $line*`
    find_ena_pair
    o="$WDIR/$line"
    run_salmon
done < $sample_list




# function verify_pair() { 
#     if [[ ${R1} == ${R2} ]]; then
#         echo "Error: R1 and R2 are the same file"
#         exit -1
#     fi
#     # check if R1 and R2 is not the same format
#     R1_R2_diff=`diff -q ${R1} ${R2} | wc -l`
#     if [[ "${R1_R2_diff}" != "1" ]]; then
#         echo "Error: R1 and R2 are not the same length (So likely not the same pair)"
#         echo "R1: ${R1}"
#         echo "R2: ${R2}"
#         exit -1
#     fi
#     # Still in progress
# }