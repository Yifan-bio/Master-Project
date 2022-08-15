# ! /bin/bash
# 9 August 2022
##########################################################################
echo "#########################################################################"
echo "RNA sequencing pseudoalignment analysis"
echo "#########################################################################"
##########################################################################
#set -o pipefail
#export LC_ALL=C

# Prerequisite:
# - salmon

# Function to do
# - add extra paramters for analysis
# - verify the salmon binary is available
# - verify the two files are paired-end

##########################################################################
#                          Adding paramters
##########################################################################

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-i <string>		Salmon Index"		
	echo -e "-l <string>		txt file with Accession IDs of required samples"
    echo -e "-f <string>		Directory with library files"
    echo -e ""
    echo -e "Optional parameters:"
    echo -e "-t <string>		Number of threads to be used (default: 4)"
    echo -e "-m <string>		Max depth to search from the library directory (-f) 
                        (default: 2; Search the directory and one subdirectory below)"
	echo -e "-o <string>		Output directory with all salmon results (default: current directory)"
    echo -e "-s <string>		Directory to salmon binary file if it is not in the PATH"
    #echo -e "-e <string>		Adding more paramters rather then default I added (default: none)"
	exit 1
}

#Initiate parameters with NULL
depth=2
WDIR="."
salmon="salmon"
extra_args=""
threads=4
mode="paired-end"

while getopts ":i:l:m:o:f:s:m:t:" op; do
	case $op in
		i) index=${OPTARG} ;;
		l) sample_list=${OPTARG} ;;
        m) depth=${OPTARG} ;;
		o) WDIR=${OPTARG} ;;
        f) dir=${OPTARG} ;;
        s) salmon=${OPTARG} ;;
        m) mode=${OPTARG} ;;
        #e) extra_args=${OPTARG} ;;
        t) threads=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# check necessary parameters
if [[ -z $index ]] || [[ -z $sample_list ]] || [[ -z $dir ]]; then
	echo -e "Guess which 3 is required"
	usage
	exit -1
fi

# check for single end mode
if [[ $mode == "single-end" ]]; then  
    if [[ -n $R2 ]]; then
        echo -e "R2 file is provided but single-end mode is selected"
        usage
        exit -1
    fi
fi

#Get absolute file path, so users can use relative/absolute as they like.
[[ ${index} != "" ]] && Index=`realpath ${index}`
[[ ${sample_list} != "" ]] && GenomeFasta=`realpath ${sample_list}`
[[ ${dir} != "" ]] && cDNA_file=`realpath ${dir}`
[[ ${WDIR} != "" ]] && WDIR=`realpath ${WDIR}`

# make directories if not exist and enter working directory.
[[ ! -d ${WDIR} ]] && mkdir -p ${WDIR}
cd ${WDIR}

##########################################################################
#                               Functions
##########################################################################

# ============================================================
# Primary functions
# ============================================================

# Separate the read1 and read2 files for the sample
function find_ena_pair() {
    for file in $fastq;
    do
        if [[ $file == *"_R1.f"* ]] || [[ $file == *"_2.f"* ]]; then
            R1=$file
            echo "Read 1 file: $R1"
        fi
        if [[ $file == *"_R2.f"* ]] || [[ $file == *"_2.f"* ]]; then
            R2=$file
            echo "Read 2 file: $R2"
        fi
    done
}

# Running paired end analysis
function paired_end() {
    echo "The full command of salmon been executed is as follows:
    $salmon quant -i $index -p $threads -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --seqBias --recoverOrphans -o $o"
    $salmon quant -i $index -l A -1 $R1 -2 $R2 -p 8 --validateMappings --gcBias --seqBias --recoverOrphans -o $o
}

# Running single end analysis
function single_end() {
    echo "The full command of salmon been executed is as follows:
    $salmon quant -i $index -l A -r $R1 -p 8 --validateMappings --seqBias -o $o"
    $salmon quant -i $index -l A -r $R1 -p 8 --validateMappings --seqBias -o $o
}

# function to block the logs
function log_block() {
    echo "#########################################################################"
    echo "#########################################################################"
    echo "#########################################################################"
}

# ============================================================
# Secondary functions
# ============================================================

function pair() {
    while IFS= read -r line; do
        # Separating the file into read1 and read2
        fastq=`find $dir -maxdepth $depth -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz" \) -print | grep -i $line*`
        find_ena_pair
        o="$WDIR/$line"
        paired_end
        log_block
    done < $sample_list
}

function single() {
    while IFS= read -r line; do
        # Separating the file into read1 and read2
        fastq=`find $dir -maxdepth $depth -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz" \) -print | grep -i $line*`
        o="$WDIR/$line"
        single_end
        log_block
    done < $sample_list
}

##########################################################################
#                               Main
##########################################################################

if [[ $mode == "paired-end" ]]; then
    echo "The mode is paired-end"
    pair
elif [[ $mode == "single-end" ]]; then
    echo "The mode is single-end"
    single
else
    echo "The mode is not specified correctly"
    usage
fi

##########################################################################
#                               End
##########################################################################
# The following is a proposed function

# Verify the pairs
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

# Verify the presence of salmon and it function