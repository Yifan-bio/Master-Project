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
    echo -e "-f <string>		Directory with library files"
    echo -e "-m <string>        PE for pair-end; SE for single-end reads (default: PE)"
    echo -e ""
    echo -e "Optional parameters:"
    echo -e "-e <string>        Extension type (default: fq.gz)"
    echo -e "-t <string>		Number of threads to be used (default: 4)"
	echo -e "-o <string>		Output directory with all salmon results (default: current directory)"
    echo -e "-s <string>		Directory to salmon binary file if it is not in the PATH"
    echo -e 'default command: $salmon quant -i index -l A -1 R1 -2 R2 -p threads --validateMappings --gcBias --seqBias --recoverOrphans -o output'
	exit 1
}

#Initiate parameters with NULL
WDIR="."
salmon="salmon"
extra_args=""
threads=4
extension="fq.gz"
mode="PE"

while getopts ":i:f:m:e:o:s:t:" op; do
	case $op in
		i) index=${OPTARG} ;;
        f) dir=${OPTARG} ;;
        m) mode=${OPTARG} ;;
        e) extension=${OPTARG} ;;
		t) threads=${OPTARG} ;;
        o) WDIR=${OPTARG} ;;
        s) salmon=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# check necessary parameters
if [[ -z $index ]] || [[ -z $dir ]]; then
	echo -e "Either the index or file dir is wrong"
	usage
	exit -1
fi

#Get absolute file path, so users can use relative/absolute as they like.
[[ ${index} != "" ]] && index=`realpath ${index}`
[[ ${dir} != "" ]] && dir=`realpath ${dir}`
[[ ${WDIR} != "" ]] && WDIR=`realpath ${WDIR}`

# make directories if not exist and enter working directory.
[[ ! -d ${WDIR} ]] && mkdir -p ${WDIR}
cd ${WDIR}

##########################################################################
#                               Main
##########################################################################

# Paired-end salmon quantifcation
if [[ ${mode} == "PE" ]]; then
    
    # PE salmon for fq.gz files
    if [[ ${extension} == "fq.gz" ]]; then
        # Detecting the files name exist in directory
        for i in $(ls ${dir}/*.fq*.gz | sed 's/[1-2].fq.gz//' | uniq); do 
        # Testing if the two files exist
        if [[ -z ${dir}/${i}1.fq.gz ]] || [[ -z ${dir}/${i}2.fq.gz ]]; then
            echo -e "Read 1 and Read 2 file cannot be detected for $i"
            usage
            exit -1
        fi
        # Output folder names
        o="$WDIR/$(basename $i)"
        # Running salmon
        "$salmon quant -i $index -l A -1 ${i}1.fq.gz -2 ${i}2.fq.gz -p $threads --validateMappings --gcBias --seqBias --recoverOrphans -o $o"
        done

    # PE salmon for fastq.gz files
    elif [[ ${extension} == "fastq.gz" ]]; then
        # Detecting the files name exist in directory
        for i in $(ls ${dir}/*.fastq*.gz | sed 's/[1-2].fastq.gz//' | uniq); do 
        # Testing if the two files exist
        if [[ -z ${dir}/${i}1.fastq.gz ]] || [[ -z ${dir}/${i}2.fastq.gz ]]; then
            echo -e "Read 1 and Read 2 file cannot be detected for $i"
            usage
            exit -1
        fi
        # Output folder names
        o="$WDIR/$(basename $i)"
        # Running salmon
        "$salmon quant -i $index -l A -1 ${i}1.fastq.gz -2 ${i}2.fastq.gz -p $threads --validateMappings --gcBias --seqBias --recoverOrphans -o $o"
        done

    # Reporting error in file extension
    else
        echo -e "Only accept file extension for fastq.gz and fq.gz"
        exit -1
    fi

# Single-end salmon quantification
elif [[ ${mode} == "SE" ]]; then
    # Testing if without extension can still run properly
    # Detecting the files exist in directory
    for i in $(ls ${dir}/*.f*q.gz | sed 's/.f[^~]*q.gz//'); do 
        # Output folder names
        o="$WDIR/$(basename $i)"
        # Running salmon
        $salmon quant -i $index -l A -r $i.fastq.gz -p $threads --validateMappings -o $o
    done

# Reporting error in read type selection
else
    echo -e "Enter either <PE> or <SE>"
    exit -1

fi
##########################################################################
#                               End
##########################################################################