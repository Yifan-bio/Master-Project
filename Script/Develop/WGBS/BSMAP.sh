# ! /bin/bash
# 29 August 2023
##########################################################################
echo "#########################################################################"
echo "BSMAP analysis"
echo "#########################################################################"
##########################################################################
#set -o pipefail
#export LC_ALL=C

# Prerequisite:
# - BSMAP
# - samtools
# - MethylDackel


##########################################################################
#                          Adding paramters
##########################################################################

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-i <string>		BSMAP Index"		
	echo -e "-r1 <string>		Read 1 file of paired-end read / OR single end file"
    echo -e "-r2 <string>		Read 2 file of paired-end read"
    echo -e "-p <string>        Name of the sample to be recorded as"
    echo -e ""
    echo -e "Optional parameters:"
    echo -e "-t <string>		Number of threads to be used (default: 6)
                    It is more ideal to provide a multiple of 3"
    echo -e "-m <string>		mode to run"
	echo -e "-o <string>		Output directory with all salmon results (default: current directory)"
    echo -e "-s <string>		Directory to salmon binary file if it is not in the PATH"
	exit 1
}

#Initiate parameters with NULL
WDIR="."
extra_args=""
threads=6
mode="paired-end"

while getopts ":i:r1:r2:o:p:s:m:t:" op; do
	case $op in
		i) index=${OPTARG} ;;
		r1) read1=${OPTARG} ;;
        r2) read2=${OPTARG} ;;
		o) WDIR=${OPTARG} ;;
        p) prefix=${OPTARG} ;;
        s) salmon=${OPTARG} ;;
        m) mode=${OPTARG} ;;
        #e) extra_args=${OPTARG} ;;
        t) threads=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# check necessary parameters
if [[ -z $index ]] || [[ -z $read1 ]] || [[ -z $prefix ]]; then
	echo -e "Guess which 3 is required"
	usage
	exit -1
fi

# check for single end mode
if [[ -z $read2 ]]; then  
    if [[ -n $R2 ]]; then
        echo -e "This is running as single end mode"
        mode="single-end"
    fi
fi

#Get absolute file path, so users can use relative/absolute as they like.
[[ ${index} != "" ]] && Index=`realpath ${index}`
[[ ${} != "" ]] && GenomeFasta=`realpath ${sample_list}`
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



Alignment() {
    bsbolt Align -F1 $r1 -F2 $r2 -DB $i -O $WDIR/$prefix -OT -t 
}


trim_galore --paired

bsmap -a $R1 -b $R2 -d $ref -o $output.sam
samtools view -bS $output.sam -o output.bam
samtools sort output.bam
samtools view -b chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18  chr19 chr20 chr21 chr22 chrX chrY
picard Markduplicates I O M
MethylDackel
