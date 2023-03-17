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
# - BSbolt
# - samtools
# - MethylDackel


##########################################################################
#                          Adding paramters
##########################################################################

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-i <string>		BSbolt Index"		
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

# Calculating the number of threads for Align
thread1=threads

##########################################################################
#                               Functions
##########################################################################

# ============================================================
# Primary functions
# ============================================================

Alignment() {
    bsbolt Align -F1 $r1 -F2 $r2 -DB $i -O $WDIR/$prefix -OT -t 
}
