# ! /bin/bash
# 2 August 2022

##########################################################################

echo "ATAC-seq analysis"

##########################################################################
set -ue
set -o pipefail
export LC_ALL=C

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
    echo -e "This script is used to analyze ATAC-seq data. It will trim the reads, align the reads, remove mitochondrial reads and duplicates, and call peaks."
    echo -e "The input files are paired-end reads in fastq format, and the output files are in bam format."
    echo -e "Default required software: Cutadapt, bowtie2, samtools, picard, HMMRATAC"
    echo -e ""
    echo -e "Required Input paramters: Read1, Read2, Index, GenomeFasta, readlength"
	echo -e "-r1 <string>		    Read1 file"
	echo -e "-r2 <string>		    Read2 file"
	echo -e "-i <string>		    Index used for bowtie2"		
	echo -e "-readlength <string>	    Read length of ATAC library"
	echo -e "-g <string>		    Genome files used to make bowtie2 index file"	

    echo -e ""
    echo -e "Customise input paramters: prefix, blacklist, WDIR, mode, peak_caller, adapter sequence"
 	echo -e "-b <string>		    Blacklist file from ENCODE project (boyle-lab)"       
    echo -e "-t <string>                 Number of threads (default: 4)"    
	echo -e "-p <string>		    Prefix for all the output files"
    echo -e "-a <string>                 Adapter sequence for trimming (CTGTCTCTTATA)"
	echo -e "-o <string>		    Output directory with all result files"
	echo -e "-m <string>		    Mode for analysis; there is strict and lenient (default: strict). Strict will remove chrM,duplicate,low-quality,multimap and improper mapped. While lenient will only run remove chrM and deduplicate"
	echo -e "-peak-caller <string>	    Peak caller selection between Genrich and HMMRATAC (default: HMMRATAC)"
	exit 1
}

#Initiate parameters with NULL
blacklist=""
GenomeFasta=""
prefix="test"
mode="strict"
WDIR="."
readlength="75"
peak_caller="HMMRATAC"
adapter="CTGTCTCTTATA"
threads=4

while getopts ":r1:r2:i:l:g:b:o:p:a:t:m:" op; do
	case $op in
        # Required parameters
		r1) R1=${OPTARG} ;;
		r2) R2=${OPTARG} ;;
		i) Index=${OPTARG} ;;
		l) readlength=${OPTARG} ;;        
		g) GenomeFasta=${OPTARG} ;;
		
        # Customise parameters
        b) blacklist=${OPTARG} ;;
		o) WDIR=${OPTARG} ;;
        p) prefix=${OPTARG} ;;
        a) adapter=${OPTARG} ;;
        t) threads=${OPTARG} ;;
		m) mode=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# check necessary parameters
if [[ -z $R1 ]] || [[ -z $R2 ]] || [[ -z $Index ]] || [[ -z $GenomeFasta ]] || [[ -z $readlength ]]; then
    echo -e "No index, no reads, no analysis"
    usage
    exit -1
fi

# For readlength, if it is longer than 100, then set it to 100
if [[ $readlength > 100 ]]; then
    $readlength=100
fi

#Get absolute file path, so users can use relative/absolute as they like.
[[ ${R1} != "" ]] && R1=`realpath ${R1}`
[[ ${R2} != "" ]] && R2=`realpath ${R2}`
[[ ${Index} != "" ]] && Index=`realpath ${Index}`
[[ ${GenomeFasta} != "" ]] && GenomeFasta=`realpath ${GenomeFasta}`
[[ ${blacklist} != "" ]] && blacklist=`realpath ${blacklist}`
[[ ${WDIR} != "" ]] && WDIR=`realpath ${WDIR}`

# make directories if not exist and enter working directory.
[[ ! -d ${WDIR} ]] && mkdir -p ${WDIR}
cd ${WDIR}


##############################################################
# Pipeline function setup starts here
##############################################################

# 1. Trimming the reads
function trim_QC() {
    # Trim the reads
    trim_dir=${WDIR}/trim
    [[ ! -d ${trim_dir} ]] && mkdir -p ${trim_dir}
    Cutadapt -j $threads -a $adapter -A $adapter  -q 20 -O 1 --minimum-length 20:20 --pair-filter=any -o ${trim_dir}/${prefix}_trim_R1.fq.gz -p ${trim_dir}/${prefix}_trim_R2.fq.gz $R1 $R2
    Read1=${trim_dir}/${prefix}_trim_R1.fq.gz
    Read2=${trim_dir}/${prefix}_trim_R2.fq.gz
}

# 2. Alignment
function Align() {
    bowtie2 --very-sensitive --end-to-end -p $threads --dovetail --no-mixed -X 2000 -t -x $Index -1 $Read1 -2 $Read2 > ${prefix}.bowtie2.log | samtools sort -@ $threads -O bam -o ${prefix}.bam
}

# 3. remove mitochondrial reads and duplicates
function QC_removal() {

    samtools view -@ $threads -h ${prefix}.bam | grep -v chrM | samtools sort -@ $threads -O bam -o ${prefix}.rmChrM.bam
    picard MarkDuplicates Input=${prefix}.rmChrM.bam Output=${prefix}.rmChrM_dedup.bam METRICS_FILE=${prefix}.dedup.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
    
    if [[ $mode == "strict" ]]; then
        samtools view -h -b -q 30 -@ $threads -F 1804 -f 2 -o ${prefix}.rmChrM_dedup_QC.bam ${prefix}.rmChrM_dedup.bam
        samtools index ${prefix}.rmChrM_dedup_QC.bam
    elif [[ $mode == "lenient" ]]; then
        samtools index ${prefix}.rmChrM_dedup.bam
    fi
}

# 4. Peak calling
function peak_calling() {
	if [[ $peak_caller == "Genrich" ]]; then
        if [[ $mode == "strict" ]]; then
            Genrich -t ${prefix}.rmChrM_dedup_QC.bam -o ${prefix}.Genrich.narrowPeak -j -y -r -v
        elif [[ $mode == "lenient" ]]; then
            Genrich -t ${prefix}.rmChrM_dedup.bam -o ${prefix}.Genrich.narrowPeak -j -y -r -v
        fi
    elif [[ $peak_caller == "HMMRATAC" ]]; then
        if [[ $mode == "strict" ]]; then
            HMMRATAC -b ${prefix}.rmChrM_dedup_QC.bam -i ${prefix}.rmChrM_dedup_QC.bam.bai -g $GenomeFasta -o ${prefix}.hmmratac -e $blacklist -m $readlength,200,400,600 --window 5000000
        elif [[ $mode == "lenient" ]]; then
            HMMRATAC -b ${prefix}.rmChrM_dedup.bam -i ${prefix}.rmChrM_dedup.bam.bai -g $GenomeFasta -o ${prefix}.hmmratac -e $blacklist -m $readlength,200,400,600 --window 5000000
        fi
    fi
}

###############################################################
# Pipeline main body start here
###############################################################

trim_QC
Align
QC_removal
peak_calling

exit 0