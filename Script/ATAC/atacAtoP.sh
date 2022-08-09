# ! /bin/bash
# 2 August 2022

##########################################################################
echo "#########################################################################"
echo "ATAC-seq analysis"
echo "#########################################################################"
##########################################################################
set -ue
set -o pipefail
export LC_ALL=C

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-r1 <string>		Read1 file"
	echo -e "-r2 <string>		Read2 file"
	echo -e "-i <string>		Index used for bowtie2"		
	echo -e "-g <string>		Genome files used to make bowtie2 index file"	
	echo -e "-p <string>		Prefix for all the output files"
	echo -e "-b <string>		Blacklist file from ENCODE project (boyle-lab)"
	echo -e "-o <string>		Output directory with all result files"
	echo -e "-m <string>		Mode for analysis; there is strict and lenient (default: strict). Strict will remove chrM,duplicate,low-quality,multimap and improper mapped. While lenient will only run remove chrM and deduplicate"
	exit 1
}

while getopts ":r1:r2:i:g:b:o:p:m:" op; do
	case $op in
		r1) R1=${OPTARG} ;;
		r2) R2=${OPTARG} ;;
		i) Index=${OPTARG} ;;
		g) GenomeFasta=${OPTARG} ;;
		b) blacklist=${OPTARG} ;;
		o) WDIR=${OPTARG} ;;
        p) prefix=${OPTARG} ;;
		m) mode=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

#Initiate parameters with NULL
blacklist=""
GenomeFasta=""
prefix="test"
mode="strict"
WDIR="."

# check necessary parameters
if [[ -z $R1 ]] || [[ -z $R2 ]] || [[ -z $Index ]]; then
	echo -e "No index, no reads, no analysis"
	usage
	exit -1
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
    trim_galore --paired --retain_unpaired --output_dir $trim_dir $R1 $R2
    # Changing the reads to trimmed
    Read1=${trim_dir}/*_val_1.f*
    Read2=${trim_dir}/*_val_2.f*
    # QC the trimmed reads
    fastqc -j 2 --outdir $trim_dir $Read1 $Read2
}

# 2. Alignment
function Align() {
    bowtie2 --very-sensitive --end-to-end -p 4 --dovetail --no-mixed -X 2000 -t -x $Index -1 $Read1 -2 $Read2 > ${prefix}.bowtie2.log | samtools sort -@ 4 -O bam -o ${prefix}.bam
}

# 3. remove mitochondrial reads and duplicates
function essential_removal() {
    samtools view -@ 4 -h ${prefix}.bam | grep -v chrM | samtools sort -@ 4 -O bam -o ${prefix}.rmChrM.bam
    picard MarkDuplicates Input=${prefix}.rmChrM.bam Output=${prefix}.rmChrM_dedup.bam METRICS_FILE=${prefix}.dedup.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
	samtools index ${prefix}.rmChrM_dedup.bam
}

# 4. remove low quality, unmapped, multimapped and impropoerly paired reads
function optimal_removal() {
	samtools view -@ 4 -h ${prefix}.bam | grep -v chrM | samtools sort -@ 4 -O bam -o ${prefix}.rmChrM.bam
    picard MarkDuplicates Input=${prefix}.rmChrM.bam Output=${prefix}.rmChrM_dedup.bam METRICS_FILE=${prefix}.dedup.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
    samtools view -h -b -q 30 -@ 4 -F 1804 -f 2 -o ${prefix}.final.bam ${prefix}.rmChrM_dedup.bam
	samtools index ${prefix}.final.bam
}

# 5. Peak calling
function peak_calling() {
	HMMRATAC -b ${prefix}.final.bam -i ${prefix}.final.bam.bai -g $GenomeFasta -o ${prefix}.hmmratac -e $blacklist -m 75,200,400,600 --window 5000000
}

###############################################################
# Pipeline main body start here
###############################################################

trim_QC
Align
if [[ $mode == "strict" ]]; then
	essential_removal
	peak_calling 
elif [[ $mode == "lenient" ]]; then
	optimal_removal
fi
peak_calling