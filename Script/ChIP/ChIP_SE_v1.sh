# ! /bin/bash
# 11 November 2022

##########################################################################
echo "#########################################################################"
echo "ChIP-seq analysis"
echo "#########################################################################"
##########################################################################

# Parameters that can added in later stage
# - blacklist region



##########################################################################

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-i <string>		Index used for bowtie2"	
	echo -e "-r <string>        Single-end read file directory"     
 	echo -e "-m <string>		Mode for analysis; there is <trim> and <no-trim> (default: no-trim)"   
	echo -e "-o <string>		Output directory with all result files (default: currect directory)"
	echo -e "-t <numeric>		Number of threads (default: 4)"
	exit 1
}

#Initiate parameters with NULL
mode="no-trim"
WDIR="."
threads=4

while getopts ":r:i:o:m:" op; do
	case $op in
		r) read=${OPTARG} ;;
		i) Index=${OPTARG} ;;
		o) WDIR=${OPTARG} ;;
		m) mode=${OPTARG} ;;
		t) threads=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))


# check necessary parameters
if [[ -z $read ]] || [[ -z $Index ]]; then
	echo -e "No index or read"
	usage
	exit -1
fi

#Get absolute file path, so users can use relative/absolute as they like.
[[ ${Index} != "" ]] && Index=`realpath ${Index}`
[[ ${read} != "" ]] && read=`realpath ${read}`
[[ ${WDIR} != "" ]] && WDIR=`realpath ${WDIR}`

# make directories if not exist and enter working directory.
[[ ! -d ${WDIR} ]] && mkdir -p ${WDIR}
cd ${WDIR}

# Making directory if trimming is selected
if [[ ${mode} == "trim" ]]; then
	trim_dir=${WDIR}/trim	# Creating a file to obtain the trimmed reads
    [[ ! -d ${trim_dir} ]] && mkdir -p ${trim_dir}
fi

##############################################################
# Pipeline function setup starts here
##############################################################

# Single end read pipeline functions

# Function for running end-to-end alignment on trimmed datasets
function Single_trim_align() {
    
	# 1. Running trimming
	# Creating a file to obtain the trimmed reads
    trim_dir=${WDIR}/trim
    [[ ! -d ${trim_dir} ]] && mkdir -p ${trim_dir}
	# Single-end trimming the reads
    trim_galore --fastqc --fastqc_args "--outdir $trim_dir" --output_dir $trim_dir $r
    # Changing the read to trimmed
    r=${trim_dir}/*_trimmed.f*q.gz

	# 2. Running end-to-end alignment
	bowtie2 --very-sensitive -p $threads -t -x $Index -U $r > ${prefix}.bowtie2.log | samtools sort -@ 4 -O bam -o ${prefix}.bam
}

# 2. Alignment
function Single_Align_local() {
	# Local alignment does not perform quality trim (Some case not even adapter trim)
    bowtie2 --very-sensitive-local -p 4 -t -x $Index -U $r > ${prefix}.bowtie2.log | samtools sort -@ 4 -O bam -o ${prefix}.bam
}

# 3. remove duplicates, Or also low quality and multimapper?
function quality_removal() {
    picard MarkDuplicates Input=${prefix}.bam Output=${prefix}.dup.bam METRICS_FILE=${prefix}.dup.txt
    samtools view -h -b -q 30 -@ 4 -F 3328 -o ${prefix}.final.bam ${prefix}.dup.bam
	samtools index ${prefix}.final.bam
}
# Seems like blacklist isnt essential

# 5. Peak calling
function peak_calling() {
	macs3 callpeak -t ${prefix}.final.bam -g genome -n NAME -B -g hs --outdir /PATH/DIR 
}

###############################################################
# Pipeline main body start here
###############################################################

for i in $(ls ${dir}/*.f*q.gz | sed 's/.f[^~]*q.gz//'); do 
	
	# if loop to run alignment
	# Running local alignment for untrimmed reads
	if [[ ${mode} == "no-trim" ]]; then
		
		bowtie2 --very-sensitive-local -p $threads -t -x $Index -U $r > $WDIR/$(basename $i).bowtie2.log | samtools sort -@ 4 -O bam -o $WDIR/$(basename $i).bam

	# Running trim and end-to-end alignment
	elif [[ ${mode} == "trim" ]]; then
		
		# 1. Running trimming
    	trim_galore --fastqc --fastqc_args "--outdir $trim_dir" --output_dir $trim_dir $r
	    # Changing the read to trimmed
    	r=${trim_dir}/$(basename $i)_trimmed.f*q.gz

		# 2. Running end-to-end alignment
		bowtie2 --very-sensitive -p $threads -t -x $Index -U $r > $WDIR/$(basename $i).bowtie2.log | samtools sort -@ 4 -O bam -o $WDIR/$(basename $i).bam

	else
		echo -e "Enter either <trim> or <no-trim> for the mode"
    	exit -1
	fi

	# After alignment, run marking duplication
	picard MarkDuplicates Input=$WDIR/$(basename $i).bam Output=$WDIR/$(basename $i).dup.bam METRICS_FILE=$WDIR/$(basename $i).dup.txt

	# Now it filtering time
	sambamba view -h -t 8 -f bam -o ./SRR111951.final.bam -F [XS] == null and not unmapped and not duplicate ./SRR111951.dup.bam
	
done
for i in $(ls ${dir}/*.f*q.gz | sed 's/.f[^~]*q.gz//'); do echo $i; done


################################################################
# Additional codes ready to be added in when need to
################################################################

# 1. Trimming the reads
# function trim_QC() {
#     # Trim the reads
#     trim_dir=${WDIR}/trim
#     [[ ! -d ${trim_dir} ]] && mkdir -p ${trim_dir}
#     trim_galore --paired --retain_unpaired --output_dir $trim_dir $R1 $R2
#     # Changing the reads to trimmed
#     Read1=${trim_dir}/*_val_1.f*
#     Read2=${trim_dir}/*_val_2.f*
#     # QC the trimmed reads
#     fastqc -j 2 --outdir $trim_dir $Read1 $Read2
# }

# 2. Alignment
# function Pair_Align_end() {
#     bowtie2 --very-sensitive -p 4 --no-mixed -X 500 -t -x $Index -1 $Read1 -2 $Read2 > ${prefix}.bowtie2.log | samtools sort -@ 4 -O bam -o ${prefix}.bam
# }
# function Pair_Align_local() {
#     bowtie2 --very-sensitive-local -p 4 --no-mixed -X 500 -t -x $Index -1 $Read1 -2 $Read2 > ${prefix}.bowtie2.log | samtools sort -@ 4 -O bam -o ${prefix}.bam
# }


