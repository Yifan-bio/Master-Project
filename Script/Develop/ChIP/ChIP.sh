# ! /bin/bash
# 15 September 2023

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
    	trim_galore --output_dir $trim_dir $r
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
	samtools view -h -b -q 20 -@ 4 -F 3328 -o ${prefix}.final.bam ${prefix}.dup.bam
	samtools index ${prefix}.final.bam

    # Peak calling
    macs3 callpeak -t ${prefix}.final.bam -g genome -n NAME -B -g hs --outdir /PATH/DIR

done

exit 0
