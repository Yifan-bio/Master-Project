#!/bin/bash
##########################################################################
#ATAC-seq pipeline

# Modified based on the code written by YenLab@github for Tn5 insertion

##########################################################################
set -ue
set -o pipefail
export LC_ALL=C

# help message
usage() {
	echo -e "This script is used to run the full atac-seq analysis by a single run."
    echo -e "However, thos os only a prototype."
	echo -e "The following are the available options"
	echo -e "-r <string>		If start from raw fastq file, please provide Read1 (_R1.fastq.gz) with this parameter."
	echo -e "-b <string>		If start from mapped bam file, please provide with this parameter and you can miss the -r/-m/-i parameter."
	echo -e "-T <string>		Indicate which trim package you wish to use (trimmomatic/fastp)[default:trimmomatic]"		
	echo -e "-m <string>		Indicate which mapper you wish to use (Bowtie2/BWA)[default:Bowtie2]"	
	echo -e "-i <string>		Location of mapping index corresponding to the mapper you provided (-m)."
	echo -e "-g <string>		Genome size for each chromosome."
	echo -e "-f <string>		Genome reference fasta file."
	echo -e "-t <string>		Tallymer mappability file directory."
	echo -e "-l <string>		[optional] Blacklist regions."
	echo -e "-p <string>		Thresholds for parallelly run this pipeline."
	echo -e "-o <string>		Working directory, all output files will be generated here."
	exit 1
}

# This script assumes all of the packages are in the PATH
# AT beginning these code will be NULL so allows a check for each case
BamFile=""
Blacklist=""
Index=""
Read_1=""
Read_2=""
Genome_fasta=""
threads=1

while getopts ":r1:r2:b:m:i:g:f:t:l:p:o:" op; do
	case $op in
		r1) Read_1=${OPTARG} ;;
		r2) Read_2=${OPTARG} ;;
        b) BamFile=${OPTARG} ;;
		i) Index=${OPTARG} ;;
		f) Genome_fasta=${OPTARG} ;;
		l) blacklist=${OPTARG} ;;
		p) threads=${OPTARG} ;;
		o) WDIR=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# verify all needed parameter
# if [[ -z $GenomeSize ]] || [[ -z $GenomeFasta ]] || [[ -z $Tallymer ]]; then
# 	echo -e "Necessary parameters missing (GenomeSize, GenomeFasta, Tallymer), please specify clearly..."
# 	usage
# 	exit -1
# fi

#Get absolute file path, so users can use relative/absolute as they like.
[[ ${Read_1} != "" ]] && Read_1=`realpath ${Read_1}`
[[ ${Read_2} != "" ]] && Read_2=`realpath ${Read_2}`
[[ ${BamFile} != "" ]] && BamFile=`realpath ${BamFile}`
[[ ${Index} != "" ]] && MapperIndex=`realpath ${Index}`
[[ ${GenomeSize} != "" ]] && GenomeSize=`realpath ${GenomeSize}`
[[ ${Genome_fasta} != "" ]] && Genome_fasta=`realpath ${Genome_fasta}`
[[ ${Blacklist} != "" ]] && Blacklist=`realpath ${Blacklist}`
[[ ${WDIR} != "" ]] && WDIR=`realpath ${WDIR}`

# make directories if not exist and enter working directory.
[[ ! -d ${WDIR} ]] && mkdir -p ${WDIR}
cd ${WDIR}

#=========================================================================================
# Pipeline main body
#=========================================================================================

# 1. Trim reads
trim_galore --paired --retain_unpaired --output_dir $WDIR/trim_galore $Read_1 $Read_2

if [[ -n ${FqFile} ]]; then
	# Get file prefixes
	RawFqPath=`dirname ${FqFile}`
	SamplePrefix=`basename ${FqFile} _R1.fastq.gz`
	echo -e "\n\n\nYou provided fastq file (${FqFile}), we will work from reads quality check ......"
	echo "------------------------------->>> Start processing ${SamplePrefix} <<<-------------------------------"
	
	echo -e "\nStep1. Checking quality of raw data ......"
	if [[ ! -f ${RawFqPath}/${SamplePrefix}_R2.fastq.gz ]];then
		echo "Sorry, BiasFreeATAC can't find fastq Read2 file (${RawFqPath}/${SamplePrefix}_R2.fastq.gz), please check!!!"
		exit -1
	fi
	fastqc --threads ${THREASHOLD} --quiet --outdir ./ ${RawFqPath}/${SamplePrefix}_R1.fastq.gz ${RawFqPath}/${SamplePrefix}_R2.fastq.gz 
	
	echo -e "\nStep2. Trimming Tn5 adapters ......"

	if [[ ${Trimmer} == "trimmomatic" ]];then
		trimmomatic PE -threads ${THREASHOLD} -phred33 \
		${RawFqPath}/${SamplePrefix}_R1.fastq.gz ${RawFqPath}/${SamplePrefix}_R2.fastq.gz \
		${SamplePrefix}_R1_trimmed.fastq ${SamplePrefix}_R1_unpaired.fastq ${SamplePrefix}_R2_trimmed.fastq ${SamplePrefix}_R2_unpaired.fastq \
		ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:1:true MINLEN:20 

		rm -f ${SamplePrefix}_R1_unpaired.fastq ${SamplePrefix}_R2_unpaired.fastq NexteraPE-PE.fa
		pigz -p ${THREASHOLD} ${SamplePrefix}_R1_trimmed.fastq ${SamplePrefix}_R2_trimmed.fastq

	elif [[ ${Trimmer} == "fastp" ]];then
		#fastp has the feature to remove ployG that frequently occur in Nova sequencer
		fastp --in1 ${RawFqPath}/${SamplePrefix}_R1.fastq.gz --in2 ${RawFqPath}/${SamplePrefix}_R2.fastq.gz \
		--out1 ${SamplePrefix}_R1_trimmed.fastq.gz --out2 ${SamplePrefix}_R2_trimmed.fastq.gz \
		--thread ${THREASHOLD} --trim_poly_g 10 --adapter_fasta NexteraPE-PE.fa --length_required 20 
		rm -f NexteraPE-PE.fa
	fi
	
	echo -e "\nStep3. Checking quality of trimmed data ......"
	fastqc --threads ${THREASHOLD} --quiet --outdir ./ ${SamplePrefix}_R1_trimmed.fastq.gz ${SamplePrefix}_R2_trimmed.fastq.gz
	
	echo -e "\nStep4. Aligning trimmed reads to reference genome and do initial filtering ......"
	# Remove read unmapped (4), mate unmapped (8), not primary alignment (256), read fails platform/vendor quality checks (512)
	# Remove low MAPQ reads (q < 30) and obtain sorted BAM file
	if [[ ${MAPPER} == "BWA" ]];then
		bwa mem ${MapperIndex} -t ${THREASHOLD} \
		<(zcat ${SamplePrefix}_R1_trimmed.fastq.gz) <(zcat ${SamplePrefix}_R2_trimmed.fastq.gz) \
		| samtools view -@ ${THREASHOLD} -F 780 -q 30 -b - \
		| samtools sort -@ ${THREASHOLD} - > ${SamplePrefix}.filtered.bam
	
	elif [[ ${MAPPER} == "Bowtie2" ]];then
		bowtie2 --threads ${THREASHOLD} --end-to-end --no-mixed --maxins 2000 --met 1 --met-file ${SamplePrefix}.alignMetrics.txt -q \
		-x ${MapperIndex} -1 ${SamplePrefix}_R1_trimmed.fastq.gz -2 ${SamplePrefix}_R2_trimmed.fastq.gz \
		| samtools view -@ ${THREASHOLD} -F 780 -q 30 -b - \
		| samtools sort -@ ${THREASHOLD} - > ${SamplePrefix}.filtered.bam
	fi
	
	samtools flagstat -@ ${THREASHOLD} ${SamplePrefix}.filtered.bam > ${SamplePrefix}.filtered.flagstat.qc
	
	echo -e "\nStep5. Getting Insert Size distribution ......"
	picard CollectInsertSizeMetrics --METRIC_ACCUMULATION_LEVEL ALL_READS --INPUT ${SamplePrefix}.filtered.bam \
	--OUTPUT ${SamplePrefix}.filtered.insertSizes.txt --Histogram_FILE ${SamplePrefix}.filtered.insertSizes.pdf
	
	echo -e "\nStep6. Marking PCR duplicates ......" 
	picard MarkDuplicates -INPUT ${SamplePrefix}.filtered.bam --OUTPUT ${SamplePrefix}.filtered.dupmark.bam \
	--METRICS_FILE ${SamplePrefix}.filtered.deduplicate.qc --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --REMOVE_DUPLICATES false
	samtools index -@ ${THREASHOLD} ${SamplePrefix}.filtered.dupmark.bam
	
	echo -e "\nStep7. Removing duplicates, Index final position sorted BAM file ......"
	# remove read is PCR or optical duplicate (1024)
	samtools view -@ ${THREASHOLD} -F 1024 -b ${SamplePrefix}.filtered.dupmark.bam > ${SamplePrefix}.filtered.dedup.bam
	samtools index -@ ${THREASHOLD} ${SamplePrefix}.filtered.dedup.bam
	samtools flagstat -@ ${THREASHOLD} ${SamplePrefix}.filtered.dedup.bam > ${SamplePrefix}.filtered.dedup.flagstat.qc
	
	rm -f ${SamplePrefix}.filtered.bam* ${SamplePrefix}.filtered.dupmark.bam*
	
	echo -e "\nStep8. Shift the reads in bam file to get Tn5 insertion sites..."
	#Adjust 4 or 5bp for the Tn5 binding site and get cut sites bed/bam file
	bedtools bamtobed -i ${SamplePrefix}.filtered.dedup.bam \
	| awk -v FS="\t" -v OFS="\t" '{ if ($6 == "+") {print $1,$2+4,$2+4+1,$4,$5,$6} else if ($6 == "-") {print $1,$3-5-1,$3-5,$4,$5,$6} }' \
	> ${SamplePrefix}.InsertSites.bed

	bedtools bedtobam -i ${SamplePrefix}.InsertSites.bed -g ${GenomeSize} | samtools sort -@ ${THREASHOLD} > ${SamplePrefix}.InsertSites.bam
	samtools index -@ ${THREASHOLD} ${SamplePrefix}.InsertSites.bam

	# Move preprocessing files
	[ ! -d ./metrics ] && mkdir -p ./metrics
	mv *zip *html *insertSizes.txt *insertSizes.pdf *alignMetrics.txt *.flagstat.qc *deduplicate.qc ./metrics
fi

if [[ -n ${BamFile} ]]; then
	BAM=${BamFile}
	echo -e "\n\n\nYou provided bam file (${BamFile}), we will direct work from bias correction ......"
	SamplePrefix=`basename ${BamFile} .bam`
	echo "------------------------------->>> Start processing ${SamplePrefix} <<<-------------------------------"
else
	BAM=${SamplePrefix}.filtered.dedup.bam
fi

echo -e "\nStep9. Correcting Tn5 bias using nucleotide dependency information ......"
#The mask scheme is slightly adapted from (Martins et al., 2017). After correction, we shifted the read position 
#to count Tn5 insertion center following (Buenrostro et al., 2013)
#The read length was used to calculate genome-wide mappability, 36 generally perform well (Derrien et al., 2012).
#kmer-size is counted from 0, so the 19mer mask give 18

kmer_mask=XNXXXCXXNNXNNNXXNNX
ReadLen=36
faspre=`basename ${GenomeFasta} .fa`

if [[ ! -f ${Tallymer}/${faspre}.tal_36.gtTxt.gz ]];then
	#If the tallymer directory (-t) doesn't contian tallymer files, BiasFreeATAC use long time to create these files (5 hours for mouse!).
	echo "BiasFreeATAC can't find tallymer file (${Tallymer}/${faspre}.tal_36.gtTxt.gz), but will build one de novo!!!"
	#Returned corrected ATAC-seq signals
	seqOutBias ${GenomeFasta} ${BAM} --skip-bed --read-size=${ReadLen} --strand-specific --custom-shift=4,-5 \
	--kmer-size=18 --plus-offset=5 --minus-offset=5 --kmer-mask ${kmer_mask} \
	--bw=${SamplePrefix}_corrected.bigWig 
	
	#Returned uncorrected ATAC-seq signals
	seqOutBias ${GenomeFasta} ${BAM} --skip-bed --read-size=${ReadLen} --strand-specific --custom-shift=4,-5 \
	--kmer-size=18 --plus-offset=5 --minus-offset=5 --no-scale \
	--bw=${SamplePrefix}_uncorrected.bigWig

	rm -f ${faspre}.sft.* ${faspre}.tal_.*

	#Move tallymer files is tallymer directory specified
	if [[ ${Tallymer} != "./" ]];then
		mv ${faspre}.tal_36.gtTxt.gz ${faspre}_36.18.5.5.tbl ${Tallymer}
	fi

elif [[ -f ${Tallymer}/${faspre}.tal_36.gtTxt.gz ]];then
	#If you provide tallymer directory (-t), BiasFreeATAC create soft link for Tallymer and SeqTable files in working directory.
	echo "BiasFreeATAC find tallymer files (${Tallymer}/${faspre}.tal_36.gtTxt.gz), will use this for bias correction..."

	if [[ ${Tallymer} != "./" ]];then
		ln -sf ${Tallymer}/${faspre}.tal_36.gtTxt.gz .
		ln -sf ${Tallymer}/${faspre}_36.18.5.5.tbl .
	fi

	#Returned corrected ATAC-seq signals
	seqOutBias ${GenomeFasta} ${BAM} --skip-bed --read-size=${ReadLen} --strand-specific --custom-shift=4,-5 \
	--kmer-size=18 --plus-offset=5 --minus-offset=5 --kmer-mask ${kmer_mask} \
	--bw=${SamplePrefix}_corrected.bigWig --tallymer=${faspre}.tal_36.gtTxt.gz
	
	#Returned uncorrected ATAC-seq signals
	seqOutBias ${GenomeFasta} ${BAM} --skip-bed --read-size=${ReadLen} --strand-specific --custom-shift=4,-5 \
	--kmer-size=18 --plus-offset=5 --minus-offset=5 --no-scale \
	--bw=${SamplePrefix}_uncorrected.bigWig --tallymer=${faspre}.tal_36.gtTxt.gz

	rm -f ${faspre}.tal_36.gtTxt.gz ${faspre}_36.18.5.5.tbl
fi

bigWigToBedGraph ${SamplePrefix}_uncorrected.bigWig /dev/stdout \
| grep "^chr" | sort --parallel=${THREASHOLD} -k1,1 -k2,2n > ${SamplePrefix}_uncorrected.bedGraph

bigWigToBedGraph ${SamplePrefix}_corrected.bigWig /dev/stdout \
| grep "^chr" | sort --parallel=${THREASHOLD} -k1,1 -k2,2n > ${SamplePrefix}_corrected.bedGraph

#If blacklist is provided, BiasFreeATAC will exclude these regions for downstream analysis
if [[ -n ${blacklist} ]]; then
	bedtools intersect -a ${SamplePrefix}_uncorrected.bedGraph -b ${blacklist} -v -wa > ${SamplePrefix}_uncorrected.bedGraph_tmp && \
	mv ${SamplePrefix}_uncorrected.bedGraph_tmp ${SamplePrefix}_uncorrected.bedGraph

	bedtools intersect -a ${SamplePrefix}_corrected.bedGraph -b ${blacklist} -v -wa > ${SamplePrefix}_corrected.bedGraph_tmp && \
	mv ${SamplePrefix}_corrected.bedGraph_tmp ${SamplePrefix}_corrected.bedGraph
fi

echo -e "\nStep10. Peakcalling using uncorrected and corrected signals ......"
GENOME_SIZE_NUM=`awk '{SUM+=$2}END{print SUM}' ${GenomeSize}`

macs2 callpeak --broad --treatment ${SamplePrefix}_uncorrected.bedGraph --format BED --gsize ${GENOME_SIZE_NUM} --qvalue 0.01 --broad-cutoff 0.01 \
--outdir ${SamplePrefix}_uncorrected_peaks --name ${SamplePrefix}_uncorrected --bdg --nomodel --max-gap 100 --shift -100 --extsize 200 

macs2 callpeak --broad --treatment ${SamplePrefix}_corrected.bedGraph --format BED --gsize ${GENOME_SIZE_NUM} --qvalue 0.01 --broad-cutoff 0.01 \
--outdir ${SamplePrefix}_corrected_peaks --name ${SamplePrefix}_corrected --bdg --nomodel --max-gap 100 --shift -100 --extsize 200 

uncorrectedBroadPeak=${SamplePrefix}_uncorrected_peaks/${SamplePrefix}_uncorrected_peaks.broadPeak
correctedBroadPeak=${SamplePrefix}_corrected_peaks/${SamplePrefix}_corrected_peaks.broadPeak

bedtools intersect -a ${uncorrectedBroadPeak} -b ${correctedBroadPeak} -wa -v > ${SamplePrefix}_uncorrected_specific_peaks.bed
bedtools intersect -a ${uncorrectedBroadPeak} -b ${correctedBroadPeak} -wa    > ${SamplePrefix}_uncorrected_shared_peaks.bed
bedtools intersect -b ${uncorrectedBroadPeak} -a ${correctedBroadPeak} -v -wa > ${SamplePrefix}_corrected_specific_peaks.bed
bedtools intersect -b ${uncorrectedBroadPeak} -a ${correctedBroadPeak} -wa    > ${SamplePrefix}_corrected_share_peaks.bed

echo "------------------------------->>> Done processing ${SamplePrefix} <<<-------------------------------"