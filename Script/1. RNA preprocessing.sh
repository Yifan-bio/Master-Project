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
    echo -e ""
    echo -e "Optional parameters:"
    echo -e "-e <string>        Extension type (default: fq.gz)"
    echo -e "-t <string>		Number of threads to be used (default: 4)"
	echo -e "-o <string>		Output directory with all salmon results (default: current directory)"
    echo -e "-s <string>		Directory to salmon binary file if it is not in the PATH"
    echo -e 'default command: $salmon quant -i index -l A -1 R1 -2 R2 -p threads --validateMappings --gcBias --seqBias --recoverOrphans -o output'
	exit 1
}

# Adjustable default paramters 
DIR="."
salmon="salmon"
extra_args=""
threads=4
extension="fq.gz"

# Adding in the paramters from command line
while getopts ":i:f:e:o:s:t:" op; do
	case $op in
		i) index=${OPTARG} ;;
        f) dir=${OPTARG} ;;
        e) extension=${OPTARG} ;;
		t) threads=${OPTARG} ;;
        o) WDIR=${OPTARG} ;;
        s) salmon=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# Verify if the essential directory are added (index and read files)
if [[ -z $index ]] || [[ -z $dir ]]; then
	echo -e "Either the index or file dir is wrong"
	usage
	exit -1
fi

# Get absolute file path, so users can use relative/absolute as they like.
[[ ${index} != "" ]] && index=`realpath ${index}`
[[ ${dir} != "" ]] && dir=`realpath ${dir}`
[[ ${WDIR} != "" ]] && WDIR=`realpath ${WDIR}`

# If the output directory is specified and does not exist then create this folder.
[[ ! -d ${WDIR} ]] && mkdir -p ${WDIR}

##########################################################################
#                               Main
##########################################################################

# There is two type of fastq file extension which is fq or fastq, two if loops accounts both.

if [[ ${extension} == "fq.gz" ]]; then
    
    # This is for paired end file, this will extract all fq file names (following ENA format)
    # and get the pairs. Then it will input one sample pair each time.
    for i in $(ls ${dir}/*.fq*.gz | sed 's/[1-2].fq.gz//' | uniq); do 
    
        # Testing if the two both Read 1 and Read 2 exist in ENA format
        if [[ -z ${dir}/${i}1.fq.gz ]] || [[ -z ${dir}/${i}2.fq.gz ]]; then
            echo -e "Read 1 and Read 2 file cannot be detected for $i"
            usage
            exit -1
        fi
    
        # Setting the output folder for salmon on this sample
        o="$WDIR/$(basename $i)"

        # Running salmon alignment. We have included several paramters to correct bias.
        $salmon quant -i $index -l A -1 ${i}1.fq.gz -2 ${i}2.fq.gz -p $threads --validateMappings --gcBias --seqBias --recoverOrphans -o $o
    done

elif [[ ${extension} == "fastq.gz" ]]; then
    for i in $(ls ${dir}/*.fastq*.gz | sed 's/[1-2].fastq.gz//' | uniq); do 
    
        # Testing if the two files exist
        if [[ -z ${dir}/${i}1.fastq.gz ]] || [[ -z ${dir}/${i}2.fastq.gz ]]; then
            echo -e "Read 1 and Read 2 file cannot be detected for $i"
            usage
            exit -1
        fi
    
        # Output file
        o="$WDIR/$(basename $i)"

        # Reporting info
        $salmon quant -i $index -l A -1 ${i}1.fastq.gz -2 ${i}2.fastq.gz -p $threads --validateMappings --gcBias --seqBias --recoverOrphans -o $o
    done

else
    echo -e "Only accept file extension for fastq.gz and fq.gz"
  	exit -1
fi

##########################################################################
#                               End
##########################################################################