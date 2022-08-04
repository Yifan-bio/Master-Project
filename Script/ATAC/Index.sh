# ! /bin/bash
# 4 August 2022
#!/bin/bash
##########################################################################
echo "#########################################################################"
echo "Indexing genome"
echo "#########################################################################"
##########################################################################
set -ue
set -o pipefail
export LC_ALL=C

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""	
	echo -e "-g <string>		Genome files used to make bowtie2 index file"	
	echo -e "-b <string>		Blacklist file from ENCODE project (boyle-lab)"
	echo -e "-o <string>		Output directory with all result files"
	exit 1
}

#Initiate parameters with NULL
GenomeFasta=""

while getopts ":r1:r2:i:g:b:o:p:" op; do
	case $op in
		r1) R1=${OPTARG} ;;
		r2) R2=${OPTARG} ;;
		i) Index=${OPTARG} ;;
		g) GenomeFasta=${OPTARG} ;;
		b) blacklist=${OPTARG} ;;
		o) WDIR=${OPTARG} ;;
        p) prefix=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

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