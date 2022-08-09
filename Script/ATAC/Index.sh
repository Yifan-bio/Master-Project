# ! /bin/bash
# 4 August 2022
#!/bin/bash
##########################################################################
echo "#########################################################################"
echo "Indexing genome"
echo "#########################################################################"
##########################################################################

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""	
	echo -e "-g <string>		Genome files used to make bowtie2 index file"	
	exit 1
}

#Initiate parameters with NULL
GenomeFasta=""

while getopts ":g:b:" op; do
	case $op in
		g) GenomeFasta=${OPTARG} ;;
		b) blacklist=${OPTARG} ;;
        o) WDIR=${OPTARG} ;;
        \?) usage ;;
	esac
done
shift $((OPTIND-1))

# check necessary parameters
if [[ -z $g ]] || [[ -z $b ]]; then
	echo -e "No index, no reads, no analysis"
	usage
	exit -1
fi

#Get absolute file path, so users can use relative/absolute as they like.
[[ ${GenomeFasta} != "" ]] && GenomeFasta=`realpath ${GenomeFasta}`
[[ ${blacklist} != "" ]] && blacklist=`realpath ${blacklist}`
[[ ${WDIR} != "" ]] && WDIR=`realpath ${WDIR}`

# make directories if not exist and enter working directory.
[[ ! -d ${WDIR} ]] && mkdir -p ${WDIR}
cd ${WDIR}

##########################################################################
echo "#########################################################################"
echo "Indexing genome"
echo "#########################################################################"
##########################################################################

#Indexing genome
bowtie2 index -p 8 -x ${WDIR}/genome -f ${GenomeFasta}

# Getting the genome size
genome_size=`awk '{sum+=$3-$2}END{print sum}' ${GenomeFasta}`
