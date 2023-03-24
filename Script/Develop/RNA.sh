# ! /bin/bash
# 9 August 2022
##########################################################################
echo "#########################################################################"
echo "RNA sequencing pseudoalignment analysis"
echo "#########################################################################"
##########################################################################

# Prerequisite:
# - salmon

##########################################################################
#                          Adding paramters
##########################################################################

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-i <string>		Salmon Index"		
    echo -e "-f <string>		Directory with library files"
    echo -e "-t <string>		Number of threads to be used (default: 4)"
	echo -e "-o <string>		Output directory with all salmon results (default: current directory)"
    echo -e "-s <string>		Directory to salmon binary file if it is not in the PATH"
    echo -e "-v <string>        Package version of salmon (defaults: 1.10.0)"
	exit 1
}

#Initiate parameters with NULL
WDIR="."
salmon="salmon"
threads=4
version="1.10.0"

while getopts ":i:f:t:o:s:v:" op; do
	case $op in
		i) index=${OPTARG} ;;
        f) dir=${OPTARG} ;;
		t) threads=${OPTARG} ;;
        o) WDIR=${OPTARG} ;;
        s) salmon=${OPTARG} ;;
        v) version=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# check necessary parameters
if [[ -z $index ]]; then
	echo -e "The index file for salmon is missing"
	usage
	exit -1
fi

if [[ -z $dir ]]; then
	echo -e "The input directory for salmon is missing"
	usage
	exit -1
fi

# Get absolute file path, so users can use relative/absolute as they like.
[[ ${index} != "" ]] && index=`realpath ${index}`
[[ ${dir} != "" ]] && dir=`realpath ${dir}`
[[ ${WDIR} != "" ]] && WDIR=`realpath ${WDIR}`

# make directories if not exist and enter working directory.
[[ ! -d ${WDIR} ]] && mkdir -p ${WDIR}
cd ${WDIR}

# A function to verify the package version of salmon is above the requirements
function verify_package_version {
    package_version=$(salmon --version | cut -d " " -f 2)
    if [[ "$package_version" > $version ]] || [[ "$package_version" == "$version" ]]; then
        echo "Salmon package version is $version, fit the requirement"
        # Execute the rest of the script here
    else
        echo "Salmon package version is less than $version"
        exit 1
    fi
}

##########################################################################
#                               Main
##########################################################################

verify_package_version

for i in $(ls ${dir}/*.fastq.gz | cut -d "_" -f 1-4 | uniq); do 
    
    # Testing if the two files exist
    if [[ -z ${i}_R1_*.fastq.gz ]] || [[ -z ${i}_R2_*.fastq.gz ]]; then
        echo -e "Read 1 and Read 2 file cannot be detected for $i"
        usage
        exit -1
    fi

    # Output file
    o="$WDIR/$(basename $i)"

    # Reporting info
    echo "$salmon quant -i $index -l A -1 ${i}_R1_001.fastq.gz -2 ${i}_R2_001.fastq.gz -p $threads --validateMappings --gcBias --seqBias --recoverOrphans -o $o"

done

##########################################################################
#                               End
##########################################################################