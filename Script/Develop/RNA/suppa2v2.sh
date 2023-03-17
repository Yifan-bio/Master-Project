# ! /bin/bash
# 24 September 2022
##########################################################################
echo "#########################################################################"
echo "SUPPA2 analysis"
echo "#########################################################################"
##########################################################################

# Prerequisite files:
# - salmon output
# - gtf file used for salmon index

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-g <string>		gtf file used for salmon index"		
    echo -e "-f <string>		Salmon output directory (The files with all salmon result as folder)"
	echo -e "-o <string>		Output directory for SUPPA2 (default: current directory)"
	echo -e ""
	echo -e "The default uses salmon output which TPM is recorded in the 4th column of quant.sh"
	echo -e ""
	echo -e ""	
    echo -e ""
	echo -e ""
	exit 1
}

##########################################################################
#                          Adding paramters
##########################################################################

#Initiate parameters with NULL
while getopts ":g:f:o:" op; do
	case $op in
		g) gtf=${OPTARG} ;;
        f) input_file=${OPTARG} ;;
        o) workdir=${OPTARG} ;;
        #e) extra_args=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# Checking the input paramter is provided
if [[ -z $gtf ]] || [[ -z $input_file ]] || [[ -z $workdir ]]; then
	echo -e "Guess which 3 is required"
	usage
	exit -1
fi

[[ ${gtf} != "" ]] && gtf=`realpath ${gtf}`
[[ ${input_file} != "" ]] && input_file=`realpath ${input_file}`
[[ ${workdir} != "" ]] && workdir=`realpath ${workdir}`

# Variables used to ease setting up
num=0 # Helps the isolate the first column
list="" # used to record the files need to be combined
name_list=""
##########################################################################
#                          Primary function
##########################################################################

function SUPPA2_preparation() {
    
    # Creating the input file for SUPPA2
    for i in $input_file/*/quant.sf;
        do
        # Getting the transcript name from salmon output as row names for SUPPA2
        if [[ $num == 0 ]]; then
            cat $i | awk -v file_name="" 'BEGIN{print skip} {if (NR!=1) {print $1}}' > $workdir/TXname.txt
            num=$(($num + 1))
        fi
    
        # Getting the name of the file to use as sample column name
        name=`realpath ${i}`
        name=$(basename $(dirname $name))
    
        # Extracting the TPM result from each replicate using the file name as column name
        cat $i | awk -v file_name="$name" 'BEGIN{print file_name} {if (NR!=1) {print $4}}' > $workdir/$name.txt
    
        # Keeping record of the location so can be easily combined
        file_dir=`realpath $workdir/$name.txt`
        list="$list $file_dir"
        name_list="$name_list $name"
    done

    # Combining files and remove intermediate files
    paste $workdir/TXname.txt $list > $workdir/suppa_input.txt
    rm $workdir/TXname.txt $list
}

##########################################################################
#                          Running the function
##########################################################################

SUPPA2_preparation
input="$workdir/suppa_input.txt"

# Running SUPPA to calucalte the AS events on input gtf file
suppa.py generateEvents -i $gtf -o $workdir/AS_events -e SE SS MX RI FL -f ioe
# Putting all into one file
awk 'FNR==1 && NR!=1 { while (/^seqname/) getline; } 1 {print}' $workdir/*ioe > $workdir/All_events.ioe

# Run SUPPA for getting the psi values of the events:
suppa.py psiPerEvent -i $workdir/All_events.ioe -e $input -o $workdir/psi_events

# Run SUPPA for obtaining the Differential splicing analysis
# Split the PSI and TPM files between the 2 conditions
~/scripts/split_file.R ~/tra2/Salmon/quantification/iso_tpm.txt SRR1513329,SRR1513330,SRR1513331 SRR1513332,SRR1513333,SRR1513334 ~/tra2/Salmon/quantification/CTRL.tpm ~/tra2/Salmon/quantification/KD.tpm 
	~/scripts/split_file.R ~/tra2/SUPPA/events.psi SRR1513329,SRR1513330,SRR1513331 SRR1513332,SRR1513333,SRR1513334 ~/tra2/SUPPA/CTRL.psi ~/tra2/SUPPA/KD.psi
	#5.2: Run SUPPA
	suppa.py diffSplice -m empirical -i ./gencodev40_events_all_events.ioe -e ./CTRL.tpm ./Late.tpm -pa -gc -p ./CTRL.psi ./Later.psi -o ./CTRL_Later
