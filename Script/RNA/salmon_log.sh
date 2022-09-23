#! /bin/bash
# 23 Sep 2022
######################################################################################
#
# This script is designed to extract the log files from salmon to get all the mapping information in one run
#
######################################################################################

dir="."
workdir="."
num=0
list=""

for i in $dir/*/aux_info/meta_info.json;
    do
    
    # Getting the name of the file
    if [[ $num == 0 ]]; then
        cat $i | sed -n '40,46p' | awk 'BEGIN{print "meta_info"} {print $1}' | tr '"' ' ' > $workdir/n.txt
        num=$(($num + 1))
    fi
    
    # Getting the name of the file to make it the name for the file record
    name=`realpath ${i}`
    name=$(basename $(dirname $(dirname $name)))
    
    # Extracting the results
    cat $i | sed -n '40,46p' | awk -v file_name="$name" 'BEGIN{print file_name}{print $2}' | tr ',' ' ' > $workdir/$name.txt
    
    # Keeping record of the location so can be easily combined
    file_dir=`realpath $workdir/$name.txt`
    list="$list $file_dir"

done

paste $workdir/n.txt $list > salmon_meta_info.csv
rm $workdir/n.txt $list

