
dir="."
workdir="."

# Variables to verify some stuff
num=0 # Helps the isolate the first column
list="" # used to record the files need to be combined

for i in $dir/*/quant.sf;
    do
    
    # Getting the name of the file
    if [[ $num == 0 ]]; then
        cat $i | awk '{print $1}' > $workdir/TXname.txt
        num=$(($num + 1))
    fi
    
    # Getting the name of the file to make it the name for the file record
    name=`realpath ${i}`
    name=$(basename $(dirname $name))
    
    # Extracting the results
    cat $i | awk -v file_name="$name" 'BEGIN{print file_name} {if (NR!=1) {print $4}}' > $workdir/$name.txt
    
    # Keeping record of the location so can be easily combined
    file_dir=`realpath $workdir/$name.txt`
    list="$list $file_dir"

done

# Combine files and and then remove temp files
paste $workdir/TXname.txt $list > suppatest.txt
rm $workdir/TXname.txt $list

