
column=4
WDIR="."

file_list=""
x=`ls /mnt/e/salmon/Input/*/quant.sf | cut --delimiter ' ' --fields 1 | readpath`
echo $x
awk 'NR!=1 {print $1}' $x > $WDIR/suppa.txt

for file in /mnt/e/salmon/Input/*/quant.sf
do
    echo $file
    # Getting the folder name
    file_name=`dirname $file`
    file_name=`basename $file_name | cut -f1`
    echo $file_name
    awk -F "\t" -v HEADER=$file_name 'BEGIN{print HEADER}; NR!=1 {print $4}' $file > $WDIR/$file_name.txt
    # add absolute path of file_name.txt to the file_list variable

    file_list="$file_list $WDIR/$file_name.txt"
done