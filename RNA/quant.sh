#! /bin/bash

# Find file contains DRR in the file name and output file that contains R1 in the file name
# function find_R1() {
#     find . -name "*DRR*R1*"
# }

i="SRR3214318"


# find files containing SRR1338 within the file name in current directory and it subdirectories
function find_srr() {
    echo "find_srr"
    for file in `find . -name "*$i*"`
    do
        echo "fastqc $file"
        #fastqc $file
    done
}
find_srr

# find file contains *SRR*R1 and make the directory of the file as $i
function find_paired() {
    echo "find_paired"
    for file in `find . -name "*$i*"`
    do
        echo "fastqc $file"
        #fastqc $file
    done
}