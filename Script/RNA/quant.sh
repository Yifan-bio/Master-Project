# ! /bin/bash

o="./24hr"
index="SRR8932929"
# input="./run.txt"

# # find files containing $line (an argument) within the file name in current directory and it subdirectories

# # function separating files based on if there is R1 in file name or R2 in file name

# function separate_reads() {
#     # if file contain R1, use the file as variable $R1
#     if [[ $file == *"R1"* ]]; then
#         R1=$file
#         echo "Read 1 file: $file"
#         continue
#     fi
#     # if the file contain R2,then echo the file name
#     if [[ $file =~ "R2" ]]; then
#         R2=$file
#         echo "Read 2 file: $file"
#     fi 
# }


# # find files containing $i within the file name in current directory and it subdirectories
# function find_srr() {
#     for file in `find . -name "*.fastq.gz"`;
#     do
#         echo $file
#         if [[ $file =~ $1 ]]; then
#             separate_reads $file
#         fi
#     done
# }

# # 

# while read line; do
#     echo "Now running paired-end salmon on $line"
#     # convert $line to SRR number
#     srr=${line##*}
#     echo "das $srr"
#     find_srr $srr
#     echo "salmon -i $index --seqBias --gcBias -validateMappings --recoverOrphans -l A -p 8 -1 $R1 -2 $R2 -o $o"
# done < $input



# function all_fq() {
#     find . -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fq.gz"
# }

# for file in $(all_fq);
# do
# echo "salmon -i $index --seqBias --gcBias -validateMappings --recoverOrphans -l A -p 8 -1 $ -2 $R2 -o $o"
# done

for i in $(ls ./*.fastq.gz | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sort -u)