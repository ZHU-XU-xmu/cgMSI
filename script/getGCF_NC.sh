#!/bin/bash

path="/data2/zhuxu_data/Nano_cgMLST/genomes/287"
cd "/data2/zhuxu_data/Nano_cgMLST/genomes/287"
files=$(ls *.fna )
cd /data2/zhuxu_data/Nano_cgMLST/
for file in $files
do 
    echo $file
    seqkit seq -n $path/$file >temp2.txt
    file=$(echo ${file} | cut -d '_' -f 1,2)
    file=$(echo ${file} | cut -d '.' -f 1)
    
    echo -e "${file}\t\c" >> /data2/zhuxu_data/Nano_cgMLST/library/Pa_NC_GCF.tsv
    cat temp2.txt | while read line
    do 
      NC_name=$(echo ${line} | cut -d ' ' -f 1 )
      echo ${NC_name}
      echo -e ${NC_name}"\t\c" >> /data2/zhuxu_data/Nano_cgMLST/library/Pa_NC_GCF.tsv
    
    done 
    echo "" >> /data2/zhuxu_data/Nano_cgMLST/library/Pa_NC_GCF.tsv

done
rm temp2.txt