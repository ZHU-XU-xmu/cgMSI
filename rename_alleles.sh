#!/bin/bash
path="/data2/zhuxu_data/Nano_cgMLST/alleles/Staphylococcus_aureus_cgMLST_alleles"
filenames=$(ls $path)
for file in $filenames
do
    #echo "$file"
    filename=$(echo $file |cut -f1 -d "." )
    #echo "$filename"
    #cat $path/${file}|seqkit replace -p .+ -r "${filename}_{nr}"  >>./Kl_alleles.fasta
   # if [ ${filename} = "KP1_RS00015" ];then
   cat $path/${file}|seqkit replace -p ^ -r "${filename//_/}_"  >>./alleles/Sa_alleles.fasta
    
   # fi
done

