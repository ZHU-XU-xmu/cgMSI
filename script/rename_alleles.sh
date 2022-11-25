#!/bin/bash
path="/data2/zhuxu_data/alleles/Enterococcus_faecalis_alleles"
filenames=$(ls $path)
for file in $filenames
do
    filename=$(echo $file |cut -f1 -d "." )
   cat $path/${file}|seqkit replace -p ^ -r "${filename//_/}_"  >>../alleles/Ef_alleles.fasta

done

