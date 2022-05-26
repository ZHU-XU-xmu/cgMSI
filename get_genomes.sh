#!/bin/bash

# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
filename="txid_1639_ftp.txt"
awk -F "\t" '$7==1639&&$12=="Complete Genome"{print $20}' assembly_summary_refseq.txt > $filename

mkdir ./genomes/1639
echo "$filename"
cat $filename | while read line
do
	echo $line
	assembly_name=${line##*/}
	echo $assembly_name
  #cp -v "/data1/zhuxu_data/snipe_nanopore/snipe_genome/${assembly_name}_genomic.fna" "/data1/zhuxu_data/MLST_nanopore/genomes/${assembly_name}_genomic.fna"
	#wget -P ./snipe_genome $line/${assembly_name}_genomic.fna.gz -c
   axel -q -n 15 $line/${assembly_name}_genomic.fna.gz -o /data2/zhuxu_data/Nano_cgMLST/genomes/1639/

done



