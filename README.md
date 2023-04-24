# cgMSI
## cgMSI: pathogen detection within species from nanopore metagenomic sequencing data
### Introduction
cgMSI (Core Genome Metagenome Strain Identify), a tool to detect pathogen from nanopore metagenomic data within species at low abundance. 
 cgMSI consists of two core modules:
- The cgMSI LIB module will create or update the library accroding to the file provided by the user. 
- The cgMSI MAP module will  identify the strain and estimate the abundance.


### Support and Contact
For any issues or concerns, please contact us at zhuxu@stu.xmu.edu.cn
 

### Pathogenic Species Supported
Species name | Number of loci
-|-
Klebsiella pneumoniae | 2358
Escherichia coli| 2531
Enterococcus faecalis | 1972
Listeria monocytogenes | 1701
Pseudomonas aeruginosa | 3867
Staphylococcus aureus| 1861
Salmonella enterica | 3002

### Software Dependencies
It is recommended to create a new conda environment:
```
conda create -n python37 python=3.7

# Activate this environment:
conda activate python37
```
       • numpy (v1.15.0)
            conda install -c conda-forge numpy
       • pandas (v0.24.2)
            conda install -c conda-forge pandas
       • minimap2 (v2.24)
            conda install -c bioconda minimap2
       • pysam (v0.15.3)
            conda install -c bioconda pysam 
       • seqkit (v2.0.0)
            conda install -c bioconda seqkit 
            
### Test
We have downloaded Klebsiella pneumoniae core gene and some reference genome in dir ./test and  added an example to show how to use cgMSI.
Firstly, generate related library by cgMSI LIB module.

  ```
  cd ../cgMSI/
  tar -zxvf ./test/Kp_alleles.tar.gz
  python cgMSI.py LIB -species Kp -genomesDir ./test/testRef/ -allelePath ./test/Kp_alleles.fasta -alleleTable ./library/Klebsiella_pneumoniae_cgMLST_count.tsv -t 12 -outPutDir ./test/library/
  ```
 Next, detection pathogen strain by cgMSI MAP module.
   ```
  python cgMSI.py MAP -species Kp -t 12  -genomesDir ./test/testRef/ -allelePath ./test/Kp_alleles.fasta -sampleFile ./test/test_01X.fna -alleleTablePath ./library/Klebsiella_pneumoniae_cgMLST_count.tsv  -genomeAlleleMatrix ./library/Kp.tsv -outPutDir ./test
  ```
  The result can be found at dir ./test.
  
### Manual
First of all, we should:
- change directory (cd) to cgMSI folder
- cd into cgMSI directory 
  ```
  cd ../cgMSI
  python cgMSI.py -h
  ```
#### LIB
We need the database of strains, which can be downloaded from NCBI. Also you can add your own genomes to the folder. First you need to make sure that 
genomes belonging to the same species are in one folder, different species are in different folders.  The allele table  and specise alleles file were download from https://www.cgmlst.org/ncs . We download 7 specises' allele tables in ./library that can use directly.The target specise allele file you
can download from the website and merge all loci into a fasta file.
create a new library for a species:
```
python cgMSI.py LIB -genomesDir genomeDIR -allelePath species_alleles.fasta -alleleTable speciesAlleleTable -species speciesName -t threadNumber 

Required arguments:

-genomesDir,              string                    Target species Reference Genome Directory Full Path 

-allelePath,              string                    alleles fasta file,can be download 

-alleleTable,             string                    path of the target specise allele table 

-species,                 string                    species name with No whitespace(if Escherichia coli ,like Ec) for distinguish different species


Optional arguments:

-t,                        int                      Number of threads to use by aligner (bowtie2) if different from default (8)
                    
```

add a genome to a existed species library:
```
python cgMSI.py LIB -addGenome -genomesDir genomeDIR -allelePath species_alleles.fasta -alleleTable speciesAlleleTable -species speciesName -genomeName addGenomeName -genomeFile addGenomeFastaFile -t threadNumber 

Required arguments:

-genomesDir,              string                    directory Full Path of target species Reference Genome  

-allelePath,              string                    alleles fasta file,can be download 

-alleleTable,             string                    path of the target specise allele table 

-species,                 string                    species name with No whitespace(if Escherichia coli ,like Ec) for distinguish different species

-genomeName               string                    the name of the genome added into the library

-genomeFile               string                    full path of the added genome fasta file

Optional arguments:

-t,                       int                       Number of threads to use by aligner (bowtie2) if different from default (8)
     
```

#### MAP
First you need to make sure that the LIB module is finished. MAP module will use library generated previously with LIB module.

call MAP module help for details

```
python cgMSI.py MAP -h

python cgMSI.py MAP  -genomesDir genomeDIR -allelePath species_alleles.fasta -alleleTable speciesAlleleTable -species speciesName -sampleFile sampleFile -outPutDir outPutDir -t threadNumber 

Required arguments:

-genomesDir,              string                    directory Full Path of target species Reference Genome  

-allelePath,              string                    alleles fasta file,can be download 

-alleleTable,             string                    full path of the target specise allele table 

-species,                 string                    species name with No whitespace(if Escherichia coli ,like Ec) for distinguish different species

-genomeName               string                    the name of the genome added into the library

-genomeFile               string                    full path of the added genome fasta file

-sampleFile               string                    full path of sample file(fasta or fastq)

-outPutDir                string                    the dir of the predict result

Optional arguments:
-t,                       int                       Number of threads to use by aligner (bowtie2) if different from default (8)

```

### Test
We have downloaded Klebsiella pneumoniae core gene and some reference genome in dir ./test and  added an example to show how to use cgMSI.
Firstly, generate related library by cgMSI LIB module.

  ```
  tar -zxvf Kp_alleles.tar.gz
  python cgMSI.py LIB -species Kp -genomesDir ./test/testRef/ -allelePath ./test/Kp_alleles.fasta -alleleTable ./library/Klebsiella_pneumoniae_cgMLST_count.tsv -t 12 -outPutDir ./test/library/
  ```
 Next, detection pathogen strain by cgMSI MAP module.
   ```
  python cgMSI.py MAP -species Kp -t 12  -genomesDir ./test/testRef/ -allelePath ./test/Kp_alleles.fasta -sampleFile ./test/test_01X.fna -alleleTablePath ./library/Klebsiella_pneumoniae_cgMLST_count.tsv  -genomeAlleleMatrix ./library/Kp.tsv -outPutDir ./test
  ```
  The result can be found at dir ./test.
