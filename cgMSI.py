#!/usr/bin/env python
# Initial author: Solaiappan Manimaran
# Wrapper file for the following modules:
# patholib: generates host/target genome libraries from ncbi nt database for given taxon IDs
# pathomap: aligns reads to host/target database independent of read type using Bowtie2
# pathoid: reassigns ambiguous reads to the correct genome using statistical models
# pathoreport: Writes sam files to xml format

#usage information: cgMSI.py -h

#       cgMSI 1.0 - Predicts strain of genomes in unassembled nanopore metagenomic data
#       Copyright (C) 2022  YuLab - Xiamen University
#
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os, sys
import pysam

import argparse
# from snipe2.patholib import pathoLib
from LIB import genomeLocus
from LIB import allelesPrepare
from MAP import mapLocus
from MAP import mcSample

from _Version import VERSION


snipedir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,snipedir)

# ===========================================================
# main ()
parser = argparse.ArgumentParser(description="cgMSI")

# create the top-level parser
parser.add_argument('--version', action='version', version=VERSION)
# parser.add_argument('--verbose', action='store_const', dest='verbose',
#     required=False, const=True, default=False, help='Prints verbose text while running')
subparsers = parser.add_subparsers(dest='subcommand', help='Select one of the following sub-commands')

#########################################
# create the parser for the "Alleles" command
# parser_a = subparsers.add_parser('Alleles', help='Alleles files prepare for cgMSI')
# parser_a.add_argument('-allelesDir', action='store', dest='allelesDir', required=True,
#     help='species alleles file dir (One dir per species)')
#
# parser_a.add_argument('-species', action='store', dest='species', required=True,
#     help='species name No whitespace(if Escherichia coli ,like Ec)')

#########################################
# create the parser for the "LIB" command
parser_a = subparsers.add_parser('LIB', help='cgMSI Library creation Module')

parser_a.add_argument('-addGenome', action='store_const', dest='addGenome', required=False,const=True,
    default=False, help='add a genome to a Already built library' )

parser_a.add_argument('-genomesDir', action='store', dest='genomesDir', required=True,
    help='species reference genome dir (One dir per species)')

parser_a.add_argument('-allelePath', action='store', dest='allelePath', required=True,
    help='alleles fasta file')

parser_a.add_argument('-alleleTablePath', action='store', dest='alleleTablePath', required=True,
    help='alleles information table')

parser_a.add_argument('-species', action='store', dest='species', required=True,
    help='species name No whitespace(if Escherichia coli ,like Ec)')

parser_a.add_argument('-t', dest='threads', required=False, default=1, help='minimap2 threads numbers' )

parser_a.add_argument('-genomeName', action='store', dest='genomeName', required=False, default="",
    help='the genome name to add (use only add one genome)' )

parser_a.add_argument('-genomeFile', action='store', dest='genomeFile', required=False, default="",
    help='the genome fasta file absolute path (use only add one genome)')

parser_a.add_argument('-outPutDir', dest='outPutDir', required=False, help='output path')
#########################################
# create the parser for the "MAP" command
parser_a = subparsers.add_parser('MAP', help='cgMSI Map Module')
# species name
parser_a.add_argument('-species', action='store', dest='species', required=True,
    help='species name without whitespace(if Escherichia coli ,like Ec)')

# species reference dir
parser_a.add_argument('-genomesDir', action='store', dest='genomesDir', required=True,
    help='species reference genome dir (One dir per species)')

# alleles file(fasta)
parser_a.add_argument('-allelePath', action='store', dest='allelePath', required=True,
    help='alleles fasta file')

# sample file，fasta or fastq formats
parser_a.add_argument('-sampleFile', action='store', dest='sampleFile', required=True,
    help='sample file path ,fasta or fastq file')

# allele information
parser_a.add_argument('-alleleTablePath', action='store', dest='alleleTablePath', required=True,
    help='alleles information table')

# genome allele matrix
parser_a.add_argument('-genomeAlleleMatrix', action='store', dest='genomeAlleleMatrix', required=True,
    help='the alleles table of genomes at each locus')


# thread number
parser_a.add_argument('-t', dest='threads', required=False, default=1, help='minimap2 threads numbers')

# output dir
parser_a.add_argument('-outPutDir', dest='outPutDir', required=False, default="./", help='output path')

def main():
    # parse some argument lists
    inputArgs = parser.parse_args()
    # if (inputArgs.subcommand=='Alleles'):
    #     allelesDir = inputArgs.allelesDir
    #     species = inputArgs.species
    #     allelesPrepare.rename_merge_file(allelesDir,species)

    if (inputArgs.subcommand=='LIB'):
        if inputArgs.addGenome == False:
            libOptions = genomeLocus.LibOptions()
            libOptions.species = inputArgs.species
            libOptions.genomesDir =  os.path.abspath(inputArgs.genomesDir)
            libOptions.allelePath =  os.path.abspath(inputArgs.allelePath)
            libOptions.alleleTablePath =  os.path.abspath(inputArgs.alleleTablePath)
            libOptions.threads = inputArgs.threads
            libOptions.outPutDir = os.path.abspath(inputArgs.outPutDir)
            # generate GCF to NC table
            GCF2NCTsvPath = genomeLocus.getGCF2NC(libOptions)
            print("Generate %s GCF number to NC number Table at %s" %(libOptions.species, GCF2NCTsvPath))
            # identity allele numbers at every locus  for each genome
            LocusTsvPath = genomeLocus.speciesGenomeLocusAnalyze(libOptions)
            print("Generate %s genome Locus Table at %s" %(libOptions.species, LocusTsvPath))

        else:
            # add a genome to the library
            libOptions = genomeLocus.LibOptions()
            libOptions.species = inputArgs.species
            libOptions.genomesDir =  os.path.abspath(inputArgs.genomesDir)
            libOptions.allelePath =  os.path.abspath(inputArgs.allelePath)
            libOptions.alleleTablePath =  os.path.abspath(inputArgs.alleleTable)
            libOptions.threads = inputArgs.threads
            libOptions.genomeName = inputArgs.genomeName
            libOptions.genomeFile = inputArgs.genomeFile
            libOptions.outPutDir = os.path.abspath(inputArgs.outPutDir)
            if libOptions.genomeFile == "" or libOptions.genomeName=="":
                print("Add a genome to library need genomeName and genomeFile")
            else:
                genomeLocus.addGenome(libOptions)


    if (inputArgs.subcommand=='MAP'):


        mapOptions = mapLocus.MapOptions()
        mapOptions.species = inputArgs.species
        mapOptions.genomesDir =  os.path.abspath(inputArgs.genomesDir)
        mapOptions.alleleTablePath = os.path.abspath(inputArgs.alleleTablePath)
        mapOptions.genomeLocusPath = inputArgs.genomeAlleleMatrix
        mapOptions.sampleFile =  os.path.abspath(inputArgs.sampleFile)
        mapOptions.allelePath =  os.path.abspath(inputArgs.allelePath)
        mapOptions.threads = inputArgs.threads
        mapOptions.outPutDir =  os.path.abspath(inputArgs.outPutDir)

        #fastq file convert ro fasta file
        if  mapOptions.sampleFile.split(".")[-1]=="fastq":
            cmd = "seqkit fq2fa %s -o %s" %(mapOptions.sampleFile, mapOptions.sampleFile.split(".")[0]+".fasta")
            os.system(cmd)

        #species allele pool  map to sample
        mapOptions.samFile = mapOptions.outPutDir + '/'+inputArgs.sampleFile.split("/")[-1].split(".")[0] + ".sam"
        cmd = "minimap2 -ax map-ont -t 10 -p 0.6 -N 20 -K 2G --sam-hit-only -c %s %s >%s "%(mapOptions.sampleFile,mapOptions.allelePath,mapOptions.samFile)
        os.system(cmd)

        # length of sample reads
        bf = pysam.Fastafile(mapOptions.sampleFile)
        lengths = bf.lengths
        readsNum = len(lengths)
        lengthSum = sum(lengths)
        avgLength = lengthSum / readsNum

        #core genome length
        genome = os.listdir(mapOptions.genomesDir)[0]
        bf = pysam.Fastafile(mapOptions.genomesDir+"/"+genome)
        genomeLength = max(bf.lengths)

        # MC method estimate abundance
        p = mcSample.readsCount(mapOptions.sampleFile, genomeLength, mapOptions.alleleTablePath)

        #analyse
        L = mapLocus.Reads_sam_analy(mapOptions)
        flag, result, count,maxGenomeLocus = mapLocus.selectReadsMiminmap(L, mapOptions)
        estCov = count / p * avgLength / genomeLength
        print(estCov,count)


        #generate result file
        if os.path.exists("%s/%s_%s.txt" % (mapOptions.outPutDir, inputArgs.sampleFile.split("/")[-1].split(".")[0],mapOptions.species)):
            os.remove("%s/%s_%s.txt" % (mapOptions.outPutDir, inputArgs.sampleFile.split("/")[-1].split(".")[0],mapOptions.species))
        with open("%s/%s_%s.txt" % (mapOptions.outPutDir, inputArgs.sampleFile.split("/")[-1].split(".")[0],mapOptions.species), "a") as f:
            if flag:
                f.write("Target pathogen (%s)  were detected\n"%mapOptions.species)
                c = 1
                for gcf, row in result.iterrows():
                    if c <= 1:
                        f.write("Strain:%s\n"%gcf)
                        c = c + 1
                f.write("Estimate coverage：%s \n"%estCov)
            else:
                f.write("No target pathogen (%s)  were detected"%mapOptions.species)


if __name__ == "__main__":
    main()
