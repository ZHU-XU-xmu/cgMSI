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
parser_a = subparsers.add_parser('Alleles', help='Alleles files prepare for cgMSI')
parser_a.add_argument('-allelesDir', action='store', dest='allelesDir', required=True,
    help='species alleles file dir (One dir per species)')

parser_a.add_argument('-species', action='store', dest='species', required=True,
    help='species name No whitespace(if Escherichia coli ,like Ec)')

#########################################
# create the parser for the "LIB" command
parser_a = subparsers.add_parser('LIB', help='cgMSI Library creation Module')
parser_a.add_argument('-genomesDir', action='store', dest='genomesDir', required=True,
    help='species reference genome dir (One dir per species)')

parser_a.add_argument('-allelePath', action='store', dest='allelePath', required=True,
    help='alleles fasta file')

parser_a.add_argument('-alleleTablePath', action='store', dest='alleleTablePath', required=True,
    help='alleles information table')

parser_a.add_argument('-species', action='store', dest='species', required=True,
    help='species name No whitespace(if Escherichia coli ,like Ec)')

parser_a.add_argument('-t', action='store_const', dest='threads', required=False, const=True,
    default=1, help='minimap2 threads numbers' )

parser_a.add_argument('-addGenome', action='store_const', dest='addGenome', required=False, type=bool,
    default=False, help='add a genome to a Already built library' )

parser_a.add_argument('-genomeName', action='store', dest='genomeName', required=False, default="",
    help='the genome name to add (use only add one genome)' )

parser_a.add_argument('-genomeFile', action='store', dest='genomeFile', required=False, default="",
    help='the genome fasta file absolute path (use only add one genome)')

#########################################
# create the parser for the "MAP" command
parser_a = subparsers.add_parser('MAP', help='cgMSI Map Module')
# 物种名
parser_a.add_argument('-species', action='store', dest='species', required=True,
    help='species name No whitespace(if Escherichia coli ,like Ec)')

# 基因组目录文件
parser_a.add_argument('-genomesDir', action='store', dest='genomesDir', required=True,
    help='species reference genome dir (One dir per species)')

# alleles文件
parser_a.add_argument('-allelePath', action='store', dest='allelePath', required=True,
    help='alleles fasta file')

# sample文件，fasta或fastq格式
parser_a.add_argument('-sampleFile', action='store', dest='sampleFile', required=True,
    help='sample file path ,fasta or fastq file')

# alleleTable文件
parser_a.add_argument('-alleleTablePath', action='store', dest='alleleTablePath', required=True,
    help='alleles information table')

# 基因组分型结果文件
parser_a.add_argument('-genomeLocusPath', action='store', dest='genomeLocusPath', required=True,
    help='the alleles table of genomes at each locus')

# 基因组之间的距离
parser_a.add_argument('-genomeDistancePath', action='store', dest='genomeDistancePath', required=True,
    help='the locus distance of each genome')

# 比对线程数目
parser_a.add_argument('-t', action='store_const', dest='threads', required=False, const=True,
    default=1, help='minimap2 threads numbers')

# 结果输出目录
parser_a.add_argument('-outPutDir', action='store_const', dest='outPutDir', required=False, default="../result",
    help='output path')

def main():
    # parse some argument lists
    inputArgs = parser.parse_args()
    if (inputArgs.subcommand=='Alleles'):
        allelesDir = inputArgs.allelesDir
        species = inputArgs.species
        allelesPrepare.rename_merge_file(allelesDir,species)

    if (inputArgs.subcommand=='LIB'):

        if inputArgs.addGenome == False:
            libOptions = genomeLocus.LibOptions()
            libOptions.species = inputArgs.species
            libOptions.genomesDir = inputArgs.genomesDir
            libOptions.allelePath = inputArgs.allelePath
            libOptions.alleleTablePath = inputArgs.alleleTablePath
            libOptions.threads = inputArgs.threads

            # identity allele numbers at every locus  for each genome
            LocusTsvPath = genomeLocus.speciesGenomeLocusAnalyze(libOptions)
            #compute distance between genomes
            print("Generate %s genome Locus Table at %s" %(libOptions.species, LocusTsvPath))
            genomeDistancePath= genomeLocus.genomesLocusDistance(LocusTsvPath, libOptions.species)
            print("Generate %s genome Distance Table at %s" %(libOptions.species, genomeDistancePath))

        else:
    #         加一个基因组到库中
            libOptions = genomeLocus.LibOptions()
            libOptions.minimapPath = inputArgs.minimapPath
            libOptions.species = inputArgs.species
            libOptions.genomesDir = inputArgs.genomesDir
            libOptions.allelePath = inputArgs.allelePath
            libOptions.alleleTablePath = inputArgs.alleleTable
            libOptions.threads = inputArgs.threads
            libOptions.genomeName = inputArgs.genomeName
            libOptions.genomeFile = inputArgs.genomeFile
            if libOptions.genomeFile == "" or libOptions.genomeName=="":
                print("Add a genome to library need genomeName and genomeFile")
            else:
                genomeLocus.addGenome(libOptions)


    if (inputArgs.subcommand=='MAP'):
        mapOptions = mapLocus.MapOptions()
        mapOptions.species = inputArgs.species
        mapOptions.genomesDir = inputArgs.genomesDir
        mapOptions.alleleTablePath = inputArgs.alleleTablePath
        mapOptions.genomeLocusPath = inputArgs.genomeLocusPath
        mapOptions.genomeDistancePath = inputArgs.genomeDistancePath
        mapOptions.sampleFile = inputArgs.sampleFile
        mapOptions.allelePath = inputArgs.allelePath
        mapOptions.threads = inputArgs.threads
        mapOptions.outPutDir = inputArgs.outPutDir

        #如果样本是fastq文件，转成fasta文件(原来是fasta文件，则不会变)
        cmd = "seqkit fq2fa %s -o %s" %(mapOptions.sampleFile, mapOptions.samFile.split(".")[0]+".fasta")
        os.system(cmd)

        # 比对
        mapOptions.samFile = mapOptions.outPutDir + inputArgs.sampleFile.split("/")[-1].split(".")[0] + ".sam"
        cmd = "minimap2 -ax map-ont -t %s -p 0.6 -N 20 -K 2G --sam-hit-only -c %s %s >%s "%(mapOptions.threads,mapOptions.sampleFile,mapOptions.allelePath,mapOptions.samFile)
        os.system(cmd)

        #样本reads 平均长度
        bf = pysam.Fastafile(mapOptions.sampleFile)
        lengths = bf.lengths
        readsNum = len(lengths)
        lengthSum = sum(lengths)
        avgLength = lengthSum / readsNum

        #core genome length
        genome = os.listdir(mapOptions.genomesDir)[0]
        bf = pysam.Fastafile(genome)
        genomeLength = max(bf.lengths)

        # MC 采样估计丰度
        p = mcSample.readsCount(mapOptions.sampleFile, genomeLength, mapOptions.alleleTablePath)

        #分析
        L = mapLocus.Reads_sam_analy(mapOptions)
        result,maxgenomeLocusRatio,count = mapLocus.selectReadsMiminmap(L, mapOptions)
        estCov = count / p * avgLength / genomeLength
        resultName = fastaFile.split("/")[-1].split(".")[0] + "_result.txt"
        with open("%s/%s" % (mapOptions.outPutDir,resultName ), "a+") as f:
            if maxgenomeLocusRatio < 0.15:
                f.write("no %s strain found in sample!"%mapOptions.species)
            else:
                f.write("Strain \t Coverage\n")

                # save the max score strain
                c = 1
                for gcf, row in result2.iterrows():
                    if c <= 1:
                        f.write("%s\t%s" % (gcf, estCov))
                    c = c + 1

if __name__ == "__main__":
    main()
