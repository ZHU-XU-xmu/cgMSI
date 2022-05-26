#!/usr/bin/env python
# Initial author: Solaiappan Manimaran
# Wrapper file for the following modules:
# patholib: generates host/target genome libraries from ncbi nt database for given taxon IDs
# pathomap: aligns reads to host/target database independent of read type using Bowtie2
# pathoid: reassigns ambiguous reads to the correct genome using statistical models
# pathoreport: Writes sam files to xml format

#usage information: snipe.py -h

#       snipe 1.0 - Predicts strains of genomes in unassembled Nextgen seq data
#       Copyright (C) 2020  YuLab - Xiamen University
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
snipedir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,snipedir)

import argparse
# from snipe2.patholib import pathoLib
from LIB import genomeLocus
from time import time
from _Version import VERSION
# ===========================================================
# main ()
parser = argparse.ArgumentParser(description="Snipe")

# create the top-level parser
parser.add_argument('--version', action='version', version=VERSION)
# parser.add_argument('--verbose', action='store_const', dest='verbose',
#     required=False, const=True, default=False, help='Prints verbose text while running')
subparsers = parser.add_subparsers(dest='subcommand', help='Select one of the following sub-commands')



# create the parser for the "LIB" command
parser_a = subparsers.add_parser('LIB', help='Nano_cgMLST Library creation Module')
parser_a.add_argument('-genomesDir', action='store', dest='genomesDir', required=True,
    help='species reference genome dir (One dir per species)')

parser_a.add_argument('-allelePath', action='store', dest='allelePath', required=True,
    help='alleles fasta file')

parser_a.add_argument('-alleleTable', action='store', dest='alleleTable', required=True,
    help='alleles information table')

parser_a.add_argument('-species', action='store', dest='species', required=True,
    help='species name No whitespace(if Escherichia coli ,like Ec)')

parser_a.add_argument('-minimap2Path', action='store', dest='minimapPath',required=True,
    help='minimap2 absolute path in your environment')

parser_a.add_argument('-t', action='store_const', dest='threads', required=False, const=True,
    default=1, help='minimap2 threads numbers' )

parser_a.add_argument('-addGenome', action='store_const', dest='addGenome', required=False, type=bool,
    default=False, help='add a genome to a Already built library' )

parser_a.add_argument('-genomeName', action='store', dest='genomeName', required=False, default="",
    help='the genome name to add (use only add one genome)' )

parser_a.add_argument('-genomeFile', action='store', dest='genomeFile', required=False, default="",
    help='the genome fasta file absolute path (use only add one genome)')

def main():
    # parse some argument lists
    inputArgs = parser.parse_args()

    #### PathoID modules ####

    start = time()

    if (inputArgs.subcommand=='LIB'):
        ################################################$
        #append taxon id in the front of sequence header
        ################################################$

        if inputArgs.addGenome == False:
            libOptions = genomeLocus.LibOptions()
            libOptions.minimapPath = inputArgs.minimapPath
            libOptions.species = inputArgs.species
            libOptions.genomesDir = inputArgs.genomesDir
            libOptions.allelePath = inputArgs.allelePath
            libOptions.alleleTable = inputArgs.alleleTable
            libOptions.threads = inputArgs.threads
            # identity allele numbers at every locus  for each genome
            genomeLocus.speciesGenomeLocusAnalyze(libOptions)
            #compute distance between genomes
            LocusTsvPath = './library/%s.tsv'%(libOptions.species)
            genomeLocus.genomesLocusDistance(LocusTsvPath, libOptions.species)

        else:
    #         加一个基因组到库中
            libOptions = genomeLocus.LibOptions()
            libOptions.minimapPath = inputArgs.minimapPath
            libOptions.species = inputArgs.species
            libOptions.genomesDir = inputArgs.genomesDir
            libOptions.allelePath = inputArgs.allelePath
            libOptions.alleleTable = inputArgs.alleleTable
            libOptions.threads = inputArgs.threads
            libOptions.genomeName = inputArgs.genomeName
            libOptions.genomeFile = inputArgs.genomeFile
            if libOptions.genomeFile == "" or libOptions.genomeName=="":
                print("add a genome to library need genomeName and genomeFile")
            else:
                genomeLocus.addGenome(libOptions)



    if (inputArgs.subcommand=='MAP'):
        pathoMapOptions = PathoMapA.PathoMapOptions()
        pathoMapOptions.verbose = inputArgs.verbose
        pathoMapOptions.outDir = inputArgs.map_outdir
        pathoMapOptions.indexDir = inputArgs.map_indexdir
        pathoMapOptions.outAlignFile = inputArgs.map_outalign
        pathoMapOptions.inReadFile = inputArgs.map_inputread
        pathoMapOptions.inReadFilePair1 = inputArgs.map_inputread1
        pathoMapOptions.inReadFilePair2 = inputArgs.map_inputread2
        pathoMapOptions.targetAlignParameters = inputArgs.map_targetalignparams
        pathoMapOptions.filterAlignParameters = inputArgs.map_filteralignparams
        if (len(inputArgs.map_targetref)>0):
            pathoMapOptions.targetRefFiles = inputArgs.map_targetref.split(",")
        if (len(inputArgs.map_filterref)>0):
            pathoMapOptions.filterRefFiles = inputArgs.map_filterref.split(",")
        if (len(inputArgs.map_targetindex)>0):
            pathoMapOptions.targetIndexPrefixes = inputArgs.map_targetindex.split(",")
        if (len(inputArgs.map_filterindex)>0):
            pathoMapOptions.filterIndexPrefixes = inputArgs.map_filterindex.split(",")
        if (len(inputArgs.map_targetalign)>0):
            pathoMapOptions.targetAlignFiles = inputArgs.map_targetalign.split(",")
        if (len(inputArgs.map_filteralign)>0):
            pathoMapOptions.filterAlignFiles = inputArgs.map_filteralign.split(",")
        pathoMapOptions.btHome = inputArgs.map_bthome
        pathoMapOptions.numThreads = inputArgs.map_numthreads
        pathoMapOptions.exp_tag = inputArgs.map_exp_tag + "-"
        PathoMapA.processPathoMap(pathoMapOptions)

    if (inputArgs.subcommand=='ID'):
        pathoIdOptions = PathoID.PathoIdOptions(inputArgs.id_ali_file)
        pathoIdOptions.ali_format = inputArgs.id_ali_format
        pathoIdOptions.verbose = inputArgs.verbose
        pathoIdOptions.out_matrix_flag = inputArgs.id_out_matrix
        pathoIdOptions.score_cutoff = inputArgs.id_score_cutoff
        pathoIdOptions.exp_tag = inputArgs.id_exp_tag
        pathoIdOptions.outdir = inputArgs.id_outdir
        pathoIdOptions.emEpsilon = inputArgs.id_emEpsilon
        pathoIdOptions.maxIter = inputArgs.id_maxIter
        pathoIdOptions.piPrior = inputArgs.id_piPrior
        pathoIdOptions.thetaPrior = inputArgs.id_thetaPrior
        pathoIdOptions.noalign = inputArgs.id_noalign
        pathoIdOptions.noCutOff = inputArgs.id_nocutoff
        PathoID.snipe_reassign(pathoIdOptions)

     ########## snipeRec model ############
    if inputArgs.subcommand=='REC':
        PathoRecOptions = PathoRec.PathoRecOptions()
        if len(inputArgs.ssr_reference)>0:
            PathoRecOptions.path_core = inputArgs.ssr_reference
        PathoRecOptions.read1 = inputArgs.rec_inputread1
        PathoRecOptions.read2 = inputArgs.rec_inputread2
        if len(inputArgs.id_rep_file)>0:
            PathoRecOptions.id_file = inputArgs.id_rep_file
        if len(inputArgs.dict_target)>0:
            PathoRecOptions.dict_ta = inputArgs.dict_target
        if len(inputArgs.dict_template)>0:
            PathoRecOptions.dict_te = inputArgs.dict_template
        if len(inputArgs.rec_outDir)>0:
            PathoRecOptions.out_path = inputArgs.rec_outDir
        PathoRecOptions.threads = inputArgs.rec_numthreads
        PathoRecOptions.l_value = inputArgs.likelihood_Ratio_L
        PathoRecOptions.p0_value = inputArgs.probability_p0
        PathoRecOptions.p1_value = inputArgs.probability_p1
        PathoRecOptions.exp_tag = inputArgs.rec_exp_tag
        PathoRec.rec(PathoRecOptions)

    if (inputArgs.subcommand=='REP'):
        pathoReportOptions = PathoReportA.PathoReportOptions(inputArgs.rep_ali_file)
        pathoReportOptions.verbose = inputArgs.verbose
        pathoReportOptions.contigFlag = inputArgs.rep_contig_flag
        pathoReportOptions.outDir = inputArgs.rep_outdir
        pathoReportOptions.samtoolsHome = inputArgs.rep_samtoolshome
        pathoReportOptions.noCutOff = inputArgs.rep_nocutoff
        mysqlConf=(inputArgs.rep_dbhost,inputArgs.rep_dbport,inputArgs.rep_dbuser,
            inputArgs.rep_dbpasswd,inputArgs.rep_db)
        pathoReportOptions.mysqlConf = mysqlConf
        PathoReportA.processPathoReport(pathoReportOptions)

    if (inputArgs.subcommand=='QC'):
        qcargs = sys.argv[2:]
        pathoqcdir = snipedir + os.path.sep + 'snipe' + os.path.sep + 'pathoqc'
        pathoqcfile = pathoqcdir + os.path.sep + 'pathoqc.py'
        if os.path.exists(pathoqcfile):
            cmd = sys.executable
            cmd += " " + pathoqcfile + " "
            cmd += " ".join(qcargs)
            print(cmd)
            os.system(cmd)
        else:
            print("PathoQC (" + pathoqcfile + ") not found. Please download pathoqc_vXXX.tar.gz and "
            "install it ("+pathoqcdir+") from http://sourceforge.net/projects/pathoscope/")

    elapsed = time() - start;
    if inputArgs.verbose:
        print("Total Elapsed Time: %d" % (elapsed))

if __name__ == "__main__":
    main()
