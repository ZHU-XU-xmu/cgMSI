import os
import sys
import pysam
import pandas as pd

class LibOptions:
    # species name
    species = "Ec"
    # genomes Dir
    genomesDir = "../"
    # allelefile file path
    allelePath = "../alleles/Ec_alleles.fasta"
    # alleleTable
    alleleTablePath = "../library/Escherichia_coli_alleles_count.tsv"
    #genome locus table path
    genomeLocusPath = ""
    #minimap2 threads
    threads = 12
    #add Genome's name
    genomeName = ""
    #add Genome file name
    genomeFile = ""
    # library dir
    outPutDir = os.path.abspath(os.path.dirname(os.path.dirname(__file__))+"/library")


def loadLocus(path):
    # 得到每个位点的名称，返回一个list
    Locus = pd.read_csv(path,sep='\t')
    speciesLocus = Locus["Locus"].tolist()
    return speciesLocus

def genomeLocusAnalyze(samfile,alleleTable):
    # 统计基因组wgMLST 分型结果
    # 输入：等位基因集与样本reads 比对sam文件
    # 输出：基因组的分型结果，字典

    speciesLocus= loadLocus(alleleTable)
    # 存放一个基因组的cgMlst的分型结果
    dict_locus = {}
    total_length = 0
    for i in speciesLocus:
        dict_locus[i] = 0

    bf = pysam.AlignmentFile(samfile,'r')

    for r in bf:
        if r.has_tag("NM") and 'S' not in str(r.cigarstring) and 'H' not in str(r.cigarstring):
            if r.get_tag("NM")==0:
                query_locus = r.query_name.split('_')[0]
                if len(r.query_name.split('_'))>2:
                    for item in r.query_name.split('_')[1:-1]:
                        query_locus += "_"+item
                query_alleles = r.query_name.split('_')[-1]
                # print(query_locus,'    ',query_alleles)
                dict_locus[query_locus]=query_alleles
                # print(query_locus,query_alleles)
                total_length +=r.query_length

    bf.close()
    return dict_locus

def addGenome(LibOption):
    # 在基因组分型表中加入一个基因组
    # 参数：genomeName:基因名称
    #      genome: 基因组序列
    # path = "/data1/zhuxu_data/MLST_nanopore/Kl.tsv"
    genomesDir = LibOption.genomesDir
    LocusDir = LibOption.outPutDir
    species = LibOption.species
    allelePath = LibOption.allelePath
    alleleTable = LibOption.alleleTablePath
    minimapPath = LibOption.minimapPath
    genomeName = LibOption.genomeName
    genomeFile = LibOption.genomeFile
    # genomeLocusPath = LibOption.genomeLocusPath
    threads = LibOption.threads
    if LibOption.minimapPath != "":
        minimap = minimapPath + "/" + "minimap2"
    else:
        minimap = "minimap2"

    samFile = "../alleles_to_%s.sam" % genomeName
    cmd = "%s -ax map-ont -t %d -c -K 2G --sam-hit-only  --secondary=no  %s/%s  %s -o %s" % (
    minimap, threads, genomesDir, genomeFile, allelePath, samFile)
    os.system(cmd)
    genomeLocu = genomeLocusAnalyze(samFile,alleleTable)
    genomeLocus = pd.read_csv('%s/%s.tsv' % (LocusDir, species), sep='\t')
    genomeLocus[genomeName] = genomeLocu.values()
    genomeLocus.to_csv('%s/%s.tsv' % (LocusDir, species), sep='\t', index=None)
    genomesLocusDistance('%s/%s.tsv' % (LocusDir, species), species)
    os.remove(samFile)


def getGCF2NC(liboption):
    genomesDir = liboption.genomesDir
    species = liboption.species
    outPutDir = liboption.outPutDir

    genomes = os.listdir(genomesDir)
    if os.path.exists("%s/%s_GCF_NC.tsv"%(outPutDir,species)):
        os.remove("%s/%s_GCF_NC.tsv"%(outPutDir,species))
    with open("%s/%s_GCF_NC.tsv"%(outPutDir,species), "a+") as file:
        file.write("GCF NAME\n")
        for genome in genomes:
            if genome.split(".")[-1]=='fna' or genome.split(".")[-1] =="fasta":
                absPath = os.path.join(genomesDir, genome)
                cmd = "seqkit seq -n %s >%s/temp.txt" % (absPath,outPutDir)
                os.system(cmd)
                genomeName = genome.split(".")[0]
                file.write(genomeName)

                with open("%s/temp.txt"%outPutDir,'r') as f:
                    for line in f:
                        file.write("\t")
                        file.write(line.split(" ")[0].split("\n")[0])
                file.write("\n")
        os.remove("%s/temp.txt"%outPutDir)
    return "%s/%s_GCF_NC.tsv"%(outPutDir,species)


def speciesGenomeLocusAnalyze(libOption):

    genomesDir= os.path.abspath(libOption.genomesDir)
    LocusDir = os.path.abspath(libOption.outPutDir)
    species = libOption.species
    allelePath = os.path.abspath(libOption.allelePath)
    alleleTable = os.path.abspath(libOption.alleleTablePath)

    # GCF2NC_Path = getGCF2NC(genomesDir,species)
    # print("Generate GCF to NC file !")
    genomeFile = os.listdir(genomesDir)
    i = 1
    # Locus = pd.DataFrame()
    if os.path.exists('%s/%s.tsv'%(LocusDir,species)):
        Locus = pd.read_csv('%s/%s.tsv'%(LocusDir,species),sep='\t')
    else:
        Locus = pd.DataFrame()
    speciesLocus = loadLocus(alleleTable)
    Locus["Locus"]=speciesLocus
    for file in genomeFile:
        # each genome mapped with allele pool
        if i > 20000:
            #for debug
            break
        else:
            i = i +1
            # genome GCF name
            fileName = file.split('.')[0]
            if fileName not in Locus.columns.values.tolist():
                samFile = "../alleles_to_%s.sam"%fileName
                if os.path.isfile(samFile)==False:
                    cmd = "minimap2 -ax map-ont -t %s -c -K 2G --sam-hit-only  --secondary=no  %s/%s  %s -o %s"%(libOption.threads,genomesDir,file,allelePath, samFile)
                    print(cmd)
                    os.system(cmd)
                genomeLocus = genomeLocusAnalyze(samFile,alleleTable)
                Locus['%s'%fileName] = genomeLocus.values()
                os.remove(samFile)
            else:
                print(fileName)

        Locus.to_csv('%s/%s.tsv'%(LocusDir,species),sep='\t',index=False)
    return '%s/%s.tsv' % (LocusDir, species)

if __name__ == '__main__':
    print("start")




