import os
import pysam
import pandas as pd

class LibOptions:
    # 物种名
    species = "Ec"
    # 基因组目录
    genomesDir = "../"
    # allelefile文件的位置
    allelePath = "../alleles/Ec_alleles.fasta"
    # alleleTable
    alleleTablePath = "../library/Escherichia_coli_alleles_count.tsv"
    #基因组的分型结果文件
    genomeLocusPath = ""
    #minimap2比对线程数
    threads = 12
    #addGenome时，加入基因的名字
    genomeName = ""
    #addGenome时，加入基因的fna文件
    genomeFile = ""
    # library目录
    LocusDir = "../library"


def loadLocus(path):
    # 得到每个位点的名称，返回一个list
    Locus = pd.read_csv(path,sep='\t')
    speciesLocus = Locus["Locus"].tolist()
    return speciesLocus

def genomesLocusDistance(path,species):
    # 计算基因组之间的分型距离
    # 输入：基因组分型tsv文件
    # 输出：基因组与基因组分型距离矩阵
    genomesLocus = pd.read_csv(path,sep='\t')
    colunmsName = genomesLocus.columns.tolist()     #列名，list
    LocusName = genomesLocus[colunmsName[0]]        #位点名，list
    LocusDistance = pd.DataFrame()                  #记录结果
    LocusDistance['GCF name'] = colunmsName[1:]
    outPutDir = os.path.dirname(path)
    outPutName = "%s_Locus_Distance"%species
    for targetGenomeName in colunmsName[1:]:
        targetGenomeLocus = genomesLocus[targetGenomeName].tolist()
        targetList = []   #当前基因组与其它基因组的分型距离

        for backGenomeName in colunmsName[1:]:
            backGenomeLocus = genomesLocus[backGenomeName].tolist()
            distanceCount = 0  # 两个基因组分型距离
            for i in range(len(LocusName)):
                if targetGenomeLocus[i]!=backGenomeLocus[i]:
                    distanceCount = distanceCount +1
                i = i +1
            targetList.append(distanceCount)
        LocusDistance[targetGenomeName] = targetList
    LocusDistance.to_csv('%s/%s.tsv'%(outPutDir,outPutName),sep='\t',index=False)

    return '%s/%s.tsv'%(outPutDir,outPutName)
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
                query_alleles = r.query_name.split('_')[1]
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
    LocusDir = LibOption.LocusDir
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


def getGCF2NC(genomeDir,species):

    genomes = os.listdir(genomeDir)
    if os.path.exists("../library/%s_GCF_NC.tsv"%species):
        os.remove("../library/%s_GCF_NC.tsv"%species)
    with open("../library/%s_GCF_NC.tsv"%species, "a+") as file:
        file.write("GCF NAME\n")
        for genome in genomes:
            absPath = os.path.join(genomeDir, genome)
            cmd = "seqkit seq -n %s >./temp.txt" % absPath
            os.system(cmd)
            genomeName = genome.split(".")[0]
            file.write(genomeName)

            with open("../library/temp.txt") as f:
                for line in f:
                    file.write("\t")
                    file.write(line.split(" ")[0].split("\n")[0])
            file.write("\n")
        os.remove("../library/temp.txt")
    return "../library/%s_test.tsv"%species

def speciesGenomeLocusAnalyze(LibOption):
    # 计算物种文件夹下所有基因组wgMLST分型
    # 输入：genomesPath 物种基因组所在文件夹
    #      LocusPath       物种等位基因集,fasta文件
    # 输出：tsv文件，记录所有基因组分型
    genomesDir= LibOption.genomesDir
    LocusDir = LibOption.LocusDir
    species = LibOption.species
    allelePath = LibOption.allelePath
    alleleTable = LibOption.alleleTablePath

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
        # 每个基因组分别与等位基因集比较
        if i > 20000:
            #调试用，控制计算基因组的数量
            break
        else:
            i = i +1
            # 取基因组GCF名
            fileName = file.split('.')[0]
            if fileName not in Locus.columns.values.tolist():
                samFile = "../alleles_to_%s.sam"%fileName
                if os.path.isfile(samFile)==False:
                    cmd = "minimap2 -ax map-ont -t %d -c -K 2G --sam-hit-only  --secondary=no  %s/%s  %s -o %s"%(LibOption.threads,genomesDir,file,allelePath, samFile)
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


