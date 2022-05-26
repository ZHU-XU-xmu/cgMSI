import pysam
import math
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd
import os
import subprocess
import argparse
from pandas import DataFrame
from openpyxl import load_workbook
dict_AS_target =[]
dict_NM_target = []
dict_ref_length_target = []
dict_AS_back =[]
dict_NM_back = []
dict_ref_length_back = []
total_number =0
target_number = 0
back_number = 0
flag4_number = 0

class MapOptions:

    species = "Ec"
    samFile  = ""
    genomesDir = "../genomes/562"
    LocusDir = "../library"
    fastaFile = ""
    genomeLocusPath = ""
    genomeDistancePath = ""
    alleleTable = "../library/Escherichia_coli_alleles_count.tsv"
    outPutDir = ""
    minimapPath = ""
    threads = 12


def loadLocus(path):
    # 得到每个位点的名称，返回一个list
    Locus = pd.read_csv(path,sep='\t')
    speciesLocus = Locus["Locus"].tolist()
    return speciesLocus

def sam_ananly(samfile):

    Kl_locus = "/data1/zhuxu_data/MLST_nanopore/Kl_locus"
    f = open('%s' % (Kl_locus), 'r')
    c = f.read()
    Kl_locus = eval(c)
    f.close()
    # dict_locus = {}
    dict_locus = [0 for i in range(len(Kl_locus))]
    # for i in Kl_locus:
    #     dict_locus[i] = 0

    bf = pysam.AlignmentFile(samfile,'r')



    for r in bf:

        if r.has_tag("NM") and 'S' not in str(r.cigarstring) and 'H' not in str(r.cigarstring):
            if r.get_tag("NM")==0:
                query_locus = r.query_name.split('_')[0]+'_'+r.query_name.split('_')[1]
                query_alleles = r.query_name.split('_')[-1]
                # print(query_locus,'    ',query_alleles)
                dict_locus[Kl_locus.index(query_locus)]=query_alleles

                # print(r.get_tag("NM"),r.is_supplementary,r.infer_read_length(),r.infer_query_length(),r.query_alignment_length,r.query_name,r.query_length,r.reference_name,r.reference_length,r.flag)
    bf.close()
    return dict_locus


def Reads_sam_ananly(samfile):

    Kl_locus = "/data1/zhuxu_data/MLST_nanopore/Kl_locus"
    f = open('%s' % (Kl_locus), 'r')
    c = f.read()
    Kl_locus = eval(c)
    f.close()
    # print(Kl_locus)
    dict_locus = [0 for i in range(len(Kl_locus))]
    # for i in Kl_locus:
    #     dict_locus[i] = 0

    bf = pysam.AlignmentFile(samfile,'r')


    R = {}
    #R[readsName] = [等位基因,得分],记录每个locus 比对上的reads的最大得分和对应的等位基因
    L = {}
    maxLocus = 0
    maxAS = {}
    query_locus = ""
    query_alleles = 0
    for r in bf:
        # start_S, end_S = cigar_S_detection(str(r.cigarstring))

        # if start_S > 10 and end_S > 10:
        #     continue


        if r.has_tag("NM"):
            query_locus_now = r.query_name.split('_')[0]+'_'+r.query_name.split('_')[1]
            query_alleles_now = r.query_name.split('_')[-1]
            query_AS_now = r.get_tag("AS")
            query_reference = r.reference_name
            # 位点未更新
            if query_locus_now == query_locus :
                # 记录当前reads在这个位点上的最大得分maxLocus
                if maxLocus < query_AS_now:
                    maxLocus = query_AS_now
                if query_reference in R:
                    # 如果当前等位基因的得分比之前等位基因得分高，则更新对应的等位基因
                    if query_AS_now > R[query_reference][1]:
                        R[query_reference][0] = query_alleles_now
                        R[query_reference][1] = query_AS_now


                # reads没有出现过，加入R中
                else:

                    R[query_reference] = [query_alleles_now,query_AS_now]

            else:
                # 位点更新时
                if query_alleles != 0:
                    #对每个位点对应的等位基因匹配到的reads数目进行统计
                    L[query_locus] = R
                    counter = {} #统计每个等位基因比对上的次数
                    maxAS = {} #统计每个等位基因比对得分最大值（在多条reads上的最大值）
                    # 分别对当前位点的每条reads统计
                    for key in R:
                        if R[key][0] not in counter :
                            if int(maxLocus)-int(R[key][1])<300:
                                counter[R[key][0]] = 1
                                maxAS[R[key][0]] = R[key][1]

                        else:
                            if int(maxLocus) - int(R[key][1]) < 300:
                                counter[R[key][0]] = counter[R[key][0]] + 1
                                if int(maxAS[R[key][0]]) < int(R[key][1]):
                                    maxAS[R[key][0]] = R[key][1]

                    x = sorted(counter.items(),key=lambda s : int(s[1]),reverse=True)
                    maxCount = x[0][1]
                    maxLocusAS = maxAS[x[0][0]]
                    maxAlleles = x[0][0]

                    for i in x :
                        if i[1]==maxCount:
                            if int(maxAS[i[0]]) > int(maxLocusAS):
                                maxLocusAS = maxAS[i[0]]
                                maxAlleles = i[0]
                    dict_locus[Kl_locus.index(query_locus)] = maxAlleles


                    R = {}
                    maxLocus = query_AS_now

                query_locus = query_locus_now
                query_alleles = query_alleles_now
                query_AS = query_AS_now
                R[query_reference] = [query_alleles,query_AS]
            # print(query_locus,'    ',query_alleles)



                # print(r.get_tag("NM"),r.is_supplementary,r.infer_read_length(),r.infer_query_length(),r.query_alignment_length,r.query_name,r.query_length,r.reference_name,r.reference_length,r.flag)
        else:
            print("----------")


    bf.close()
    return dict_locus


def Locus_EM(L,maxIter, emEpsilon,LocusNumber,LocusNames):


    loucusCount = len(LocusNumber)
    pi = []
    for key in L:
        for read in L[key]:
            readSum = 0
            # 位点中每个reads
            Tu = L[key][read]
            for t in range(len(Tu)):
                readSum +=Tu[t][1]
            for i in range(len(Tu)):
                # 每个reads比对上的等位基因
                L[key][read][i] = (Tu[i][0], Tu[i][1]/readSum)

    for i in range(loucusCount):
        locuNumber = LocusNumber[i]
        pi.append([1. /locuNumber  for _ in range(locuNumber)])
    pi_old = []
    k = 0 #位点编号
    for key in LocusNames:
        pi_old = [i for i in pi[k]]
        if key not in L:
            k = k +1
            continue
        # print(key, LocusNumber[k],L[key] )
        # if key  == "b0824":
        #     print(pi[k],L[key])
        for j in range(maxIter):

            # E step
            pi_old = [i for i in pi[k]]
            lReadSum = [0 for _ in range(LocusNumber[k])]
            for read in L[key]:
                Tu = L[key][read]
                loc = []
                for l in Tu:
                    loc.append(int(l[0])-1)
                ltemp = [q[1]*pi[k][int(q[0])-1] for q in Tu]
                lsum = sum(ltemp)
                m = 0
                lnorm = []

                for t in range(LocusNumber[k]):
                    if m<len(Tu):
                        if t in loc:
                            lnorm.append(1.* ltemp[m]/lsum)
                            lReadSum[t] += 1.* ltemp[m]/lsum
                            m += 1
                        else:
                            lnorm.append(0)

                    else:
                        lnorm.append(0)

            #M step
            lReadLocusSum = sum(lReadSum)
            if lReadLocusSum == 0:
                # pi[k]=0

                break
            for t in range(LocusNumber[k]):
                pi[k][t] = 1. * lReadSum[t] / lReadLocusSum

            #判断是否继续迭代
            cutoff = 0.0
            for t in range(len(pi[k])):
                cutoff += abs(pi_old[t] - pi[k][t])

            # print("[%d]%g" % (k, cutoff))
            if (cutoff <= emEpsilon ):
                break

        k += 1


    return pi


def rescaleAS(L,maxAS):
    for key in L:
        # 每个位点
        maxScore = maxAS[key]
        if(maxScore==0):
            print(key)
        scalingFactor = 100.0 / (maxScore)
        for read in L[key]:
            # 位点中每个reads
            Tu = L[key][read]
            for i in range(len(Tu)):
                #每个reads比对上的等位基因
                L[key][read][i] = (Tu[i][0],math.exp(Tu[i][1]*scalingFactor))
    return L

#筛选后的reads与筛选后的基因组比对
def selectReadsMiminmap(pi,L,MapOptions):
    # 基因组分型结果文件
    genomeLocusPath = MapOptions.genomeLocusPath
    genomeDistancePath = MapOptions.genomeDistancePath
    outPutDir = MapOptions.outPutDir
    fastaFile = MapOptions.fastaFile
    # 要去掉fasta文件中的注释才能筛出文件来
    fastaFileRename = fastaFile.split(".")[0]+"_rename.fasta"
    genomesDir = MapOptions.genomesDir
    species = MapOptions.species
    # 基因组分型距离文件
    genomeNameDict = {"Ec": "B1109", "Ef": "B537", "Sa": "B41012", "Se": "B4242", "Pa": "B3509", "Lm": "B33116"}

    select_reads_name = fastaFile.split("/")[-1].split(".")[0] + "_selected.txt"
    select_reads_file = fastaFile.split("/")[-1].split(".")[0] + "_selected" + ".fasta"
    select_ref_file = fastaFile.split("/")[-1].split(".")[0] + "_selected" + "_ref.fna"
    selectMapResultFile = fastaFile.split("/")[-1].split(".")[0] + "_selected" + ".sam"

    genomeLocus = pd.read_csv('%s' % (genomeLocusPath), sep='\t')
    genomeDistance = pd.read_csv('%s' % (genomeDistancePath), index_col='GCF name', sep='\t')
    genomeName = genomeLocus.columns.tolist()
    LocusName = genomeLocus[genomeName[0]]
    Count = {}

    targetName = genomeNameDict[species]
    # 样本和参考基因组一一比较
    for genomename in genomeName[1:]:
        count = 0
        ref_locus = genomeLocus[genomename].tolist()
        for j in range(len(ref_locus)):
            if int(ref_locus[j]) != 0:
                count += pi[j][int(ref_locus[j]) - 1]

        target_ref_distance = genomeDistance.loc[targetName, genomename]
        Count[genomename] = [count, target_ref_distance]
    result = pd.DataFrame.from_dict(Count, orient='index',
                                    columns=['genome score', 'target and genome distance'])
    result = result.sort_values(by=['genome score'], ascending=False)

    #选择排名靠前的五个基因组，第五名后和第五名相同得分，则也选择
    fifthScore = result['genome score'][4]
    cut = 5
    while(result['genome score'][cut]==fifthScore):
        cut += 1
    Ref = result.index.tolist()[:cut]
    print(Ref)
    for k in Ref :
        print(Count[k])
    #选择的基因组放入一个fasta文件中
    if os.path.exists("%s/%s"%(outPutDir, select_ref_file)):
        os.remove("%s/%s"%(outPutDir, select_ref_file))
    for i in range(len(Ref)):
        refname = Ref[i]
        cmd = "cat %s/%s*.fna >> %s/%s" % (genomesDir, refname, outPutDir, select_ref_file)
        os.system(cmd)


    #选择的reads重新写入一个fasta文件中

    R = []
    for i in range(len(Ref)):
        k = 0
        ref_locus = genomeLocus[Ref[i]].tolist()
        for key in LocusName:
            if key not in L:
                k += 1
                continue
            if int(ref_locus[k]) != 0 and pi[k][int(ref_locus[k])-1] != 0 and len(L[key]) :
                # 判断每条reads是否比对到了候选菌株对应的等位基因上（等位基因的得分是否在前五）
                for r in L[key].keys():
                    for i in range(len(L[key][r])):
                        if int(L[key][r][i][0]) == int(ref_locus[k]):
                            R.append(r)
                            break

            k += 1

    #列表去重，去除重复的reads
    Rset = set(R)
    R = list(Rset)

    #生成筛选后的fasta文件
    if os.path.exists("%s/%s" % (outPutDir, select_reads_name)):
        os.remove("%s/%s" % (outPutDir, select_reads_name))

    if os.path.exists("%s/%s" % (outPutDir, select_reads_file)):
        os.remove("%s/%s" % (outPutDir , select_reads_file))
    with open("%s/%s" % (outPutDir , select_reads_name), "a+") as f:
        for key in R :
                newline = key + "\n"
                f.write(newline)
    cmd ="/home/zhuxu/miniconda3/bin/seqkit replace -p \"\s.+\"  %s > %s" %(fastaFile,fastaFileRename)
    os.system(cmd)
    cmd = "/home/zhuxu/miniconda3/bin/seqkit grep -n -f %s/%s %s -o %s/%s" % (outPutDir , select_reads_name, fastaFileRename, outPutDir , select_reads_file)
    os.system(cmd)

    #
    if os.path.exists("%s/%s" % (outPutDir,selectMapResultFile)):
        os.remove("%s/%s" % (outPutDir,selectMapResultFile))
    cmd = "/home/zhuxu/miniconda3/bin/minimap2 -ax map-ont -t 10 -c -K 2G --sam-hit-only  -N 10 %s/%s %s/%s > %s/%s" % \
    (outPutDir, select_ref_file,outPutDir ,select_reads_file,outPutDir,selectMapResultFile)
    os.system(cmd)
    samplePath = outPutDir+'/'+selectMapResultFile
    mapScore,count = selectedSamAnaly(species,samplePath,Ref)
    Score = sorted(mapScore.items(),key=lambda x: x[1],reverse=True)
    #将两个分数都进行缩放和归一化

    MlstScaleFactor = 100/Count[Ref[0]][0]
    if Score[0][1] ==0:
        MapScaleFactor = 0
    else:
        MapScaleFactor = 100/Score[0][1]
    MlstSum = 0
    MapSum = 0
    for k in Ref:
        Count[k][0] = math.exp(Count[k][0] * MlstScaleFactor)
        MlstSum  = MlstSum + Count[k][0]
        mapScore[k] = math.exp( mapScore[k] * MapScaleFactor )
        MapSum = MapSum + mapScore[k]
    result = {}
    for k in Ref:
        mapScore[k] = mapScore[k] / MapSum
        Count[k][0] =Count[k][0] / MlstSum
        print(k,Count[k][0])
        result[k] = [(Count[k][0] * 0.2 + mapScore[k] * 0.8), Count[k][1]]

    result = pd.DataFrame.from_dict(result, orient='index',
                                    columns=['genome score', 'target and genome distance'])
    result = result.sort_values(by=['genome score'], ascending=False)
    # print(Count)
    print(mapScore)
    print(result)



    return result,count



def selectedSamAnaly(species,samPath,Ref):
    #分析选择的reads比对后的结果
    genomeRefPath = "/data2/zhuxu_data/Nano_cgMLST/library/%s_NC_GCF.tsv"%species
    bf = pysam.AlignmentFile(samPath, 'r')
    score = {}
    count = 0
    R = []
    l = []
    G = {}
    GCF_NC_dict = getGCF_NC_Number(genomeRefPath)
    for r in bf:
        if r.has_tag("AS"):
            if int(r.get_tag('AS'))/int(r.query_alignment_length) > 1.1:
                # 统计比对到Kl 基因组的reads里有多少错误的，这些错误的reads得分如何，
                # 正确的reads得分如何，如何正确区分两种reads
                reads_name = r.query_name
                ref_name = r.reference_name
                # print(r.query_name,r.get_tag("AS"))
                if r.query_name not in R:
                    R.append(r.query_name)
                    count += 1
                    l.append(r.query_length)

                if ref_name in G:
                    G[ref_name][0].append(reads_name)
                    G[ref_name][1].append(r.get_tag('AS'))
                    G[ref_name][2].append(r.query_alignment_length)
                    G[ref_name][3] += int(r.get_tag('AS'))

                else:
                    G[ref_name] = [[reads_name], [r.get_tag('AS')], [r.query_length], int(r.get_tag('AS'))]
    bf.close()
    for key in Ref:
        score[key] = 0
        for NC in G:
            if NC in GCF_NC_dict[key]:
                score[key] += G[NC][3]
    print(count)



    return score, count


def Reads_sam_analy3(MapOptions):
    # EM算法计算可能菌株
    alleleTable = MapOptions.alleleTable
    samFile = MapOptions.samFile
    # speciesLocus = loadLocus(alleleTable)
    # dict_locus = [0 for i in range(len(speciesLocus))]
    bf = pysam.AlignmentFile(samFile, 'r')

    # lk 表示样本中大于li+ 2l2的reads比例

    # locusPath = "/data2/zhuxu_data/cgMLST/Klebsiella_pneumoniae_cgMLST_count.tsv"

    LocusFile = pd.read_csv('%s' % (alleleTable), sep='\t')
    LocusNumber = LocusFile['Allele Count']
    LocusNames = LocusFile['Locus']
    # LocusLength = LocusFile['Length']
    # 1.得到每个reads在每个等位基因上的得分,同时记录最大得分，比对上的位点

    R_test = []

    R = {}
    Temp = {}
    #R[readsName] = [等位基因,得分],记录每个locus 比对上的reads的最大得分和对应的等位基因
    L = {}
    maxLocusAS = 0
    maxAS = {}
    query_locus = ""
    query_alleles = 0
    for r in bf:
        if r.has_tag("NM") and int(r.get_tag('AS'))/int(r.query_alignment_length) > 1:
            query_locus_now = r.query_name.split('_')[0]
            query_alleles_now = r.query_name.split('_')[1]
            query_AS_now = r.get_tag("AS")
            query_reference = r.reference_name

            # 位点未更新
            if query_locus_now == query_locus :
                # 记录当前reads在这个位点上的最大得分maxLocus
                if maxLocusAS < query_AS_now:
                    maxLocusAS = query_AS_now
                if query_reference in R:
                    R[query_reference][0].append((query_alleles_now,query_AS_now))
                    R[query_reference][1] += 1
                    if query_AS_now > R[query_reference][2]:
                        R[query_reference][2] = query_AS_now
                # reads没有出现过，加入R中
                else:
                    R[query_reference] = [[(query_alleles_now,query_AS_now)],1,query_AS_now]

            else:
                # 位点更新时
                if query_alleles != 0:
                    #对每个位点对应的等位基因匹配到的reads数目进行统计
                    # L[query_locus] = R
                    #统计每个位点比对得分最大值（在多条reads上的最大值）
                    maxAS[query_locus] = maxLocusAS
                    # 对R中的reads进行筛选
                    # 1）一个reads如果比对上的等位基因数目如果少于等位基因数*60%个，则认为reads不是肺克reads
                    # 2）reads在一个位点的最高得分小于位点所有reads最高得分的80%，则认为reads质量过低或非肺克
                    # 3）对每个reads ,保留前5名得分对应等位基因
                    for read in R:

                        if R[read][1] < 0.1 * LocusFile.loc[LocusFile['Locus'] == query_locus, 'Allele Count'].item() or \
                                R[read][1] < 5 or R[read][2] < 0.75 * int(maxLocusAS):
                            continue
                        # 对比对得分排序
                        x = sorted(R[read][0], key=lambda s: int(s[1]), reverse=True)
                        # 保留前五名
                        R_test.append(read)
                        Temp[read] = [x[0],x[1],x[2],x[3],x[4]]
                        # Temp[read] = [x[0], x[1], x[2]]
                    L[query_locus] = Temp
                    # rescaleAS(L,maxAS)
                    Temp = {}
                    R = {}
                    maxLocusAS = query_AS_now

                query_locus = query_locus_now
                query_alleles = query_alleles_now
                query_AS = query_AS_now
                R[query_reference] = [[(query_alleles_now, query_AS_now)], 1, query_AS_now]
                maxLocusAS = query_AS_now
    R_test = set(R_test)
    # 2.对得分进行缩放
    L = rescaleAS(L,maxAS)
    # 3.EM算法得到每个位点的中等位基因的概率
    maxIter = 200
    emEpsilon = 0.00001
    pi = Locus_EM(L,maxIter, emEpsilon,LocusNumber,LocusNames)
    # 4.计算每个基因组的概率
    bf.close()






    bf.close()
    return pi , L

def Reads_sam_ananly2(samfile):

    Kl_locus = "/data1/zhuxu_data/MLST_nanopore/Kl_locus"
    f = open('%s' % (Kl_locus), 'r')
    c = f.read()
    Kl_locus = eval(c)
    f.close()
    # print(Kl_locus)

    dict_locus = [0 for i in range(len(Kl_locus))]

    bf = pysam.AlignmentFile(samfile,'r')

    L = {}
    # A记录每个位点比对的情况，每个位点包含多个等位基因的比对
    A = {}
    #R[query_reference] = [query_alleles,query_AS],记录每个locus 比对上的reads的最大得分和对应的等位基因
    # 当前比对的位点记录，初始值为""
    query_locus = ""
    query_alleles = 0
    for r in bf:

        if r.has_tag("NM"):
            query_locus_now = r.query_name.split('_')[0]+'_'+r.query_name.split('_')[1]
            query_alleles_now = r.query_name.split('_')[-1]
            query_AS_now = r.get_tag("AS")

            if query_locus_now == query_locus :
                # 如果比对的位点没有更新
                if query_alleles_now not in A:
                    A[query_alleles_now] = query_AS_now

                else:
                    A[query_alleles_now] = A[query_alleles_now]+query_AS_now

            else:
                # 新的比对位点
                if query_locus != '':
                    # 不是第一个比对的位点
                    # 对位点的等位基因比对上的reads数目进行排序
                    x = sorted(A.items(),key=lambda s : int(s[1]), reverse=True)
                    # 取比对数目最多的等位基因做为该位点的等位基因
                    dict_locus[Kl_locus.index(query_locus)] = x[0][0]
                    A ={}
                query_locus = query_locus_now
                A[query_alleles_now] = query_AS_now
                # print(r.get_tag("NM"),r.is_supplementary,r.infer_read_length(),r.infer_query_length(),r.query_alignment_length,r.query_name,r.query_length,r.reference_name,r.reference_length,r.flag)
        else:
            print("----------")


    bf.close()
    return dict_locus

# 阴性样本判断
def Exist_Judge(p1, x, n):
    # x表示排名第一的菌株匹配的位点数
    # n表示总的位点数目
    # p1表示H1（目标物种存在）的条件下一个位点被匹配的概率
    # p0表示H0(目标物种不存在)的条件下一个位点匹配的概率
    p0 = 0.05
    l = 1
    G = math.pow(p0/p1, x) * math.pow((1-p0)/(1-p1), n-x)
    result = 1 / (1 + l * G)
    return result

def genomeZoroLocusCount(target_locus):
    # 统计目标菌株基因组位点为0的数目
    count = 0
    for i in target_locus:
        if i == 0:
            count = count + 1
    print("target locus 0 count:", count)


def getGCF_NC_Number(genomeRefPath):
    #由参考基因组获取GCF对应NC表
    GCF_NC = pd.read_csv(genomeRefPath,sep='\t',index_col="GCF NAME")
    GCF_NC = GCF_NC.to_dict(orient ='index')
    GCF_NC_dict = {}
    # print(GCF_NC)
    for key in GCF_NC:
        NC_list = set(GCF_NC[key].values())
        NC_list = list(NC_list)
        GCF_NC_dict[key] = NC_list
    # print(GCF_NC_dict)
    return GCF_NC_dict


def Zymo():

    # report文件目录
    report_path = "/data2/zhuxu_data/Nano_cgMLST"
    report_name = "ZymoReport4.tsv"
    coverage_list = ['01', '03', '05', '1', '3', '5']
    species_list = ["Ec", "Ef", "Sa", "Se", "Pa", "Lm"]
    # species_list = ["Lm"]
    # coverage_list = ["5"]
    genomeNameDict = {"Ec": "B1109", "Ef": "B537", "Sa": "B41012", "Se":"B4242", "Pa": "B3509", "Lm": "B33116"}
    genomeMCDict = {"Ec": 0.545, "Ef": 0.313, "Sa": 0.389, "Se":0.620, "Pa": 0.695, "Lm": 0.383}
    alleleTabelDict = {"Ec": "../library/Escherichia_coli_alleles_count.tsv", "Ef": "../library/Enterococcus_faecalis_alleles_count.tsv",
                       "Sa": "../library/Staphylococcus_aureus_alleles_count.tsv", "Se": "../library/Salmonella_enterica_alleles_count.tsv",
                       "Pa": "../library/Pseudomonas_aeruginosa_alleles_count.tsv", "Lm": "../library/Listeria_monocytogenes_alleles_count.tsv"}
    genomeSizeDict = {"Ec": 4875000, "Ef": 2845000, "Sa": 2730000, "Se":4760000, "Pa": 6792000, "Lm": 2992000}
    avgLength = 4398.3
    Count = {}

    if os.path.exists('%s/%s'%(report_path ,report_name)):
        report = pd.read_csv('%s/%s'%(report_path ,report_name),sep='\t')
    else:
        report = pd.DataFrame(columns=("sampleID","coverage", "best hit", "best hit score", "best hit distance", "second hit",
                                       "second hit score", "second hit distance", "third hit", "third hit score",
                                       "third hit distance", "fourth hit", "fourth hit score", "fourth hit distance",
                                       "fifth hit", "fifth hit score", "fifth hit distance","coverage estimation"))

    index = 288
    for i in range(9,11):
        for species in species_list:


            genomeName = genomeNameDict[species]
            for coverage in coverage_list:
                mapOptions = MapOptions()
                mapOptions.species = species
                mapOptions.alleleTable = alleleTabelDict[species]
                mapOptions.minimapPath = "/home/zhuxu/miniconda3/bin"
                mapOptions.samFile = "/data2/zhuxu_data/Nano_cgMLST/Zymo/%s/%s_%sX_%s.sam" % (species, species,coverage,i)
                mapOptions.fastaFile = "/data2/zhuxu_data/Nano_cgMLST/Zymo/%s/%s_%sX_%s.fasta" % (species, species,coverage,i)
                mapOptions.genomesDir = "/data2/zhuxu_data/Nano_cgMLST/genomes/%s"%species
                mapOptions.genomeLocusPath = "../library/%s.tsv"%species
                mapOptions.genomeDistancePath = "../library/%s_Locus_Distance.tsv"%species
                mapOptions.outPutDir = "/data2/zhuxu_data/Nano_cgMLST/Zymo/result"

                genomeLocus = pd.read_csv("%s"%mapOptions.genomeLocusPath,sep="\t")
                target_locus = genomeLocus[genomeName]
                pi , L = Reads_sam_analy3(mapOptions)
                result,count = selectReadsMiminmap(pi,L,mapOptions)
                sample = "%s_%s_%s" % (mapOptions.species,coverage,i)
                print(sample)
                Count[sample] = count
                sampleResult = [sample]
                sampleResult.append(str(coverage)+"X")
                c = 1
                for gcf, row in result.iterrows():
                    if c <=5 :
                        sampleResult.append(gcf)
                        sampleResult.append(row["genome score"])
                        sampleResult.append(row["target and genome distance"])
                    c = c + 1
                # MC estimate coverage
                estCov = count /genomeMCDict[species] * avgLength /genomeSizeDict[species]
                sampleResult.append(estCov)
                # 结果记录
                report.loc[index] = sampleResult
                index += 1
                print(index)
                report.to_csv('%s/%s' % (report_path, report_name), sep='\t', index=False)
                os.remove(mapOptions.samFile)

Zymo()