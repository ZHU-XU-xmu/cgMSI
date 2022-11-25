import pysam
import math
import pandas as pd
import os
from scipy import stats
import time


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
    # 物种名
    species = ""
    # 比对结果文件
    samFile = ""
    # fasta或fastq文件
    sampleFile = ""
    # alleleTable文件
    alleleTablePath = "../library/Escherichia_coli_alleles_count.tsv"
    # 基因组分型结果文件
    genomeLocusPath = ""
    # 基因组之间的距离
    genomeDistancePath = ""
    # 基因组目录文件
    genomesDir = ""
    # 结果输出目录
    outPutDir = ""
    # 比对线程数目
    threads = 10
    # library目录
    LocusDir = "../library"


def loadLocus(path):
    # 得到每个位点的名称，返回一个list
    Locus = pd.read_csv(path,sep='\t')
    speciesLocus = Locus["Locus"].tolist()
    return speciesLocus


def selectReadsMiminmap(L, MapOptions):
    # 筛选的reads与筛选后的基因组比对
    genomeLocusPath = MapOptions.genomeLocusPath
    outPutDir = MapOptions.outPutDir
    fastaFile = MapOptions.sampleFile

    genomesDir = MapOptions.genomesDir
    species = MapOptions.species
    alleleTable = MapOptions.alleleTablePath
    genomeDistancePath = MapOptions.genomeDistancePath
    # 基因组分型距离文件

    select_reads_name = fastaFile.split("/")[-1].split(".")[0] + "_selected.txt"
    select_reads_file = fastaFile.split("/")[-1].split(".")[0] + "_selected" + ".fasta"
    select_ref_file = fastaFile.split("/")[-1].split(".")[0] + "_selected" + "_ref.fna"
    selectMapResultFile = fastaFile.split("/")[-1].split(".")[0] + "_selected" + ".sam"

    # genomeDistance = pd.read_csv('%s' % (genomeDistancePath), index_col='GCF name', sep='\t')
    LocusFile = pd.read_csv('%s' % (alleleTable), sep='\t')
    LocusNumber = LocusFile['Allele Count']
    genomeLocus = pd.read_csv('%s' % (genomeLocusPath), sep='\t')
    genomeName = genomeLocus.columns.tolist()
    LocusName = genomeLocus[genomeName[0]]

    Count = {}
    Sum = []
    loucusCount = len(LocusNumber)
    for i in range(loucusCount):
        locuNumber = LocusNumber[i]
        Sum.append([0 for _ in range(locuNumber)])

    # 样本和参考基因组一一比较
    for key in L:
        for read in L[key]:
            # 位点中每个reads
            Tu = L[key][read]
            for t in range(len(Tu)):
                index = int(LocusName[LocusName.values == key].index[0])
                Sum[index][int(Tu[t][0]) - 1] += Tu[t][1]
    predictLocusNum = {}
    for genomename in genomeName[1:]:
        count = 0
        number = 0
        ref_locus = genomeLocus[genomename].tolist()
        for j in range(len(ref_locus)):
            # 去除等位基因为0的情况
            if int(ref_locus[j]) != 0:
                count += Sum[j][int(ref_locus[j]) - 1]
                if Sum[j][int(ref_locus[j]) - 1]!=0:
                    number = number + 1
        predictLocusNum[genomename] = number

        Count[genomename] = [count]

    predictLocusNum = pd.DataFrame.from_dict(predictLocusNum, orient='index',
                                    columns=['genome locus num'])
    predictLocusNum = predictLocusNum.sort_values(by=['genome locus num'], ascending=False)

    # maxGenome = predictLocusNum.index[0]
    maxgenomeLocus = predictLocusNum["genome locus num"][0]

    result = pd.DataFrame.from_dict(Count, orient='index',
                                    columns=['genome score'])
    result = result.sort_values(by=['genome score'], ascending=False)
    #选择排名靠前的十个基因组，第十一名后和第十名相同得分，则也确认为候选菌株
    fifthScore = result['genome score'][9]
    cut = 10
    while(result['genome score'][cut]==fifthScore and cut <20 ):
        cut += 1
    Ref = result.index.tolist()[:cut]
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
            if int(ref_locus[k]) != 0 and Sum[k][int(ref_locus[k]) - 1] != 0 and len(L[key]):
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
        os.remove("%s/%s" % (outPutDir, select_reads_file))
    with open("%s/%s" % (outPutDir, select_reads_name), "a+") as f:
        for key in R:
                newline = key + "\n"
                f.write(newline)

    # 要去掉fasta文件中的注释才能筛出文件来
    fastaFileRename = fastaFile.split(".")[0]+"_rename.fasta"
    cmd ="seqkit replace -p \"\s.+\"  %s > %s" %(fastaFile,fastaFileRename)
    os.system(cmd)
    cmd = "seqkit grep -n -f %s/%s %s -o %s/%s" % (outPutDir , select_reads_name, fastaFileRename, outPutDir , select_reads_file)
    os.system(cmd)

    #筛选出的reads和候选菌株比对
    if os.path.exists("%s/%s" % (outPutDir,selectMapResultFile)):
        os.remove("%s/%s" % (outPutDir,selectMapResultFile))
    cmd = "minimap2 -ax map-ont -t 10 -c -K 2G --sam-hit-only  -N 10 %s/%s %s/%s > %s/%s" % \
    (outPutDir, select_ref_file,outPutDir ,select_reads_file,outPutDir,selectMapResultFile)
    os.system(cmd)
    samplePath = outPutDir+'/'+selectMapResultFile
    mapScore,count = selectedSamAnaly(samplePath,Ref)
    Score = sorted(mapScore.items(),key=lambda x: x[1],reverse=True)


    #将分数进行缩放和归一化
    MapScaleFactor = 100/Score[0][1]
    MlstSum = 0
    MapSum = 0
    for k in Ref:
        mapScore[k] = math.exp( mapScore[k] * MapScaleFactor )
        MapSum = MapSum + mapScore[k]
        MlstSum = MlstSum + Count[k][0]


    result2 = {}  # 作后验时
    for k in Ref:
        mapScore[k] = mapScore[k] / MapSum
        Count[k][0] = Count[k][0] / MlstSum
        result2[k] = [(Count[k][0] * mapScore[k]), Count[k][0]]

    result2 = pd.DataFrame.from_dict(result2, orient='index',
                                     columns=['genome score', 'target and genome distance'])
    result2 = result2.sort_values(by=['genome score'], ascending=False)

    return result2,maxgenomeLocus/locuNumber,count


def selectedSamAnaly(samPath,Ref):
    #分析选择的reads比对后的结果
    # genomeRefPath = "../library/%s_NC_GCF.tsv"%species
    genomeRefPath = "/data2/zhuxu_data/cgMLST/NC_GCF.tsv"
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
                reads_name = r.query_name
                ref_name = r.reference_name
                if r.query_name not in R:
                    R.append(r.query_name)
                    l.append(r.query_length)
                    count += 1

                if ref_name in G:
                    #一个基因组上一次reads 可能有多个比对结果（基因的复制）
                    if reads_name not in G[ref_name][0]:
                        G[ref_name][0].append(reads_name)
                        G[ref_name][1].append(r.get_tag('AS'))
                        G[ref_name][2].append(r.query_alignment_length)
                        G[ref_name][3] += int(r.get_tag('AS'))
                else:
                    G[ref_name] = [[reads_name], [r.get_tag('AS')], [r.query_length], int(r.get_tag('AS'))]
    bf.close()

    #GCF和NC编号可能是一对多的关系，将同一GCF下的NC得分累加
    for key in Ref:
        score[key] = 0
        for NC in G:
            if NC in GCF_NC_dict[key]:
                score[key] += G[NC][3]

    return score, count



def Reads_sam_analy(MapOptions):
    e=0.1

    alleleTable = MapOptions.alleleTablePath
    samFile = MapOptions.samFile
    bf = pysam.AlignmentFile(samFile, 'r')
    LocusFile = pd.read_csv('%s' % (alleleTable), sep='\t')
    # 1.得到每个reads在每个等位基因上的得分,同时记录最大得分，比对上的位点
    R = {}
    Temp = {}
    L = {}
    minLocusNM = 999999999
    maxNM = {}
    query_locus = ""
    query_alleles = 0

    for r in bf:
        if r.has_tag("NM"):
            # and int(r.get_tag('AS')) / int(r.query_alignment_length) > 1
            query_locus_now = r.query_name.split('_')[0] + '_' + r.query_name.split('_')[1]
            query_alleles_now = r.query_name.split('_')[2]
            query_NM_now = r.get_tag("NM")
            # print(query_NM_now)
            query_reference = r.reference_name
            query_length_now = r.query_length

            # 位点未更新
            if query_locus_now == query_locus:
                # 记录当前reads在这个位点上的最大得分maxLocus
                if minLocusNM > query_NM_now:
                    minLocusNM = query_NM_now
                if query_reference in R:
                    R[query_reference][0].append((query_alleles_now, query_NM_now, query_length))
                    R[query_reference][1] += 1
                    if query_NM_now < R[query_reference][2]:
                        R[query_reference][2] = query_NM_now
                # reads没有出现过，加入R中
                else:
                    # Pe = Comb_Pro(query_length_now, query_NM_now, e)
                    # print(Pe)
                    R[query_reference] = [[(query_alleles_now, query_NM_now, query_length)], 1, query_NM_now]

            else:
                # 位点更新时
                if query_alleles != 0:
                    #对每个位点对应的等位基因匹配到的reads数目进行统计
                    # L[query_locus] = R
                    #统计每个位点比对得分最大值（在多条reads上的最大值）
                    maxNM[query_locus] = minLocusNM
                    # 对R中的reads进行筛选
                    # 1）一个reads如果比对上的等位基因数目如果少于等位基因数*60%个，则认为reads不是肺克reads
                    # 2）reads在一个位点的最高得分小于位点所有reads最高得分的80%，则认为reads质量过低或非肺克
                    # 3）对每个reads ,保留前5名得分对应等位基因
                    for read in R:
                        if R[read][1] < 0.1 * LocusFile.loc[LocusFile['Locus'] == query_locus, 'Allele Count'].item() or R[read][1] < 10 or R[read][2] > 1.25 * int(minLocusNM):
                            continue
                        # 对比对得分排序
                        x = sorted(R[read][0], key=lambda s: int(s[1]), reverse=False)
                        # 保留前五名
                        temp = []
                        for i in range(5):
                            temp.append((x[i][0],Comb_Pro(int(x[i][2]),int(x[i][1]),e),x[i][1]))
                        Temp[read] = temp

                    L[query_locus] = Temp
                    Temp = {}
                    R = {}
                    minLocusNM = query_NM_now

                query_locus = query_locus_now
                query_alleles = query_alleles_now
                query_NM = query_NM_now
                R[query_reference] = [[(query_alleles_now, query_NM_now,query_length_now)], 1, query_NM_now]
                minLocusNM = query_NM_now
                query_length = query_length_now
    bf.close()

    return L



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




if __name__ == "__main__":
    ANI_test8()
    # ANI_test10()
    # ANI_test13()
    # ANI_test14()
    time_performance()


    # decoy_test1()
    # ANI_Distance("GCF_018449065","GCF_000364385")
