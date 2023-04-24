import pysam
import math
import pandas as pd
import os
from scipy import stats



class MapOptions:
    # species name
    species = ""
    # sam file path
    samFile = ""
    # fasta or fastq file
    sampleFile = ""
    # alleleTable file
    alleleTablePath = ""
    # genome Locus path
    genomeLocusPath = ""
    # reference genomes dir
    genomesDir = ""
    # out put dir
    outPutDir = ""
    # thread
    threads = 10
    # library dir
    LocusDir = "../library"
    # pathogen existing Judge parameter
    e = 0.008


def loadLocus(path):
    # Get the name of each locus，return a list
    Locus = pd.read_csv(path,sep='\t')
    speciesLocus = Locus["Locus"].tolist()
    return speciesLocus

def genomeZoroLocusCount(target_locus):
    # The number of zero genomic loci of the target strain was counted.
    count = 0
    for i in target_locus:
        if i == 0:
            count = count + 1
    print("target locus 0 count:", count)


def getGCF_NC_Number(genomeRefPath):
    #get the relationship of genome GCF num and NC number（得到基因组的GCF编号与NC编号的关系）
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


def rescaleAS(L,maxAS):
    for key in L:
        maxScore = maxAS[key]
        if(maxScore==0):
            print(key)
        scalingFactor = 100.0 / (maxScore)
        for read in L[key]:
            Tu = L[key][read]
            for i in range(len(Tu)):
                L[key][read][i] = (Tu[i][0],math.exp(Tu[i][1]*scalingFactor))
    return L


def Reads_sam_analy(MapOptions):
    e=0.1

    alleleTable = MapOptions.alleleTablePath
    samFile = MapOptions.samFile
    bf = pysam.AlignmentFile(samFile, 'r')
    LocusFile = pd.read_csv('%s' % (alleleTable), sep='\t')

    # Get the score of each reads on each allele, and at the same time record the maximum score, the loci mapped (得到每个reads在每个等位基因上的得分,同时记录最大得分，比对上的位点)
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

            # locus is not update
            if query_locus_now == query_locus:
                if minLocusNM > query_NM_now:
                    minLocusNM = query_NM_now
                if query_reference in R:
                    R[query_reference][0].append((query_alleles_now, query_NM_now, query_length))
                    R[query_reference][1] += 1
                    if query_NM_now < R[query_reference][2]:
                        R[query_reference][2] = query_NM_now
                #  a new read
                else:
                    R[query_reference] = [[(query_alleles_now, query_NM_now, query_length)], 1, query_NM_now]

            else:
                # a new locus
                if query_alleles != 0:
                    maxNM[query_locus] = minLocusNM
                    for read in R:
                        if R[read][1] < 0.2 * LocusFile.loc[LocusFile['Locus'] == query_locus, 'Allele Count'].item() or R[read][1] < 10 or R[read][2] > 1.25 * int(minLocusNM):
                            continue
                        # sort by mapped score
                        x = sorted(R[read][0], key=lambda s: int(s[1]), reverse=False)
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


def selectReadsMiminmap(L, MapOptions):
    # Screened reads vs. screened genomes
    genomeLocusPath = MapOptions.genomeLocusPath
    outPutDir = MapOptions.outPutDir
    fastaFile = MapOptions.sampleFile
    genomesDir = MapOptions.genomesDir
    species = MapOptions.species
    alleleTable = MapOptions.alleleTablePath
    e = MapOptions.e

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

    # Counting the results for each genome
    for key in L:
        for read in L[key]:
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
            if int(ref_locus[j]) != 0:
                count += Sum[j][int(ref_locus[j]) - 1]
                if Sum[j][int(ref_locus[j]) - 1]!=0:
                    number = number + 1
        predictLocusNum[genomename] = number

        Count[genomename] = [count]
    # print(predictLocusNum)

    # Count the maximum number of matched loci
    predictLocusNum = pd.DataFrame.from_dict(predictLocusNum, orient='index',
                                    columns=['genome locus num'])
    predictLocusNum = predictLocusNum.sort_values(by=['genome locus num'], ascending=False)

    # maxGenome = predictLocusNum.index[0]
    maxgenomeLocus = predictLocusNum["genome locus num"][0]

    if not existJudge(maxgenomeLocus,e,len(LocusFile)):
        return False,0,0,0

    result = pd.DataFrame.from_dict(Count, orient='index',
                                    columns=['genome score'])
    result = result.sort_values(by=['genome score'], ascending=False)
    # Select the top K genomes, and the same score after the K will also be selected
    fifthScore = result['genome score'][9]
    cut = 10
    while(result['genome score'][cut]==fifthScore and cut <20 ):
        cut += 1
    Ref = result.index.tolist()[:cut]
    # Selected genome into a fasta file
    if os.path.exists("%s/%s"%(outPutDir, select_ref_file)):
        os.remove("%s/%s"%(outPutDir, select_ref_file))

    for i in range(len(Ref)):
        refname = Ref[i]
        cmd = "cat %s/%s*.fna >> %s/%s" % (genomesDir, refname, outPutDir, select_ref_file)
        os.system(cmd)

    # Selected reads are written to a fasta file
    R = []
    for i in range(len(Ref)):
        k = 0
        ref_locus = genomeLocus[Ref[i]].tolist()
        for key in LocusName:
            if key not in L:
                k += 1
                continue
            if int(ref_locus[k]) != 0 and Sum[k][int(ref_locus[k]) - 1] != 0 and len(L[key]):
                for r in L[key].keys():
                    for i in range(len(L[key][r])):
                        if int(L[key][r][i][0]) == int(ref_locus[k]):
                            R.append(r)
                            break
            k += 1
    Rset = set(R)
    R = list(Rset)
    #get th fasta file of selected reads
    if os.path.exists("%s/%s" % (outPutDir, select_reads_name)):
        os.remove("%s/%s" % (outPutDir, select_reads_name))
    if os.path.exists("%s/%s" % (outPutDir, select_reads_file)):
        os.remove("%s/%s" % (outPutDir, select_reads_file))
    with open("%s/%s" % (outPutDir, select_reads_name), "a+") as f:
        for key in R:
                newline = key + "\n"
                f.write(newline)

    # Remove comments from fasta files
    fastaFileRename = fastaFile.split(".")[0]+"_rename.fasta"
    cmd ="seqkit replace -p \"\s.+\"  %s > %s" %(fastaFile,fastaFileRename)
    os.system(cmd)
    cmd = "seqkit grep -n -f %s/%s %s -o %s/%s" % (outPutDir , select_reads_name, fastaFileRename, outPutDir , select_reads_file)
    os.system(cmd)

    #selected reads mapped to candidate strains
    if os.path.exists("%s/%s" % (outPutDir,selectMapResultFile)):
        os.remove("%s/%s" % (outPutDir,selectMapResultFile))
    cmd = "minimap2 -ax map-ont -t 10 -c -K 2G --sam-hit-only  -N 10 %s/%s %s/%s > %s/%s" %  (outPutDir, select_ref_file,outPutDir,select_reads_file,outPutDir,selectMapResultFile)
    os.system(cmd)

    samplePath = outPutDir+'/'+selectMapResultFile
    mapScore,count = selectedSamAnaly(species,samplePath,Ref)
    Score = sorted(mapScore.items(),key=lambda x: x[1],reverse=True)

    MapScaleFactor = 100/Score[0][1]
    MlstSum = 0
    MapSum = 0
    for k in Ref:
        mapScore[k] = math.exp( mapScore[k] * MapScaleFactor )
        MapSum = MapSum + mapScore[k]
        MlstSum = MlstSum + Count[k][0]

    result = {}
    for k in Ref:
        mapScore[k] = mapScore[k] / MapSum
        Count[k][0] = Count[k][0] / MlstSum
        result[k] = [(Count[k][0] * mapScore[k]), Count[k][0]]

    result = pd.DataFrame.from_dict(result, orient='index',
                                     columns=['genome score', 'target and genome distance'])
    result = result.sort_values(by=['genome score'], ascending=False)

    if os.path.exists("%s/%s"%(outPutDir , select_reads_file)):
        os.remove("%s/%s"%(outPutDir , select_reads_file))
    if os.path.exists("%s/%s" % (outPutDir, selectMapResultFile)):
        os.remove("%s/%s" % (outPutDir, selectMapResultFile))
    if os.path.exists("%s/%s" % (outPutDir, select_reads_name)):
        os.remove("%s/%s" % (outPutDir, select_reads_name))
    if os.path.exists("%s/%s" % (outPutDir, select_ref_file)):
        os.remove("%s/%s" % (outPutDir, select_ref_file))
    if os.path.exists(fastaFileRename):
        os.remove(fastaFileRename)


    return True,result,count,maxgenomeLocus


def selectedSamAnaly(species,samPath,Ref):
    #Analyze the results after mapping the selected reads to candidate strains

    genomeRefPath = "./library/%s_NC_GCF.tsv"%species
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
                    count += 1
                    l.append(r.query_length)

                if ref_name in G:
                    if reads_name not in G[ref_name][0]:
                        G[ref_name][0].append(reads_name)
                        G[ref_name][1].append(r.get_tag('AS'))
                        G[ref_name][2].append(r.query_alignment_length)
                        G[ref_name][3] += int(r.get_tag('AS'))
                else:
                    G[ref_name] = [[reads_name], [r.get_tag('AS')], [r.query_length], int(r.get_tag('AS'))]
    bf.close()

    # A GCF numbet may correspond to many NC number
    for key in Ref:
        score[key] = 0
        for NC in G:
            if NC in GCF_NC_dict[key]:
                score[key] += G[NC][3]

    return score, count


def Comb_Pro(n, k, e):
    probability = stats.binom(n,e)
    return probability.pmf(k)

# Judging the presence of target pathogen species
def existJudge(maxGenomeLocus, e,locusNumber):
    if e*locusNumber <= maxGenomeLocus:
        return True
    else:
        return False


if __name__ == "__main__":
    print("start analyse")
