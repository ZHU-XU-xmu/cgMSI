import pysam
import pandas as pd
import random
import os
import csv

def readsCount(sampleFile,GenomeLength,alleleTable):

    GenomeLength = int(GenomeLength)
    #1. generate genome gene section
    LocusFile = pd.read_csv('%s' % (alleleTable), sep='\t')
    LocusLength = LocusFile['Length']
    G = []
    preEnd = 0
    LocusStart = LocusFile['Start']
    for s in range(len(LocusStart)):
        # add gap section
        temp0 = [preEnd,LocusStart[s]-1,0]
        G.append(temp0)
        #add gene section
        locusEnd = LocusStart[s] + LocusLength[s]
        temp1 = [LocusStart[s],locusEnd,1]
        G.append(temp1)
        preEnd = locusEnd + 1
    # print(G)
    G.append([preEnd,GenomeLength,0])
    # 2. get sample length distribution
    bf = pysam.Fastafile(sampleFile)
    lengths = bf.lengths
    readsNum = len(lengths)
    countT = 0
    successNum = 0
    while(countT<1000000):

        # 3. generate random location in the genome
        # r = L[t,l] s:random reads start, e : reads end
        s = random.randint(0,GenomeLength)
        e = s + lengths[countT%readsNum]
        if e > GenomeLength:
            continue
        countT += 1

        # 4. binary searchs r.s,r.e section in G
        flag = 0
        left = 0
        right = len(G)
        while not flag :
            mid = int(left + (right - left)/2)
            if s > G[mid][1]:
                left = mid + 1
            elif s < G[mid][0]:
                right = mid -1
            else:
                flag = 1
        sSection = mid
        eSection = sSection
        while e > G[eSection][0]:
            if e > G[eSection][1]:
                eSection += 1
            else:
                break

        # 5. count complete gene section
        if eSection - sSection > 2:
            successNum += 1
        elif eSection - sSection == 2:
            sflag = G[sSection][2]
            if not sflag:
                successNum += 1
    print(successNum/countT)
    return successNum/countT

'''
sampleFile = "/data2/zhuxu_data/Zymo/Ec/Ec_5X_1.fasta"
species = "Ec"
GenomeLength = "4875000"
alleleTable = "/data2/zhuxu_data/cgMSI/library/Escherichia_coli_alleles_count.tsv"
readsLengthCount(sampleFile,GenomeLength,alleleTable)

sampleFile = "/data2/zhuxu_data/Zymo/Ef/Ef_5X_1.fasta"
species = "Ef"
GenomeLength = "2845000"
alleleTable = "/data2/zhuxu_data/cgMSI/library/Enterococcus_faecalis_alleles_count.tsv"
readsLengthCount(sampleFile,GenomeLength,alleleTable)


sampleFile = "/data2/zhuxu_data/Zymo/Lm/Lm_5X_1.fasta"
species = "Lm"
GenomeLength = "2992000"
alleleTable = "/data2/zhuxu_data/cgMSI/library/Listeria_monocytogenes_alleles_count.tsv"
readsLengthCount(sampleFile,GenomeLength,alleleTable)


sampleFile = "/data2/zhuxu_data/Zymo/Sa/Sa_5X_1.fasta"
species = "Sa"
GenomeLength = "2730000"
alleleTable = "/data2/zhuxu_data/cgMSI/library/Staphylococcus_aureus_alleles_count.tsv"
readsLengthCount(sampleFile,GenomeLength,alleleTable)

sampleFile = "/data2/zhuxu_data/Zymo/Se/Se_5X_1.fasta"
species = "Se"
GenomeLength = "4760000"
alleleTable = "/data2/zhuxu_data/cgMSI/library/Salmonella_enterica_alleles_count.tsv"
readsLengthCount(sampleFile,GenomeLength,alleleTable)

sampleFile = "/data2/zhuxu_data/Zymo/Pa/Pa_5X_1.fasta"
species = "Pa"
GenomeLength = "6792000"
alleleTable = "/data2/zhuxu_data/cgMSI/library/Pseudomonas_aeruginosa_alleles_count.tsv"
readsLengthCount(sampleFile,GenomeLength,alleleTable)
'''


