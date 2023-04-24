import pandas as pd

def rename_Locus(path):

    Locus = pd.read_csv(path,sep='\t')
    # print(Locus["Locus"])
    LocusName = Locus["Locus"].tolist()
    LocusName_uni = []
    for locus in LocusName:

        LocusName_uni.append(locus.replace("_",""))
    Locus["Locus"] = LocusName_uni
    # print(Locus["Locus"])
    Locus.to_csv(path,sep='\t', index=None)
path = "/data2/zhuxu_data/cgMSI/library/Enterococcus_faecalis_alleles_count.tsv"
rename_Locus(path)
