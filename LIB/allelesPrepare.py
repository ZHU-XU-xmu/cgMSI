import os



def rename_merge_file(alllelesDir,species):

    loci = os.listdir(alllelesDir)
    if os.path.exists("%s/%s_alleles.fasta" % (os.path.dirname(alllelesDir), species)):
        os.remove("%s/%s_alleles.fasta" % (os.path.dirname(alllelesDir), species))

    for locus in loci:
        absPath = os.path.join(alllelesDir, locus)
        cmd = "cat %s|seqkit replace -p ^ -r \"${filename//_/}_\"  >> %s/%s_alleles.fasta" % (absPath, os.path.dirname(alllelesDir), species)
        os.system(cmd)
