def extract_geneID(gff3):
    """
    This function extracts gene ids from gff3 annotation file
    """
    res = open(gff3,'r')
    gene_list = []
    geneFile = '/data/shangzhong/CHO2Human/geneID.txt'
    output = open(geneFile,'w')
    for line in res:
        try:
            index = line.index('GeneID:')
            index = index + 7
        except:
            continue
        geneid = ''
        while line[index] !=';' and line[index] != ',':
                geneid = geneid + line[index]
                index = index + 1
                
        gene_list.append(geneid)
    gene_list = list(set(gene_list))    
    for item in gene_list:
        output.write(item + '\n')
    output.close()

gff3 = '/home/shangzhong/chok1.gff3'
extract_geneID(gff3)