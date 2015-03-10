import subprocess,os

def extract_geneID(outputFile,gff3):
    """
    This function extracts gene ids from gff3 annotation file
    """
    res = open(gff3,'r')
    gene_list = []
    output = open(outputFile,'w')
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

gff3 = '/opt/genome/cho/chok1.gff3'

def extract_gff4cufflink_symbol2ID(outputFile,gffFile):
    """
    This function extracts gene id,gene symbol,chromosome and stores these inoformation
    into a file with 3 colums. 
    
    * outputFile: outputFile name
    * gffFile: gffFile name
    return a file with 3 columns: [id,symbol,chromosome]
    """
    handle = open(gffFile,'r')
    f = open('inter.txt','w')
    for line in handle:
        item = line.split('\t')
        chm = item[0]
        des = item[-1]
        # get gene id
        try:
            index = des.index('GeneID:')
            index = index + 7
            ids = ''
            while des[index] !=';' and des[index] != ',':
                    ids = ids + des[index]
                    index = index + 1
        except:
            continue
        # get gene symbol
        try:
            index = des.index('gene=')
            index = index + 5
            symbol = ''
            while des[index] !=';' and des[index] != ',' and des[index]!='\n':
                symbol = symbol + des[index]
                index = index + 1
#             inter = des[index+5:]
#             index = inter.index(';')
#             symbol = inter[:index]
        except:
            continue
        f.write(ids+'\t'+symbol+'\t'+chm+'\n')
    f.close()
    cmd = ('sort inter.txt | uniq > {output}').format(output=outputFile)
    subprocess.call(cmd,shell=True)
    os.remove('inter.txt')
#extract_gff4cufflink_symbol2ID('/data/shangzhong/test.txt',gff3)          