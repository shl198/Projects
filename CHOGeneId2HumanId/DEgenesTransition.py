import pandas as pd
import os,sys
import ipdb
from natsort import natsorted
import subprocess


#source = "/data/shangzhong/DE/GlycoKnock"
source = sys.argv[1]
os.chdir(source)
deseq_files = [f for f in os.listdir(source) if f.endswith('csv')]
mappingFile = '/data/shangzhong/CHO2Mouse/finalMergeWithmRNA.final.txt'
res = open(mappingFile,'r')
dic = {}
for line in res:
    item = line[:-1].split('\t')
    if ';' in item[1]:
        split = item[1].split(';')
    else:
        split = item[1].split(',')
    dic[item[0]] = split
res.close()

for f in deseq_files:
    outFile = f[:-4] + '_mouse.txt'
    output = open(outFile,'w')
    res = open(f,'r')
    firstline = res.readline()[:-1].split(',')
    firstline[0] = 'ID'
    for i in range(1,len(firstline)):
        firstline[i] = firstline[i].strip("'" '"')
    output.write('\t'.join(firstline)+'\n')
    for line in res:
        item = line[:-1].split(',')
        item[0] = item[0].strip("'" '"')
        try:
            mapids = dic[str(item[0])]
        except:
            print item[0],'is not in database'
            continue
        for gene in mapids:
            output.write(gene + '\t' + '\t'.join(item[1:]) + '\n')
    output.close()
    # remove duplicates
    df = pd.read_csv(outFile,sep='\t',header=0)
    genes = df['ID'].tolist()
    res_df = pd.DataFrame()
    for gene in genes:
        gene_df = df[df['ID'].values==gene]
        shape = gene_df.shape
        if shape[0] == 1:
            res_df = pd.concat([res_df,gene_df])
        else:
            index = shape[0]/2
            res_df = pd.concat([res_df,gene_df.iloc[[index],:]])
    finalresult = outFile[:-3] + 'uniq.xlsx'    
    res_df.to_excel(finalresult,index=False)
    
    
