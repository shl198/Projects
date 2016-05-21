import pandas as pd
import os
import ipdb
import numpy
from Bio import SeqIO
import subprocess
from Modules.f00_Message import Message
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib as mpl
from multiprocessing import Process
mpl.style.use('ggplot')
import pysam
from Bio.Seq import Seq

#===============================================================================
#      merge cufflink and htseq for pgsa
#===============================================================================
# from Modules.f05_IDConvert import addProduct2CufflinkResultWithNCBIAnnotation
# htseq = '/data/shangzhong/DE/pgsa/htseq/htseq.name.csv'
# cufflink = '/data/shangzhong/DE/pgsa/cufflinks/sample_FPKM.csv.txt.csv'
# 
# # remove duplicate gene id + gene symbol in cufflinks results
# cuff_df = pd.read_csv(cufflink,sep='\t',header=0)
# cuff_sort = cuff_df.sort(['gene_short_name','Entrez_GeneID'])
# cuff_rmdup = cuff_sort.drop_duplicates(['gene_short_name','Entrez_GeneID'])
# cuff_res = cuff_sort.groupby(['gene_short_name','Entrez_GeneID'],as_index=False).sum()
# cuff_res.to_csv('/data/shangzhong/inter.csv',index=False,sep='\t')
# cuff_rmdupFile = addProduct2CufflinkResultWithNCBIAnnotation('/data/shangzhong/Database/141028chok1gene_info.txt','/data/shangzhong/inter.csv')
# cuffres_df = pd.read_csv(cuff_rmdupFile,sep='\t',header=0)
# # merge cufflinks and htseq
# htseq = pd.read_csv(htseq,header=0,sep='\t')
# htseq['GeneID'] = htseq['GeneID'].astype(str)
# htseq = htseq.rename(columns={'GeneID':'Entrez_GeneID','GeneSymbol':'gene_short_name'})
# res = pd.merge(cuffres_df,htseq,how='outer',on=['gene_short_name','Entrez_GeneID'])
# res.to_csv('/data/shangzhong/pgsa.csv',index=False,sep='\t')


# change name
# path = '/data/shangzhong/DE/fpkm/change_name'
# os.chdir(path)
# folders = [f for f in os.listdir(path) if os.path.isdir(f)]
# folders = natsorted(folders)
# for f in folders:
#     item =  f.split('_')
#     new_name = '_'.join(['base',item[1],item[2]])
#     os.rename(f,new_name)



# import pybedtools as pybed
# cdsFile = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam/db/01_pr_cds.txt'
# df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
# df.to_csv('inter.bed',sep='\t',header=None,index=False)
# a = pybed.BedTool('inter.bed')
# c = a.merge(c='4,5,6',o='distinct,distinct,distinct')
# os.remove('inter.bed')
# merge_df = pd.read_csv(c.fn,sep='\t',header=None,names=['chr','start','end','geneid','praccess','strand'])
# cri = merge_df['strand'].map(lambda x: ',' in x)
# overlapped_region = merge_df[cri] 
# print 'done'



#===============================================================================
#             get coverge plot for a single gene
#===============================================================================
def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]
"""
fn = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam/04_gene_total_count/s05_cov_geneCount.txt'
df = pd.read_csv(fn,sep='\t',header=0)
df = df.sort(['count'])
gene = df['GeneID'].tolist()
gene.reverse()
pr_fn = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam/03_pr_pos_cov/s05_cov_prpos.txt'
tr_fn = '/data/shangzhong/RibosomeProfiling/TotalRNA_align/03_tr_pos_cov/s02_cov_trpos.txt'
handle1 = open(pr_fn,'r')
handle2 = open(tr_fn,'r')
g = gene[24]  # 7
color = 'black'
print g
for line in handle1:
    item = line[:-1].split('\t')
    if item[0] == g:
        cov = item[2:]
        cov = [int(p) for p in cov]
        cover = cov
        n = 19#len(cover)
        plt.figure()
        plt.bar(range(n),cover[-19:],width=1.0,facecolor=color, edgecolor=color)
        plt.ylim(0,1000)
        plt.xlabel(g+'_pr')
        break
for line in handle2:
    item = line[:-1].split('\t')
    if item[0] == g:
        cov = item[3:]
        cov = [int(p) for p in cov]
        cover=cov
        cover = chunks(cov,3)
        cover = [sum(p) for p in cover]
        n = len(cover)
        plt.figure()
        plt.bar(range(312),cover[-312:],width=1.0,facecolor=color, edgecolor=color)
        plt.ylim(0,450)
        plt.xlabel(g+'_tr')
        break
plt.show()
"""

# from Modules.f03_samtools import *
# path = '/data/shangzhong/DetectVirus/Non_Autism'
# os.chdir(path)
# bamfiles = [f for f in os.listdir(path) if f.endswith('.bam')] 
# merge_bam(bamfiles,'/data/shangzhong/DetectVirus/virus.bam')
# import shutil
# 
# 


# fn = '/data/shangzhong/Pacbio/fa/CHOS/CHOS_pacbio.fa'
# length = []
# for record in SeqIO.parse(open(fn,'r'),'fasta'):
#     n = len(record.seq)
#     length.append(n)
# 
# 
# from collections import Counter
# c = Counter(length)
# print c
# handle = open('/data/shangzhong/Pacbio/Results/01_reads_length.txt','w')
# handle.write('length'+'\t'+ 'count'+'\n')
# for key in c:
#     handle.write(str(key)+'\t'+str(c[key])+'\n')




            






    








