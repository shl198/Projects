from __future__ import division
import os
from natsort import natsorted
import pandas as pd
from f02_RiboDataModule import *
from multiprocessing import Process

#===============================================================================
#                         1. get coverage file for mRNA seqs
#===============================================================================
# Run the 8th part in the 02_RiboSignalpCoverPipeline
# path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align'
# os.chdir(path)
# bamFiles = [f for f in os.listdir(path) if f.endswith('.sort.bam')]
# cov5end = pos5Coverage(bamFiles,batch=6)
# # move to a folder
# folder = '01_cov5end'
# if not os.path.exists(folder):
#     os.mkdir(folder)
# for f in cov5end:
#     cmd = ('mv {input} {out}').format(input=f,out=folder)
#     subprocess.call(cmd.split(' '))
#===============================================================================
#                         2. get mRNA coverage for all positions of genes
#===============================================================================
def wrap_covAllPos(target_path,exonCovFile,exonBedFile,geneIDs,chr_len_file):
    resFile = RNAcovAllPos(exonCovFile,exonBedFile,geneIDs,chr_len_file)
    cmd = ('mv {old} {new}').format(old=resFile,new=target_path)
    subprocess.call(cmd.split(' '))
"""
# 1).read chromosome length
chr_len_file = '/data/shangzhong/RibosomeProfiling/Database/combined_ChrLen.txt'
# 2).sp genes and no sp genes
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
df = pd.read_csv(gene_class_file,sep='\t',header=0)
gene_sp = df['gene_sp'].dropna().astype(str).tolist()
gene_no_sp = df['gene_no_sp'].dropna().astype(str).tolist()
genes = gene_sp + gene_no_sp
# 3).read cds position file
exonBedFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_cho_ab_pr_sp.gene.exon.bed'
# 4) read coverage file
path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align/01_cov5end'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('5endCov.txt')]
coverFiles = natsorted(coverFiles)
# 5) define target path
target_path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align/02_allGenePos_cov'
if not os.path.exists(target_path):
    os.mkdir(target_path)
# get position
#filename = covAllPos(coverFiles[0],exonBedFile,genes,chr_len_file,up,down)

proc = []
for f in coverFiles:
    p1 = Process(target=wrap_covAllPos,args=(target_path,f,exonBedFile,genes,chr_len_file))    
    p1.start()
    proc.append(p1)
for p in proc:
    p.join()
"""
#===============================================================================
#                         3. Calculate translation efficiency
#===============================================================================
"""
# 1. total count of reads in riboseq
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_cov5end'
os.chdir(path)
coverFiles = [os.path.join(path,f) for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
chrPosCovFiles = natsorted(coverFiles)
# 2. coverage for count in riboseq
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/04_Gene_raw_cov'
os.chdir(path)
coverFiles = [os.path.join(path,f) for f in os.listdir(path) if f.endswith('geneRawCount.txt')]
geneCovFiles = natsorted(coverFiles)
# 3. run for ribo seq data
outFile = '/data/shangzhong/RibosomeProfiling/cho_pr/13_pr_gene_rpkm.txt'
mergeGeneExpress(outFile,geneCovFiles,chrPosCovFiles,'rpkm')

# 4. total count of count in rna
path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align/01_cov5end'
os.chdir(path)
coverFiles = [os.path.join(path,f) for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
chrPosCovFiles = natsorted(coverFiles)
# 5. coverage for count in riboseq
path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align/02_allGenePos_cov'
os.chdir(path)
coverFiles = [os.path.join(path,f) for f in os.listdir(path) if f.endswith('geneRawCount.txt')]
geneCovFiles = natsorted(coverFiles)
# 6. run for ribo seq data
outFile = '/data/shangzhong/RibosomeProfiling/cho_pr/13_rna_gene_rpkm.txt'
mergeGeneExpress(outFile,geneCovFiles,chrPosCovFiles,'rpkm')
"""
#===============================================================================
#                         4. extract off diagnal genes transfer to mouse genes
#===============================================================================

res_df = pd.DataFrame()
# 1. read ribo count
riboFile = '/data/shangzhong/RibosomeProfiling/cho_pr/13_pr_gene_rpkm.txt'
ribo_df = pd.read_csv(riboFile,header=0,sep='\t',index_col=0)
res_df['day3_ribo'] = (ribo_df.iloc[:,0:3]).mean(axis=1)
res_df['day6_ribo'] = (ribo_df.iloc[:,3:6]).mean(axis=1)
# 2. read rna count
rnaFile = '/data/shangzhong/RibosomeProfiling/cho_pr/13_rna_gene_rpkm.txt'
rna_df = pd.read_csv(rnaFile,header=0,sep='\t',index_col=0)
res_df['day3_rna'] = (rna_df.iloc[:,0:3]).mean(axis=1)
res_df['day6_rna'] = (rna_df.iloc[:,3:6]).mean(axis=1)
# 3. calculate translation efficiency
# res_df = res_df + 0.001
cri = (res_df['day3_ribo']>1) & (res_df['day6_ribo']>1) & (res_df['day3_rna']>1) & (res_df['day6_rna']>1)
res_df = res_df[cri]
res_df['ribo_change'] = (res_df['day6_ribo']/res_df['day3_ribo']).apply(np.log2)
res_df['mrna_change'] = (res_df['day6_rna']/res_df['day3_rna']).apply(np.log2)
# get the diagnal genes
cri = ((res_df['ribo_change']>1) & (res_df['mrna_change']<-0)) | ((res_df['ribo_change']>0) & (res_df['mrna_change']<-1))
riboUp = res_df[cri]
up_genes = riboUp.index.tolist()
print 'up_genes number',len(up_genes)
cri = ((res_df['ribo_change']<-1) & (res_df['mrna_change']>0)) | ((res_df['ribo_change']<-0) & (res_df['mrna_change']>1))
riboDown = res_df[cri]
down_genes = riboDown.index.tolist()
print 'down_genes number',len(down_genes)
# convert to mouse genes
# build dictionary
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

df1 = pd.DataFrame({'riboUp':up_genes})
df2 = pd.DataFrame({'riboDown':down_genes})
res_df = pd.concat([df1,df2],axis=1)
res_df.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/15_TransEffiOffDiag.csv',index=False,sep='\t')

up_mouse = [];down_mouse=[]
for g in up_genes:
    try:
        up_mouse.extend(dic[g])
    except:
        print g,'not in dictionary'

for g in down_genes:
    try:
        down_mouse.extend(dic[g])
    except:
        print g,'not in dictionary'

# output
df1 = pd.DataFrame({'riboUp':up_mouse})
df2 = pd.DataFrame({'riboDown':down_mouse})
res_df = pd.concat([df1,df2],axis=1)
res_df.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/15_TransEffiOffDiag_mouse.csv',index=False,sep='\t')
print df1.shape
print df2.shape





