import pandas as pd
import gffutils
import re
from natsort import natsorted
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('ggplot')


def tr2pr(tr,ensembl_trpr,ncbi_trpr):
    if tr.startswith('ENS'):
        if '_' in tr:
            tr = '_'.join(tr.split('_')[:-1])
        if tr in ensembl_trpr:
            return ensembl_trpr[tr]
        else:
            return tr
    else:
        if '.' in tr:
            tr = tr.split('.')[0]
        if tr in ncbi_trpr:
            return ncbi_trpr[tr]
        else:
            return tr

def parse_mutation(mutation):
    trim = mutation[2:]
    if trim == '?':
        return ['?']*3
    res = re.findall('[^\d_]+|\d+',trim)
    ref = res[0]
    pos = res[1]
    alt = res[2]
    if 'del' in res[-1]:
        ref = res[-1][3:]
        alt = '.'
    return [pos,ref,alt]

def remove_version(access):
    access= str(access)
    if '.' in access:
        res = access.split('.')[0]
    else:
        res = access
    return res
gene2ref = '/data/shangzhong/Database/human/141026gene2refseq.human.txt'
mutFile = '/data/shangzhong/Dream/mutations.csv'
"""
# 1. build transcript to protein dictionary
# ensembl annotation dictionary
gffFile = '/data/shangzhong/Dream/Homo_sapiens.GRCh37.82.gff3'
handle = open(gffFile,'r')
tr_pr = {}  # {tr:pr}
for line in handle:
    if line.startswith('#'):
        continue
    item = line[:-1].split('\t')
    if item[2]=='CDS':
        anno = item[-1]
        index = anno.index(';')
        pr = anno[7:index]  # len('ID=CDS:')
        anno = anno[index+1:]
        index = anno.index(';')
        tr = anno[len('Parent=transcript:'):index]
        if tr not in tr_pr:
            tr_pr[tr] = pr

# ncbi dictionary
gene2ref = '/data/shangzhong/Database/human/141026gene2refseq.human.txt'
g_ref_df = pd.read_csv(gene2ref,sep='\t',header=None,names=['tr','pr'])
g_ref_df = g_ref_df.drop_duplicates()
g_ref_df['new_tr'] = g_ref_df['tr'].map(lambda x: remove_version(x))
g_ref_df['new_pr'] = g_ref_df['pr'].map(lambda x: remove_version(x))
ncbi_trpr = g_ref_df.set_index('new_tr')['new_pr'].to_dict()
# 2. get the input for provean
mutFile = '/data/shangzhong/Dream/mutations.csv'
mut_df = pd.read_csv(mutFile,header=0,usecols=[1,14],names=['transcript','mut'])
mut_df['protein'] = mut_df['transcript'].map(lambda x: tr2pr(x,tr_pr,ncbi_trpr))
mut_df['mutation'] = mut_df['mut'].map(lambda row: parse_mutation(row))
out = open('/data/shangzhong/Dream/protein_mutation.csv','w')
for row in mut_df.itertuples(index=False):
    out.write(row[2]+','+','.join(row[3])+'\n')
out.close()
"""

# # 3. merge provean with the origional
# path = '/data/shangzhong/Dream'
# os.chdir(path)
# files = [f for f in os.listdir(path) if f.endswith('.tsv')]
# files = natsorted(files)
# dfs = []
# for f in files:
#     df = pd.read_csv(f,header=None,sep='\t',usecols=[2,6,7],names=['protein','provean_score','effect'],skiprows=[0])
#     dfs.append(df)
# merge_df = pd.concat(dfs,ignore_index=True)
#     
# mut_df = pd.read_csv(mutFile,header=0)
# res = pd.concat([mut_df,merge_df],axis=1)
# res.to_csv('/data/shangzhong/Dream/anno_mut.csv',index=False)



# # Score = 3 * ln(X * e^(M-u))
# exp_file = '/data/shangzhong/Dream/Molecular_Exp.csv'
# mut_file = '/data/shangzhong/Dream/Molecular_Mut.csv'
# exp_df = pd.read_csv(exp_file,header=0,index_col=0)
# print exp_df
# gene1 = exp_df.loc[:,'22RV1']
# gene1 = gene1[gene1>0]
# gene1.hist()
# print gene1
# plt.show()
# #mut_file = pd.read_csv(mut_file,header=0,index_col=0)
# #score = 3 *  ((exp_df.loc[:,:] * (mut_file['22RV1']+2.5)).apply(np.exp)).apply(np.log2)
# print 'done'




#===============================================================================
#                 plot correlation between cnv data and expression data
#===============================================================================
# sample = '22RV1'
# f1 = '/data/shangzhong/Dream/gex.csv'
# df1 = pd.read_csv(f1,header=0)
# df1 = df1[['symbol',sample]]
# dic1 = df1.set_index('symbol')[sample].to_dict()
# df1 = pd.DataFrame.from_dict(dic1,orient='index')
# df1.columns = [sample+'_1']
# 
# f2 = '/data/shangzhong/Dream/CNV_all.csv'
# df2 = pd.read_csv(f2,header=0)
# df2 = df2[['symbol',sample]]
# dic2 = df2.set_index('symbol')[sample].to_dict()
# df2 = pd.DataFrame.from_dict(dic2,orient='index')
# df2.columns = [sample + '_2']
# 
# df = pd.concat([df1,df2],axis=1)
# df = df.dropna()
# plt.scatter(df['22RV1_1'],df['22RV1_2'])
# plt.show()
# print 'done'


#===============================================================================
#             run codra
#===============================================================================
import subprocess
path = '/data/shangzhong/Dream/CORDA'
os.chdir('/data/shangzhong/Dream/CORDA')
num = range(72,83,2)
num = num[:-1]
cmd = ''
for n in num:
    start=n+1
    if n ==80:
        end=83
    else:
        end = n+2
    cmd = cmd + ('matlab -r \"createDREAMmodels({start},{end})\" & ').format(start=start,end=end)
subprocess.call(cmd,shell=True)


# # # merge the results
# path = '/data/shangzhong/Dream/CORDA'
# os.chdir(path)
# files = [f for f in os.listdir(path) if f.startswith('MetabolicCapabilities')]
# files = natsorted(files)
# res_df = pd.read_csv(files[0],header=0)
# for f in files[1:]:
#     df = pd.read_csv(f,header=0)
#     res_df = res_df + df
# res_df.to_csv(path + '/merge_metaCap.csv',index=False)




#===============================================================================
#         merge CNV and mutation
#===============================================================================
# cnvFile = '/data/shangzhong/Dream/Molecular_CNV.csv'
# cnv_df = pd.read_csv(cnvFile,header=0,index_col=0)
# cnv_df[cnv_df<0] = 0
# 
# mutFile = '/data/shangzhong/Dream/Molecular_Mut.csv'
# mut_df = pd.read_csv(mutFile,header=0,index_col=0)
# mut_df[mut_df==-100] = 3
# 
# res_df = cnv_df + mut_df
# res_df.to_csv('/data/shangzhong/Dream/cnv_mut.csv')


#===============================================================================
#                     imputate the lipinski score
#===============================================================================
# fn = '/data/shangzhong/Dream/Drug_info_release.csv'
# df = pd.read_csv(fn,header=0)
# lipinski = []
# for row in df.itertuples():
#     res = 0
#     if row[3] > 10: res = res + 1
#     if row[4] > 5: res = res + 1
#     if row[5] > 5: res = res + 1
#     if row[-1] > 500: res = res + 1
#     lipinski.append(res)
# df['Lipinski'] = lipinski
# df.to_csv('/data/shangzhong/Dream/Drug_info.csv',index=False)
    







