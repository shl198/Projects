from __future__ import division
import os
import pandas as pd
from natsort import natsorted
from f02_RiboDataModule import *
from Modules.f05_IDConvert import addGeneIDorNameForDESeqResult

signalP_path = '/data/shangzhong/RibosomeProfiling/signalP_part'
ribo_bam_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam'
ribo_gene_count_path = ribo_bam_path + '/04_gene_total_count'

rna_bam_path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align'
rna_gene_count_path = rna_bam_path + '/01_gene_count'
#===============================================================================
#                         1. generate count data (rpm) for day3 and day6 proteomaps
#===============================================================================
# =========== convert cho genes to mouse genes ========================
def splitMouseIDs(gene_id):
    if ';' in gene_id:
        return gene_id.split(';')
    else:
        return gene_id.split(',')
 
def cho2mouseID(dic,gene):
    if gene in dic:
        return dic[gene]
    else:
        return '0'
import time
start_time = time.time()
 
def cho_2mouse_count(gene_count_df,cho_geneIDs,cho_mus_dict):
    """
    This function change the count dataframe with cho gene id to mouse id
    * gene_count_df: df. 
    * cho_geneIDs: list. A list of cho_geneIDs
    * cho_mus_dict: dictionary. {choid:[mouse ids]}
    return mouse_df
    """
    n = 0
    columns = gene_count_df.columns.tolist()
    mouse_df = pd.DataFrame(columns=columns)
    for gene in cho_geneIDs:
        if gene in cho_mus_dict:
            mus_gene = cho_mus_dict[gene]
        else:
            n = n + 1
            print gene,'is not in the cho_mus_dictionary'
            continue
        unit_count = (gene_count_df.loc[gene]/float(len(mus_gene))).tolist()
        for g in mus_gene:
            if g not in mouse_df.index:
                mouse_df.loc[g] =  unit_count
            else:
                mouse_df.loc[g] = mouse_df.loc[g] + unit_count
    return mouse_df
# 1). generate dictionary {cho:mouse}
cho2musFile = '/data/shangzhong/CHO2Mouse/finalMergeWithmRNA.final.txt'
cho2mus_df = pd.read_csv(cho2musFile,header=None,sep='\t',names=['cho','mus'])
cho2mus_df = cho2mus_df.astype(str)
cho2mus_df['mus'] = cho2mus_df['mus'].map(lambda x: splitMouseIDs(x))
cho_mus_dict = cho2mus_df.set_index('cho')['mus'].to_dict()
cho_mus_dict['heavychain'] = ['heavychain']
cho_mus_dict['lightchain'] = ['lightchain']
cho_mus_dict['NeoRKanR'] = ['NeoRKanR']
# 2). read rpm data data
gene_countFile = signalP_path + '/08_ribo_rna_rpm.csv'
gene_count_df = pd.read_csv(gene_countFile,sep='\t',header=0,index_col=0,names=['geneid','ribo_day3','ribo_day6','rna_day3','rna_day6'])
cho_geneIDs = gene_count_df.index.astype(str).tolist()
   
mouse_df = cho_2mouse_count(gene_count_df,cho_geneIDs,cho_mus_dict)
mouse_df[['ribo_day3']].to_csv(signalP_path+'/09_day3_rpm.csv',sep='\t',header=None)
mouse_df[['ribo_day6']].to_csv(signalP_path+'/09_day6_rpm.csv',sep='\t',header=None)
mouse_df.to_csv(signalP_path+'/10_mouse_ribo_rna_rpm.csv',sep='\t')
 
# 3) build replicate count dataframe
os.chdir(ribo_gene_count_path)
ribo_count_files = [f for f in os.listdir(ribo_gene_count_path) if f.endswith('Count.txt')]
ribo_count_files = natsorted(ribo_count_files)
dfs = []
for f in ribo_count_files:
    df = pd.read_csv(f,header=0,sep='\t',index_col=0,usecols=[0,1],names=['geneid','ribo_'+f[:3]])
    dfs.append(df)
   
os.chdir(rna_gene_count_path)
rna_count_files = [f for f in os.listdir(rna_gene_count_path) if f.endswith('Count.txt')]
rna_count_files = natsorted(rna_count_files)
for f in rna_count_files:
    df = pd.read_csv(f,header=0,sep='\t',index_col=0,usecols=[0,1],names=['geneid','rna_'+f[:3]])
    dfs.append(df)
ribo_rna_count_df = pd.concat(dfs,axis=1,join='inner')
 
mouse_df = cho_2mouse_count(ribo_rna_count_df,cho_geneIDs,cho_mus_dict)
mouse_df.to_csv(signalP_path+'/10_mouse_ribo_rna_count.csv',sep='\t')
#----------- calculate each gene percentage for each replicate ------------
total = ribo_rna_count_df.sum().tolist()
percent_df = ribo_rna_count_df.div(total)
columns = ['ribo_day3_rep1','ribo_day3_rep2','ribo_day3_rep3',
              'ribo_day6_rep1','ribo_day6_rep2','ribo_day6_rep3',
              'rna_day3_rep1','rna_day3_rep2','rna_day3_rep3',
              'rna_day6_rep1','rna_day6_rep2','rna_day6_rep3']
percent_df.columns = columns
percent_df['ribo_day3_mean'] = percent_df.iloc[:,0:3].mean(axis=1)
percent_df['ribo_day6_mean'] = percent_df.iloc[:,3:6].mean(axis=1)
percent_df['rna_day3_mean'] = percent_df.iloc[:,6:9].mean(axis=1)
percent_df['rna_day6_mean'] = percent_df.iloc[:,9:12].mean(axis=1)
 
percent_df['ribo_day3_std'] = percent_df.iloc[:,0:3].std(axis=1)
percent_df['ribo_day6_std'] = percent_df.iloc[:,3:6].std(axis=1)
percent_df['rna_day3_std'] = percent_df.iloc[:,6:9].std(axis=1)
percent_df['rna_day6_std'] = percent_df.iloc[:,9:12].std(axis=1)
outFile = signalP_path + '/12_cho_gene_percent_rep.csv'
percent_df.to_csv(outFile,sep='\t')
# add gene symbol to the file
MapFile = '/data/shangzhong/Database/cho/gff_chok1_ID_symbol.txt'
addGeneIDorNameForDESeqResult(signalP_path + '/12_cho_gene_percent_rep.csv',MapFile,addType='gene_symbol',IDVersion='no')

######## the following is faster, but needs more code.
# ribo_day3_dict = {}; ribo_day6_dict={}; rna_day3_dict = {}; rna_day6_dict = {}
# for gene in cho_geneIDs:
#     if gene in cho_mus_dict:
#         mus_gene = cho_mus_dict[gene]
#     else:
#         n = n + 1
#         print gene,'is not in the cho_mus_dictionary'
#         continue
#          
#     per_count3 = gene_count_df.loc[gene,'ribo_day3']/len(mus_gene)  # evenly distribute to all the mus_gene
#     per_count6 = gene_count_df.loc[gene,'ribo_day6']/len(mus_gene)
#     rna_count3 = gene_count_df.loc[gene,'rna_day3']/len(mus_gene)
#     rna_count6 = gene_count_df.loc[gene,'rna_day6']/len(mus_gene)
#     for mus in mus_gene:
#         if mus not in ribo_day3_dict:
#             ribo_day3_dict[mus] = per_count3
#             ribo_day6_dict[mus] = per_count6
#             rna_day3_dict[mus] = rna_count3
#             rna_day6_dict[mus] = rna_count6
#         else:
#             ribo_day3_dict[mus] = ribo_day3_dict[mus] + per_count3
#             ribo_day6_dict[mus] = ribo_day6_dict[mus] + per_count6
#             rna_day3_dict[mus] = rna_day3_dict[mus] + rna_count3
#             rna_day6_dict[mus] = rna_day6_dict[mus] + rna_count6
# ribo_day3_df = pd.DataFrame(ribo_day3_dict.items(),columns=['gene','ribo_day3'])
# ribo_day6_df = pd.DataFrame(ribo_day6_dict.items(),columns=['gene','ribo_day6'])
# rna_day3_df = pd.DataFrame(rna_day3_dict.items(),columns=['gene','rna_day3'])
# rna_day6_df = pd.DataFrame(rna_day6_dict.items(),columns=['gene','rna_day6'])
# ribo_day3_df = ribo_day3_df.set_index(['gene'])
# ribo_day6_df = ribo_day6_df.set_index(['gene'])
# rna_day3_df = rna_day3_df.set_index(['gene'])
# rna_day6_df = rna_day6_df.set_index(['gene'])
#   
# res_df = pd.concat([ribo_day3_df,ribo_day6_df,rna_day3_df,rna_day6_df],axis=1)
# ribo_day3_df.to_csv(signalP_path+'/09_day3_rpm.csv',sep='\t',header=None)
# ribo_day6_df.to_csv(signalP_path+'/09_day6_rpm.csv',sep='\t',header=None)
# res_df.to_csv(signalP_path+'/10_mouse_ribo_rna_rpm.csv',sep='\t')
# print n
# #===============================================================================
# #                         2. function protein percentage
# #===============================================================================
# #--------------- 1. read tree file ------------------
# treeFile = '/data/shangzhong/RibosomeProfiling/figures/Proteomaps/KO_gene_hierarchy_2014.tms'
# handle = open(treeFile,'r')
# fun_kegg_dic = {}
# func = ''
# for line in handle:
#     item = line.rstrip().split('\t')
#     if len(item) == 1 or len(item) == 3:
#         continue
#     if len(item) == 2:
#         func = item[1]
#     if len(item) == 4:
#         if func in fun_kegg_dic:
#             fun_kegg_dic[func].append(item[3])
#         else:  # function not in the dictionary
#             fun_kegg_dic[func] = []
# #--------------- 2. read kegg gene mapping file ------------------
# kegg_gene_file = '/data/shangzhong/RibosomeProfiling/figures/Proteomaps/mmu_mapping.csv'
# k_g_df = pd.read_csv(kegg_gene_file,sep='\t',header=None,names=['id','genename','kegg','fullname'])
# k_g_df['id'] = k_g_df['id'].astype(str)
# k_g_dic = {k:list(v) for k,v in k_g_df.groupby('kegg')['id']}
# #--------------- 3. build function:gene dictionary ---------------
# func_gene_dic = {}
# n = 0
# m = 0
# for func in fun_kegg_dic:
#     func_gene_dic[func] = []
#     keggs = fun_kegg_dic[func]
#     for k in keggs:
#         try:
#             genes = k_g_dic[k]
#             func_gene_dic[func].extend(genes)
#         except:
#             n = n + 1
#             print k,'not in the fun_kegg_dic'
# func_gene_dic['heavychain'] = ['heavychain']
# func_gene_dic['lightchain'] = ['lightchain']
# func_gene_dic['NeoRKanR'] = ['NeoRKanR']
# db_genes = []
# for key in func_gene_dic:
#     db_genes.extend(func_gene_dic[key])
# #---------------- 4. calculate percentage for each category -------
# ribo_rna_rpm_file = signalP_path +'/10_mouse_ribo_rna_rpm.csv'
# ribo_rna_rpm_df = pd.read_csv(ribo_rna_rpm_file,sep='\t',header=0,index_col=0)
# ribo_rna_rpm_df = ribo_rna_rpm_df[ribo_rna_rpm_df.index.isin(db_genes)]
# total = ribo_rna_rpm_df.sum().tolist()
# outFile = signalP_path + '/11_proteomap_percentage.csv'
# handle = open(outFile,'w')
# handle.write('\t'.join(['Category','ribo_day3','ribo_day6','rna_day3','rna_day6'])+'\n')
# for func in func_gene_dic:
#     genes = func_gene_dic[func]
#     sub_df = ribo_rna_rpm_df[ribo_rna_rpm_df.index.isin(genes)]
#     if sub_df.empty:
#         continue
#     percent = ((sub_df.sum()).div(total)).tolist()
#     percent = [str(p) for p in percent]
#     handle.write(func+'\t'+'\t'.join(percent)+'\n')
# handle.close()
# df = pd.read_csv(outFile,sep='\t',header=0,index_col=0)
# df = df.sort(['ribo_day3'],ascending=False)
# df.to_csv(outFile,sep='\t')
# #----------------- 5. calculate the percentage for each replicate ---------
# ribo_rna_count_file = signalP_path + '/10_mouse_ribo_rna_count.csv'
# ribo_rna_count_df = pd.read_csv(ribo_rna_count_file,sep='\t',header=0,index_col=0)
# ribo_rna_count_df = ribo_rna_count_df[ribo_rna_count_df.index.isin(db_genes)]
# total = ribo_rna_count_df.sum().tolist()
# gene_percent_df = ribo_rna_count_df.div(total)
# columns = ['ribo_day3_rep1','ribo_day3_rep2','ribo_day3_rep3',
#               'ribo_day6_rep1','ribo_day6_rep2','ribo_day6_rep3',
#               'rna_day3_rep1','rna_day3_rep2','rna_day3_rep3',
#               'rna_day6_rep1','rna_day6_rep2','rna_day6_rep3']
# gene_percent_df.columns = columns
# 
# outFile = signalP_path + '/11_proteomap_percentage_rep.csv' 
# handle = open(outFile,'w')
# handle.write('Category'+'\t'+'\t'.join(columns)+'\n')
# for func in func_gene_dic:
#     genes = func_gene_dic[func]
#     sub_df = ribo_rna_count_df[ribo_rna_count_df.index.isin(genes)]
#     if sub_df.empty:
#         continue
#     percent = ((sub_df.sum()).div(total)).tolist()
#     percent = [str(p) for p in percent]
#     handle.write(func+'\t'+'\t'.join(percent)+'\n')
# handle.close()
# df = pd.read_csv(outFile,sep='\t',header=0,index_col=0)
# df = df.sort([columns[0]],ascending=False)
#  
# df['ribo_day3_std'] = df.iloc[:,0:3].std(axis=1)
# df['ribo_day6_std'] = df.iloc[:,3:6].std(axis=1)
# df['rna_day3_std'] = df.iloc[:,6:9].std(axis=1)
# df['rna_day6_std'] = df.iloc[:,9:12].std(axis=1)
# df.to_csv(outFile,sep='\t')
# 
# gene_percent_df['ribo_day3_std'] = gene_percent_df.iloc[:,0:3].std(axis=1)
# gene_percent_df['ribo_day6_std'] = gene_percent_df.iloc[:,3:6].std(axis=1)
# gene_percent_df['rna_day3_std'] = gene_percent_df.iloc[:,6:9].std(axis=1)
# gene_percent_df['rna_day6_std'] = gene_percent_df.iloc[:,9:12].std(axis=1)
# gene_percent_df.to_csv(signalP_path + '/12_mouse_gene_percent_rep.csv',sep='\t')


#### rethink from here
#===============================================================================
#                         2. custamize protein annotation for proteomaps
#===============================================================================
"""
# 1. read the id file
human_disease = '/data/shangzhong/RibosomeProfiling/figures/Proteomaps/KO_human_disease.txt'
human_disease_df = pd.read_csv(human_disease,sep='\t',header=0,names=['path_name','kegg_id'])
name_id_dict = {k:list(v) for k,v in human_disease_df.groupby('kegg_id')['path_name']}
kegg = human_disease_df['kegg_id'].tolist()
# 2. read mouse id mapping
mmu = '/data/shangzhong/RibosomeProfiling/figures/Proteomaps/mmu_mapping.csv'
mmu_df = pd.read_csv(mmu,sep='\t',header=None,names=['GeneID','GeneSymbol','kegg_id','gene_name'])
# 3. generate the annotation table
criteria = mmu_df['kegg_id'].map(lambda x: x in kegg)
anno_df = mmu_df[criteria]
anno_df = anno_df.reset_index()
anno_df['Pathway'] = anno_df['kegg_id'].apply(lambda x: name_id_dict[x][0])
anno_df['Org'] = pd.Series(['mmu'] * (anno_df.shape[0]))
anno_df = anno_df[['Org','GeneID','GeneSymbol','kegg_id','Pathway','gene_name']]
anno_df.to_csv(mmu[:-3]+'reanno.csv',sep='\t',index=False)
"""

#===============================================================================
#                         3. check percentage of reads that map to the whole genome are represented proteomaps
#===============================================================================
"""
# 1) get total count for each sample
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_cov5end'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
coverFiles = natsorted(coverFiles)
totalCount = []
for f in coverFiles:
    df = pd.read_csv(f,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
    total = df['coverage'].sum()
    totalCount.append(total)
# 2) calculate percentage of reads that are included by the genes mappable in the cho2mouse id map
raw_count_file = '/data/shangzhong/RibosomeProfiling/cho_pr/12_pr_gene_rawCount.txt'
raw_count_df = pd.read_csv(raw_count_file,sep='\t',header=0)
# map file
cho2musFile = '/data/shangzhong/CHO2Mouse/finalMergeWithmRNA.final.txt'
cho2mus_df = pd.read_csv(cho2musFile,header=None,sep='\t',names=['cho','mus'])
cho2mus_df = cho2mus_df.astype(str)
mappableIDs = cho2mus_df['cho'].tolist()
print totalCount
cri = raw_count_df['GeneID'].map(lambda x: x in mappableIDs)
filtered_genes = raw_count_df[cri]
sum_count = filtered_genes.iloc[:,1:].sum(axis=0).tolist()
percent = [c/total for c,total in zip(sum_count,totalCount)]
print 'All genes mapping count percent is:',percent

cri = raw_count_df['GeneID'].map(lambda x: x in ['heavychain','lightchain','NeoRKanR'])
filtered_genes = raw_count_df[cri]
sum_count = filtered_genes.iloc[:,1:].sum(axis=0).tolist()
percent = [c/total for c,total in zip(sum_count,totalCount)]
print 'Antibody mapping count percent is:',percent
"""

#===============================================================================
#                         4. get how much percent of protein coding genes mapping reads are representaed by proteomaps 
#===============================================================================
def mappableID(gene,cho_mus_dict,mouseMapIDs):
    if str(gene) in cho_mus_dict:
        return set(cho_mus_dict[str(gene)]).intersection(mouseMapIDs) != set()
    else:
        return False
"""
# 1) get the cho 2 mouse dictionary
cho2musFile = '/data/shangzhong/CHO2Mouse/finalMergeWithmRNA.final.txt'
cho2mus_df = pd.read_csv(cho2musFile,header=None,sep='\t',names=['cho','mus'])
cho2mus_df = cho2mus_df.astype(str)
cho2mus_df['mus'] = cho2mus_df['mus'].map(lambda x: splitMouseIDs(x))
cho_mus_dict = cho2mus_df.set_index('cho')['mus'].to_dict()
# 2) get the mouse mappable ids in proteomaps
mmuMapFile = '/data/shangzhong/RibosomeProfiling/figures/Proteomaps/mmu_mapping.csv'
mmuMap_df = pd.read_csv(mmuMapFile,sep='\t',header=None,usecols=[0],names=['GeneID'])
mouseMapIDs = mmuMap_df['GeneID'].astype(str).tolist()
# 3) read gene coverage file
gene_rawCount_file = '/data/shangzhong/RibosomeProfiling/cho_pr/12_pr_gene_rawCount.txt'
gene_raw_count_df = pd.read_csv(gene_rawCount_file,sep='\t',header=0,index_col=0)
antibody_df = gene_raw_count_df[gene_raw_count_df.index.isin(['heavychain','lightchain','NeoRKanR'])]  # 3 antibody raw count
total = gene_raw_count_df.sum().tolist()  # total number of reads that map to protein coding genes
cri = gene_raw_count_df.index.map(lambda x: mappableID(x,cho_mus_dict,mouseMapIDs))
mappable_df = gene_raw_count_df[cri].append(antibody_df) # genes in cho that can be mapped to proteomaps
count_sum = mappable_df.sum().tolist()  # total count of the reads that are represented by proteomaps
percent = [m/n for m,n in zip(count_sum,total)]
print total
print count_sum
print percent

antibody_df.div(count_sum)
"""





