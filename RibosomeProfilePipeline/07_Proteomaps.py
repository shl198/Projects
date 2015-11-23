from __future__ import division
import os
import pandas as pd
from natsort import natsorted
from f02_RiboDataModule import *


#===============================================================================
#                         1. generate count data (rpm) for day3 and day6 proteomaps
#===============================================================================

def riboGeneCount(posFiles,totalCounts,countType='rawCount'):
    """
    This function does the count 
    """
    gene_count_df = pd.DataFrame()
    for f,total in zip(posFiles,totalCounts):
        geneIDs = []
        handle = open(f,'r')
        counts = []
        for line in handle:
            item = line[:-1].split('\t')
            gene = item[0]
            geneIDs.append(gene)
            cov = item[1:]
            for i in range(len(cov)):
                if cov[i]=='-': cov[i] = 0
                else: cov[i] = int(cov[i])
            if countType == 'rawCount':
                count = sum(cov)
            elif countType == 'rpm':
                count = sum(cov)/total*(10**6)
            elif countType == 'rpkm':
                count = sum(cov)/total/len(cov)*(10**9)
            counts.append(count)        
        gene_count_df[f[:3]] = pd.Series(counts)
    gene_count_df.insert(0,'GeneID',geneIDs)
    outFile = '/data/shangzhong/RibosomeProfiling/cho_pr/12_pr_gene_' + countType + '.txt'
    gene_count_df.to_csv(outFile,sep='\t',index=False)
"""
# 1. read coverage file and get total count
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_cov5end'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
coverFiles = natsorted(coverFiles)
totalCount = []
for f in coverFiles:
    df = pd.read_csv(f,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
    total = df['coverage'].sum()
    totalCount.append(total)
# 2. get count for each gene
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/04_Gene_raw_cov'
os.chdir(path)
posFiles = [f for f in os.listdir(path) if f.endswith('sort.geneAllposCov.txt')]
posFiles = natsorted(posFiles)

riboGeneCount(posFiles,totalCount,countType='rawCount')
riboGeneCount(posFiles,totalCount,countType='rpm')
"""
# =========== convert cho genes to mouse genes ========================
def splitMouseIDs(gene_id):
    if ';' in gene_id:
        return gene_id.split(';')
    else:
        return gene_id.split(',')
"""
# 1). generate dictionary {cho:mouse}
cho2musFile = '/data/shangzhong/CHO2Mouse/finalMergeWithmRNA.final.txt'
cho2mus_df = pd.read_csv(cho2musFile,header=None,sep='\t',names=['cho','mus'])
cho2mus_df = cho2mus_df.astype(str)
cho2mus_df['mus'] = cho2mus_df['mus'].map(lambda x: splitMouseIDs(x))
cho_mus_dict = cho2mus_df.set_index('cho')['mus'].to_dict()
# 2). read count data
gene_countFile = '/data/shangzhong/RibosomeProfiling/cho_pr/12_pr_gene_rpm.txt'
gene_count_df = pd.read_csv(gene_countFile,sep='\t',header=0,index_col=0)
gene_count_df = gene_count_df.fillna(0)
cho_geneIDs = gene_count_df.index.astype(str).tolist()
count_df = pd.DataFrame()
count_df['day3'] = gene_count_df.iloc[:,0:3].mean(axis=1)
count_df['day6'] = gene_count_df.iloc[:,3:6].mean(axis=1)
# 3). define transfered dictionary
day3_dict = {}; day6_dict={}
n = 0
for gene in cho_geneIDs:
    try:
        mus_gene = cho_mus_dict[gene]
    except:
        n = n + 1
        print gene,'is not in the cho_mus_dictionary'
        #continue
        mus_gene = [gene]  # this would output the cho gene ids that don't map to proteomaps
        
    per_count3 = count_df.loc[gene,'day3']/len(mus_gene)  # evenly distribute to all the mus_gene
    per_count6 = count_df.loc[gene,'day6']/len(mus_gene)
    for mus in mus_gene:
        if mus not in day3_dict:
            day3_dict[mus] = per_count3
            day6_dict[mus] = per_count6
        else:
            day3_dict[mus] = day3_dict[mus] + per_count3
            day6_dict[mus] = day6_dict[mus] + per_count6
day3_df = pd.DataFrame(day3_dict.items(),columns=['gene','count'])
day6_df = pd.DataFrame(day6_dict.items(),columns=['gene','count'])
day3_df.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/13_day3_rpm_mouse.csv',sep='\t',index=False,header=None)
day6_df.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/13_day6_rpm_mouse.csv',sep='\t',index=False,header=None)
print n
"""
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





