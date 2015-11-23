import numpy as np
import scipy.stats as sp_stats
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('ggplot')
import os
import pandas as pd
from natsort import natsorted
from f02_RiboDataModule import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from scipy import stats
#===============================================================================
#            plot gene cds length stats
#===============================================================================
# df = pd.read_csv('/data/shangzhong/RibosomeProfiling/cho_pr/10_CDS_len_stats_count.csv',sep='\t',header=0)
# df['index'] = pd.Series(df.index)
# matrix = df.as_matrix()
# df = df.apply(lambda x: np.log2(x))
# 
# gene_no_sp, = plt.plot(df['index'],df['gene_no_sp'],label='gene_no_sp')
# gene_sp, = plt.plot(df['index'],df['gene_sp'],label='gene_sp')
# plt.xlabel('log2(position)')
# plt.ylabel('log2(number of genes)')
# plt.title('gene cds length stats')
# plt.legend(handles=[gene_no_sp,gene_sp])
# plt.savefig('/data/shangzhong/RibosomeProfiling/cho_pr/CDS_len_stats_count.pdf')
# plt.show()

def getLowCountGenes(coverFiles,appendix,count):
    """
    This function gets gene list with counts lower than the threshold count
    
    * coverFiles: list. A list of files with row as genes, columns as percentile. values are count
    * appendix: str. Name which will be append to filename in coverFiles, and then be the column name
    * count: int. Footprints less than count will be extracted
    """
    genes_count = pd.DataFrame()
    for f in coverFiles:
        df = pd.read_csv(f,sep='\t',header=0,index_col=0)  # column: percent. row: genes.
        df['sum'] = df.sum(axis=1)
        genes_count[f[:3]+appendix] = df['sum']
    low_genes_df = genes_count[genes_count.min(axis=1)<count]
    print low_genes_df.shape
    genes = low_genes_df.index.tolist()
    return genes_count,low_genes_df,genes

def getMeanMedianStd(coverFiles,appendix,drop_genes):
    """
    This function calculates the mean value, median, std of the files provided. Each file should be
    in the format: row is gene, column is percentile, values are count of ribo footprint.
    
    * coverFiles:list. A list of file of replicates
    * appendix: str. Name which will be append to filename in coverFiles, and then be the column name
    * drop_genes:list. A list of genes that will be excluded from the dataframe
    """
    mean = pd.DataFrame()
    median = pd.DataFrame()
    std = pd.DataFrame()
    for f in coverFiles:
        df = pd.read_csv(f,sep='\t',header=0,index_col=0)
        df = df.drop(drop_genes)
        df = df.transpose()
        df = df/df.sum()
        mean[f[:3]+appendix+'_mean'] = df.mean(axis=1)
        median[f[:3]+appendix+'_median'] = df.median(axis=1)
        std[f[:3]+appendix+'_std'] = df.std(axis=1)
    return [mean,median,std]
#===============================================================================
#             plot the percentage coverage
#===============================================================================
"""
path = '/data/shangzhong/RibosomeProfiling/Ribo_align'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('sort.5endCov_sp_percent5end.txt')]
coverFiles = natsorted(coverFiles)
# =====1. analysis sp genes
# 1). get counts for genes with lower than 100 count in day3 and day6
min_count = 100
day3_genes_count,day3_low_genes_df,day3_low_genes = getLowCountGenes(coverFiles[0:3],'_sp_count',min_count)
day6_genes_count,day6_low_genes_df,day6_low_genes = getLowCountGenes(coverFiles[3:6],'_sp_count',min_count)
# 2). get percent coverage for day3 sample with sp
[mean3,median3,std3] = getMeanMedianStd(coverFiles[0:3],'_sp',day3_low_genes)
[mean6,median6,std6] = getMeanMedianStd(coverFiles[3:6],'_sp',day6_low_genes)
res_sp_df = pd.DataFrame()
res_sp_df['day3_sp_mean'] = mean3.mean(axis=1)
res_sp_df['day3_sp_std'] = ((std3.iloc[:,0]**2+std3.iloc[:,1]**2+std3.iloc[:,2]**2)/3).apply(np.sqrt)
res_sp_df['day6_sp_mean'] = mean6.mean(axis=1)
res_sp_df['day6_sp_std'] = ((std6.iloc[:,0]**2+std6.iloc[:,1]**2+std6.iloc[:,2]**2)/3).apply(np.sqrt)

diff_df = pd.concat([day3_low_genes_df,day6_low_genes_df],axis=1)
diff_df = diff_df[diff_df.isnull().any(axis=1)]
diff_df.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/11_sp_Diff3_6_low_genes.csv',sep='\t')
diff_df = pd.concat([day3_genes_count,day6_genes_count],axis=1)
diff_df.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/11_sp_day3_6_genes.csv',sep='\t')
overlap_genes = list(set(day3_low_genes).intersection(day6_low_genes))
print len(overlap_genes)
# 3) plot
res_sp_df.plot(title='sp genes')

# ===== 2. analysis no sp genes
coverFiles = [f for f in os.listdir(path) if f.endswith('sort.5endCov_no_sp_percent5end.txt')]
coverFiles = natsorted(coverFiles)
day3_genes_count,day3_low_genes_df,day3_low_genes = getLowCountGenes(coverFiles[0:3],'_no_sp_count',min_count)
day6_genes_count,day6_low_genes_df,day6_low_genes = getLowCountGenes(coverFiles[3:6],'_no_sp_count',min_count)
# 2). get percent coverage for day3 sample with sp
[mean3,median3,std3] = getMeanMedianStd(coverFiles[0:3],'_sp',day3_low_genes)
[mean6,median6,std6] = getMeanMedianStd(coverFiles[3:6],'_sp',day6_low_genes)
res_sp_df = pd.DataFrame()
res_sp_df['day3_no_sp_mean'] = mean3.mean(axis=1)
res_sp_df['day3_no_sp_std'] = ((std3.iloc[:,0]**2+std3.iloc[:,1]**2+std3.iloc[:,2]**2)/3).apply(np.sqrt)
res_sp_df['day6_no_sp_mean'] = mean6.mean(axis=1)
res_sp_df['day6_no_sp_std'] = ((std6.iloc[:,0]**2+std6.iloc[:,1]**2+std6.iloc[:,2]**2)/3).apply(np.sqrt)

diff_df = pd.concat([day3_low_genes_df,day6_low_genes_df],axis=1)
diff_df = diff_df[diff_df.isnull().any(axis=1)]
diff_df.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/11_no_sp_Diff3_6_low_genes.csv',sep='\t')
diff_df = pd.concat([day3_genes_count,day6_genes_count],axis=1)
diff_df.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/11_no_sp_day3_6_genes.csv',sep='\t')
overlap_genes = list(set(day3_low_genes).intersection(day6_low_genes))
print len(overlap_genes)
# 3) plot
res_sp_df.plot(title='no sp genes')
plt.show()
"""

#===============================================================================
#         plot sp length distribution
#===============================================================================
# fig,axs = plt.subplots(nrows=1,ncols=2)
# filename = '/data/shangzhong/RibosomeProfiling/cho_pr/03_choPrRefseq_sp.gene.txt'
# df = pd.read_csv(filename,sep='\t',header=0)
# sp_df = df[df['?']=='Y']
# max_len = sp_df['pos.1'].max()
# min_len = sp_df['pos.1'].min()
# ax = sp_df['pos.1'].plot(ax=axs[0],kind='hist',bins=max_len-min_len+1,title='length of signal peptide')
# ax.set_xlabel('length')
# ax.set_ylabel('number of proteins')
# 
# gene_sp_df = sp_df[['GeneID','pos.1']].drop_duplicates()
# gene_sp_df.shape
# max_len = gene_sp_df['pos.1'].max()
# min_len = gene_sp_df['pos.1'].min()
# ax1 = gene_sp_df['pos.1'].plot(ax=axs[1],kind='hist',bins=max_len-min_len+1,title='length of signal peptide')
# ax1.set_xlabel('length')
# ax1.set_ylabel('number of genes')
# fig.savefig('/data/shangzhong/RibosomeProfiling/figures/04_sp_length_distribution.pdf')
# plt.show()



#===============================================================================
#         plot coverage around end of signal peptide on nt level
#===============================================================================
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
# 2. read position file
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/03_SpEnd_CodonCov'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('5endNearSPendCov.txt')]
coverFiles = natsorted(coverFiles)
# 3. plot
calType = 'total'
up = 10; down=86
window_up = 10; window_down=86
x = np.array(range(-window_up,window_down))
colors = ['#FF4500','#008B8B']
ylabel = 'total count (RPM)'

day3_sp_end = groupCovForPosRep(coverFiles[0:3],totalCount[0:3],up,down,calType)
day6_sp_end = groupCovForPosRep(coverFiles[3:6],totalCount[3:6],up,down,calType)

# chose position ranges to plot
day3_sp_end = day3_sp_end.loc[str(0):str(window_up + window_down-1)]  #### this is used for generating subset of the plot
day6_sp_end = day6_sp_end.loc[str(0):str(window_up + window_down-1)]
sp_end = pd.DataFrame();sp_end_err = pd.DataFrame()
sp_end['day3'] = day3_sp_end['mean']
sp_end['day6'] = day6_sp_end['mean']
sp_end_err['day3'] = day3_sp_end['std']
sp_end_err['day6'] = day6_sp_end['std']
# plot the end site
f, ax = plt.subplots(2, sharex=True)
for i in range(2):
    ax[i].bar(x,sp_end.iloc[:,i],color=colors[i],align='center',alpha=0.69,yerr=[tuple([0]*(sp_end.shape[0])),tuple(sp_end_err.iloc[:,i].values)])
    #ax[i].set_ylim([0,4])
    ax[i].set_xlim([-window_up,window_down])
    ax[i].set_title(sp_end.columns.values[i])
ax[1].set_xlabel('distance from end of signal peptide')
#ax[0].set_ylim([0,0.6]);ax[1].set_ylim([0,0.8])
f.text(0.06, 0.5, ylabel, ha='center', va='center', rotation='vertical',fontsize=12)
plt.xticks(x,range(-window_up,window_down),rotation=90)
plt.suptitle('coverage near end of signal peptide',fontsize=14,color='black')
# # define ylim
# ax[0].set_ylim([0,0.5]);ax[1].set_ylim([0,1]);ax[2].set_ylim([0,0.5]);ax[3].set_ylim([0,0.4])

# ax[0].set_title('day3_Antibody');ax[1].set_title('day3_NeoR');ax[2].set_title('day6_Antibody');ax[3].set_title('day6_NeoR')
# plt.suptitle('Stop sites of antibody')
#plt.savefig('/data/shangzhong/RibosomeProfiling/figures/05_CodonCovSpEnd_total.svg',dpi=300)
plt.show()
"""
#===============================================================================
#             plot clustering at TSS and TSE sites
#===============================================================================
def mergeRepNormPosCov(files,gene_sp,totalCount,up,down,center):
    """
    This function first normalize read count at each position for each file, then extract the interested
    window of sequences, at last take average of all samples.
    
    * files: list. A list of file names.
    * gene_sp: list. A list of genes with signal peptide.
    * totalCount: list. A list of total count for each sample.
    * up: int. # of nucleotides upstream of center site.
    * down: int. # of nts downstream of center site.
    * center: int. center position
    """
    df_list = [];gene_list = []
    for f,total in zip(files,totalCount):
        df = normGeneCovWindow(f,gene_sp,total,up,down,center)
        genes  = df.index.tolist()
        df_list.append(df);gene_list.append(genes)
    gene_intersect = list(set.intersection(*map(set,gene_list)))
    for i in range(len(df_list)):
        df_list[i] = df_list[i].loc[gene_intersect]
    mean_df = sum(df_list)/len(df_list)
    
    return mean_df

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
# 2) read file
up = 50;down=50;center=0
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/04_allGenePos_cov'
os.chdir(path)
posCovFile = [f for f in os.listdir(path) if f.endswith('geneAllposCov.txt')]
# 3).sp genes and no sp genes
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
df = pd.read_csv(gene_class_file,sep='\t',header=0)
gene_sp = df['gene_sp'].dropna().astype(str).tolist()
# 4) merge replicates
day3_df = mergeRepNormPosCov(posCovFile[0:3],gene_sp,totalCount[0:3],up,down,center)
day6_df = mergeRepNormPosCov(posCovFile[3:6],gene_sp,totalCount[3:6],up,down,center)
day3_df.to_csv('/data/shangzhong/RibosomeProfiling/Ribo_align/05_TSS_TSE_heatmap/day3_TSS.txt',sep='\t')
day6_df.to_csv('/data/shangzhong/RibosomeProfiling/Ribo_align/05_TSS_TSE_heatmap/day6_TSS.txt',sep='\t')
"""

#===============================================================================
#             plot translation efficiency for each gene
#===============================================================================
"""
# 0). get the genes with sp and without sp
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/05_cho_ab_pr_sp.gene.classify.txt'
gene_class_df = pd.read_csv(gene_class_file,header=0,sep='\t')
gene_class_df = gene_class_df.astype(str)
gene_sp = gene_class_df['with_sp'].tolist()
gene_sp = list(set(gene_sp))
gene_no_sp = gene_class_df['no_sp'].tolist()
gene_no_sp = list(set(gene_no_sp))
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
cri = (res_df['day3_ribo']>0) & (res_df['day6_ribo']>0) & (res_df['day3_rna']>0) & (res_df['day6_rna']>0)
res_df = res_df[cri]
res_df['ribo_change'] = (res_df['day6_ribo']/res_df['day3_ribo']).apply(np.log2)
res_df['mrna_change'] = (res_df['day6_rna']/res_df['day3_rna']).apply(np.log2)
# 4. plot
plt.plot(res_df['mrna_change'],res_df['ribo_change'],'.',color='#6699CC')

# gene_nosp_df = res_df[res_df.index.isin(gene_no_sp)]
# plt.plot(gene_nosp_df['ribo_change'],gene_nosp_df['mrna_change'],'.')
# 
# gene_sp_df = res_df[res_df.index.isin(gene_sp)]
# plt.plot(gene_sp_df['ribo_change'],gene_sp_df['mrna_change'],'.',color='b',alpha=1)

heavy_df = res_df[res_df.index.isin(['heavychain'])]
plt.plot(heavy_df['mrna_change'],heavy_df['ribo_change'],'.',color='red',alpha=1)

light_df = res_df[res_df.index.isin(['lightchain'])]
plt.plot(light_df['mrna_change'],light_df['ribo_change'],'.',color='yellow',alpha=1)

neo_df = res_df[res_df.index.isin(['NeoRKanR'])]
plt.plot(neo_df['mrna_change'],neo_df['ribo_change'],'.',color='black',alpha=1)

plt.xlabel('log2 mRNA change')
plt.ylabel('log2 rpf change')
plt.title('Translation change Day6 VS Day3')
plt.axhline(y=0,color = 'k')
plt.axvline(x=0,color='k')
plt.show()
"""
#===============================================================================
#             plot cumulative codon fractions
#===============================================================================
# This is on the codon level
"""
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/05_codon_cov'
os.chdir(path)
codonCovFiles = [f for f in os.listdir(path) if f.endswith('.txt')]
codonCovFiles = natsorted(codonCovFiles)

AllPos = []
for f in codonCovFiles:
    handle = open(f,'r')
    for line in handle:
        norm = normCodonCov(line,by='median')
        if norm == '': continue
        pos = [np.log2(float(p)+0.0001) for p in norm[15:-10]]
#         pos_norm = []
#         for p in pos:
#             if p == 0.0: pos_norm.append(0.0001)
#             else: pos_norm.append(p)
#         pos_norm = [p + 0.0001 for p in pos if p == 0.0]
#         pos = list(np.log2(pos_norm))
        AllPos.extend(pos)
"""
# plt.hist(AllPos,bins=1000,normed=True,alpha=0.7,color='#4169E1',log=True)
# #plt.xlim([-10,15])
# #plt.ylim([0,1])
# plt.axvline(x=np.log10(25))
# plt.xlabel('Ribosome Footprints VS median')
# plt.ylabel('Frequency of codon')
# plt.title('Stalling sites')
"""
plt.hist(AllPos, bins=800000, normed=True,histtype='stepfilled', cumulative=-1,alpha=0.7,color='#4169E1')
plt.axvline(x=25)
plt.yscale('log')
plt.xlabel('Ribosome Footprints VS median')
plt.ylabel('Cumulative Fraction of codon')
plt.title('Stalling sites')
plt.xlim([0,100])
plt.text(27,0.1,'ribosomal pausing',fontsize=12,color='#4169E1')
plt.arrow(25,0.05,30,0,shape='full',head_width=0.02,head_length=3,fc='#4169E1',ec='#4169E1')
#plt.savefig('/data/shangzhong/RibosomeProfiling/figures/10_CodonPausingCumulativePlot.png',dpi=300)
"""
# plt.show()


"""
f,ax = plt.subplots(2,sharex=True)
ax[0].bar(range(len(day3_nt_dic)),day3_nt_dic.values(),align='center',color=c_colors)
ax[0].set_ylabel('day3',color='black',fontsize=16)
ax[0].set_title('codon frequencies in day3 and day6',fontsize=24)
ax[1].bar(range(len(day6_nt_dic)),day6_nt_dic.values(),align='center',color=c_colors)
ax[1].set_ylabel('day6',color='black',fontsize=16)
plt.xticks(range(len(day6_nt_dic)),day6_nt_dic.keys(),rotation=90)
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color=colors[0], label='significant')
blue_patch = mpatches.Patch(color=colors[1], label='non significant')
ax[0].legend(handles=[red_patch,blue_patch])

f,ax = plt.subplots(2,sharex=True)
ax[0].bar(range(len(day3_aa_dic)),day3_aa_dic.values(),align='center',color=a_colors)
ax[0].set_ylabel('day3',color='black',fontsize=16)
ax[0].set_title('AA frequencies in day3 and day6',fontsize=24)
ax[1].bar(range(len(day6_aa_dic)),day6_aa_dic.values(),align='center',color=a_colors)
ax[1].set_ylabel('day6',color='black',fontsize=16)
plt.xticks(range(len(day6_aa_dic)),day6_aa_dic.keys())
red_patch = mpatches.Patch(color=colors[0], label='significant')
blue_patch = mpatches.Patch(color=colors[1], label='non significant')
ax[0].legend(handles=[red_patch,blue_patch])
plt.show()
"""
#===============================================================================
#             plot cds, intron,5utr,3utr
#===============================================================================
"""
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/12_cds_utr_intron'
os.chdir(path)
ribo_class_files = [f for f in os.listdir(path) if f.endswith('.txt')]
ribo_class_files = natsorted(ribo_class_files)
class_count_table = pd.DataFrame()
for f in ribo_class_files:
    df = pd.read_csv(f,sep='\t',header=0,index_col=0)
    sum_count = df.sum().tolist()
    class_count_table[f[:3]] = pd.Series(sum_count)
class_count_table.index = pd.Series(['cds','5UTR','3UTR','intron'])
total = class_count_table.sum().tolist()
class_count_table = class_count_table.div(total)
print 'done'
class_count_table['day3'] = class_count_table.iloc[:,0:3].mean(axis=1)
class_count_table['day6'] = class_count_table.iloc[:,3:6].mean(axis=1)
error = pd.DataFrame()
error['day3'] = class_count_table.iloc[:,0:3].std(axis=1)
error['day6'] = class_count_table.iloc[:,3:6].std(axis=1)

ax = class_count_table.loc[:,['day3','day6']].plot(kind='bar',align='center',yerr=error,alpha=0.69)
ax.set_ylabel('percentage')
ax.set_xlabel('regions in genes')
ax.set_title('Ribo seq distribution among different regions of genes')
plt.savefig('/data/shangzhong/RibosomeProfiling/figures/13_ribo_distribution_cds_utr_intron.svg')
plt.show()
"""
