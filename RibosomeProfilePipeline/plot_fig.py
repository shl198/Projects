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
#                         plot coverage around TSS and TSE
#===============================================================================
def groupCovForPosRep(files,totalCounts,up,down,calType='total'):
    """
    This function groups the replicates of samples' coverage around TSS and TSE sites.
    return two columns: ['mean'/'total'/'median'/'geometricmean,'std']
    
    * files: list. A list of replicate files
    * totalCounts: list. A list of int with same length of files. Stores the total mapped reads in files.
    * calType: str. Calculation type.
    """
    res_df = pd.DataFrame()
    mean = pd.DataFrame()
    std = pd.DataFrame()
    for f,total in zip(files,totalCounts):
        df = pd.read_csv(f,sep='\t',header=0,index_col=0,low_memory=False)
        # remove the antibody
        try:
            df = df.drop('heavychain');df=df.drop('lightchain')
            #df = df.loc[['heavychain','lightchain']]
        except:
            df = df.drop('NeoRKanR')
            #df = df.loc[['NeoRKanR']]
            
        # filter by median of gene
        df['median'] = df.median(axis=1)
        #df = df[df['median']>=0.5]
        del df['median']
        up_df = df.iloc[:,0:up]
        down_df = df.iloc[:,up:up+down]
        try:
            up_df = up_df.replace('-',np.nan).dropna().astype('int').T
            down_df = down_df.replace('-',np.nan).dropna().astype('int').T
        except:
            pass
        up_df = up_df/total*(10**6)   # row: position. col: gene
        down_df = down_df/total*(10**6)
        if calType == 'total':
            up_df[calType] = up_df.sum(axis=1);down_df[calType] = down_df.sum(axis=1)
        if calType == 'mean':
            up_df[calType] = up_df.mean(axis=1);down_df[calType] = down_df.mean(axis=1)
        if calType == 'median':
            up_df[calType] = up_df.median(axis=1);down_df[calType] = down_df.median(axis=1)
        if calType == 'geoMean':
            up_df = up_df.replace([0],[1]);down_df=down_df.replace([0],[1])
            up_df[calType] = sp_stats.gmean(up_df.values,axis=1)
            down_df[calType] = sp_stats.gmean(down_df.values,axis=1)
        df = df.T
        df[calType]= pd.concat([up_df[calType],down_df[calType]])
        up_df['std']=up_df.std(axis=1);down_df['std']=down_df.std(axis=1)
        df['std'] = pd.concat([up_df['std'],down_df['std']])
        mean[f+calType] = df[calType]
        std[f+'std'] = df['std']
    res_df['mean'] = mean.mean(axis=1)
    res_df['std'] = ((std**2).sum(axis=1)/std.shape[1]).apply(np.sqrt)
    return res_df
"""
# read files
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_cov5end'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
coverFiles = natsorted(coverFiles)
totalCount = []
up = 50; down=60
for f in coverFiles:
    df = pd.read_csv(f,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
    total = df['coverage'].sum()
    totalCount.append(total)
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/35/02_TSS_TSE_cov'
os.chdir(path)
start_sp_file = [f for f in os.listdir(path) if f.startswith('start_sp_')]; start_sp_file = natsorted(start_sp_file)
end_sp_file = [f for f in os.listdir(path) if f.startswith('end_sp_')]; end_sp_file = natsorted(end_sp_file)
start_nosp_file = [f for f in os.listdir(path) if f.startswith('start_nosp_')]; start_nosp_file = natsorted(start_nosp_file)
end_nosp_file = [f for f in os.listdir(path) if f.startswith('end_nosp_')]; end_nosp_file = natsorted(end_nosp_file)
# calculate coverage at each position
calType = 'total'
ylabel = 'total count (RPM)'
day3_TSS_sp = groupCovForPosRep(start_sp_file[0:3],totalCount[0:3],up,down,calType)
day3_TSE_sp = groupCovForPosRep(end_sp_file[0:3],totalCount[0:3],up,down,calType)
day6_TSS_sp = groupCovForPosRep(start_sp_file[3:6],totalCount[3:6],up,down,calType)
day6_TSE_sp = groupCovForPosRep(end_sp_file[3:6],totalCount[3:6],up,down,calType)
day3_TSS_nosp = groupCovForPosRep(start_nosp_file[0:3],totalCount[0:3],up,down,calType)
day3_TSE_nosp = groupCovForPosRep(end_nosp_file[0:3],totalCount[0:3],up,down,calType)
day6_TSS_nosp = groupCovForPosRep(start_nosp_file[3:6],totalCount[3:6],up,down,calType)
day6_TSE_nosp = groupCovForPosRep(end_nosp_file[3:6],totalCount[3:6],up,down,calType)

TSS = pd.DataFrame();TSE=pd.DataFrame();TSS_err=pd.DataFrame();TSE_err=pd.DataFrame()
TSS['day3_sp']=day3_TSS_sp['mean'];TSS['day3_nosp']=day3_TSS_nosp['mean']
TSS['day6_sp']=day6_TSS_sp['mean'];TSS['day6_nosp']=day6_TSS_nosp['mean']
TSS_err['day3_sp']=day3_TSS_sp['std'];TSS_err['day3_nosp']=day3_TSS_nosp['std']
TSS_err['day6_sp']=day6_TSS_sp['std'];TSS_err['day6_nosp']=day6_TSS_nosp['std']

TSE['day3_sp']=day3_TSE_sp['mean'];TSE['day3_nosp']=day3_TSE_nosp['mean']
TSE['day6_sp']=day6_TSE_sp['mean'];TSE['day6_nosp']=day6_TSE_nosp['mean']
TSE_err['day3_sp']=day3_TSE_sp['std'];TSE_err['day3_nosp']=day3_TSE_nosp['std']
TSE_err['day6_sp']=day6_TSE_sp['std'];TSE_err['day6_nosp']=day6_TSE_nosp['std']

#==================== plot start site =============================
x = np.array(range(-up,down))
f, ax = plt.subplots(4, sharex=True)
colors = ['#FF4500','#008B8B','#6495ED','#808080']

for i in range(4):
    ax[i].bar(x,TSS.iloc[:,i],color=colors[i],align='center')#,yerr=[tuple([0]*(TSS.shape[0])),tuple(TSS_err.iloc[:,i].values)])
    #ax[i].set_ylim([0,4])
    ax[i].set_xlim([-50,60])
    ax[i].set_title(TSS.columns.values[i])
ax[3].set_xlabel('distance from TSS')
f.text(0.06, 0.5, ylabel, ha='center', va='center', rotation='vertical',fontsize=12)
plt.xticks(x,range(-up,down),rotation=90 )
plt.suptitle('Start sites without antibody (35 nt)',fontsize=12)
# # define ylim
# ax[0].set_ylim([0,1]);ax[1].set_ylim([0,1.5]);ax[2].set_ylim([0,1]);ax[3].set_ylim([0,1])

# ax[0].set_title('day3_Antibody');ax[1].set_title('day3_NeoR');ax[2].set_title('day6_Antibody');ax[3].set_title('day6_NeoR')
# plt.suptitle('Start sites of antibody')
# =================== plot stop site ==============================
f, ax = plt.subplots(4, sharex=True)
for i in range(4):
    ax[i].bar(x,TSE.iloc[:,i],color=colors[i],align='center')#,yerr=[tuple([0]*(TSE.shape[0])),tuple(TSE_err.iloc[:,i].values)])
    #ax[i].set_ylim([0,4])
    ax[i].set_xlim([-50,60])
    ax[i].set_title(TSE.columns.values[i])
ax[3].set_xlabel('distance from TSE')
f.text(0.06, 0.5, ylabel, ha='center', va='center', rotation='vertical',fontsize=12)
plt.xticks(x,range(-up,down),rotation=90)
plt.suptitle('Stop sites without antibody (35 nt)',fontsize=12)
# # define ylim
# ax[0].set_ylim([0,0.5]);ax[1].set_ylim([0,1]);ax[2].set_ylim([0,0.5]);ax[3].set_ylim([0,0.4])

# ax[0].set_title('day3_Antibody');ax[1].set_title('day3_NeoR');ax[2].set_title('day6_Antibody');ax[3].set_title('day6_NeoR')
# plt.suptitle('Stop sites of antibody')
plt.show()
"""
#===============================================================================
#         plot coverage around end of signal peptide
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
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/03_SpEnd_cov'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('5endNearSPendCov.txt')]
coverFiles = natsorted(coverFiles)
# 3. plot
calType = 'mean'
up=50;down=50
x = np.array(range(-up,down))
colors = ['#FF4500','#008B8B']
ylabel = 'mean count (RPM)'

day3_sp_end = groupCovForRep(coverFiles[0:3],totalCount[0:3],up,down,calType)
day6_sp_end = groupCovForRep(coverFiles[3:6],totalCount[3:6],up,down,calType)
sp_end = pd.DataFrame();sp_end_err = pd.DataFrame()
sp_end['day3'] = day3_sp_end['mean']
sp_end['day6'] = day6_sp_end['mean']
sp_end_err['day3'] = day3_sp_end['std']
sp_end_err['day6'] = day6_sp_end['std']
# plot the end site
f, ax = plt.subplots(2, sharex=True)
for i in range(2):
    ax[i].bar(x,sp_end.iloc[:,i],color=colors[i],align='center',yerr=[tuple([0]*(sp_end.shape[0])),tuple(sp_end_err.iloc[:,i].values)])
    #ax[i].set_ylim([0,4])
    ax[i].set_xlim([-up,down])
    ax[i].set_title(sp_end.columns.values[i])
ax[1].set_xlabel('distance from end of signal peptide')
ax[0].set_ylim([0,0.6]);ax[1].set_ylim([0,0.8])
f.text(0.06, 0.5, ylabel, ha='center', va='center', rotation='vertical',fontsize=12)
plt.xticks(x,range(-up,down),rotation=90)
plt.suptitle('coverage near end of signal peptide')
# # define ylim
# ax[0].set_ylim([0,0.5]);ax[1].set_ylim([0,1]);ax[2].set_ylim([0,0.5]);ax[3].set_ylim([0,0.4])

# ax[0].set_title('day3_Antibody');ax[1].set_title('day3_NeoR');ax[2].set_title('day6_Antibody');ax[3].set_title('day6_NeoR')
# plt.suptitle('Stop sites of antibody')
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
res_df = pd.DataFrame()
# 1. read ribo count
riboFile = '/data/shangzhong/RibosomeProfiling/cho_pr/12_pr_gene_rpkm.txt'
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
plt.plot(res_df['ribo_change'],res_df['mrna_change'],'.')
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
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/09_codon_cov'
os.chdir(path)
codonCovFiles = [f for f in os.listdir(path) if f.endswith('.txt')]
codonCovFiles = natsorted(codonCovFiles)

AllPos = []
for f in codonCovFiles:
    handle = open(f,'r')
    for line in handle:
        item = line[:-1].split('\t')
        posCov = item[2:]
        posCov = [float(p) for p in posCov]
        if len(posCov)<25: continue
        if np.median(np.array(posCov))<1: continue
        pos = codonStallSites(posCov,by='median')
        if pos == 0:
            continue
        AllPos.extend(pos)

plt.hist(AllPos, bins=800000, normed=True,histtype='stepfilled', cumulative=-1,alpha=0.7,color='#4169E1')
plt.axvline(x=25)
plt.yscale('log')
plt.xlabel('Ribosome Footprints VS median')
plt.ylabel('Cumulative Fraction of codon')
plt.title('codon level')
plt.xlim([0,100])
plt.text(27,0.1,'ribosomal pausing',fontsize=12,color='#4169E1')
plt.arrow(25,0.05,30,0,shape='full',head_width=0.02,head_length=3,fc='#4169E1',ec='#4169E1')
plt.show()
"""

#===============================================================================
#                         plot AA and codon frequencys
#===============================================================================
"""
freqFile1 = '/data/shangzhong/RibosomeProfiling/cho_pr/17_day3_codon_frequency.txt'
freqFile2 = '/data/shangzhong/RibosomeProfiling/cho_pr/17_day6_codon_frequency.txt'
day3_dic,day3_aa_dic = getcodon_AAFreq(freqFile1)
day6_dic,day6_aa_dic = getcodon_AAFreq(freqFile2)
del day6_dic['TNN']; del day6_dic['NNN']
del day6_aa_dic['X']
print day3_dic;print day6_dic
print day3_aa_dic;print day6_aa_dic
print len(day3_dic);print len(day6_dic)
print len(day3_aa_dic);print len(day6_aa_dic)

c_color = ['AAA','AAC','ACG','AGA','CCG','CGG','GCA','GTC','TGG']
c_colors = []
for key in day6_dic.keys():
    if key in c_color:
        c_colors.append('r')
    else:
        c_colors.append('b')
a_color = 'W'
a_colors = []

for key in day6_aa_dic.keys():
    if key == 'W':
        a_colors.append('r')
    else:
        a_colors.append('b')

f,ax = plt.subplots(2,sharex=True)
ax[0].bar(range(len(day3_dic)),day3_dic.values(),align='center',color=c_colors)
ax[0].set_ylabel('day3')
ax[0].set_title('codon frequencies in day3 and day6')
ax[1].bar(range(len(day6_dic)),day6_dic.values(),align='center',color=c_colors)
ax[1].set_ylabel('day6')
plt.xticks(range(len(day6_dic)),day6_dic.keys(),rotation=90)
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='red', label='significant')
blue_patch = mpatches.Patch(color='b', label='non significant')
ax[0].legend(handles=[red_patch,blue_patch])

f,ax = plt.subplots(2,sharex=True)
ax[0].bar(range(len(day3_aa_dic)),day3_aa_dic.values(),align='center',color=a_colors)
ax[0].set_ylabel('day3')
ax[0].set_title('AA frequencies in day3 and day6')
ax[1].bar(range(len(day6_aa_dic)),day6_aa_dic.values(),align='center',color=a_colors)
ax[1].set_ylabel('day6')
plt.xticks(range(len(day6_aa_dic)),day6_aa_dic.keys())

red_patch = mpatches.Patch(color='red', label='significant')
blue_patch = mpatches.Patch(color='b', label='non significant')
ax[0].legend(handles=[red_patch,blue_patch])
plt.show()
"""










