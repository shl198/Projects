import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('ggplot')
import os
import pandas as pd
from natsort import natsorted

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
    This function calculates the meean value, median, std of the files provided. Each file should be
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
#plt.show()


#===============================================================================
#             plot translational efficiency
#===============================================================================
ribo_file = '/data/shangzhong/
ribo_df = pd.read_csv()






print 'done'









