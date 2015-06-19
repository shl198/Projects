from __future__ import division
import subprocess,os,sys
import pandas as pd
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from f02_RiboDataModule import *
mpl.style.use('ggplot')
sys.path.append('/data/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer



# 2. read cds file
cdsPosFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_combined.CDS.bed'
cds_df = pd.read_csv(cdsPosFile,header=None,sep='\t',names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
cds_df['GeneID'] = cds_df['GeneID'].astype(str)
cds_df['length'] = cds_df['cds_end'] - cds_df['cds_start'] 
# 3. get gene ids
geneIDs = cds_df['GeneID'].tolist()
geneIDs = list(set(geneIDs))
# 4. get the genes with sp and without sp
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/05_choPrRefseq_sp.gene.sp_classify.txt'
gene_class_df = pd.read_csv(gene_class_file,header=0,sep='\t')
gene_class_df = gene_class_df.astype(str)
gene_sp = gene_class_df['with_sp'].tolist()
gene_sp = list(set(gene_sp))
gene_no_sp = gene_class_df['no_sp'].tolist()
gene_no_sp = list(set(gene_no_sp))
try:
    gene_sp.remove('-');gene_no_sp.remove('-')
except:
    pass
print '# of genes with sp:',len(gene_sp)
print '# of genes without ps:',len(gene_no_sp)
# 5. exclude the genes with multiple scaffolds
gene_multi_file = '/data/shangzhong/RibosomeProfiling/cho_pr/04_choPrRefseq_sp.gene.multichr.txt'
gene_multi_df = pd.read_csv(gene_multi_file,header=None,sep='\t',engine='python',names=['GeneID','chr1','chr2','chr3','chr4','chr5'])
gene_multi_df['GeneID'] = gene_multi_df['GeneID'].astype(str)
gene_multi_list = gene_multi_df['GeneID'].tolist()
gene_multi_list = list(set(gene_multi_list))
gene_sp = [g for g in gene_sp if g not in gene_multi_list]
gene_no_sp = [g for g in gene_no_sp if g not in gene_multi_list]
print '# of genes with sp after remove multi_chrom:',len(gene_sp)
print '# of genes without sp after remove multi_chrom:',len(gene_no_sp)
# 6. exclude the genes with overlapped CDS
gene_cds_ovlap_file = '/data/shangzhong/RibosomeProfiling/cho_pr/07_overlap_CDS.txt'
gene_cds_ovlap_df = pd.read_csv(gene_cds_ovlap_file,header=0,sep='\t',names=['Gene1','Gene2','Strand'])
gene_cds_ovlap_df = gene_cds_ovlap_df.astype(str)
gene1_list = gene_cds_ovlap_df['Gene1'].tolist()
gene2_list = gene_cds_ovlap_df['Gene2'].tolist()
gene_cds_ovlap_list = list(set(gene1_list + gene2_list))
print '# of overlapped genes:',len(gene_cds_ovlap_list)
gene_sp = [g for g in gene_sp if g not in gene_cds_ovlap_list]
gene_no_sp = [g for g in gene_no_sp if g not in gene_cds_ovlap_list]
print '# of genes with sp after remove overlap genes:',len(gene_sp)
print '# of genes without sp after remove overlap genes:',len(gene_no_sp)



#===============================================================================
#             get length stats
#===============================================================================
# df = pd.DataFrame()
# res = geneCDSStats(cdsPosFile,gene_no_sp)
# df['gene_no_sp'] = pd.Series(res)
# res = geneCDSStats(cdsPosFile,gene_sp)
# df['gene_sp'] = pd.Series(res)
# df = df.fillna(0)
# df.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/10_CDS_len_stats_count.csv',sep='\t',index=False)
#===============================================================================
#     # 7.1 read coverage data
#===============================================================================
# gene_sp = gene_no_sp = ['100750439']
# coverFiles = '/data/shangzhong/RibosomeProfiling/Ribo_align/s04.sort.5endCov.txt'
"""
coverFiles = sys.argv[1]
if coverFiles.endswith('sort.txt'):
    exon_cov_df = pd.read_csv(coverFiles,sep='\t',header=None,dtype={'GeneID':str},names=['Chr','ex_start','ex_end','GeneID','None','Strand','pos','coverage'])
    #gene_sp = gene_no_sp = ['100750439','103164226']
    # 1). get percent coverage for sp genes
    gene_sp_cov_df = percent_cov(exon_cov_df,cds_df,gene_sp)
    gene_sp_cov_df.to_csv(coverFiles[:-4]+'_sp_percent.txt',sep='\t')
    print 'sp genes done'
    # 2). get percent coverage for non sp genes
    gene_no_sp_cov_df = percent_cov(exon_cov_df,cds_df,gene_no_sp)
    gene_no_sp_cov_df.to_csv(coverFiles[:-4]+'_no_sp_percent.txt')
    print 'no sp genes done'
#===============================================================================
#     7.2 read coverage data 5' end
#===============================================================================
elif coverFiles.endswith('sort.5endCov.txt'):
    exon_cov_df = pd.read_csv(coverFiles,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
    # 1) get percent coverage for sp genes
    gene_sp_cov_df = percent5endCov(exon_cov_df,cds_df,gene_sp)
    gene_sp_cov_df.to_csv(coverFiles[:-4]+'_sp_percent5end.txt',sep='\t')
    print 'sp genes done'
    # 2) get percent coverage for no sp genes
    gene_no_sp_cov_df = percent5endCov(exon_cov_df,cds_df,gene_no_sp)
    gene_no_sp_cov_df.to_csv(coverFiles[:-4]+'_no_sp_percent5end.txt',sep='\t')
    print 'no sp genes done'
"""
#===============================================================================
#         9. plot
#===============================================================================
# filepath = '/data/shangzhong/RibosomeProfiling/Ribo_align'
# os.chdir(filepath)
# cal_type = 'mean'
# coverFiles = [f for f in os.listdir(filepath) if f.endswith('sort_sp_percent.txt')]
# coverFiles = natsorted(coverFiles)
# day3_sp_df = merge_percent_cov(coverFiles[:3],'day3_sp',cal_type)
# day6_sp_df = merge_percent_cov(coverFiles[3:6],'day6_sp',cal_type)
# coverFiles = [f for f in os.listdir(filepath) if f.endswith('sort_no_sp_percent.txt')]
# coverFiles = natsorted(coverFiles)
# day3_sp_df = merge_percent_cov(coverFiles[:3],'day3__no_sp',cal_type)
# day6_no_sp_df = merge_percent_cov(coverFiles[:3],'day6_no_sp',cal_type)
# # merge
# df = pd.concat([day3_sp_df,day6_sp_df],axis=1)
# # plot
# ax = df.plot(title='Average coverage at each percentage')
# ax.set_xlabel('Gene (%)')
# ax.set_ylabel('Coverage (%)')
# fig = ax.get_figure()
# fig.savefig('/data/shangzhong/RibosomeProfiling/sp_cov.pdf')
# plt.show()ls


# filepath = '/data/shangzhong/RibosomeProfiling/Ribo_align'
# os.chdir(filepath)
# coverFiles = [f for f in os.listdir(filepath) if f.endswith('sort_sp_percent.txt')]
# coverFiles = natsorted(coverFiles)
# # define dataframe
# mean_df = pd.DataFrame()
# median_df = pd.DataFrame()
# std_df = pd.DataFrame()
# # get mean,median


#===============================================================================
#             calculate coverage around TSS and TSE sites
#===============================================================================
"""
filepath = '/data/shangzhong/RibosomeProfiling/Ribo_align'
os.chdir(filepath)
coverFiles = [f for f in os.listdir(filepath) if f.endswith('sort.5endCov.txt')]
coverFiles = natsorted(coverFiles)
res_sp = pd.DataFrame()
res_no_sp = pd.DataFrame()
TSS_up = 50;TSS_down=100
TSE_up = 50;TSE_down= 30
for f in coverFiles:
    exon_cov_df = pd.read_csv(f,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
    # gene sp
    start_sum,stop_sum = covNearTSS_TSE(exon_cov_df,cds_df,gene_sp,TSS_up,TSS_down,TSE_up,TSE_down)
    df = pd.DataFrame({f[:-16]+'start':start_sum})
    res_sp = pd.concat([res_sp,df],axis=1)
    df = pd.DataFrame({f[:-16]+'stop':stop_sum})
    res_sp = pd.concat([res_sp,df],axis=1)
    # gene no sp
    start_sum,stop_sum = covNearTSS_TSE(exon_cov_df,cds_df,gene_no_sp,TSS_up,TSS_down,TSE_up,TSE_down)
    df = pd.DataFrame({f[:-16]+'start':start_sum})
    res_no_sp = pd.concat([res_no_sp,df],axis=1)
    df = pd.DataFrame({f[:-16]+'stop':stop_sum})
    res_no_sp = pd.concat([res_no_sp,df],axis=1)
res_sp.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/08_covTSS_TSE_sp.csv',sep='\t',index=False)
res_no_sp.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/08_covTSS_TSE_no_sp.csv',sep='\t',index=False)
# plot the figures 

covTS_sp_df = pd.read_csv('/data/shangzhong/RibosomeProfiling/cho_pr/08_covTSS_TSE_sp.csv',sep='\t',header=0)
covTS_no_sp_df = pd.read_csv('/data/shangzhong/RibosomeProfiling/cho_pr/08_covTSS_TSE_no_sp.csv',sep='\t',header=0)
up = 50; down = 100
x = range(-up,down)
# plot sp
df = covTS_sp_df[['s04.start','s05.start','s06.start','s10.start','s11.start','s12.start']]
ax = df.plot(kind='bar',subplots=True,title='Start Sites')
ax[5].set_xticklabels(x)
ax[5].set_xlabel('distance from TSS')
ax[2].set_ylabel('total count')
fig = ax[1].get_figure()
fig.savefig('/data/shangzhong/RibosomeProfiling/cho_pr/covTSS_sp.pdf')

df = covTS_sp_df[['s04.stop','s05.stop','s06.stop','s10.stop','s11.stop','s12.stop']]
ax = df.plot(kind='bar',subplots=True,title='Stop Sites')
ax[5].set_xticklabels(x)
ax[5].set_xlabel('distance from TSE')
ax[2].set_ylabel('total count')
fig = ax[1].get_figure()
fig.savefig('/data/shangzhong/RibosomeProfiling/cho_pr/covTSE_sp.pdf',format='pdf')
# plot no sp
df = covTS_no_sp_df[['s04.start','s05.start','s06.start','s10.start','s11.start','s12.start']]
ax = df.plot(kind='bar',subplots=True,title='Start Sites')
ax[5].set_xticklabels(x)
ax[5].set_xlabel('distance from TSS')
ax[2].set_ylabel('total count')
fig = ax[1].get_figure()
fig.savefig('/data/shangzhong/RibosomeProfiling/cho_pr/covTSS_no_sp.pdf')

df = covTS_sp_df[['s04.stop','s05.stop','s06.stop','s10.stop','s11.stop','s12.stop']]
ax = df.plot(kind='bar',subplots=True,title='Stop Sites')
ax[5].set_xticklabels(x)
ax[5].set_xlabel('distance from TSE')
ax[2].set_ylabel('total count')
fig = ax[1].get_figure()
fig.savefig('/data/shangzhong/RibosomeProfiling/cho_pr/covTSE_no_sp.pdf')
plt.show()
"""
    












