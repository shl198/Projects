import pysam
import os,subprocess,sys
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from natsort import natsorted
from f02_RiboDataModule import *
from multiprocessing import Process,Pool
from Bio import SeqIO
from scipy import stats

import pdb
#=============== parameters =======================
#------ need to define manually ----------
ribo_offset_file = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam/02_TSS_TSE_cov/ribo_offset.txt'
bam_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam'
db_path = '/data/shangzhong/RibosomeProfiling/Database'
rna_bam_path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align'
#----- don't need to define manually -----
general_path = bam_path + '/general_files'
fig_path = bam_path + '/figures'

fwd_rev_path = bam_path + '/01_cov'
pr_cov_path = bam_path + '/03_pr_pos_cov'    # store results from part 1
gene_count_path = bam_path + '/04_gene_total_count' 
stall_path = bam_path + '/05_stall_sites'
codon_AA_freq_path = bam_path + '/06_codon_AA_freq'
utr3_cov_path = bam_path + '/07_utr3_cov'
pr_ntcov_path = bam_path + '/08_pr_ntpos_cov'
frame_cov_path = bam_path + '/09_frame_cov'

all_id_file = db_path + '/combined_AllIDS.txt'  # part 8
cdsFile = db_path + '/01_pr_cds.txt' # part 8
exnFile = db_path + '/01_pr_rna.txt' # part 8
# rna
rna_count_path = rna_bam_path + '/01_gene_count'

def chunk(l,n):
    n = max(1,n)
    res = [l[i:i+n] for i in range(0,len(l),n)]
    return res
#===============================================================================
#                     1. get position coverage for each A site for all protiens. Results stored in folder 03_pr_pos_cov
#===============================================================================
# os.chdir(fwd_rev_path)
# covFiles = [f for f in os.listdir(fwd_rev_path) if f.endswith('cov.txt')]
# covFiles = natsorted(covFiles)
# cdsFile = db_path + '/01_pr_cds.txt'
# #gene_pr_cov(covFiles[1],cdsFile,ribo_offset_file,pr_cov_path+'/old')
# #------------- pr_pos_cov ------------
# proc = [Process(target=gene_pr_cov,args=(covFile,cdsFile,ribo_offset_file,pr_cov_path,)) for covFile in covFiles]
# for p in proc:
#     p.start()
# # for p in proc:
# #     p.join()
# #------------- pr_ntpos_cov ------------
# proc = [Process(target=gene_pr_cov,args=(covFile,cdsFile,ribo_offset_file,pr_ntcov_path,'pos',[],'nt',)) for covFile in covFiles]
# for p in proc:
#     p.start()
# for p in proc:
#     p.join()
# 
# #===============================================================================
# #                     2. gene total ribosome count for each gene (restuls in 04_gene_total_count)
# #===============================================================================
# os.chdir(fwd_rev_path)
# covFiles = [f for f in os.listdir(fwd_rev_path) if f.endswith('cov.txt')]
# covFiles = natsorted(covFiles)
# cdsFile = db_path + '/01_pr_cds.txt'
# #gene_pr_cov(covFiles[0],cdsFile,ribo_offset_file,gene_count_path,'gene',['heavychain'])
# proc = [Process(target=gene_pr_cov,args=(covFile,cdsFile,ribo_offset_file,gene_count_path,'gene',)) for covFile in covFiles]
# for p in proc:
#     p.start()
# for p in proc:
#     p.join()
# 
# #===============================================================================
# #                     3. Detect Stalling sites (results in 05_stall_sites)
# #===============================================================================
# #--------------- 1) read cds file -----------------------------
# cdsFile = db_path + '/01_pr_cds.txt'
# cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False,names=['chr','start','end','geneid','access','strand'])
# #--------------- 2) read position coverage file ---------------
# pr_cov_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam/03_pr_pos_cov'
# os.chdir(pr_cov_path)
# pr_cov_files = [f for f in os.listdir(pr_cov_path) if f.endswith('pos.txt')]
# pr_cov_files = natsorted(pr_cov_files)
# #CodonStallSites(pr_cov_files[0],cds_df,stall_path,15,10)
# proc = [Process(target=CodonStallSites,args=(f,cds_df,stall_path,15,10,)) for f in pr_cov_files]
# for p in proc:
#     p.start()
# for p in proc:
#     p.join()

# #===============================================================================
# #                        4. list all the start sites and stop sites in the chromosome for every protein
# #===============================================================================
# cdsFile = db_path + '/01_pr_cds.txt'
# cds_pr_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False,names=['chr','start','end','geneid','access','strand'])
# cds_obj = trpr(cds_pr_df)
# chr_pr_start_end_sites = cds_obj.all_pr_start_end_pos()
# if not os.path.exists(general_path): os.mkdir(general_path)
# os.chdir(general_path)
# chr_pr_start_end_sites.to_csv(general_path + '/01_chr_pr_start_end.txt',index=False,sep='\t')

# #===============================================================================
# #                        5. codon frequencies (resultsin 06_codon_AA_freq)
# #===============================================================================
# #---------- 1. find the protein sequence inconsistancy between refseq and ref genome --------------
# cdsFile = db_path + '/01_pr_cds.txt'
# prFaFile = db_path + '/combined_pr.fa'
# refFaFile = db_path + '/combined.fa'
# pr_seq_file = general_path + '/02_prInconsistant.txt'
# pr_id_file = general_path + '/02_prInconsistantID.txt'
# refseq_reffa_inconsist_pr(cdsFile,refFaFile,prFaFile,pr_id_file,pr_seq_file)
# 
# #--------- 2. codon frequencies -----------------------------
# # 1. read gene with all proteins file
# cdsFile = db_path + '/01_pr_cds.txt'
# # 2. read genome file
# refFaFile = db_path + '/combined.fa'
# refDNA_dic = SeqIO.index(refFaFile,'fasta')
# # 3). pr inconsistant id
# prInconID = general_path + '/02_prInconsistantID.txt'
# # 4) list chromosome protein start and end sites
# chr_pr_start_end_file = general_path + '/01_chr_pr_start_end.txt'
# # 5) calculate frequencies
# os.chdir(stall_path)
# stallFiles = [os.path.join(stall_path,f) for f in os.listdir(stall_path) if f.endswith('stallsites.txt')]
# stallFiles = natsorted(stallFiles)
# target_path = codon_AA_freq_path
# if not os.path.exists(target_path):
#     os.mkdir(target_path)
# os.chdir(target_path)
# for f in stallFiles:
#     StallSeq([f],cdsFile,refDNA_dic,prInconID,chr_pr_start_end_file,f)
#   
# target_path = codon_AA_freq_path + '/day'
# if not os.path.exists(target_path):
#     os.mkdir(target_path)
# os.chdir(target_path)
# StallSeq(stallFiles[0:3],cdsFile,refDNA_dic,prInconID,chr_pr_start_end_file,target_path + '/day3.txt') 
# StallSeq(stallFiles[3:6],cdsFile,refDNA_dic,prInconID,chr_pr_start_end_file,target_path + '/day6.txt')
#===============================================================================
#                        6. plot AA and codon frequencys (figures in folder figures)
#===============================================================================
#=================== change nt to AA and calculate the t-test===============
"""
# 1. get the significant AA and codons
path = codon_AA_freq_path
os.chdir(path)
codon_files = [f for f in os.listdir(path) if f.endswith('codon_freq.txt')]
codon_files = natsorted(codon_files)
AA = ['A','C','D','E','F','G','I','H','K','M','L','Q','P','S','R','T','W','V','Y']
bases = ['A','C','G','T']
codons = [a + b + c for a in bases for b in bases for c in bases]
codon_df = pd.DataFrame()
aa_df = pd.DataFrame()
for f in codon_files:
    codon_dic,aa_dic = getcodon_AAFreq(f)
    codon = []; aa = []
    for c in codons:
        if c in codon_dic:
            codon.append(codon_dic[c])
        else:
            codon.append(0)
    codon_df[f[:3]] = pd.Series(codon)
    for a in AA:
        if a in aa_dic:
            aa.append(aa_dic[a])
        else:
            aa.append(0)
    aa_df[f[:3]] = pd.Series(aa)
codon_df.index=codons
aa_df.index=AA
# t-test
codon_df['t_score'] = pd.Series(stats.ttest_ind(codon_df.iloc[:,0:3],codon_df.iloc[:,3:6],axis=1)[0],index=codon_df.index)
codon_df['p_value'] = pd.Series(stats.ttest_ind(codon_df.iloc[:,0:3],codon_df.iloc[:,3:6],axis=1)[1],index=codon_df.index)
print codon_df[codon_df['p_value'].values<0.05]
aa_df['t_score'] = pd.Series(stats.ttest_ind(aa_df.iloc[:,0:3],aa_df.iloc[:,3:6],axis=1)[0],index=aa_df.index)
aa_df['p_value'] = pd.Series(stats.ttest_ind(aa_df.iloc[:,0:3],aa_df.iloc[:,3:6],axis=1)[1],index=aa_df.index)
print aa_df[aa_df['p_value'].values<0.05]
# # wilcoxin rank sum test
# codon_df['t_score'] = pd.Series(stats.wilcoxon(codon_df.iloc[:,0:3],codon_df.iloc[:,3:6])[0],index=codon_df.index)
# codon_df['p_value'] = pd.Series(stats.wilcoxon(codon_df.iloc[:,0:3],codon_df.iloc[:,3:6])[1],index=codon_df.index)
# print codon_df[codon_df['p_value'].values<0.05]
# aa_df['t_score'] = pd.Series(stats.wilcoxon(aa_df.iloc[:,0:3],aa_df.iloc[:,3:6])[0],index=aa_df.index)
# aa_df['p_value'] = pd.Series(stats.wilcoxon(aa_df.iloc[:,0:3],aa_df.iloc[:,3:6])[1],index=aa_df.index)
# print aa_df[aa_df['p_value'].values<0.05]
# 2. plot the frequencies,get sequence percentage
freqFile1 = codon_AA_freq_path + '/day/day3codon_freq.txt'
freqFile2 = codon_AA_freq_path + '/day/day6codon_freq.txt'
day3_nt_dic,day3_aa_dic = getcodon_AAFreq(freqFile1)
day6_nt_dic,day6_aa_dic = getcodon_AAFreq(freqFile2)
print day3_nt_dic;print day6_nt_dic
print day3_aa_dic;print day6_aa_dic
print len(day3_nt_dic);print len(day6_nt_dic)
print len(day3_aa_dic);print len(day6_aa_dic)
# 3. define the colors
colors = ['#FF4500','#808080','#6495ED']
c_color = codon_df[codon_df['p_value'].values<0.05].index.tolist() # list of codons that are significant different.
c3_colors = []; c6_colors = []
codon = day6_nt_dic.keys()
codon = [str(Seq(aa,generic_dna).translate()) + '_'+str(aa) for aa in codon]  # generate codon_AA format
keys = sorted(codon)
for key in keys:
    if key[2:] in c_color:
        c3_colors.append(colors[0])
        c6_colors.append(colors[0])
    else:
        c3_colors.append(colors[1])
        c6_colors.append(colors[2])
a_color = aa_df[aa_df['p_value'].values<0.05].index.tolist()  # list of AA that are significant different.
a3_colors = []; a6_colors = []

keys = sorted(day6_aa_dic.keys())
for key in keys:
    if key in a_color:
        a3_colors.append(colors[0])
        a6_colors.append(colors[0])
    else:
        a3_colors.append(colors[1])
        a6_colors.append(colors[2])
        
# 4. plot day3 and day6 together
# plot codon frequency
codon = sorted(day3_nt_dic.keys())
codon = [str(c) for c in codon]
day3 = [];day6 = []
for c in codon:
    day3.append(day3_nt_dic[c])
    day6.append(day6_nt_dic[c])
df = pd.DataFrame()
df['day3'] = pd.Series(day3)
df['day6'] = pd.Series(day6)
codon = [str(Seq(aa,generic_dna).translate()) + '_'+str(aa) for aa in codon]
df.index = pd.Series(codon)
df = df.sort_index()
ax = df.plot(kind='bar',align='center',color=[c3_colors,c6_colors],legend=False,figsize=(14.5,8))
ax.set_title('codon frequencies in day3 and day6',fontsize=24)
ax.set_ylabel('percentage',color='black',fontsize=16)
ax.set_xlabel('Codon',color='black',fontsize=16)
import matplotlib.patches as mpatches
d3_patch = mpatches.Patch(color=colors[1], label='day3')
d6_patch = mpatches.Patch(color=colors[2], label='day6')
sig_patch = mpatches.Patch(color=colors[0], label='significant different')
ax.legend(handles=[d3_patch,d6_patch,sig_patch])
plt.savefig(fig_path + '/04_codon_freq.png')
plt.savefig(fig_path + '/04_codon_freq.svg')
# plot AA frequency
aa = sorted(day3_aa_dic.keys())
aa = [str(a) for a in aa]
day3 = [];day6 = []
for a in aa:
    day3.append(day3_aa_dic[a])
    day6.append(day6_aa_dic[a])
df = pd.DataFrame()
df['day3'] = pd.Series(day3)
df['day6'] = pd.Series(day6)
df.index = pd.Series(aa)
ax = df.plot(kind='bar',align='center',color=[a3_colors,a6_colors],figsize=(14.5,8))
ax.set_title('AA frequencies in day3 and day6',fontsize=24)
ax.set_ylabel('percentage',color='black',fontsize=16)
ax.set_xlabel('Amino Acids',color='black',fontsize=16)
plt.savefig(fig_path + '/04_AA_freq.png')
plt.savefig(fig_path + '/04_AA_freq.svg')
#plt.show()
"""
#===============================================================================
#                         7. plot translation and transcription change for each gene
#===============================================================================
# #-------------- 1). get total riboseq count -------------------
# os.chdir(bam_path)
# bamFiles = [f for f in os.listdir(bam_path) if f.endswith('.bam')]
# totalCount = []
# for bam in bamFiles:
#     handle = pysam.AlignmentFile(bam,'rb')
#     totalCount.append(handle.mapped)
# #-------------- 2). get ribo rpkm of each sample -------------------
# os.chdir(gene_count_path)
# ribo_count_files = [f for f in os.listdir(gene_count_path) if f.endswith('geneCount.txt')]
# ribo_count_files = natsorted(ribo_count_files)
# ribo_rpkm = pd.DataFrame()
# for f,total in zip(ribo_count_files,totalCount):
#     df = pd.read_csv(f,sep='\t',names=['GeneID','count','length'],header=0)
#     df[f] = df['count']/df['length']/float(total)*(10**9)   # calculate rpkm
#     ribo_rpkm[f]= df[f]
# ribo_rpkm.index=df['GeneID']
# ribo_rpkm['ribo_day3'] = ribo_rpkm.iloc[:,0:3].mean(axis=1)
# ribo_rpkm['ribo_day6'] = ribo_rpkm.iloc[:,3:6].mean(axis=1)
# ribo_rpkm = ribo_rpkm[(ribo_rpkm['ribo_day3']>0) & (ribo_rpkm['ribo_day6']>0)]
# #-------------- 3). get total rna count ---------------------
# os.chdir(rna_bam_path)
# bamFiles = [f for f in os.listdir(rna_bam_path) if f.endswith('.bam')]
# totalCount = []
# for bam in bamFiles:
#     handle = pysam.AlignmentFile(bam,'rb')
#     totalCount.append(handle.mapped)
# #-------------- 4). get rna fpkm for each sample ---------------------
# os.chdir(rna_count_path)
# rna_count_files = [f for f in os.listdir(rna_count_path) if f.endswith('Count.txt')]
# rna_count_files = natsorted(rna_count_files)
#  
# rna_rpkm = pd.DataFrame()
#  
# for f,total in zip(rna_count_files,totalCount):
#     df = pd.read_csv(f,sep='\t',header=0,names=['GeneID','count','length'])
#     df[f] = df['count']/df['length']/float(total)*(10**9)  # calculate rpkm
#     rna_rpkm[f] = df[f]
# rna_rpkm.index = df['GeneID']
# rna_rpkm['rna_day3'] = rna_rpkm.iloc[:,0:3].mean(axis=1)
# rna_rpkm['rna_day6'] = rna_rpkm.iloc[:,3:6].mean(axis=1)
# rna_rpkm = rna_rpkm[(rna_rpkm['rna_day3']>0) & (rna_rpkm['rna_day6']>0)]
# overlap_genes = set(ribo_rpkm.index.tolist()).intersection(rna_rpkm.index.tolist())  # overlapped genes
# ribo_rpkm = ribo_rpkm[ribo_rpkm.index.isin(overlap_genes)]
# rna_rpkm = rna_rpkm[rna_rpkm.index.isin(overlap_genes)]
# ribo_rna_df = pd.concat([ribo_rpkm[['ribo_day3','ribo_day6']],rna_rpkm[['rna_day3','rna_day6']]],axis=1)
# # get translation efficiency
# ribo_rna_df['day3_trans'] = ribo_rna_df['ribo_day3']/ribo_rna_df['rna_day3']
# ribo_rna_df['day6_trans'] = ribo_rna_df['ribo_day6']/ribo_rna_df['rna_day6']
# # get mrna and ribo fold change
# ribo_rna_df['ribo_change'] = (ribo_rna_df['ribo_day6']/ribo_rna_df['ribo_day3']).apply(np.log2)
# ribo_rna_df['mrna_change'] = (ribo_rna_df['rna_day6']/ribo_rna_df['rna_day3']).apply(np.log2)
# #---------- plot -------------
# plt.plot(ribo_rna_df['mrna_change'],ribo_rna_df['ribo_change'],'.',color='#6699CC')
#  
# heavy_df = ribo_rna_df[ribo_rna_df.index.isin(['heavychain'])]
# plt.plot(heavy_df['mrna_change'],heavy_df['ribo_change'],'.',color='red',alpha=1)
#   
# light_df = ribo_rna_df[ribo_rna_df.index.isin(['lightchain'])]
# plt.plot(light_df['mrna_change'],light_df['ribo_change'],'.',color='yellow',alpha=1)
#   
# neo_df = ribo_rna_df[ribo_rna_df.index.isin(['NeoRKanR'])]
# plt.plot(neo_df['mrna_change'],neo_df['ribo_change'],'.',color='black',alpha=1)
#   
# plt.xlabel('log2 mRNA change')
# plt.ylabel('log2 rpf change')
# plt.title('Translation change Day6 VS Day3')
# plt.axhline(y=0,color = 'k')
# plt.axvline(x=0,color='k')
# plt.savefig(fig_path + '/05_translation_efficiency.png')
# plt.savefig(fig_path + '/05_translation_efficiency.svg')
#===============================================================================
#                         8. 3'UTR of all genes
#===============================================================================
# #----------------- 1. read coverage file ------------------------------
# os.chdir(fwd_rev_path)
# covFiles = [f for f in os.listdir(fwd_rev_path) if f.endswith('cov.txt')]
# covFiles = natsorted(covFiles)
# #utr3_cov(exnFile,cdsFile,all_id_file,covFiles[0],utr3_cov_path,genes=[])
# proc = [Process(target=utr3_cov,args=(exnFile,cdsFile,all_id_file,f,utr3_cov_path,)) for f in covFiles]
# for p in proc:
#     p.start()
# for p in proc:
#     p.join()
# #===============================================================================
# #                         9. frame coverage of each frame
# #===============================================================================
# os.chdir(pr_ntcov_path)
# covFiles = [f for f in os.listdir(pr_ntcov_path) if f.endswith('txt')]
# covFiles = natsorted(covFiles)
# # genes_frame_cov(covFiles[0],frame_cov_path,genes=['heavychain'],calType='sum')
# proc = [Process(target=genes_frame_cov,args=(f,frame_cov_path)) for f in covFiles]
# for p in proc:
#     p.start()
# for p in proc:
#     p.join()


