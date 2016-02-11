from __future__ import division
import subprocess,os,sys
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
import pandas as pd
from natsort import natsorted
import numpy as np
import scipy.stats as sp_stats
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from f02_RiboDataModule import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from multiprocessing import Process
from Modules.p05_ParseGff import extractAllPr
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('ggplot')
import matplotlib.gridspec as gridspec

bam_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam'
db_path = '/data/shangzhong/RibosomeProfiling/Database'
signalP_path = '/data/shangzhong/RibosomeProfiling/signalP_part'   # part 8
rna_bam_path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align'

fig_path = signalP_path + '/figures'
fwd_rev_path = bam_path + '/01_cov'   # mapping file. part 8
gene_count_path = bam_path + '/04_gene_total_count'  # part 8
rna_gene_count_path = rna_bam_path + '/01_gene_count'
rna_cov_path = rna_bam_path + '/02_cov'
stall_site_path = bam_path + '/05_stall_sites'
frame_cov_path = bam_path + '/09_frame_cov'

exnFile = db_path +'/01_pr_rna.txt'  # part 8,
cdsFile = db_path + '/01_pr_cds.txt' # part 8,
ribo_offset_file = bam_path + '/02_TSS_TSE_cov/ribo_offset.txt'      # part 8,
refFaFile = db_path + '/combined.fa'
gffFile = db_path + '/combined.gff'
mapFile = db_path + '/combined_AllIDS.txt'
all_id_file = db_path + '/combined_AllIDS.txt'

def chunk(l,n):
    n = max(n,1)
    result = [l[i:i+n] for i in range(0,len(l),n)]
    return result
#===============================================================================
#                         0. prepare protein files for signalP, add antibody proteins.
#===============================================================================
# ref_dict = SeqIO.index(refFaFile, "fasta")
# output_handle = open(signalP_path + '/01_antibody_pr.fa', "a")
# # get pr sequence
# neo = ref_dict['neoR_contig']
# neo_pr = Seq.translate(neo.seq[154:949])
# neo.id = 'NeoRKanR_prid';neo.seq = neo_pr
# heavy = ref_dict['heavychain_contig']
# heavy_pr = Seq.translate(heavy.seq[233:1637])
# heavy.id = 'heavychain_prid';heavy.seq = heavy_pr
# light = ref_dict['lightchain_contig']
# light_pr = Seq.translate(light.seq[88:805])
# light.id = 'lightchain_prid';light.seq=light_pr
# SeqIO.write(neo, output_handle, "fasta")
# SeqIO.write(heavy,output_handle,'fasta')
# SeqIO.write(light,output_handle,'fasta')
# output_handle.close()
# #===============================================================================
# #                         1. run signalP
# #===============================================================================
# inputFa = db_path + '/combined_pr.fa'
# choPrRefseq_sp = signalP(inputFa)
# os.rename(choPrRefseq_sp,signalP_path+'/02_cho_ab_pr_sp.txt')                                              # 02_cho_ab_pr_sp.txt
# #===============================================================================
# #                         2. insert gene id and accession of mRNA
# #===============================================================================
# choPrRefseq_sp = signalP_path + '/02_cho_ab_pr_sp.txt'
# choPrRefseq_sp_gene = addGeneIDSymbolChr2signalP(choPrRefseq_sp,mapFile,organism='',mapSource='gff')
# choPrRefseq_sp_gene = changeFileName(choPrRefseq_sp_gene)                       # 03_cho_ab_pr_sp.txt
# #===============================================================================
# #                         3. get genes mapping to multiple scaffold
# #===============================================================================
# choPrRefseq_sp_gene = signalP_path + '/03_cho_ab_pr_sp.gene.txt'
# choPrRefseq_gene_multichr =  genesMap2diffChrom(choPrRefseq_sp_gene,mapFile)
# choPrRefseq_gene_multichr = changeFileName(choPrRefseq_gene_multichr)            # 04_cho_ab_pr_sp.gene.multichr.txt
# #===============================================================================
# #                         4. classify gene ids into 3 groups. with sp, without sp, with and without sp
# #===============================================================================
# choPrRefseq_sp_gene = signalP_path + '/03_cho_ab_pr_sp.gene.txt'
# choPrRefseq_sp_gene_classify = gene_sp_classify(choPrRefseq_sp_gene)
# choPrRefseq_sp_gene_classify = changeFileName(choPrRefseq_sp_gene_classify,2)    # 05_cho_ab_pr_sp.gene.classify.txt
# #===============================================================================
# #      5. get target sp genes and non sp genes, remove those genes
# #      mapping to muliple scaffolds and genes have overlapped cds.
# #===============================================================================
# 
# #------------ 1). read cds file
# cds_df = pd.read_csv(db_path+'/01_pr_cds.txt',header=0,sep='\t',low_memory=False,names=['chr','start','end','geneid','access','strand'])
# cds_obj = trpr(cds_df)
# #------------ 2). get gene ids
# geneIDs = cds_df['geneid'].tolist()
# geneIDs = list(set(geneIDs))
# #------------ 3). get the genes with sp and without sp
# gene_class_file = signalP_path + '/05_cho_ab_pr_sp.gene.classify.txt'
# gene_class_df = pd.read_csv(gene_class_file,header=0,sep='\t')
# gene_class_df = gene_class_df.astype(str)
# gene_sp = gene_class_df['with_sp'].tolist()
# gene_sp = list(set(gene_sp))
# gene_no_sp = gene_class_df['no_sp'].tolist()
# gene_no_sp = list(set(gene_no_sp))
# try:
#     gene_sp.remove('-');gene_no_sp.remove('-')
# except:
#     pass
# print '# of genes with sp:',len(gene_sp)
# print '# of genes without sp:',len(gene_no_sp)
# #------------ 4). exclude the genes with multiple scaffolds
# gene_multi_list = cds_obj.genes2multi_chr()
# gene_sp = [g for g in gene_sp if g not in gene_multi_list]
# gene_no_sp = [g for g in gene_no_sp if g not in gene_multi_list]
# print '# of genes with sp after remove multi_chrom:',len(gene_sp)
# print '# of genes without sp after remove multi_chrom:',len(gene_no_sp)
# gene_sp_df = pd.DataFrame({'gene_sp':gene_sp})
# gene_no_sp_df = pd.DataFrame({'gene_no_sp':gene_no_sp})
# df = pd.concat([gene_sp_df,gene_no_sp_df],axis=1)
# sp_nosp_gene_file = signalP_path + '/06_sp_no_sp_genes.txt'        # 06_sp_no_sp_genes.txt
# df.to_csv(sp_nosp_gene_file,sep='\t',index=False)
# #===============================================================================
# #                         6. coverage around signal peptide end site
# #===============================================================================
# def covNearSpEnd(pr_pos_cov_file,sp_len_dic,up,down,total):
#     """This function calculates coverage of each AA around the cutting site of signal peptide
#      
#     * pr_pos_cov_file: str. A tab delimited file, first column is gene id, second is pr access, others are coverage at each position
#     * sp_len_dic: dict. {praccess:length}
#     * up: int. number of nucleotide upstream of the cutting nt.
#     * down: int. number of nucleotide downstream of the cutting nt.
#     """
#     up = int(up); down = int(down)
#     sp_end_cov_df = pd.DataFrame() # row position, column protein
#     handle = open(pr_pos_cov_file,'r')
#     for line in handle:
#         item = line[:-1].split('\t')
#         pr = item[1]
#         cov =item[2:]
#         if pr not in sp_len_dic:
#             continue
#         sp_len = sp_len_dic[pr]
#         pr_cov = []
#         for index in range(sp_len-up,sp_len+down):   #  range(0,down)
#             if index < 0 or index >= len(cov):
#                 pr_cov.append(0)
#             else:
#                 pr_cov.append(cov[index])
#         sp_end_cov_df[pr]=pd.Series(pr_cov)
#     sp_end_cov_df = sp_end_cov_df.astype(int)
#     sp_end_cov_df = sp_end_cov_df/float(total)*(10**6)
#     res = pd.DataFrame()
#     res[pr_pos_cov_file.split('_')[0]] = sp_end_cov_df.sum(axis=1)
#     return res
#  
# # #-------------- 1) get total count --------------------
# os.chdir(bam_path)
# bamFiles = [f for f in os.listdir(bam_path) if f.endswith('.sort.bam')]
# bamFiles = natsorted(bamFiles)
# totalCount = []
# for bam in bamFiles:
#     handle = pysam.AlignmentFile(bam,'rb')
#     totalCount.append(handle.mapped)
# #-------------- 2) read sp genes ----------------------
# sp_nosp_gene_file = signalP_path + '/06_sp_no_sp_genes.txt'
# gene_df = pd.read_csv(sp_nosp_gene_file,sep='\t',header=0)
# gene_sp = gene_df['gene_sp'].dropna().astype(str).tolist()
# gene_sp.remove('heavychain');gene_sp.remove('lightchain')
# #-------------- 3) build sp length dictionary. {praccess:sp length} ------------
# choPrRefseq_sp_gene = signalP_path + '/03_cho_ab_pr_sp.gene.txt'
# df = pd.read_csv(choPrRefseq_sp_gene,sep='\t',header=0)
# sp_anno_df = df[df['GeneID'].isin(gene_sp)] 
# sp_len_df = sp_anno_df[['Pr_Access','posY']]
# sp_len_dic = sp_anno_df.set_index('Pr_Access')['posY'].to_dict()
# #-------------- 4) get coverage for all sp proteins around sp end----------------------------
# up = 10; down = 100
# pr_pos_cov_path = bam_path + '/03_pr_pos_cov'
# os.chdir(pr_pos_cov_path)
# pr_pos_cov_files = [f for f in os.listdir(pr_pos_cov_path) if f.endswith('txt')]
# pr_pos_cov_files = natsorted(pr_pos_cov_files)
#   
# sp_cov_list = []   # rows sum at each position, column sample
# for f,total in zip(pr_pos_cov_files,totalCount):
#     df = covNearSpEnd(f,sp_len_dic,up,down,total)
#     sp_cov_list.append(df)
# sp_cov_df = pd.concat(sp_cov_list,axis=1)
#   
# sp_cov_df['day3'] = sp_cov_df.iloc[:,0:3].mean(axis=1)
# sp_cov_df['day6'] = sp_cov_df.iloc[:,3:6].mean(axis=1)
# df = sp_cov_df[['day3','day6']]
# #-------------- 5) plot ------------------------------
# if not os.path.exists(fig_path): os.mkdir(fig_path)
# # plot the coverage for each replicates
# sp_cov_df.index = range(-up,down)
# ax = sp_cov_df.iloc[:,0:6].plot(kind='line',title='Coverage around cutting site of signal peptide',figsize=(14.5,8))
# ax.set_xlabel('Distance to end of signal peptide',color='black',fontsize=16)
# x = np.array(range(-up,down))
# plt.xticks(x,range(-up,down),rotation=90)
# plt.savefig(fig_path + '/01_cov_near_sp_sample_level.svg')
# # plot the coverage for day 3 and day6
# f,ax=plt.subplots(2,sharex=True)
# f.set_size_inches(14.5,8)
# x = np.array(range(-up,down))
# #x = np.array(range(0,down))
# for i in range(2):
#     ax[i].bar(x,df.iloc[:,i],align='center',alpha=0.89)
#     ax[i].set_xlim([-up,down])
#     ax[i].set_title(df.columns.values[i])
# ax[1].set_xlabel('Distance to end of signal peptide',color='black',fontsize=18)
# f.text(0.07, 0.5, 'total (rpm)', ha='center', va='center', rotation='vertical',fontsize=18)
# plt.xticks(x,range(-up,down),rotation=90)
# plt.suptitle('coverage around cutting site of signal peptide',fontsize=26,color='black')
# plt.savefig(fig_path+'/02_sp_end_cov.pdf')
# plt.savefig(fig_path+'/02_sp_end_cov.svg')
# #plt.show()                                   # 02_sp_end_cov.pdf is in folder figures
#===============================================================================
#                         7. get RNAseq transcriptional distribution (calculate how many percentage each class of transcripts occupies)
#                            result file 07_rna_distribution is in folder signalP_part
#===============================================================================
#------- 1. get the genes with sp and without sp --------------------
def occupancy_distribution(sp_nosp_gene_file,count_path,outFile):
    """This function calculates the occupancy distribution for riboseq or rnaseq"""
    
    gene_df = pd.read_csv(sp_nosp_gene_file,sep='\t',header=0)
    gene_sp = gene_df['gene_sp'].dropna().astype(str).tolist()
    gene_no_sp = gene_df['gene_no_sp'].dropna().astype(str).tolist()
      
    genes = gene_sp + gene_no_sp
    gene_sp.remove('heavychain');gene_sp.remove('lightchain')
    gene_no_sp.remove('NeoRKanR')
      
    #------- 2. read the count data for all replicates -------------------
    os.chdir(count_path)
    countFiles = [f for f in os.listdir(count_path) if f.endswith('Count.txt')]
    countFiles = natsorted(countFiles)
    dfs = []
    for f in countFiles:
        df = pd.read_csv(f,sep='\t',header=0,index_col=0,usecols=[0,1],names=['geneid',f])
        dfs.append(df)
    gene_df = pd.concat(dfs,axis=1)
    gene_df = gene_df[gene_df.index.isin(genes)]    # only concern the genes that are used in riboseq
    #gene_df = gene_df.iloc[:-5,:]                  # this will concern all the genes
    gene_df.loc['total'] = gene_df.sum()
    gene_df.loc['other'] = (gene_df[~gene_df.index.isin(genes+['total'])]).sum()
    gene_df.loc['sp'] = (gene_df[gene_df.index.isin(gene_sp)]).sum()
    gene_df.loc['no_sp'] = (gene_df[gene_df.index.isin(gene_no_sp)]).sum()
    # stores in signalP path
    percent_df = gene_df.div(gene_df.loc['total'])
    percent_df.loc[['sp','no_sp','heavychain','lightchain','NeoRKanR','other','total']].to_csv(outFile,sep='\t')
# sp_nosp_gene_file = signalP_path + '/06_sp_no_sp_genes.txt'
# outFile = signalP_path+'/07_rna_distribution.csv'
# occupancy_distribution(sp_nosp_gene_file,rna_gene_count_path,outFile)
#===============================================================================
#                     8. get Riboseq distribution
#===============================================================================
# sp_nosp_gene_file = signalP_path + '/06_sp_no_sp_genes.txt'
# outFile = signalP_path+'/07_ribo_distribution.csv'
# occupancy_distribution(sp_nosp_gene_file,gene_count_path,outFile)
#===============================================================================
#                 9. 3'UTR of recombinant genes
#===============================================================================
def gene_cov_at_trpt_level(totalCount,covFiles,gene,exnFile):
    """
    This function calculates the ribo seq coverage on the whole transcript. and then merge the replicates.
    * totalCount: list. A list of total count for replicated bam files.
    * covFiles: list. A list of filenames with 6 columns. # 'num','chr','end5','end3','strand','len'
    * gene: str. Target gene
    """
    exn_df = pd.read_csv(exnFile,sep='\t',header=0,low_memory=False)
    gene_cov_df = pd.DataFrame()
    for f in covFiles:
        cov_df = pd.read_csv(f,sep='\t',header=None,names=['num','chr','end5','end3','strand','len'])
        gene_df = exn_df[exn_df['GeneID'].values==gene]
        gene_obj = trpr(gene_df)
        trpr_id = list(set(gene_df['access'].tolist()))[0]
        pos = gene_obj.get_trpr_pos(trpr_id)
        gene_array_set = gene_obj.htseq_genome_array_set(False)
        cov_pos = posCoverage(gene_df,cov_df,pos,gene_array_set)
        gene_cov_df[f] = pd.Series(cov_pos)
    gene_cov_df = gene_cov_df.div(totalCount) * (10**6)
    gene_cov_df['mean'] = gene_cov_df.mean(axis=1)
    return gene_cov_df['mean']
"""
#-------------- 1) get total count --------------------
os.chdir(bam_path)
bamFiles = [f for f in os.listdir(bam_path) if f.endswith('.sort.bam')]
bamFiles = natsorted(bamFiles)
totalCount = []
for bam in bamFiles:
    handle = pysam.AlignmentFile(bam,'rb')
    totalCount.append(handle.mapped)
#-------------- 2) get coverage for heavy/lightchain -------------------
outpath = signalP_path + '/01_heavy_light_nt_cov'
os.chdir(fwd_rev_path)
covFiles = [f for f in os.listdir(fwd_rev_path) if f.endswith('cov.txt')]
covFiles = natsorted(covFiles)
heavy_df = pd.DataFrame()
heavy_df['day3'] = gene_cov_at_trpt_level(totalCount[0:3],covFiles[0:3],'heavychain',exnFile)
heavy_df['day6'] = gene_cov_at_trpt_level(totalCount[3:6],covFiles[3:6],'heavychain',exnFile)
#---------------- heavy chain -----------------
def height_limit(x,n):
    if x > n:
        return n
    else:
        return x
plt_df = heavy_df#.loc[1637:]
plt_df = plt_df + 0.001
#plt_df = plt_df.apply(np.log2)

f,ax=plt.subplots(2,sharex=True)
f.set_size_inches(14.5,8)
for i in range(2):
    ax[i].bar(plt_df.index,plt_df.iloc[:,i],align='center',color='black',edgecolor='black')
    ax[i].set_yscale('log')
ax[1].set_xlabel('Distance to start of mRNA',color='black',fontsize=18)
f.text(0.07, 0.5, 'mean (rpm)', ha='center', va='center', rotation='vertical',fontsize=18)
plt.suptitle('coverage around 3 UTR of heavychain', fontsize=26,color='black')
plt.savefig(fig_path+'/03_heavy_ribo_mrna_cov.svg',format='svg')
plt.savefig(fig_path+'/03_heavy_ribo_mrna_cov.pdf')

#-------------- light chain---------
light_df = pd.DataFrame()
light_df['day3'] = gene_cov_at_trpt_level(totalCount[0:3],covFiles[0:3],'lightchain',exnFile)
light_df['day6'] = gene_cov_at_trpt_level(totalCount[3:6],covFiles[3:6],'lightchain',exnFile)
plt_df = light_df#.loc[1312:]
plt_df = plt_df + 0.001
#plt_df = plt_df.apply(np.log2)

f,ax=plt.subplots(2,sharex=True)
f.set_size_inches(14.5,8)
for i in range(2):
    ax[i].bar(plt_df.index,plt_df.iloc[:,i],align='center',width=1,color='black',edgecolor='black')
    ax[i].set_yscale('log')
ax[1].set_xlabel('Distance to start of mRNA',color='black',fontsize=18)
f.text(0.07, 0.5, 'mean (rpm)', ha='center', va='center', rotation='vertical',fontsize=18)
plt.suptitle('coverage around 3 UTR of lightchain', fontsize=26,color='black')
plt.savefig(fig_path+'/03_light_ribo_mrna_cov.svg',format='svg')
plt.savefig(fig_path+'/03_light_ribo_mrna_cov.pdf')
plt.show()
"""
#===============================================================================
#                         10. find out one endogenous gene coverage (this part cannot be automated because there are some parameters
#                             in part 2 that I need to modify for each different gene) (for Eef1a1, heavychain,lightchain,NeoRKanR)
#===============================================================================
#============================ ribo seq plot ===============================
def plot_cds_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,genename,ylim1,ylim2,hp,dy0_times,dy1_times,outFile,shift=0,more_shift=0):
    """
    * shift: int. shift for distance between 5'end and TSS of ribosome.
    * more_shift: int. this is distance betwween transription and translation start site.
    """
    #-------------- 3) get each exn start and end pos, CDS first and last pos -----------------------
    # exon start and end
    exn_df = pd.read_csv(exnFile,sep='\t',header=0,low_memory=False)
    gene_df = exn_df[exn_df['GeneID'].values==gene]
    gene_obj = trpr(gene_df)
    exn_start,exn_end = gene_obj.get_exn_cds_start_end_pos(gene)
    ex_start = [abs(n-exn_start[0]) for n in exn_start]
    ex_end = [abs(n-exn_start[0]) for n in exn_end]
    # cds start and end
    cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
    gene_df = cds_df[cds_df['GeneID'].values==gene]
    gene_obj = trpr(gene_df)
    cds_start,cds_end = gene_obj.get_exn_cds_start_end_pos(gene)
    cds_start = abs(cds_start[0]-exn_start[0])
    cds_end = abs(cds_end[-1]-exn_start[0])
    print 'exon start is:',ex_start
    print 'exon end is:', ex_end
    print 'cds start is:',cds_start
    print 'cds end is:', cds_end
    #-------------- 3) get coverage -----------------------
    gene_df = pd.DataFrame()
    #gene_df['day3'] = merge_rep_gene_pos_cov(totalCount[0:3],covFiles[0:3],gene,exn_df,covType)
    gene_df['day6'] = merge_rep_gene_pos_cov(totalCount[3:6],covFiles[3:6],gene,exn_df,covType)
    #-------------- 4) plot day3 and day6 ------------------
    plt_df = gene_df.copy()
    if shift != 0 or more_shift != 0:
        for start, end in zip(ex_start[1:],ex_end[:-1]):
            df1 = (plt_df.loc[start-shift+1:start]).copy()
            df2 = (plt_df.loc[end-shift+1:end]).copy()
            df2.index = df1.index
            plt_df.loc[start-shift+1:start] = df1.add(df2)
            plt_df.loc[end-shift+1:end] = 0
        plt_df.index = plt_df.index+shift-more_shift
    y1ratio = (ylim1[1]-ylim1[0])/(ylim2[1]-ylim2[0]+ylim1[1]-ylim1[0])
    y2ratio = (ylim2[1]-ylim2[0])/(ylim2[1]-ylim2[0]+ylim1[1]-ylim1[0])
    
    f = plt.figure()
    gs = gridspec.GridSpec(2,1,height_ratios=[1,3])
    ax0 = f.add_subplot(gs[0]);ax1 = f.add_subplot(gs[1]);#ax3 = f.add_subplot(gs[2]);ax4 = f.add_subplot(gs[3])
    f.set_size_inches(14.5,8)
    ax0.bar(plt_df.index,plt_df[gene_df.columns[0]],align='center',color='black',edgecolor='black')
    ax1.bar(plt_df.index,plt_df[gene_df.columns[0]],align='center',color='black',edgecolor='black')
    # limit the view
    ax0.set_ylim(ylim1)
    ax1.set_ylim(ylim2)
    plt.subplots_adjust(hspace=hp)
    ax0.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax0.xaxis.tick_top()
    ax0.tick_params(labeltop='off')
    ax1.xaxis.tick_bottom()
    
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(color='k', clip_on=False) # transform=ax1.transAxes,
    xlim = ax1.get_xlim()
    dx = 0.02 * (xlim[1]-xlim[0])
    dy = 0.003 * (ylim1[1]-ylim1[0])/y1ratio
    ax0.plot((xlim[0]-dx,xlim[0]+dx), (ylim1[0]-dy*dy0_times,ylim1[0]+dy*dy0_times), **kwargs)
    ax0.plot((xlim[1]-dx,xlim[1]+dx), (ylim1[0]-dy*dy0_times,ylim1[0]+dy*dy0_times), **kwargs)
    
    dy = 0.003 * (ylim2[1]-ylim2[0])/y2ratio
    ax1.plot((xlim[0]-dx,xlim[0]+dx), (ylim2[1]-dy*dy1_times,ylim2[1]+dy*dy1_times), **kwargs)
    ax1.plot((xlim[1]-dx,xlim[1]+dx), (ylim2[1]-dy*dy1_times,ylim2[1]+dy*dy1_times), **kwargs)
    ax0.set_xlim(xlim)
    ax1.set_xlim(xlim)
    
    for s,e in zip(ex_start,ex_end):
        ax0.axvline(s-more_shift, color='red', linestyle='--',label='exon_start')
        ax0.axvline(e-more_shift, color='blue',linestyle='--',label='exon_end')
    ax0.axvline(cds_start-more_shift, color='green',linestyle='--',label='cds_start')
    ax0.axvline(cds_end-more_shift, color='yellow',linestyle='--',label='cds_end')
    handles, labels = ax0.get_legend_handles_labels()
    ax0.legend(handles[-4:], labels[-4:])
#     # NeoR plot
#     ax0.axvline(168-more_shift,color='black',linestyle='--',label='alter_cds_start')
#     ax0.axvline(326-more_shift,color='orange',linestyle='--',label='alter_cds_end')
#     ax0.axvline(327-more_shift,color='black',linestyle='--',label='alter_cds_start')
#     ax0.axvline(713-more_shift,color='orange',linestyle='--',label='alter_cds_end')
#     ax0.axvline(717-more_shift,color='black',linestyle='--',label='alter_cds_start')
#     ax0.axvline(818-more_shift,color='orange',linestyle='--',label='alter_cds_end')
#     handles,labels = ax0.get_legend_handles_labels()
#     ax0.legend(handles[-10:-4],labels[-10:-4])
#     ax1.axvline(168-more_shift,color='black',linestyle='--',label='alter_cds_start')
#     ax1.axvline(326-more_shift,color='orange',linestyle='--',label='alter_cds_end')
#     ax1.axvline(327-more_shift,color='black',linestyle='--',label='alter_cds_start')
#     ax1.axvline(713-more_shift,color='orange',linestyle='--',label='alter_cds_end')
#     ax1.axvline(717-more_shift,color='black',linestyle='--',label='alter_cds_start')
#     ax1.axvline(818-more_shift,color='orange',linestyle='--',label='alter_cds_end')
    
    for s,e in zip(ex_start,ex_end):
        ax1.axvline(s-more_shift, color='red', linestyle='--',label='exon_start')
        ax1.axvline(e-more_shift, color='blue',linestyle='--',label='exon_end')
    ax1.axvline(cds_start-more_shift, color='green',linestyle='--',label='cds_start')
    ax1.axvline(cds_end-more_shift, color='yellow',linestyle='--',label='cds_end')
    
    
    ax0.set_title('{name} position coverage at {day} (Riboseq)'.format(name=genename,day=gene_df.columns.tolist()[0]))
    plt.xlabel('gene position')
    plt.ylabel('rpm')
    plt.savefig(outFile+'.svg')
    plt.savefig(outFile+'.pdf')
    

def plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,genename,outFile,ylim=[],seqtype='RNAseq',scale='',shift=0,more_shift=0):
    # exon start and end
    exn_df = pd.read_csv(exnFile,sep='\t',header=0,low_memory=False)
    gene_df = exn_df[exn_df['GeneID'].values==gene]
    gene_obj = trpr(gene_df)
    exn_start,exn_end = gene_obj.get_exn_cds_start_end_pos(gene)
    ex_start = [abs(n-exn_start[0]) for n in exn_start]
    ex_end = [abs(n-exn_start[0]) for n in exn_end]
    # cds start and end
    cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)    
    gene_df = cds_df[cds_df['GeneID'].values==gene]
    gene_obj = trpr(gene_df)
    cds_start,cds_end = gene_obj.get_exn_cds_start_end_pos(gene)
    cds_start = abs(cds_start[0]-exn_start[0])
    cds_end = abs(cds_end[-1]-exn_start[0])
    print 'exon start is:',ex_start
    print 'exon end is:', ex_end
    print 'cds start is:',cds_start
    print 'cds end is:', cds_end
    gene_df = pd.DataFrame()
    #gene_df['day3'] = merge_rep_gene_pos_cov(totalCount[0:3],covFiles[0:3],gene,exn_df,covType)
    gene_df['day6'] = merge_rep_gene_pos_cov(totalCount[3:6],covFiles[3:6],gene,exn_df,covType)
    #-------------- 3) plot figure -----------------------------
    plt_df = gene_df.copy()
    if shift != 0 or more_shift !=0:
        for start, end in zip(ex_start[1:],ex_end[:-1]):
            df1 = (plt_df.loc[start-shift+1:start]).copy()
            df2 = (plt_df.loc[end-shift+1:end]).copy()
            df2.index = df1.index
            plt_df.loc[start-shift+1:start] = df1.add(df2)
            plt_df.loc[end-shift+1:end] = 0
        plt_df.index = plt_df.index+shift-more_shift    
    f = plt.figure()
    f.set_size_inches(14.5,8)
    ax = plt.subplot(111)
    ax.bar(plt_df.index,plt_df[gene_df.columns[0]],align='center',color='black',edgecolor='black')
    if scale == 'log':
        ax.set_yscale('log',basey=10)
    if ylim!=[]:
        ax.set_ylim(ylim)
    
    for s,e in zip(ex_start,ex_end):
        ax.axvline(s-more_shift, color='red', linestyle='--',label='exon_start')
        ax.axvline(e-more_shift, color='blue',linestyle='--',label='exon_end')
    ax.axvline(cds_start-more_shift, color='green',linestyle='--',label='cds_start')
    ax.axvline(cds_end-more_shift, color='yellow',linestyle='--',label='cds_end')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[-4:], labels[-4:])
    
    ax.set_title('{name} position coverage at {day} ({seqtype})'.format(name=genename,seqtype=seqtype,day=gene_df.columns[0]))
    plt.xlabel('gene position')
    plt.ylabel('rpm')
    plt.savefig(outFile +'.pdf')
    plt.savefig(outFile +'.svg')

# #-------------- 1) get total count --------------------
# os.chdir(bam_path)
# bamFiles = [f for f in os.listdir(bam_path) if f.endswith('.sort.bam')]
# bamFiles = natsorted(bamFiles)
# totalCount = []
# for bam in bamFiles:
#     handle = pysam.AlignmentFile(bam,'rb')
#     totalCount.append(handle.mapped)
# #-------------- 2) get genes and set all parameters -----------------------
# fn = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam/04_gene_total_count/s05_cov_geneCount.txt'
# df = pd.read_csv(fn,sep='\t',header=0)
# df = df.sort(['count'])
# genes = df['GeneID'].tolist()
# genes.reverse()
# gene1 = genes[4]
# covType = 'gene'
# #-------------- 3) get coverage files -----------------------
# os.chdir(fwd_rev_path)
# covFiles = [f for f in os.listdir(fwd_rev_path) if f.endswith('cov.txt')]
# covFiles = natsorted(covFiles)
# # #-------------- 4) plot Eef1a1 ribo ------------------
# ylim1 = [300,910]
# ylim2 = [0,75]
# hp = 0.05
# dy0_times = 25  # this is for the slip of diagnal line that connecting broken axis
# dy1_times = 1
# # outFile = fig_path + '/07_Eef1a1_ribo_day3_cov'
# # plot_cds_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene1,'Eef1a1',ylim1,ylim2,hp,dy0_times,dy1_times,outFile,0,0)
# ylim1= [70,500]
# outFile = fig_path + '/15_Eef1a1_ribo_day6_cov_align'
# plot_cds_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene1,'Eef1a1',ylim1,ylim2,hp,dy0_times,dy1_times,outFile,15,34)
# #-------------- 5) plot heavychain ribo ------------------
# ylim1 = [8000,9100]
# ylim2 = [0,3000]
# hp = 0.05
# dy0_times = 10  # this is for the slip of diagnal line that connecting broken axis
# dy1_times = 9
# gene = 'heavychain'
# # outFile = fig_path + '/07_heavychain_ribo_day3_cov'
# # plot_cds_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'heavychain',ylim1,ylim2,hp,dy0_times,dy1_times,outFile,0,0)
# ylim1=[600,1500]
# ylim2=[0,500]
# outFile = fig_path + '/15_heavychain_ribo_day6_cov_align'
# plot_cds_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'heavychain',ylim1,ylim2,hp,dy0_times,dy1_times,outFile,15,233)
# #-------------- 6) plot lightchain ribo ------------------
# ylim1 = [700,800]
# ylim2 = [0,500]
# hp = 0.05
# dy0_times = 5  # this is for the slip of diagnal line that connecting broken axis
# dy1_times = 9
# gene = 'lightchain'
# # outFile = fig_path + '/07_lightchain_ribo_day3_cov'
# # plot_cds_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'lightchain',ylim1,ylim2,hp,dy0_times,dy1_times,outFile,0,0)
# ylim1=[65,120]
# ylim2=[0,50]
# outFile = fig_path + '.;/15_lightchain_ribo_day6_cov_align'
# plot_cds_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'lightchain',ylim1,ylim2,hp,dy0_times,dy1_times,outFile,15,88)
# #-------------- 7) plot NeoRkanR ribo ------------------
# ylim1 = [1500,9000]
# ylim2 = [0,1000]
# hp = 0.05
# dy0_times = 19  # this is for the slip of diagnal line that connecting broken axis
# dy1_times = 1
# gene = 'NeoRKanR'
# # outFile = fig_path + '/07_NeoR_ribo_day3_cov'
# # plot_cds_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'NeoR',ylim1,ylim2,hp,dy0_times,dy1_times,outFile,0,0)
# ylim1=[800,2000]
# ylim2=[0,600]
# outFile = fig_path + '/15_NeoR_ribo_day6_cov_align'
# plot_cds_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'NeoR',ylim1,ylim2,hp,dy0_times,dy1_times,outFile,16,154)
# #---------------- including the alternative plot -----------------------
# gene = 'NeoRKanR'
# outFile = fig_path + '/14_NeoR_ribo_day3_cov_alter_align'
# plot_cds_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'NeoR',ylim1,ylim2,hp,dy0_times,dy1_times,outFile,16,154)

# #------------------- use log scale ----------------
# gene = gene1
# outFile = fig_path + '/09_Eef1a1_day3_ribo_log_cov'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'Eef1a1',outFile,[0.1,1000],'Riboseq','log')
# outFile = fig_path + '/10_Eef1a1_day3_ribo_log_cov_align'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'Eef1a1',outFile,[0.1,1000],'Riboseq','log',15,34)
#  
# gene = 'heavychain'
# outFile = fig_path + '/09_heavychain_day3_ribo_log_cov'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'heavychain',outFile,[0.1,10000],'Riboseq','log')
# outFile = fig_path + '/10_heavychain_day3_ribo_log_cov_align'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'heavychain',outFile,[0.1,10000],'Riboseq','log',15,233)
#  
# gene = 'lightchain'
# outFile = fig_path + '/09_lightchain_day3_ribo_log_cov'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'lightchain',outFile,[0.1,1000],'Riboseq','log')
# outFile = fig_path + '/10_lightchain_day3_ribo_log_cov_align'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'lightchain',outFile,[0.1,1000],'Riboseq','log',15,88)
#  
# gene = 'NeoRKanR'
# outFile = fig_path + '/09_NeoR_day3_ribo_log_cov'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'NeoR',outFile,[0.1,10000],'Riboseq','log')
# outFile = fig_path + '/10_NeoR_day3_ribo_log_cov_align'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'NeoR',outFile,[0.1,10000],'Riboseq','log',16,154)
# #========================= plot the rnaseq coverage ===========================
# covType = 'gene'
# #-------------- 1) get total count --------------------
# os.chdir(rna_bam_path)
# bamFiles = [f for f in os.listdir(rna_bam_path) if f.endswith('.sort.bam')]
# bamFiles = natsorted(bamFiles)
# totalCount = []
# for bam in bamFiles:
#     handle = pysam.AlignmentFile(bam,'rb')
#     totalCount.append(handle.mapped)
# #-------------- 2) read rnaseq cov files ------------------
# os.chdir(rna_cov_path)
# covFiles = [f for f in os.listdir(rna_cov_path) if f.endswith('cov.txt')]
# covFiles = natsorted(covFiles)
# #-------------- 3) plot Eef1a1 rna ------------------
# gene = '100689276'
# # outFile = fig_path + '/07_Eef1a1_day3_rna_cov'
# # plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'Eef1a1',outFile,[],'RNAseq','',0,0)
# outFile = fig_path + '/16_Eef1a1_day6_rna_cov_align'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'Eef1a1',outFile,[],'RNAseq','',0,34)
# #-------------- 4) plot heavychain rna ------------------
# gene = 'heavychain'
# # outFile = fig_path + '/07_heavychain_day3_rna_cov'
# # plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'heavychain',outFile,[],'RNAseq','',0,0)
# outFile = fig_path + '/16_heavychain_day6_rna_cov_align' 
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'heavychain',outFile,[],'RNAseq','',0,233)
# #-------------- 5) plot lightchain rna ------------------
# gene = 'lightchain'
# # outFile = fig_path + '/07_lightchain_day3_rna_cov'
# # plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'lightchain',outFile,[],'RNAseq','',0,0)
# outFile = fig_path + '/16_lightchain_day6_rna_cov_align'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'lightchain',outFile,[],'RNAseq','',0,88)
# #-------------- 6) plot NeoR rna ------------------
# gene = 'NeoRKanR'
# # outFile = fig_path + '/07_NeoR_day3_rna_cov'
# # plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'NeoR',outFile,[],'RNAseq','',0,0)
# outFile = fig_path + '/16_NeoR_day6_rna_cov_align'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'NeoR',outFile,[],'RNAseq','',0,154)
# #---------- log coverage ------
# gene = '100689276' # [0.1,10000],'Riboseq','log'
# outFile = fig_path + '/09_Eef1a1_day3_rna_log_cov_cut'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'Eef1a1',outFile,[0.1,1000],'RNAseq','log',0,0)
# outFile = fig_path + '/10_Eef1a1_day3_rna_log_cov_align_cut'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'Eef1a1',outFile,[0.1,1000],'RNAseq','log',0,34)
#   
# gene = 'heavychain'
# outFile = fig_path + '/09_heavychain_day3_rna_log_cov_cut'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'heavychain',outFile,[0.1,1000],'RNAseq','log',0,0)
# outFile = fig_path + '/10_heavychain_day3_rna_log_cov_align_cut'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'heavychain',outFile,[0.1,1000],'RNAseq','log',0,233)
#   
# gene = 'lightchain'
# outFile = fig_path + '/09_lightchain_day3_rna_log_cov_cut'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'lightchain',outFile,[0.1,1000],'RNAseq','log',0,0)
# outFile = fig_path + '/10_lightchain_day3_rna_log_cov_align_cut'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'lightchain',outFile,[0.1,1000],'RNAseq','log',0,88)
#   
# gene = 'NeoRKanR'
# outFile = fig_path + '/09_NeoR_day3_rna_log_cov_cut'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'NeoR',outFile,[0.1,1000],'RNAseq','log',0,0)
# outFile = fig_path + '/10_NeoR_day3_rna_log_cov_align_cut'
# plot_rna_cov(exnFile,cdsFile,totalCount,covFiles,covType,gene,'NeoR',outFile,[0.1,1000],'RNAseq','log',0,154)
# plt.show()
#===============================================================================
#                         11. coverage plot in sub region (for Eef1a1, heavychain,lightchain,NeoRKanR)
#===============================================================================
def sub_region_plot(cdsFile,exnFile,totalCount,covFiles,cds_pos,utr_pos,covType,gene,refFaFile,outFile1,outFile2):
    """This function plot some sub regions of some genes. Aiming to show the codon periodicity coverage of Riboseq
    
    * cdsFile: str. File name of the cds information
    * exnFile: str. File name of the exon information
    * totalCount: list. A list of integers indicating total number of reads
    * covFiles: list. A list of position covergae for each gene.
    * cds_pos: list. A list of CDS position relative to chromosome.
    * utr_pos: list. A list of utr positions.
    * covType: str. If is gene, it will include the introns.
    * gene: str. gene name
    * refFaFile: str. Fasta file that has reference sequence.
    * outFile1: str. Plot the CDS subregion plot
    * outFile2: str. Plot the 3'UTR subrergion plot
    """
    # get pr pos
    cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
    gene_df = cds_df[cds_df['GeneID'].values==gene]
    gene_obj = trpr(gene_df)
    long_pr = gene_obj.get_longest_trpr(gene)
    pr_pos = gene_obj.get_trpr_pos(long_pr)
    # get gene pos and tr pos
    exn_df = pd.read_csv(exnFile,sep='\t',header=0,low_memory=False)
    gene_df = exn_df[exn_df['GeneID'].values==gene]
    gene_obj = trpr(gene_df)
    trpr_id = list(set(gene_df['access'].tolist()))[0]
    chr_id = list(set(gene_df['chr'].tolist()))[0]
    tr_pos = gene_obj.get_trpr_pos(trpr_id)
    if covType == 'gene':
        if tr_pos[0]>tr_pos[1]:  # - strand
            gene_pos = range(tr_pos[0],tr_pos[-1]-1,-1)
        else:
            gene_pos = range(tr_pos[0],tr_pos[-1]+1)
    #------------- 4) prepare chromosome  ----------------
    refDNA_dic = SeqIO.index(refFaFile,'fasta')
    chr_seq = str(refDNA_dic[chr_id].seq)
    #-------------  4) plot the sub coverage ----------------
    gene_df = pd.DataFrame()
    gene_df['day3'] = merge_rep_gene_pos_cov(totalCount[0:3],covFiles[0:3],gene,exn_df,covType)
    plt_df = gene_df.copy()
    
    cds_chr_pos = [gene_pos[p] for p in cds_pos]
    cds_seq = [chr_seq[p-1].upper() for p in cds_chr_pos]
    frame = (pr_pos.index(cds_chr_pos[0]))%3 + 1
    print 'frame is:', frame
    cds_x = [str(p)+'_'+l for p,l in zip(cds_pos,cds_seq)]
    
    utr_chr_pos = [gene_pos[p] for p in utr_pos]
    utr_seq = [chr_seq[p-1].upper() for p in utr_chr_pos]
    utr_x = [str(p)+'_'+l for p,l in zip(utr_pos,utr_seq)]
    print 'cds sequence is:',cds_seq
    print 'utr sequence is:',utr_seq
    
    plt.figure()
    ax = plt_df.loc[cds_pos]['day3'].plot(kind='bar',title='cds region pattern')
    ax.set_xticklabels(cds_x)
    plt.savefig(outFile1) #fig_path+'/11_Eef1a1_cds_day3_cov.svg')
    plt.figure()
    ax = plt_df.loc[utr_pos]['day3'].plot(kind='bar',title='3UTR region pattern')
    ax.set_xticklabels(utr_x)
    plt.savefig(outFile2) # fig_path+'/11_Eef1a1_3UTR_day3_cov.svg')
# #-------------- 1) get total count --------------------
# os.chdir(bam_path)
# bamFiles = [f for f in os.listdir(bam_path) if f.endswith('.sort.bam')]
# bamFiles = natsorted(bamFiles)
# totalCount = []
# for bam in bamFiles:
#     handle = pysam.AlignmentFile(bam,'rb')
#     totalCount.append(handle.mapped)
# #-------------- 2) get coverage files -----------------------
# os.chdir(fwd_rev_path)
# covFiles = [f for f in os.listdir(fwd_rev_path) if f.endswith('cov.txt')]
# covFiles = natsorted(covFiles)
# #-------------- 3) parameters -------------------------------
# covType = 'gene'
# cds_pos = range(1333,1354) # these position are relative to start of gene
# utr_pos = range(2087,2108)
# gene = '100689276'
# outFile1 = fig_path+'/11_Eef1a1_cds_day3_cov.svg'
# outFile2 = fig_path+'/11_Eef1a1_3UTR_day3_cov.svg'
# sub_region_plot(cdsFile,exnFile,totalCount,covFiles,cds_pos,utr_pos,covType,gene,refFaFile,outFile1,outFile2)
# 
# cds_pos = range(779,800) # these position are relative to start of gene 617,638
# utr_pos = range(1740,1760)
# gene = 'heavychain'
# outFile1 = fig_path+'/11_heavychain_cds_day3_cov.svg'
# outFile2 = fig_path+'/11_heavychain_3UTR_day3_cov.svg'
# sub_region_plot(cdsFile,exnFile,totalCount,covFiles,cds_pos,utr_pos,covType,gene,refFaFile,outFile1,outFile2)
# 
# cds_pos = range(430,451) # these positions are relative to start of gene 617,638
# utr_pos = range(1042,1063)
# gene = 'lightchain'
# outFile1 = fig_path+'/11_lightchain_cds_day3_cov.svg'
# outFile2 = fig_path+'/11_lightchain_3UTR_day3_cov.svg'
# sub_region_plot(cdsFile,exnFile,totalCount,covFiles,cds_pos,utr_pos,covType,gene,refFaFile,outFile1,outFile2)
# 
# cds_pos = range(385,406) # these positions are relative to start of gene 617,638
# utr_pos = range(1020,1041)
# gene = 'NeoRKanR'
# outFile1 = fig_path+'/11_NeoR_cds_day3_cov.svg'
# outFile2 = fig_path+'/11_NeoR_3UTR_day3_cov.svg'
# sub_region_plot(cdsFile,exnFile,totalCount,covFiles,cds_pos,utr_pos,covType,gene,refFaFile,outFile1,outFile2)
# plt.show()
# #===============================================================================
# #                         12. count fraction for different frame to show coverage periodicity (only consider Eef1a1, heavychain,lightchain,NeoR)
# #===============================================================================
#------------------- 1. get read coverage information file ------------------
os.chdir(fwd_rev_path)
covFiles = [f for f in os.listdir(fwd_rev_path) if f.endswith('cov.txt')]
covFiles = natsorted(covFiles)
outpath = signalP_path + '/01_tr_ribo_pos_cov'
outFile1 = signalP_path + '/13_cds_frame_fraction.csv'
outFile2 = signalP_path + '/13_utr5_frame_fraction.csv'
outFile3 = signalP_path + '/13_utr3_frame_fraction.csv'
 
genes = ['100689276','heavychain','lightchain','NeoRKanR']
#------------------- 2. get coverage percentage for each frame in coding region ------------------
# get the position coverage for each transcript
# for f in covFiles:
#     gene_tr_cov(exnFile,cdsFile,all_id_file,f,outpath,genes,'yes',ribo_offset_file)
# get cds positions, whole transcript positions
cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
cds_obj = trpr(cds_df)
exn_df = pd.read_csv(exnFile,sep='\t',header=0,low_memory=False)
exn_obj = trpr(exn_df)
# get frame coverage
os.chdir(outpath)
pr_pos_cov_files = [f for f in os.listdir(outpath) if f.endswith('trpos.txt')]
pr_pos_cov_files = natsorted(pr_pos_cov_files)
# 1) cds frame calculation
frame_df = pd.DataFrame(columns=['frame1','frame2','frame3'])
utr5_df = pd.DataFrame(columns=['frame1','frame2','frame3'])
utr3_df = pd.DataFrame(columns=['frame1','frame2','frame3'])
for f in pr_pos_cov_files:
    handle = open(f,'r')
    # loop for genes
    for line in handle:
        item = line[:-1].split('\t')
        geneid = item[0]
        trid = item[1]
        prid = item[2]
        # get distance between tr and pr
        tr_pos = exn_obj.get_trpr_pos(trid)
        pr_pos = cds_obj.get_trpr_pos(prid)
        s_index = tr_pos.index(pr_pos[0])
        e_index = tr_pos.index(pr_pos[-1])
           
        cov = item[3:]
        # 1) cds coverage calculation
        cds_cov = cov[s_index:e_index+1]
        sum_count = sum([int(p) for p in cds_cov])
        # loop for each position
        frames = frame_count(cds_cov)
        frame_df.loc[geneid+'_'+f[:3]]=frames
        # 2) 5' UTR frame
        utr5_cov = cov[:s_index]
        if len(utr5_cov) == 1:
            utr5_cov = utr5_cov[1:]
        elif len(utr5_cov) == 2:
            utr5_cov = utr5_cov[2:]
        sum5_count = sum([int(p) for p in utr5_cov])
        frames = frame_count(utr5_cov)
        if sum5_count == 0:
            print geneid,'dones not have 5 utr'
            utr5_df.loc[geneid+'_'+f[:3]] = [0] * 3
        else:
            utr5_df.loc[geneid+'_'+f[:3]]=[frame/float(sum5_count) for frame in frames]
        # 3) 3' UTR frame
        utr3_cov = cov[e_index+1:]
        if len(utr3_cov) == 1:
            utr3_cov = utr3_cov[:-1]
        elif len(utr3_cov) == 2:
            utr3_cov = utr3_cov[:-2]
        sum3_count = sum([int(p) for p in utr3_cov])
        frames = frame_count(utr3_cov)
        if sum3_count == 0:
            print geneid,'does not have 3 utr'
            utr3_df.loc[geneid+'_'+f[:3]] = [0] * 3
        else:
            utr3_df.loc[geneid+'_'+f[:3]]=[frame/float(sum3_count) for frame in frames]
   
frame_df = frame_df.sort_index()
print frame_df
# frame_df.to_csv(outFile1,sep='\t')
# 
# utr5_df = utr5_df.sort_index()
# utr5_df.to_csv(outFile2,sep='\t')
#   
# utr3_df = utr3_df.sort_index()
# utr3_df.to_csv(outFile3,sep='\t')
#------------------- 3. plot the frame ----------------------------
def frame_plot(frameFile,outFile,title):
    # 1) plot cds freame
    df = pd.read_csv(frameFile,sep='\t',index_col=0,header=0)
    genes = df.index.tolist()
    # get uniuqe genes
    genes = [g.split('_')[0] for g in genes]
    genes = list(set(genes))
    mean_df = pd.DataFrame()
    std_df = pd.DataFrame()
    for g in genes:
        cri = df.index.map(lambda x: g in x)
        g_df = df[cri]
        g_df = g_df.T # column is gene id of each replicate, row is frame
        mean_df[g] = g_df.mean(axis=1)
        std_df[g] = g_df.std(axis=1)
    mean_df.plot(kind='bar',subplots=True,layout=(2,2),figsize=(10, 8),title=title,yerr=std_df,legend=False)
    plt.savefig(outFile+'.pdf')
    plt.savefig(outFile+'.svg')
# frame_plot(outFile1,fig_path+'/12_cds_frame','CDS frame percentage')
# frame_plot(outFile2,fig_path+'/12_utr5_frame','5UTR frame percentage')
# frame_plot(outFile3,fig_path+'/12_utr3_frame','3UTR frame percentage')
# plt.show()
# #------------------- 4. plot the frame for NeoR for nts that only in main CDS not in alternate CDS ----------------------------
# out_path = '/data/shangzhong/RibosomeProfiling/signalP_part/01_tr_ribo_pos_cov'
# os.chdir(out_path)
# pr_pos_cov_files = [f for f in os.listdir(out_path) if f.endswith('trpos.txt')]
# pr_pos_cov_files = natsorted(pr_pos_cov_files)
def NeoR_frame_cov(pr_pos_cov_files,start,end,add=[],calType='sum'):
    """frame coverage only for NeoR
    * pr_pos_cov_files: list. A list of files with transcript position coverage.
    * start: int. Start position of interested region.
    * end: int. End position of interested region.
    * calType: str. if sum, it plots sum coverage of each frame.
                    if median, it plits median coverage of each frame.
    """
    frame_df = pd.DataFrame(columns=['frame1','frame2','frame3'])
    for f in pr_pos_cov_files:
        handle = open(f,'r')
        for line in handle:
            item = line[:-1].split('\t')
            geneid = item[0]
            if geneid == 'NeoRKanR':
                cov = item[3:]     # get coverage
                cov = [float(p) for p in cov]
                cds_cov = cov[start:end]  # coverage of interested region
                cds_cov = add + cds_cov
                print len(cds_cov)
                sum_count = sum(cds_cov)
                sum_count = sum_count
                if calType != 'sum':
                    sum_count = 1
                frames = frame_count(cds_cov,calType)
                frame_df.loc[geneid+'_'+f[:3]]=[frame/float(sum_count) for frame in frames]
    cov_mean = frame_df.mean().tolist()
    cov_std = frame_df.std().tolist()
    print frame_df
     
    return cov_mean,cov_std
 
# x = [1,2,3]
# tick = ['frame1','frame2','frame3']
#  
# cov_mean,cov_std = NeoR_frame_cov(pr_pos_cov_files,154,167)
# f, ((ax1, ax2,ax3), (ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3)
# f.set_size_inches(18.5, 10.5)
#  
# ax1.bar(x,cov_mean,align='center',yerr=cov_std,color='#AEC7E8')#1F77B4
# ax1.set_xticks(x)
# ax1.set_xticklabels(tick)
# ax1.set_title('NeoR 155:167')
#  
# cov_mean,cov_std = NeoR_frame_cov(pr_pos_cov_files,167,326,[0])
# ax2.bar(x,cov_mean,align='center',yerr=cov_std,color='#AEC7E8')
# ax2.set_xticks(x)
# ax2.set_xticklabels(tick)
# ax2.set_title('NeoR 168:326')
#  
# cov_mean,cov_std = NeoR_frame_cov(pr_pos_cov_files,326,713,[0])
# ax3.bar(x,cov_mean,align='center',yerr=cov_std,color='#AEC7E8')
# ax3.set_xticks(x)
# ax3.set_xticklabels(tick)
# ax3.set_title('NeoR 327:713')
#  
# cov_mean,cov_std = NeoR_frame_cov(pr_pos_cov_files,716,818,[0])
# ax4.bar(x,cov_mean,align='center',yerr=cov_std,color='#AEC7E8')
# ax4.set_xticks(x)
# ax4.set_xticklabels(tick)
# ax4.set_title('NeoR 717:818')
#  
# cov_mean,cov_std = NeoR_frame_cov(pr_pos_cov_files,818,949,[0])
# ax5.bar(x,cov_mean,align='center',yerr=cov_std,color='#AEC7E8')
# ax5.set_xticks(x)
# ax5.set_xticklabels(tick)
# ax5.set_title('NeoR 819:949')
#  
# cov_mean,cov_std = NeoR_frame_cov(pr_pos_cov_files,154,326)
# ax6.bar(x,cov_mean,align='center',yerr=cov_std,color='#AEC7E8')
# ax6.set_xticks(x)
# ax6.set_xticklabels(tick)
# ax6.set_title('NeoR 155:326')
#  
# cov_mean,cov_std = NeoR_frame_cov(pr_pos_cov_files,154,713)
# ax7.bar(x,cov_mean,align='center',yerr=cov_std,color='#AEC7E8')
# ax7.set_xticks(x)
# ax7.set_xticklabels(tick)
# ax7.set_title('NeoR 155:713')
#  
# cov_mean,cov_std = NeoR_frame_cov(pr_pos_cov_files,154,818)
# ax8.bar(x,cov_mean,align='center',yerr=cov_std,color='#AEC7E8')
# ax8.set_xticks(x)
# ax8.set_xticklabels(tick)
# ax8.set_title('NeoR 155:818')
#  
# cov_mean,cov_std = NeoR_frame_cov(pr_pos_cov_files,154,949)
# ax9.bar(x,cov_mean,align='center',yerr=cov_std,color='#AEC7E8')
# ax9.set_xticks(x)
# ax9.set_xticklabels(tick)
# ax9.set_title('NeoR 155:949')
#  
# f.suptitle("frame is based on main ORF", fontsize=14)
# plt.savefig(fig_path + '/13_NeoR_periodicity.svg')
# plt.savefig(fig_path + '/13_NeoR_periodicity.pdf')
# 
# plt.show()

#===============================================================================
#                         13. RNA coverage VS ribo coverage (for rpkm or rpm) scatter plot
#===============================================================================
def rpm_condition(totalCount,count_files,covType='rpm'):
    """This function calculates the rpm for each gene and merge the replicates
    * totalCount: list. A list of total count for replicates
    * count_files: list. A list of file names for replicates
    * covType: str. rpm or rpkm 
    """
    res_df = pd.DataFrame()
    for f,total in zip(count_files,totalCount):
        df = pd.read_csv(f,sep='\t',header=0,index_col=0,names=['geneid','count','len'])
        if covType == 'rpkm':
            df['rpkm'] = (df['count'].div(df['len']))/total*(10**9)
            res_df[f] = df['rpkm']
        elif covType == 'rpm':
            df['rpm'] = df['count']/total*(10**6)
            res_df[f] = df['rpm']
    res_df.index = df.index
    res_df['mean'] = res_df.mean(axis=1)
    return res_df['mean']

def trans_effi_plot(gene_count_path,rna_gene_count_path,bam_path,covType):
    """This function calculates the rpm and rpkm for each gene of RNAseq and Riboseq at day3 and day6,
    and then plot a scatter plot of rna vs ribo.
    
    """
    #---------------- 1. read ribo count files ------------
    ribo_gene_count_files = [f for f in os.listdir(gene_count_path) if f.endswith('Count.txt')]
    ribo_gene_count_files = natsorted(ribo_gene_count_files)
    #---------------- 2. read rna count files -------------
    rna_gene_count_files = [f for f in os.listdir(rna_gene_count_path) if f.endswith('Count.txt')]
    rna_gene_count_files = natsorted(rna_gene_count_files)
    #---------------- 3. get ribo total count -------------
    os.chdir(bam_path)
    bamFiles = [f for f in os.listdir(bam_path) if f.endswith('.bam')]
    bamFiles = natsorted(bamFiles)
    ribo_totalCount = []
    for bam in bamFiles:
        handle = pysam.AlignmentFile(bam,'rb')
        ribo_totalCount.append(handle.mapped)
    #---------------- 4. get rna total count -------------
    os.chdir(rna_bam_path)
    bamFiles = [f for f in os.listdir(rna_bam_path) if f.endswith('.bam')]
    bamFiles = natsorted(bamFiles)
    rna_totalCount = []
    for bam in bamFiles:
        handle = pysam.AlignmentFile(bam,'rb')
        rna_totalCount.append(handle.mapped)
    #---------------- 5. get sp and non sp genes ---------
    gene_class_file = signalP_path + '/05_cho_ab_pr_sp.gene.classify.txt'
    gene_class_df = pd.read_csv(gene_class_file,header=0,sep='\t')
    gene_class_df = gene_class_df.astype(str)
    gene_sp = gene_class_df['with_sp'].tolist()
    gene_sp = list(set(gene_sp))
    gene_no_sp = gene_class_df['no_sp'].tolist()
    gene_no_sp = list(set(gene_no_sp))
    #---------------- 6. get ribo rpm for each condition -----------
    # ribo
    os.chdir(gene_count_path)
    ribo_df = pd.DataFrame()
    ribo_df['ribo_day3'] = rpm_condition(ribo_totalCount[0:3],ribo_gene_count_files[0:3],covType)
    ribo_df['ribo_day6'] = rpm_condition(ribo_totalCount[3:6],ribo_gene_count_files[3:6],covType)
    # rna
    os.chdir(rna_gene_count_path)
    rna_df = pd.DataFrame()
    rna_df['rna_day3'] = rpm_condition(rna_totalCount[0:3],rna_gene_count_files[0:3],covType)
    rna_df['rna_day6'] = rpm_condition(rna_totalCount[3:6],rna_gene_count_files[3:6],covType)
    #----------------- 6. remove non coding rna in rna_df -------------
    ribo_genes = ribo_df.index.tolist()
    rna_genes = rna_df.index.tolist()
    print [f for f in ribo_genes if f not in rna_genes]
    rna_df = rna_df[rna_df.index.isin(ribo_genes)]
    ribo_rna_rpm = pd.concat([ribo_df,rna_df],axis=1,join='inner')
    #----------------- 7. add translation efficiency if it is rpkm --------
    if covType == 'rpkm':
        ribo_rna_rpm['trans_eff3'] = ribo_rna_rpm['ribo_day3']/ribo_rna_rpm['rna_day3']
        ribo_rna_rpm['trans_eff6'] = ribo_rna_rpm['ribo_day6']/ribo_rna_rpm['rna_day6']
        ribo_rna_rpm.to_csv(signalP_path+'/08_ribo_rna_rpkm_trans_eff.csv',sep='\t')
        # get log2 value
        rpkm_df = ribo_rna_rpm[['trans_eff3','trans_eff6']]
        rpkm_df = rpkm_df.replace([np.inf, -np.inf], np.nan)
        rpkm_df = rpkm_df.dropna()
        rpkm_df = rpkm_df + 0.001
        rpkm_df = rpkm_df.apply(np.log2)
        # plot day 3
        ax = rpkm_df['trans_eff3'].plot(kind='hist',title='day3 translation efficiency',bins=100,color='k',alpha=0.5)
        ax.set_ylabel('Frequency')
        ax.set_xlabel('log2 translation efficiency')
        ax.set_xlim([-8,8])
        ax.set_ylim([0,1200])
          
        ax.axvline(rpkm_df.loc['heavychain','trans_eff3'], color='red', linestyle='--',label='hevaychain')
        ax.axvline(rpkm_df.loc['lightchain','trans_eff3'], color='blue', linestyle='--',label='lightchain')
        ax.axvline(rpkm_df.loc['NeoRKanR','trans_eff3'], color='green', linestyle='--',label='NeoR')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:-1], labels[:-1])
        plt.savefig(fig_path + '/05_day3_trans_effi.pdf')
        plt.savefig(fig_path + '/05_day3_trans_effi.svg')
        # plot day 6
        ax = plt.figure()
        ax = rpkm_df['trans_eff6'].plot(kind='hist',title='day6 translation efficiency',bins=100,color='k',alpha=0.5)
        ax.set_ylabel('Frequency')
        ax.set_xlabel('log2 translation efficiency')
        ax.set_xlim([-8,8])
        ax.set_ylim([0,1200])
          
        ax.axvline(rpkm_df.loc['heavychain','trans_eff3'], color='red', linestyle='--',label='hevaychain')
        ax.axvline(rpkm_df.loc['lightchain','trans_eff3'], color='blue', linestyle='--',label='lightchain')
        ax.axvline(rpkm_df.loc['NeoRKanR','trans_eff3'], color='green', linestyle='--',label='NeoR')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:-1], labels[:-1])
        plt.savefig(fig_path + '/05_day6_trans_effi.pdf')
        plt.savefig(fig_path + '/05_day6_trans_effi.svg')
    else:
        ribo_rna_rpm.to_csv(signalP_path+'/08_ribo_rna_rpm.csv',sep='\t')      
    #---------- get log value ------------
    ribo_rna_rpm = ribo_rna_rpm + 0.001
    ribo_rna_rpm = ribo_rna_rpm.apply(np.log10)
    #----------------- 8. plot the results --------------------------
    sp_color ='green'; no_sp_color='red'
    heavy_color = 'yellow'; light_color = 'blue';neoR_color = 'black'
    ax = ribo_rna_rpm[ribo_rna_rpm.index.isin(gene_no_sp)].plot(kind='scatter',x='ribo_day3',y='rna_day3',color=no_sp_color,label='no_sp',title='day3 rna VS ribo (log10)')
    ribo_rna_rpm[ribo_rna_rpm.index.isin(gene_sp)].plot(kind='scatter',x='ribo_day3',y='rna_day3',color=sp_color,label='sp',ax=ax)
    ribo_rna_rpm[ribo_rna_rpm.index.isin(['heavychain'])].plot(kind='scatter',x='ribo_day3',y='rna_day3',color=heavy_color,label='heavy chain',ax=ax)
    ribo_rna_rpm[ribo_rna_rpm.index.isin(['lightchain'])].plot(kind='scatter',x='ribo_day3',y='rna_day3',color=light_color,label='light chain',ax=ax)
    ribo_rna_rpm[ribo_rna_rpm.index.isin(['NeoRKanR'])].plot(kind='scatter',x='ribo_day3',y='rna_day3',color=neoR_color,label='neoR',ax=ax)
      
    lim = [-4,8]
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel('ribo_day3_' + covType)
    ax.set_ylabel('rna_day3_' + covType)
    ax.legend(loc='upper left')
      
    plt.savefig(fig_path + '/04_day3_rna_vs_ribo_' + covType + '.pdf')
    plt.savefig(fig_path + '/04_day3_rna_vs_ribo_' + covType + '.svg')
       
    ax = ribo_rna_rpm[ribo_rna_rpm.index.isin(gene_no_sp)].plot(kind='scatter',x='ribo_day6',y='rna_day6',color=no_sp_color,label='no_sp',title='day6 rna VS ribo (log10)')
    ribo_rna_rpm[ribo_rna_rpm.index.isin(gene_sp)].plot(kind='scatter',x='ribo_day6',y='rna_day6',color=sp_color,label='sp',ax=ax)
    ribo_rna_rpm[ribo_rna_rpm.index.isin(['heavychain'])].plot(kind='scatter',x='ribo_day6',y='rna_day6',color=heavy_color,label='heavy chain',ax=ax)
    ribo_rna_rpm[ribo_rna_rpm.index.isin(['lightchain'])].plot(kind='scatter',x='ribo_day6',y='rna_day6',color=light_color,label='light chain',ax=ax)
    ribo_rna_rpm[ribo_rna_rpm.index.isin(['NeoRKanR'])].plot(kind='scatter',x='ribo_day6',y='rna_day6',color=neoR_color,label='neoR',ax=ax)
      
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel('ribo_day6_' + covType)
    ax.set_ylabel('rna_day6_' + covType)
    ax.legend(loc='upper left')
      
    plt.savefig(fig_path + '/04_day6_rna_vs_ribo_' + covType + '.pdf')
    plt.savefig(fig_path + '/04_day6_rna_vs_ribo_' + covType + '.svg')
    #plt.show()
# trans_effi_plot(gene_count_path,rna_gene_count_path,bam_path,'rpm')
# trans_effi_plot(gene_count_path,rna_gene_count_path,bam_path,'rpkm')

#===============================================================================
#                         14. RNA coverage VS ribo coverage (only for sp genes and in percentage)
#             It plots scatter plots for day3 and day6. In the figure,each dot represents one gene, 
#             the axis values means how many perentages
#             does each gene take in the whole rna or ribo pool.
#===============================================================================
# #---------------- 1. read ribo count files ------------
# ribo_gene_count_files = [f for f in os.listdir(gene_count_path) if f.endswith('Count.txt')]
# ribo_gene_count_files = natsorted(ribo_gene_count_files)
# #---------------- 2. read rna count files -------------
# rna_gene_count_files = [f for f in os.listdir(rna_gene_count_path) if f.endswith('Count.txt')]
# rna_gene_count_files = natsorted(rna_gene_count_files)
# #---------------- 3. get sp and non sp genes ---------
# gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/05_cho_ab_pr_sp.gene.classify.txt'
# gene_class_df = pd.read_csv(gene_class_file,header=0,sep='\t')
# gene_class_df = gene_class_df.astype(str)
# gene_sp = gene_class_df['with_sp'].tolist()
# gene_sp = list(set(gene_sp))
# #gene_sp.append('NeoRKanR')
# #---------------- 4. read sp genes count ---------------
# def perc_condition(count_files,sp_genes):
#     """This function calculates the percentage for each sp genes and merge the replicates
#     * count_files: list. A list of file names for replicates
#     """
#     res_df = pd.DataFrame()
#     for f in count_files:
#         df = pd.read_csv(f,sep='\t',header=0,index_col=0,names=['geneid','count','len'])
#         df = df[df.index.isin(sp_genes)]
#         res_df[f] = df['count']
#     totalCount = res_df.sum().tolist()
#     res_df = res_df.div(totalCount)
#     res_df.index = df.index
#     res_df['mean'] = res_df.mean(axis=1)
#     return res_df['mean']
# # get ribo percentage
# os.chdir(gene_count_path)
# ribo_df = pd.DataFrame()
# ribo_df['ribo_day3'] = perc_condition(ribo_gene_count_files[0:3],gene_sp)
# ribo_df['ribo_day6'] = perc_condition(ribo_gene_count_files[3:6],gene_sp)
# # get rna percentage
# os.chdir(rna_gene_count_path)
# rna_df = pd.DataFrame()
# rna_df['rna_day3'] = perc_condition(rna_gene_count_files[0:3],gene_sp)
# rna_df['rna_day6'] = perc_condition(rna_gene_count_files[3:6],gene_sp)
# #----------------- 5. remove non coding rna in rna_df -------------
# ribo_genes = ribo_df.index.tolist()
# rna_genes = rna_df.index.tolist()
# print [f for f in ribo_genes if f not in rna_genes]
# rna_df = rna_df[rna_df.index.isin(ribo_genes)]
# ribo_rna_rpm = pd.concat([ribo_df,rna_df],axis=1,join='inner')
#  
# ribo_rna_rpm = ribo_rna_rpm + 0.0000001
# ribo_rna_rpm = ribo_rna_rpm.apply(np.log10)
# #----------------- 6. plot the percentage -------------
# sp_color ='green'
# heavy_color = 'red'; light_color = 'blue';neoR_color = ''
# ax = ribo_rna_rpm[ribo_rna_rpm.index.isin(gene_sp)].plot(kind='scatter',x='ribo_day3',y='rna_day3',color=sp_color,label='sp',title='day3 rna VS ribo (log10)')
# ribo_rna_rpm[ribo_rna_rpm.index.isin(['heavychain'])].plot(kind='scatter',x='ribo_day3',y='rna_day3',color=heavy_color,label='heavy chain',ax=ax)
# ribo_rna_rpm[ribo_rna_rpm.index.isin(['lightchain'])].plot(kind='scatter',x='ribo_day3',y='rna_day3',color=light_color,label='light chain',ax=ax)
# #ribo_rna_rpm[ribo_rna_rpm.index.isin(['NeoRKanR'])].plot(kind='scatter',x='ribo_day3',y='rna_day3',color=neoR_color,label='neoR',ax=ax)
#    
# lim = [-8,0.1]
# ax.set_xlim(lim)
# ax.set_ylim(lim)
# ax.set_xlabel('ribo_day3_percentage')
# ax.set_ylabel('rna_day3_percentage')
# ax.legend(loc='upper left')
# 
# plt.savefig(fig_path + '/06_day3_rna_vs_ribo_percent.pdf')
# plt.savefig(fig_path + '/06_day3_rna_vs_ribo_percent.svg')
#     
# ax = ribo_rna_rpm[ribo_rna_rpm.index.isin(gene_sp)].plot(kind='scatter',x='ribo_day6',y='rna_day6',color=sp_color,label='sp',title='day6 rna VS ribo (log10)')
# ribo_rna_rpm[ribo_rna_rpm.index.isin(['heavychain'])].plot(kind='scatter',x='ribo_day6',y='rna_day6',color=heavy_color,label='heavy chain',ax=ax)
# ribo_rna_rpm[ribo_rna_rpm.index.isin(['lightchain'])].plot(kind='scatter',x='ribo_day6',y='rna_day6',color=light_color,label='light chain',ax=ax)
# #ribo_rna_rpm[ribo_rna_rpm.index.isin(['NeoRKanR'])].plot(kind='scatter',x='ribo_day6',y='rna_day6',color=neoR_color,label='neoR',ax=ax)
#    
# ax.set_xlim(lim)
# ax.set_ylim(lim)
# ax.set_xlabel('ribo_day6_percentage')
# ax.set_ylabel('rna_day6_percentage')
# ax.legend(loc='upper left')
# 
# plt.savefig(fig_path + '/06_day6_rna_vs_ribo_percent.pdf')
# plt.savefig(fig_path + '/06_day6_rna_vs_ribo_percent.svg')  
# plt.show()
#===============================================================================
#                         15. merge the stallsites for recombinant genes
#===============================================================================
# os.chdir(stall_site_path)
# files = [f for f in os.listdir(stall_site_path) if f.endswith('stallsites.txt')]
# files = natsorted(files)
# genes = ['heavychain','lightchain','NeoRKanR']
# 
# dfs = []
# for f in files:
#     df = pd.read_csv(f,sep='\t',header=0,names=['Chr','GeneID','Chr_pos','Pr_Access','Pr_pos',f[:3]])
#     df = df[df['GeneID'].isin(genes)]
#     df = df.set_index(['Chr','GeneID','Chr_pos','Pr_Access','Pr_pos'])
#     dfs.append(df)
# merge_df = pd.concat(dfs,axis=1)
# merge_df.to_csv(signalP_path + '/15_recom_genes_stall_sties.csv',sep='\t')
# print merge_df

# #===============================================================================
# #                         16. plot frame coverage
# #===============================================================================
# os.chdir(frame_cov_path)
# frame_files = [f for f in os.listdir(frame_cov_path) if f.endswith('frame.txt')]
# frame_files = natsorted(frame_files)
# frame_df = []
# row = 2; col = 3
# fig,axes = plt.subplots(row,col,sharex=True,figsize=(14.5,8))
# names = []
# for f in frame_files:
#     df = pd.read_csv(f,sep='\t',header=0,index_col=0)
#     frame_df.append(df)
#     names.append(f[:3])
# for i in range(row):
#     for j in range(col):
#         df = frame_df[i*col+j]
#         df.plot(kind='box',ax=axes[i,j],title=names[i*col+j])
#         cri = df.index.map(lambda x: x in ['NeoRKanR','heavychain','lightchain'])
#         recom_df = df[cri]
#         columns = recom_df.columns
#         num = len(columns)
#         for n in range(num):
#             p1,p2,p3 = axes[i,j].plot([n+1.1],recom_df.iloc[0,n],'r.',[n+1.05],recom_df.iloc[1,n],'g.',[n+0.95],recom_df.iloc[2,n],'b.')
# fig.legend((p1,p2,p3),('NeoR','heavy','light'),'upper left')
# outFile = fig_path + '/17_frame_cov'
# plt.savefig(outFile+'.svg')
# plt.savefig(outFile+'.pdf')
# # plt.show()















########## rewrite here ###############
#===============================================================================
#                          15. plot metagene coverage for sp, non sp, antibody at codon level, 5 needs minutes.
#===============================================================================
def groupCodonCovRep(coverFiles,chrCoverFiles,window,genes='',calType='total'):
    # 1. get totalcount
    totalCount = []
    for f in chrCoverFiles:
        df = pd.read_csv(f,header=None,names=['coverage','Chr','pos'],delim_whitespace=True,low_memory=False)
        total = df['coverage'].sum()
        totalCount.append(total)
    # 2. get coverage for all genes
    n = 0
    cal_df = pd.DataFrame()
    std_df = pd.DataFrame()
    for f,total in zip(coverFiles,totalCount):
        n = n + 1
        cov_df = pd.DataFrame()
        handle = open(f,'r')
        for line in handle:
            item = line[:-1].split('\t')
            if item[0] not in genes: continue
            cov = [int(p) for p in item[2:]]
            # if not calculate total, should filter the median
            if calType != 'total':
                if np.median(np.array(cov)) < 1: continue
            if len(cov) < window:
                cov = cov + ['-']*(window-len(cov))
            else:
                cov = cov[0:window]
            cov_df[item[0]] = pd.Series(cov)
        cov_df = cov_df.T.replace('-',np.nan).dropna().astype('int').T
        cov_df = cov_df/float(total)*(10**6)   # rpm
        if calType == 'total':
            cal_df['rep_'+str(n)] = cov_df.sum(axis=1)
        elif calType == 'mean':
            cal_df['rep_'+str(n)] = cov_df.mean(axis=1)
        elif calType =='median':
            cal_df['rep_'+str(n)] = cov_df.median(axis=1)
        std_df['rep_'+str(n)] = cov_df.std(axis=1)
    res_df = pd.DataFrame()
    res_df[calType] = cal_df.mean(axis=1)  # row is the position
    res_df['std'] = cal_df.std(axis=1)
    
    return res_df
"""
# 1. read coverage file and get total count
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_cov5end'
os.chdir(path)
coverFiles = [os.path.join(path,f) for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
chrCoverFiles = natsorted(coverFiles)
# 2. read position file
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/05_codon_cov'
os.chdir(path)
coverFiles = [os.path.join(path,f) for f in os.listdir(path) if f.endswith('.txt')]
coverFiles = natsorted(coverFiles)
# 3. get genes
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/05_cho_ab_pr_sp.gene.classify.txt'
gene_class_df = pd.read_csv(gene_class_file,header=0,sep='\t')
gene_class_df = gene_class_df.astype(str)
gene_sp = gene_class_df['with_sp'].tolist()
gene_sp = list(set(gene_sp))
gene_no_sp = gene_class_df['no_sp'].tolist()
gene_no_sp = list(set(gene_no_sp))
# 4. plot
calType = 'total'
window = 100
gene_sp.remove('heavychain');gene_sp.remove('lightchain')
gene_no_sp.remove('NeoRKanR')
# day3
day3_sp = groupCodonCovRep(coverFiles[0:3],chrCoverFiles[0:3],window,gene_sp,calType)
day3_nosp = groupCodonCovRep(coverFiles[0:3],chrCoverFiles[0:3],window,gene_no_sp,calType)
day3_heavy = groupCodonCovRep(coverFiles[0:3],chrCoverFiles[0:3],window,['heavychain'],calType)
day3_light = groupCodonCovRep(coverFiles[0:3],chrCoverFiles[0:3],window,['lightchain'],calType)
day3_neoR = groupCodonCovRep(coverFiles[0:3],chrCoverFiles[0:3],window,['NeoRKanR'],calType)
# day6
day6_sp = groupCodonCovRep(coverFiles[3:6],chrCoverFiles[3:6],window,gene_sp,calType)
day6_nosp = groupCodonCovRep(coverFiles[3:6],chrCoverFiles[3:6],window,gene_no_sp,calType)
day6_heavy = groupCodonCovRep(coverFiles[3:6],chrCoverFiles[3:6],window,['heavychain'],calType)
day6_light = groupCodonCovRep(coverFiles[3:6],chrCoverFiles[3:6],window,['lightchain'],calType)
day6_neoR = groupCodonCovRep(coverFiles[3:6],chrCoverFiles[3:6],window,['NeoRKanR'],calType)

plt.subplot('511')
day3_sp_plt, = plt.plot(day3_sp.index,day3_sp[calType],color='r',label='day3',linewidth=2.0)
day6_sp_plt, = plt.plot(day6_sp.index,day6_sp[calType],color='b',label='day6',linewidth=2.0)
plt.legend(handles=[day3_sp_plt,day6_sp_plt])
plt.title('signal peptide')

plt.subplot('512')
day3_nosp_plt, = plt.plot(day3_sp.index,day3_nosp[calType],color='r',label='day3_nosp',linewidth=2.0)
day6_nosp_plt, = plt.plot(day6_sp.index,day6_nosp[calType],color='b',label='day6_nosp',linewidth=2.0)
plt.title('non signal peptide')

plt.subplot('513')
heavy3, = plt.plot(day3_heavy.index,day3_heavy[calType],color='r',label='day3_heavy',linewidth=2.0)
heavy6, = plt.plot(day6_heavy.index,day6_heavy[calType],color='b',label='day6_heavy',linewidth=2.0)
plt.title('heavychain')

plt.subplot('514')
light3, = plt.plot(day3_light.index,day3_light[calType],color='r',label='day3_light',linewidth=2.0)
light6, = plt.plot(day6_light.index,day6_light[calType],color='b',label='day6_light',linewidth=2.0)
plt.title('lightchain')

plt.subplot('515')
neo3, = plt.plot(day3_neoR.index,day3_neoR[calType],color='r',label='day3_NeoR',linewidth=2.0)
neo6, = plt.plot(day6_neoR.index,day6_neoR[calType],color='b',label='day6_NeoR',linewidth=2.0)
plt.title('NeoR')
plt.xlabel('Distance to the Tanslation start site')
plt.savefig('/data/shangzhong/RibosomeProfiling/figures/14_codon_cov_all.svg')
plt.show()
"""











#===============================================================================
#                         18. find longest transcripts and longest proteins consistancy
#===============================================================================
"""
# 1. get multiple protein encoded genes
gene_multiPr = '/data/shangzhong/RibosomeProfiling/cho_pr/11_gene_multiPrs.txt'
handle = open(gene_multiPr,'r')
genes = []
for line in handle:
    item = line[:-1].split('\t')
    genes.append(item[0])
# 2. get longest protein and transcript
# cds pr file
rna = [];prs=[];longest_test=[]  # list which will be output to files
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt'
cds_pr_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
cds_pr_df['cds_start'] = cds_pr_df['cds_start']-1
gene_pr_df = cds_pr_df[['GeneID','Pr_Access']].drop_duplicates()
gene_pr_dict = {k:list(v) for k,v in gene_pr_df.groupby('GeneID')['Pr_Access']}  # {gene:[pr1,pr2]}
# get {tr:pr} dictionary
all_id_file = '/data/shangzhong/RibosomeProfiling/Database/combined_AllIDS.txt'
all_id_df = pd.read_csv(all_id_file,sep='\t',header=0)
tr_pr_dic = all_id_df.set_index('TrAccess')['PrAccess'].to_dict()
# exon tr file
exonFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_rna_exon.txt'
exon_tr_df = pd.read_csv(exonFile,sep='\t',header=0,low_memory=False)
exon_tr_df['ex_start'] = exon_tr_df['ex_start'] -1
#####!!!!! should only include protein coding transcipts
cir = exon_tr_df['Tr_Access'].map(lambda x: x in tr_pr_dic.keys())
exon_tr_df = exon_tr_df[cir]
gene_tr_df = exon_tr_df[['GeneID','Tr_Access']].drop_duplicates()
gene_tr_dict = {k:list(v) for k,v in gene_tr_df.groupby('GeneID')['Tr_Access']}
for gene in genes:
    tr_pos,chrome,tr = posLongestPr(gene,exon_tr_df,gene_tr_dict,idType='tr')
    pr_pos,chrome,pr = posLongestPr(gene,cds_pr_df,gene_pr_dict,idType='pr')
    rna.append(tr);prs.append(pr)
    if tr_pr_dic[tr] == pr: 
        longest_test.append('Y')
    else:
        longest_test.append('N')
# results
res_df = pd.DataFrame()
res_df['GeneID'] = pd.Series(genes)
res_df['Longest_Trid'] = pd.Series(rna)
res_df['Longest_Prid'] = pd.Series(prs)
res_df['Tr_Pr_consistency'] = pd.Series(longest_test)
res_df.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/19_longest_tr_pr.txt',sep='\t',index=False)
"""

#===============================================================================
#                         19. detect 5'TOP mRNA
#===============================================================================
"""
# before this part should run part 4 in 05_compareDay3Day6 first
# 1. get all 5'TOP mRNA
mRNAfile = '/data/shangzhong/RibosomeProfiling/Database/combined_mRNA.fa'
mRNA_handle = SeqIO.index(mRNAfile, "fasta")
rna_access = []
sep = 6
for item in mRNA_handle:
    rnaSeq = mRNA_handle[item].seq
    rnaSeq = rnaSeq.upper()
    if rnaSeq.startswith('C'):  # must start with C
        first5 = rnaSeq[1:sep]
        if 'G' in first5 or 'A' in first5: # 2 to 6 must be CT
            continue
        else:
            second = rnaSeq[sep:15]
            m = 0
            for letter in second:
                if letter == 'G' or letter =='A':
                    m = m + 1
            if m/float(len(second)) > 1:
                continue
            else:
                rna_access.append(item)
            
# print '\n'.join(rna_access)
print '# of total 5 TOP mRNA',len(rna_access)
# 2. build rna_gene dictionary
allID_df = pd.read_csv(mapFile,sep='\t',header=0)   # GeneID    GeneSymbol    Chrom    TrAccess    PrAccess
rna_gene_dic = allID_df.set_index('TrAccess')['GeneID'].to_dict()
# 3. get 5'TOP mRNA genes
genes = []
for rna in rna_access:
    try:
        genes.append(rna_gene_dic[rna])
    except:
        print rna,'not in the dict'
genes = list(set(genes))
print '# of total 5 TOP gene',len(genes)
# 4. read the off diagnal genes
off_file = '/data/shangzhong/RibosomeProfiling/cho_pr/15_TransEffiOffDiag.csv'
off_df = pd.read_csv(off_file,sep='\t',header=0)
down_genes = off_df['riboDown'].dropna().tolist()
print '# of ribo down genes:',len(down_genes)
# find top and non top rna
top_genes = []; nontop_genes = []
for gene in down_genes:
    if gene in genes:
        top_genes.append(gene)
    else:
        nontop_genes.append(gene)
print '# of top genes is:',len(top_genes)
print '# of non top genes is:',len(nontop_genes)
"""
#===============================================================================
#                         20. Secondary structure analysis
#===============================================================================
# rnaSecStrucFile = '/data/shangzhong/RibosomeProfiling/cho_pr/20_rnaSecondaryStrucute.txt'
# handle = open(rnaSecStrucFile,'r')
# rna = []; energy = []; pos = []
# for line in handle:
#     if line[0].isalpha():
#         continue
#     if line.startswith('>'):
#         access = line.split('')[1:]  # accession number of mRNA
#         rna.append(access)
#         continue
#     else:
#         item = line[:-1].split('\t')
        
    

# #===============================================================================
# #                         6. get gene length stats
# #===============================================================================
# sp_nosp_gene_file = signalP_path + '/08_sp_no_sp_genes.txt'
# cdsBedFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_cho_ab_pr_sp.gene.CDS.bed'
# gene_df = pd.read_csv(sp_nosp_gene_file,sep='\t',header=0)
# gene_sp = gene_df['gene_sp'].dropna().astype(str).tolist()
# gene_no_sp = gene_df['gene_no_sp'].astype(str).tolist()
# df = pd.DataFrame()
# res = geneCDSStats(cdsBedFile,gene_no_sp)
# df['gene_no_sp'] = pd.Series(res)
# res = geneCDSStats(cdsBedFile,gene_sp)
# df['gene_sp'] = pd.Series(res)
# df = df.fillna(0)
# cds_len_stats = '/data/shangzhong/RibosomeProfiling/cho_pr/09_CDS_len_stats_count.csv'     
# df.to_csv(cds_len_stats,sep='\t',index=False)                                                # 09_CDS_len_stats_count.csv

# #===============================================================================
# #                         9. percentage coverage
# #===============================================================================
# # sp genes and no sp genes
# gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
# df = pd.read_csv(gene_class_file,sep='\t',header=0)
# gene_sp = df['gene_sp'].dropna()
# gene_no_sp = df['gene_no_sp'].dropna()
# # cds coverage
# cdsBedFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_cho_ab_pr_sp.gene.CDS.bed'
# cds_df = pd.read_csv(cdsBedFile,sep='\t',header=None,names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
# #coverFiles = sys.argv[1]
# path = '/data/shangzhong/RibosomeProfiling/Ribo_align'
# os.chdir(path)
# coverFiles = [f for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
# coverFiles = natsorted(coverFiles)
# for coverFile in coverFiles:
#     exon_cov_df = pd.read_csv(coverFile,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
#     # 1) get percent coverage for sp genes
#     gene_sp_cov_df = ribo_percent5endCov(exon_cov_df,cds_df,gene_sp)
#     gene_sp_cov_df.to_csv(coverFile[:-4]+'_sp_percent5end.txt',sep='\t')
#     print 'sp genes done'
#     # 2) get percent coverage for no sp genes
#     gene_no_sp_cov_df = ribo_percent5endCov(exon_cov_df,cds_df,gene_no_sp)
#     gene_no_sp_cov_df.to_csv(coverFile[:-4]+'_no_sp_percent5end.txt',sep='\t')
#     print 'no sp genes done'
#     print coverFile,'done'








