from __future__ import division
import subprocess,os,sys
import pandas as pd
from natsort import natsorted
import numpy as np
import scipy.stats as sp_stats
import matplotlib.pyplot as plt
import matplotlib as mpl
from f02_RiboDataModule import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from multiprocessing import Process
from Modules.p05_ParseGff import extractAllPr
mpl.style.use('ggplot')
sys.path.append('/data/shangzhong/Codes/Pipeline')
#sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
import pdb

bam_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam'
db_path = '/data/shangzhong/RibosomeProfiling/Database'
signalP_path = '/data/shangzhong/RibosomeProfiling/signalP_part'
fig_path = signalP_path + '/figures'
rna_htseq = '/data/shangzhong/RibosomeProfiling/TotalRNA_align/htseq'

refFaFile = db_path + '/combined.fa'
gffFile = db_path + '/combined.gff'
mapFile = db_path + '/combined_AllIDS.txt'
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
# plt.savefig(fig_path+'/figures/02_sp_end_cov.pdf')
# plt.savefig(fig_path+'/02_sp_end_cov.svg')
# #plt.show()

#===============================================================================
#                         8. get RNAseq transcriptional distribution
#===============================================================================
#------- 1. get the genes with sp and without sp --------------------
sp_nosp_gene_file = signalP_path + '/06_sp_no_sp_genes.txt'
gene_df = pd.read_csv(sp_nosp_gene_file,sep='\t',header=0)
gene_sp = gene_df['gene_sp'].dropna().astype(str).tolist()
gene_no_sp = gene_df['gene_no_sp'].dropna().astype(str).tolist()

genes = gene_sp + gene_no_sp
gene_sp.remove('heavychain');gene_sp.remove('lightchain')
gene_no_sp.remove('NeoRKanR')

#------- 2. read the count data for all replicates -------------------
os.chdir(rna_htseq)
countFiles = [f for f in os.listdir(rna_htseq) if f.endswith('Count.txt')]
countFiles = natsorted(countFiles)
dfs = []
for f in countFiles:
    df = pd.read_csv(f,sep='\t',header=None,index_col=0,names=['geneid',f])
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
percent_df.loc[['sp','no_sp','heavychain','lightchain','NeoRKanR','other','total']].to_csv(signalP_path+'/07_rna_distribution.csv',sep='\t')


# # # # # # # # # # # # # # # rewrite here
#===============================================================================
#                         7. plot metagene coverage for sp, non sp, antibody at codon level, 5 needs minutes.
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








