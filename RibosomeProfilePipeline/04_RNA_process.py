import os,sys
sys.path.append('/home/shangzhong/Codes/Pipeline')
from Modules.f04_htseq import htseq_count_py
from natsort import natsorted
from multiprocessing import Process
import subprocess
from Modules.f05_IDConvert import geneSymbol2EntrezID
import pandas as pd
from f02_RiboDataModule import *
import shutil
#===============================================================================
#                     define parameters
#===============================================================================
bam_path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align'
db_path = '/data/shangzhong/RibosomeProfiling/Database'
rna_htseq = bam_path + '/htseq'
count_path = bam_path + '/01_gene_count'
fwd_rev_path = bam_path + '/02_cov'
tr_pos_cov_path = bam_path + '/03_tr_pos_cov'
gene_id_symbol_file = '/data/shangzhong/Database/cho/gff_chok1_ID_symbol.txt'

exnFile = db_path + '/01_pr_rna.txt'
all_id_file = db_path + '/combined_AllIDS.txt'
cdsFile = db_path + '/01_pr_cds.txt'

max_len = 50; seq_len = 50
#===============================================================================
#                         1. get total count for each gene 
#===============================================================================
# # # # # # # # gffFile = db_path + '/combined.gff'    # for this part (use htseq instead of this)
# # # # # # # # os.chdir(bam_path)
# # # # # # # # bamFiles = [f for f in os.listdir(bam_path) if f.endswith('.bam')]
# # # # # # # # bamFiles = natsorted(bamFiles)
# # # # # # # # proc = [Process(target=htseq_count_py,args=(gffFile,bam)) for bam in bamFiles]
# # # # # # # # for p in proc:
# # # # # # # #     p.start()
# # # # # # # # for p in proc:
# # # # # # # #     p.join()
# # # # # # # # os.chdir(bam_path)
# # # # # # # # countFiles = [f for f in os.listdir(rna_htseq) if f.endswith('Count.txt')]
# # # # # # # # countFiles = natsorted(countFiles)
# # # # # # # #------------ gene id mapping ------------
# # # # # # # # geneSymbol2EntrezID(gene_id_symbol_file,count_path,rna_htseq,sym2ID='yes')
#------------ get gene length ------------- 
def gene_len(obj,gene):
    try:
        gene_length = len(obj.get_trpr_pos(gene,level='gene'))
    except:
        gene_length = -1   # gene mapping to multiple chromosomes
    return gene_length
#-------------- 1. read htseq count ------------------
if not os.path.exists(count_path): os.mkdir(count_path)
os.chdir(rna_htseq)
rna_count_files = [f for f in os.listdir(rna_htseq) if f.endswith('Count.txt')]
rna_count_files = natsorted(rna_count_files)
#-------------- 2. read exon position file ------------
exn_file = db_path + '/01_pr_rna.txt'
exn_df = pd.read_csv(exn_file,sep='\t',header=0,low_memory=False)
exn_obj = trpr(exn_df)
 
gene_length = []
for f in rna_count_files:
    df = pd.read_csv(f,sep='\t',header=None,names=['GeneID','count'])
    if gene_length == []:
        gene_length = (df['GeneID'].apply(lambda x: gene_len(exn_obj,str(x)))).tolist()
    else:
        pass
    df['length'] = pd.Series(gene_length)
    df.to_csv(count_path + '/' + f,sep='\t',index=False)

#===============================================================================
#                         2. get forward reverse coverage 
#===============================================================================
# os.chdir(bam_path)
# bamFiles = [f for f in os.listdir(bam_path) if f.endswith('.bam')]
# bamFiles = natsorted(bamFiles)
#  
# proc = [Process(target=fwd_rev_cov,args=(bam,max_len,seq_len,)) for bam in bamFiles]              # in folder 02_cov
# for p in proc:
#     p.start()
# for p in proc:
#     p.join()
# # move file into 02_cov
# fwd_rev_path = bam_path + '/02_cov'
# if not os.path.exists(fwd_rev_path): os.mkdir(fwd_rev_path)
# for f in os.listdir(bam_path):
#     if f.endswith('cov.txt'):
#         if os.path.exists(fwd_rev_path + '/' + f):
#             os.remove(fwd_rev_path + '/' + f)
#         shutil.move(f,fwd_rev_path)
#===============================================================================
#                         3. get position count of each mRNA
#===============================================================================
# #------------- read cover file -------------
# os.chdir(fwd_rev_path)
# covFiles = [f for f in os.listdir(fwd_rev_path) if f.endswith('cov.txt')]
# covFiles = natsorted(covFiles)
# # gene_tr_cov(exnFile,cdsFile,all_id_file,covFiles[0],tr_pos_cov_path)
# proc = [Process(target=gene_tr_cov,args=(exnFile,cdsFile,all_id_file,covFile,tr_pos_cov_path,)) for covFile in covFiles]
# for p in proc:
#     p.start()
# for p in proc:
#     p.join()


    







