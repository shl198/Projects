from __future__ import division
import subprocess,os,sys
import pandas as pd
from natsort import natsorted
import numpy as np
import scipy.stats as sp_stats
import matplotlib.pyplot as plt
import matplotlib as mpl
from f02_RiboDataModule import *
from Bio import SeqIO,Seq
from multiprocessing import Process
from Modules.p05_ParseGff import extractAllPr
mpl.style.use('ggplot')
sys.path.append('/data/shangzhong/Codes/Pipeline')
#sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
import pdb

mapFile = '/data/shangzhong/RibosomeProfiling/Database/combined_AllIDS.txt'
gffFile = '/data/shangzhong/RibosomeProfiling/Database/combined.gff'

#===============================================================================
#                         0. prepare protein files
#===============================================================================
# from Bio import Seq
# filename = '/data/shangzhong/RibosomeProfiling/Database/combined.fa'
# dict = SeqIO.index(filename, "fasta")
# output_handle = open('/data/shangzhong/RibosomeProfiling/cho_pr/00_antibody_pr.fa', "a")
# # get pr sequence
# neo = dict['neoR_contig']
# neo_pr = Seq.translate(neo.seq[154:949])
# neo.seq = neo_pr
# heavy = dict['heavychain_contig']
# heavy_pr = Seq.translate(heavy.seq[233:1637])
# light = dict['lightchain_contig']
# light_pr = Seq.translate(light.seq[88:805])
# SeqIO.write(neo, output_handle, "fasta")
# SeqIO.write(heavy,output_handle,'fasta')
# SeqIO.write(light,output_handle,'fasta')
# output_handle.close()
# #===============================================================================
# #                         1. run signalP
# #===============================================================================
# inputFa = '/data/shangzhong/RibosomeProfiling/cho_pr/01_cho_ab_pr.fa'
# choPrRefseq_sp = signalP(inputFa)
# #===============================================================================
# #                         2. insert gene id and accession of mRNA
# #===============================================================================
# choPrRefseq_sp = '/data/shangzhong/RibosomeProfiling/cho_pr/02_cho_ab_pr_sp.txt'
# choPrRefseq_sp_gene = addGeneIDSymbolChr2signalP(choPrRefseq_sp,mapFile,organism='',mapSource='gff')
# choPrRefseq_sp_gene = changeFileName(choPrRefseq_sp_gene)
# #===============================================================================
# #                         3. get genes mapping to multiple scaffold
# #===============================================================================
# choPrRefseq_sp_gene = '/data/shangzhong/RibosomeProfiling/cho_pr/03_cho_ab_pr_sp.gene.txt'
# choPrRefseq_gene_multichr =  genesMap2diffChrom(choPrRefseq_sp_gene,mapFile)
# choPrRefseq_gene_multichr = changeFileName(choPrRefseq_gene_multichr)
# #===============================================================================
# #                         4. classify gene ids into 3 groups. with sp, without sp, with and without sp
# #===============================================================================
# choPrRefseq_sp_gene = '/data/shangzhong/RibosomeProfiling/cho_pr/03_cho_ab_pr_sp.gene.txt'
# choPrRefseq_sp_gene_classify = gene_sp_classify(choPrRefseq_sp_gene)
# choPrRefseq_sp_gene_classify = changeFileName(choPrRefseq_sp_gene_classify,2)
# #===============================================================================
# #                         5. generate bed file for CDS
# #===============================================================================
# choPrRefseq_sp_gene = '/data/shangzhong/RibosomeProfiling/cho_pr/03_cho_ab_pr_sp.gene.txt'
# df = pd.read_csv(choPrRefseq_sp_gene,header=0,sep='\t')
# genes = df['GeneID'].astype(str).tolist()
# genes = list(set(genes))
# exonBed = choPrRefseq_sp_gene[:-4]+'.bed'
# exonBedFile = gff2Bed4Genes(gffFile,genes,'exon',exonBed)
# outBed = changeFileName(exonBedFile,3)
# cdsBed = choPrRefseq_sp_gene[:-4] + '.bed'
# cdsBedFile = gff2Bed4Genes(gffFile,genes,'CDS',cdsBed)
# outBed = changeFileName(cdsBedFile,3)
#===============================================================================
#      6. get target sp genes and non sp genes, remove those genes 
#      mapping to muliple scaffolds and genes have overlapped cds.
#===============================================================================
"""
# 1). read cds bed file
cdsBedFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_cho_ab_pr_sp.gene.CDS.bed'
cds_df = pd.read_csv(cdsBedFile,header=None,sep='\t',names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
cds_df['GeneID'] = cds_df['GeneID'].astype(str)
cds_df['length'] = cds_df['cds_end'] - cds_df['cds_start'] 
# 2). get gene ids
geneIDs = cds_df['GeneID'].tolist()
geneIDs = list(set(geneIDs))
# 3). get the genes with sp and without sp
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/05_cho_ab_pr_sp.gene.classify.txt'
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
# 4). exclude the genes with multiple scaffolds
gene_multi_file = '/data/shangzhong/RibosomeProfiling/cho_pr/04_cho_ab_pr_sp.gene.multichr.txt'
gene_multi_df = pd.read_csv(gene_multi_file,header=None,sep='\t',engine='python',names=['GeneID','chr1','chr2','chr3','chr4','chr5'])
gene_multi_df['GeneID'] = gene_multi_df['GeneID'].astype(str)
gene_multi_list = gene_multi_df['GeneID'].tolist()
gene_multi_list = list(set(gene_multi_list))
gene_sp = [g for g in gene_sp if g not in gene_multi_list]
gene_no_sp = [g for g in gene_no_sp if g not in gene_multi_list]
print '# of genes with sp after remove multi_chrom:',len(gene_sp)
print '# of genes without sp after remove multi_chrom:',len(gene_no_sp)
# 5). exclude the genes with overlapped CDS
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
gene_sp_df = pd.DataFrame({'gene_sp':gene_sp})
gene_no_sp_df = pd.DataFrame({'gene_no_sp':gene_no_sp})
df = pd.concat([gene_sp_df,gene_no_sp_df],axis=1)
sp_nosp_gene_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
df.to_csv(sp_nosp_gene_file,sep='\t',index=False)
"""
# #===============================================================================
# #                         7. get gene length stats
# #===============================================================================
# sp_nosp_gene_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
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
# df.to_csv(cds_len_stats,sep='\t',index=False)
# #===============================================================================
# #                         8. get cover data for mapping bam files
# #===============================================================================
# path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align'
# os.chdir(path)
# bamFiles = [f for f in os.listdir(path) if f.endswith('.sort.bam')]
# cov5end = pos5Coverage(bamFiles,batch=6)
# # move to a folder
# folder = '01_cov5end'
# if not os.path.exists(folder):
#     os.mkdir(folder)
# for f in cov5end:
#     cmd = ('mv {input} {out}').format(input=f,out=folder)
#     subprocess.call(cmd.split(' '))
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
#===============================================================================
#                         10. calculate coverage around TSS and TSE sites
#===============================================================================
def wrap_covNearTSS_TSE(target_path,exonCovFile,cdsBedFile,genes,chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down,gene_type='sp'):
    """
    This function is a wrapper of covNearTSS_TSE
    """
    #index = exonCovFile.rindex('/')
    f = exonCovFile#[index+1:]
    df_start,df_end = covNearTSS_TSE(exonCovFile,cdsBedFile,genes,chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down)
    if gene_type == 'sp':
        startFile = 'start_sp_' + f[:-16]+'txt';endFile='end_sp_'+f[:-16]+'txt'
    else:
        startFile = 'start_nosp_' + f[:-16]+'txt';endFile='end_nosp_'+f[:-16]+'txt'
    df_start.to_csv(startFile,sep='\t')
    df_end.to_csv(endFile,sep='\t')
    cmd = ('mv {old} {new}').format(old=startFile,new=target_path)
    subprocess.call(cmd.split(' '))
    cmd = ('mv {old} {new}').format(old=endFile,new=target_path)
    subprocess.call(cmd.split(' '))
"""
# 1).read chromosome length
chr_len_file = '/data/shangzhong/RibosomeProfiling/Database/combined_ChrLen.txt'
# 2).sp genes and no sp genes
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
df = pd.read_csv(gene_class_file,sep='\t',header=0)
gene_sp = df['gene_sp'].dropna()
gene_no_sp = df['gene_no_sp'].dropna()
# 3).read cds position file
cdsBedFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_cho_ab_pr_sp.gene.CDS.bed'
# 4). define target folder
target_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/23/02_TSS_TSE_cov'
if not os.path.exists(target_path):
    os.mkdir(target_path)
# 5).define results
res_sp = pd.DataFrame()
res_no_sp = pd.DataFrame()
TSS_up = 50;TSS_down= 60
TSE_up = 50;TSE_down= 60
# 6).parallel
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/23'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
coverFiles = natsorted(coverFiles)
#covNearTSS_TSE(coverFiles[0],cdsBedFile,['100689328'],chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down)
proc = []
for f in coverFiles:
    p1 = Process(target=wrap_covNearTSS_TSE,args=(target_path,f,cdsBedFile,gene_sp,chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down,'sp',))    
    p1.start()
    proc.append(p1)
    p2 = Process(target=wrap_covNearTSS_TSE,args=(target_path,f,cdsBedFile,gene_no_sp,chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down,'nosp',))    
    p2.start()
    proc.append(p2)
for p in proc:
    p.join()
"""
#===============================================================================
#                         11. Get bed file for protein cds. (5th column is protein access)
#===============================================================================
# gffFile = '/data/shangzhong/RibosomeProfiling/Database/combined.gff'
#outFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt'
# extractAllPr(gffFile,feature='CDS',out=outFile)
#===============================================================================
#                         12. coverage around signal peptide end site
#===============================================================================
def wrap_covNearSpEnd(target_path,exonCovFile,id_file,spFile,cdsFile,chr_len_file,gene_sp,up,down):
    """
    This function is a wrapper of covNearSpEnd
    """
    resFile = covNearSpEnd(exonCovFile,id_file,spFile,cdsFile,chr_len_file,gene_sp,up,down)
    cmd = ('mv {old} {new}').format(old=resFile,new=target_path)
    subprocess.call(cmd.split(' '))
"""
# 1).read chromosome length
chr_len_file = '/data/shangzhong/RibosomeProfiling/Database/combined_ChrLen.txt'
# 2) get {gene:chr} from all ID files
id_file = '/data/shangzhong/RibosomeProfiling/Database/combined_AllIDS.txt'
# 3) read gene_sp genes
geneFile = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
gene_df = pd.read_csv(geneFile,sep='\t',header=0)
gene_sp = gene_df['gene_sp'].dropna().astype(str).tolist()
# 4) get sp predicted results
spFile = '/data/shangzhong/RibosomeProfiling/cho_pr/03_cho_ab_pr_sp.gene.txt'
# 5) read 10_pr_cds.txt file
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt'
# 6) read coverage file
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_cov5end'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('5endCov.txt')]
up = 50; down = 50
# 7). define target folder
target_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/03_SpEnd_cov'
if not os.path.exists(target_path):
    os.mkdir(target_path)

#covNearSpEnd(coverFiles[0],id_file,spFile,cdsFile,chr_len_dict,gene_sp,up,down)
# run parrallel
proc = []
for f in coverFiles:
    p1 = Process(target=wrap_covNearSpEnd,args=(target_path,f,id_file,spFile,cdsFile,chr_len_file,gene_sp,up,down))    
    p1.start()
    proc.append(p1)
for p in proc:
    p.join()
"""
#===============================================================================
#                         13. get coverage for all positions of genes
#===============================================================================
def wrap_covAllPos(target_path,exonCovFile,cdsBedFile,geneIDs,chr_len_file,up,down):
    resFile = covAllPos(exonCovFile,cdsBedFile,geneIDs,chr_len_file,up,down)
    cmd = ('mv {old} {new}').format(old=resFile,new=target_path)
    subprocess.call(cmd.split(' '))
"""
# 1).read chromosome length
chr_len_file = '/data/shangzhong/RibosomeProfiling/Database/combined_ChrLen.txt'
# 2).sp genes and no sp genes
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
df = pd.read_csv(gene_class_file,sep='\t',header=0)
gene_sp = df['gene_sp'].dropna().astype(str).tolist()
gene_no_sp = df['gene_no_sp'].dropna().astype(str).tolist()
genes = gene_sp + gene_no_sp
# 3).read cds position file
cdsBedFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_cho_ab_pr_sp.gene.CDS.bed'
# 4) read coverage file
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_cov5end'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('5endCov.txt')]
coverFiles = natsorted(coverFiles)
# 5) define target path
target_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/04_allGenePos_cov'
if not os.path.exists(target_path):
    os.mkdir(target_path)
# get position
up = 50;down=50
#filename = covAllPos(coverFiles[0],cdsBedFile,genes,chr_len_file,up,down)
proc = []
for f in coverFiles:
    p1 = Process(target=wrap_covAllPos,args=(target_path,f,cdsBedFile,genes,chr_len_file,up,down))    
    p1.start()
    proc.append(p1)
for p in proc:
    p.join()
"""

#===============================================================================
#         calculate coverage of house keeping genes at each position
#===============================================================================
# ========== 1. analyze ribo data
"""
# 1) read cds file
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_combined.CDS.bed'
cds_df = pd.read_csv(cdsFile,sep='\t',header=None,names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
# 2) read cov file
path = '/data/shangzhong/RibosomeProfiling/Ribo_align'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('5endCov.txt')]
coverFiles = natsorted(coverFiles)
# 3) define genes
genes = ['heavychain','lightchain','100689477','100765489','100689409','100736557','100757362','100689395',
         '100689467','100754509','100756201','100769768','100757453','100750732']
genes_label = ['heavychain','lightchain','Actb','Actr5','B2m','Gapdh','Gusb','Tfrc','Pgk1','Hsp90ab1','Rplp0','Hprt1','Nono','Sdha']
# 4) loop for each file
gene_ribo_cov_df = pd.DataFrame()
for f in coverFiles:
    gene_cov = []
    cov_df = pd.read_csv(f,header=0,names=['coverage','Chr','pos'],delim_whitespace=True)
    total = cov_df['coverage'].sum()
    for gene in genes:
        gene_cds_df = cds_df[cds_df['GeneID'].values==gene]
        chrome = list(set(gene_cds_df['Chr'].tolist()))
        if len(chrome) > 1:
            assert False,'gene has multi chromosomes'
        pos = getRiboCDSpos(gene_cds_df)
        gene_cov_df = cov_df[cov_df['Chr'].values==chrome[0]]
        cov = getCDSCov(gene_cov_df,pos)
        gene_count = sum(cov)
        gene_rpkm = gene_count*(10**9)/total/len(pos)
        gene_cov.append(gene_rpkm)
    gene_ribo_cov_df[f[:3]] = pd.Series(gene_cov)
gene_ribo_cov_df.index = genes
gene_ribo_cov_df['day3'] = gene_ribo_cov_df.iloc[:,0:3].mean(axis = 1)
gene_ribo_cov_df['day6'] = gene_ribo_cov_df.iloc[:,3:6].mean(axis = 1)
gene_ribo_cov_df['day3_heavy'] = gene_ribo_cov_df.loc['heavychain','day3']/gene_ribo_cov_df['day3']
gene_ribo_cov_df['day3_light'] = gene_ribo_cov_df.loc['lightchain','day3']/gene_ribo_cov_df['day3']
gene_ribo_cov_df['day6_heavy'] = gene_ribo_cov_df.loc['heavychain','day6']/gene_ribo_cov_df['day6']
gene_ribo_cov_df['day6_light'] = gene_ribo_cov_df.loc['lightchain','day6']/gene_ribo_cov_df['day6']
ribo_df = gene_ribo_cov_df.iloc[:,-4:]
ribo_df.index = genes_label
ribo_df = ribo_df[2:]
ribo_df = np.log2(ribo_df)
ax = ribo_df.plot(kind='bar',title='ribosome log2(Ab/house keeping gene)')
ax.set_ylabel('log2 ratio')
#=========== 2. analyze total mRNA =========
# 1) read exon file
exonFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_combined.exon.bed'
exon_df = pd.read_csv(exonFile,sep='\t',header=None,names=['Chr','ex_start','ex_end','GeneID','None','Strand'])
# 2) read cov file
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/total_RNA'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('5endCov.txt')]
coverFiles = natsorted(coverFiles)
# 4) loop for each file
gene_rna_cov_df = pd.DataFrame()
for f in coverFiles:
    gene_cov = []
    cov_df = pd.read_csv(f,header=0,names=['coverage','Chr','pos'],delim_whitespace=True)
    total = cov_df['coverage'].sum()
    for gene in genes:
        gene_exon_df = exon_df[exon_df['GeneID'].values==gene]
        chrome = list(set(gene_exon_df['Chr'].tolist()))
        if len(chrome) > 1:
            assert False,'gene has multi chromosomes'
        pos = getGenepos(gene_exon_df)
        gene_cov_df = cov_df[cov_df['Chr'].values==chrome[0]]
        cov = getCDSCov(gene_cov_df,pos)
        gene_count = sum(cov)
        gene_rpkm = gene_count*(10**9)/total/len(pos)
        gene_cov.append(gene_rpkm)
    gene_rna_cov_df[f[:3]] = pd.Series(gene_cov)
gene_rna_cov_df.index = genes
gene_rna_cov_df['day3'] = gene_rna_cov_df.iloc[:,0:3].mean(axis = 1)
gene_rna_cov_df['day6'] = gene_rna_cov_df.iloc[:,3:6].mean(axis = 1)
gene_rna_cov_df['day3_heavy'] = gene_rna_cov_df.loc['heavychain','day3']/gene_rna_cov_df['day3']
gene_rna_cov_df['day3_light'] = gene_rna_cov_df.loc['lightchain','day3']/gene_rna_cov_df['day3']
gene_rna_cov_df['day6_heavy'] = gene_rna_cov_df.loc['heavychain','day6']/gene_rna_cov_df['day6']
gene_rna_cov_df['day6_light'] = gene_rna_cov_df.loc['lightchain','day6']/gene_rna_cov_df['day6']
rna_df = gene_rna_cov_df.iloc[:,-4:]
rna_df.index = genes_label
rna_df = rna_df[2:]
rna_df = np.log2(rna_df)
ax = rna_df.plot(kind='bar',title='total rna log2(Ab/house keeping gene)')
ax.set_ylabel('log2 ratio')
# 3. translation efficiency
trans_df = pd.DataFrame((gene_ribo_cov_df[['day3','day6']].values) / (gene_rna_cov_df[['day3','day6']].values),columns=['day3','day6'],index=gene_rna_cov_df.index)
trans_df['day3_heavy'] = trans_df.loc['heavychain','day3']/trans_df['day3']
trans_df['day3_light'] = trans_df.loc['lightchain','day3']/trans_df['day3']
trans_df['day6_heavy'] = trans_df.loc['heavychain','day6']/trans_df['day6']
trans_df['day6_light'] = trans_df.loc['lightchain','day6']/trans_df['day6']
trans_df = trans_df.iloc[:,-4:]
trans_df.index = genes_label
trans_df = trans_df[2:]
ax = trans_df.plot(kind='bar',title='translation efficiency (Ab/house keeping gene)')
ax.set_ylabel('ratio')
plt.show()
print 'done'
"""







