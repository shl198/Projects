from __future__ import division
import os
from natsort import natsorted
import pandas as pd
from f02_RiboDataModule import *
from multiprocessing import Process
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('ggplot')
import matplotlib.ticker as mtick


#===============================================================================
#                         3. check start and stop exon consistancy between proteins from the same gene
#===============================================================================
"""
# 1. read protein gene cds file
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt'
cds_pr_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
# 2. read gene multi proteins file
gene_multiPr = '/data/shangzhong/RibosomeProfiling/cho_pr/11_gene_multiPrs.txt'
handle = open(gene_multiPr,'r')
n = 0
for line in handle:
    fir_exon = pd.DataFrame(columns=['Chr','cds_start','cds_end','GeneID','Pr_Access','Strand']);lst_exon=fir_exon.copy()
    item = line[:-1].split('\t')
    gene = item[0]
    prs = item[1:]
    for pr in prs:
        pr_df = cds_pr_df[cds_pr_df['Pr_Access'].values==pr]
        Strand = pr_df['Strand'].tolist()
        if '-' in Strand:
            fir_exon = fir_exon.append(pr_df.tail(1))
            lst_exon = lst_exon.append(pr_df.iloc[0])
        else:
            lst_exon = lst_exon.append(pr_df.tail(1))
            fir_exon = fir_exon.append(pr_df.iloc[0])
    fir_exon = fir_exon.drop('Pr_Access',axis=1).drop_duplicates()
    lst_exon = lst_exon.drop('Pr_Access',axis=1).drop_duplicates()
#     if fir_exon.shape[0]>1:
#         print gene,'has different start sites'
#         num_start=num_start+1
#     if lst_exon.shape[0]>1:
#         print gene,'has different stop sites'
#         num_stop=num_stop+1
    if (fir_exon.shape[0]>1) or (lst_exon.shape[0]>1):
        n = n + 1
print n
"""
#===============================================================================
#                         7. Ribo seq occupancy distribution mapped from human
#===============================================================================
def splitMouseIDs(gene_id):
    if ';' in gene_id:
        return gene_id.split(';')
    else:
        return gene_id.split(',')
"""
#------------------------------ Human predicted secreted proteins mapped genes --------------------------------
# 1. build human {gene ensembl: gene id} dictionary
id_file = '/data/shangzhong/Database/human/20150522gene2ensembl.human.txt'
id_df = pd.read_csv(id_file,sep='\t',header=None,usecols=[1,2],names=['GeneID','ensembl'])
id_df = id_df.drop_duplicates()
ensembl_id_dic = {k:list(v) for k,v in id_df.groupby('ensembl')['GeneID']}
# 2. get human predicted secreted genes
homo_pr_sp_file = '/data/shangzhong/RibosomeProfiling/cho_pr/18_homo_pr_sp_predict.csv'
homo_pr_sp_df = pd.read_csv(homo_pr_sp_file,sep='\t',header=0,usecols=[2],names=['ensembl'])
homo_pr_sp_genes = homo_pr_sp_df['ensembl'].tolist() 
human_genes = []
for g in homo_pr_sp_genes:
    try:
        human_genes.extend(ensembl_id_dic[g])
    except:
        print g,'not in the {ensembl:geneid} dictionary'
        continue
human_genes = [str(g) for g in human_genes]
# 3. transfer human genes to cho genes
cho2homoFile = '/data/shangzhong/CHO2Human/finalMergeWithmRNA.final.txt'
cho2homo_df = pd.read_csv(cho2homoFile,header=None,sep='\t',names=['cho','homo'])
cho2homo_df = cho2homo_df.astype(str)
cho2homo_df['homo'] = cho2homo_df['homo'].map(lambda x: splitMouseIDs(x))
cho_homo_dict = cho2homo_df.set_index('cho')['homo'].to_dict()
cho_genes = []
for g in cho_homo_dict:
    if set(cho_homo_dict[str(g)]).intersection(human_genes) != set():
        cho_genes.append(g)
"""
"""
# #-------------------------------- test sp genes ---------------------
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
df = pd.read_csv(gene_class_file,sep='\t',header=0)
gene_sp = df['gene_sp'].dropna().astype(str).tolist() 
cho_genes = gene_sp
#---------------------------------------------------------------------------------
# 4. read the gene count file
count_res = pd.DataFrame()
rawCountFile = '/data/shangzhong/RibosomeProfiling/cho_pr/12_pr_gene_rawCount.txt'
raw_df = pd.read_csv(rawCountFile,sep='\t',header=0,index_col=0)
total_count = raw_df.sum().tolist()
# sp genes count
sp_gene_df = raw_df[raw_df.index.isin(cho_genes)]
sp_sum = sp_gene_df.sum().tolist()
count_res['sp'] = pd.Series(sp_sum)
# antibody genes count
ab_gene_df = raw_df[raw_df.index.isin(['heavychain','lightchain','NeoRKanR'])]
count_res.index = ab_gene_df.T.index
count_res = pd.concat([count_res,ab_gene_df.T],axis=1)
# non sp percentage
nosp_df = raw_df[~raw_df.index.isin(cho_genes+['heavychain','lightchain','NeoRKanR'])]
nosp_sum = nosp_df.sum().tolist()
count_res['no_sp'] = pd.Series(nosp_sum,index=ab_gene_df.T.index)
count_res['total'] = count_res.sum(axis=1)
res = count_res.div(count_res['total'],axis='index').T
res['day3'] = res.iloc[:,0:3].mean(axis=1)
res['day6'] = res.iloc[:,3:6].mean(axis=1)
error = pd.DataFrame()
error['day3'] = res.iloc[:,0:3].std(axis=1)
error['day6'] = res.iloc[:,3:6].std(axis=1)

ax = res.loc[:'no_sp','day3':'day6'].plot(kind='bar',title='ribo occupancy distribution',alpha=0.69,yerr=error.loc[:'no_sp','day3':'day6'])
ax.set_xlabel('protein type',fontsize=12,color='black')
ax.set_ylabel('percentage',fontsize=12,color='black')
ax.set_xticklabels(res.index[:-1],rotation=0)
vals = ax.get_yticks()
ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in vals])
plt.savefig('/data/shangzhong/RibosomeProfiling/figures/12_ribo_distribution.svg',dpi=300)
plt.show()
#res.to_csv('/data/shangzhong/ribo_occupancy.csv',sep='\t')
"""
#===============================================================================
#                         8. cds 5'utr 3'utr intron
#===============================================================================
def RiboClassDistribution(gene_class_file,mapFile,cdsFile,exonFile,chr_len_file,coverFile,target):
    # 1. get all proteins. 2. get correspond transcripts. 3. 
    # 1. get all genes
    df = pd.read_csv(gene_class_file,sep='\t',header=0)
    gene_sp = df['gene_sp'].dropna().astype(str).tolist()
    gene_no_sp = df['gene_no_sp'].dropna().astype(str).tolist()
    genes = gene_sp + gene_no_sp
    # 2. get {gene:[proteins]} dictionary
    map_df = pd.read_csv(mapFile,sep='\t',header=0)
    gene_pr_dic = {k:list(v) for k,v in map_df.groupby('GeneID')['PrAccess']}
    # 3. build {protine:transcript dictionary}
    pr_tr_dic = map_df.set_index('PrAccess')['TrAccess'].to_dict()
    # 4. read proteins cds file
    cds_pr_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
    cds_pr_df['cds_start'] = cds_pr_df['cds_start'] -1
    # 5. read transcripts exon file
    exon_tr_df = pd.read_csv(exonFile,sep='\t',header=0,low_memory=False)
    exon_tr_df['ex_start'] = exon_tr_df['ex_start'] - 1
    # 6. read chromosome length file
    df = pd.read_csv(chr_len_file,sep='\t',header=0)
    chr_len_dict = df.set_index('Chr')['Length'].to_dict()
    # 7. read cover files
    f = coverFile
    cov_df = pd.read_csv(f,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
    # 8. define target path
    if not os.path.exists(target):
        os.mkdir(target)
    os.chdir(target)
    # 9. loop for each gene
    count_df = pd.DataFrame()
    cds_count=[];utr5_count=[];utr3_count=[];all_count=[]
    for g in genes:
        proteins = gene_pr_dic[g]
        all_5utr=[]; all_3utr=[];cds=[];all_part=[]
        for p in proteins:
        # 1) get the protein positions and transcript positions
            tr = pr_tr_dic[p] # transcript
            pr_df = cds_pr_df[cds_pr_df['Pr_Access'].values==p]
            tr_df = exon_tr_df[exon_tr_df['Tr_Access'].values==tr]
            pr_pos = getGenepos(pr_df,feature_type='cds')
            tr_pos = getGenepos(tr_df,feature_type='exon')
            if pr_pos[0]>pr_pos[1]:
                pr_pos = range(pr_pos[0]+15,pr_pos[0],-1) + pr_pos[:-15]
                all_pos = range(max(tr_pos[0],pr_pos[0]+15),tr_pos[-1]-1,-1) # the max function make sure the 1st base should include the first 15 nts upstream of TSS
            else:
                pr_pos = range(pr_pos[0]-15,pr_pos[0]) + pr_pos[:-15]
                all_pos = range(min(tr_pos[0],pr_pos[0]-15),tr_pos[-1]+1)  # the min function has the sam function with max
            if pr_pos[0]>pr_pos[1]:
                UTR5 = range(tr_pos[0],pr_pos[0],-1)
                UTR3 = range(pr_pos[-1]-1,tr_pos[-28],-1)
                tr_pos = range(tr_pos[0]+15,tr_pos[0],-1) + tr_pos
            else:
                UTR5 = range(tr_pos[0],pr_pos[0])
                UTR3 = range(pr_pos[-1]+1,tr_pos[-28])
                tr_pos = range(tr_pos[0]-15,tr_pos[0]) + tr_pos
            all_5utr.extend(UTR5)
            all_3utr.extend(UTR3)
            cds.extend(pr_pos)
            dist = min(abs(all_pos[-1]-pr_pos[-1]),28)
            all_part.extend(all_pos[:-dist])
            # remove intersection
            all_5utr = set(all_5utr)-set(all_5utr).intersection(cds)
            all_3utr = set(all_3utr)-set(all_3utr).intersection(cds)
            utr5_3 = set(all_5utr).intersection(all_3utr)
            all_5utr = list(all_5utr-utr5_3)
            all_3utr = list(all_3utr-utr5_3)
            all_part = list(set(all_part)-utr5_3)
            #intr.extend([n for n in all_part if (n not in all_5utr) & (n not in all_3utr) & (n not in cds)])
            all_5utr = sorted(list(set(all_5utr)))  # position 
            all_3utr = sorted(list(set(all_3utr)))
            cds = sorted(list(set(cds)))
            all_part = sorted(list(set(all_part)))
            #intr = sorted(list(set(intr)))
        # 2) get the coverage
        chrom = list(set(pr_df['Chr'].tolist()))[0]
        gene_cov_df = cov_df[cov_df['Chr'].values==chrom]
        cds_cov = getCDSCov(gene_cov_df,cds,chr_len_dict[chrom])
        for i in range(len(cds_cov)):
            if cds_cov[i]=='-': cds_cov[i]= 0
            else: cds_cov[i] = int(cds_cov[i])
        all_cov = getCDSCov(gene_cov_df,all_part,chr_len_dict[chrom])
        for i in range(len(all_cov)):
            if all_cov[i]=='-': all_cov[i]= 0
            else: all_cov[i] = int(all_cov[i])
        utr5_cov = getCDSCov(gene_cov_df,all_5utr,chr_len_dict[chrom])
        utr3_cov = getCDSCov(gene_cov_df,all_3utr,chr_len_dict[chrom])
        #intr_cov = getCDSCov(gene_cov_df,intr,chr_len_dict[chrom])
        # 3) get count data
        cds_count.append(sum(cds_cov))
        utr5_count.append(sum(utr5_cov))
        utr3_count.append(sum(utr3_cov))
        all_count.append(sum(all_cov)-sum(cds_cov)-sum(utr5_cov)-sum(utr3_cov))
        #intron.append(sum(intr_cov))
    count_df['GeneID'] = pd.Series(genes)
    count_df['cds'] = pd.Series(cds_count)
    count_df['5UTR'] = pd.Series(utr5_count)
    count_df['3UTR'] = pd.Series(utr3_count)
    #count_df['intron'] = pd.Series(intron)
    count_df['intron'] = pd.Series(all_count)
    count_df.to_csv(f[:3]+'_distribution.txt',sep='\t',index=False)
"""
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
mapFile = '/data/shangzhong/RibosomeProfiling/Database/combined_AllIDS.txt'
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt'
exonFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_rna_exon.txt'
chr_len_file = '/data/shangzhong/RibosomeProfiling/Database/combined_ChrLen.txt'
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_cov5end'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
coverFiles = natsorted(coverFiles)
target = '/data/shangzhong/RibosomeProfiling/Ribo_align/12_cds_utr_intron'

#RiboClassDistribution(gene_class_file,mapFile,cdsFile,exonFile,chr_len_file,coverFiles[0],target)

proc = []
for f in coverFiles:
    p1 = Process(target=RiboClassDistribution,args=(gene_class_file,mapFile,cdsFile,exonFile,chr_len_file,f,target))    
    p1.start()
    proc.append(p1)
for p in proc:
    p.join()
"""


