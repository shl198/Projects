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
#===============================================================================
#                         1. detect stalling sites
#===============================================================================
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
up = 50; down=50
cdsBedFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_cho_ab_pr_sp.gene.CDS.bed'
# 2).sp genes and no sp genes
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
df = pd.read_csv(gene_class_file,sep='\t',header=0)
gene_sp = df['gene_sp'].dropna().astype(str).tolist() # here gene_sp is used to add the _sp suffix for signal peptide genes

path = '/data/shangzhong/RibosomeProfiling/Ribo_align/04_allGenePos_cov'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('geneAllposCov.txt')]
coverFiles = natsorted(coverFiles)
# 3) run stall sites in parallel
#StallSites(coverFiles[0],cdsBedFile,gene_sp,total,up,down,by='median')
by='median'
proc = []
for f,total in zip(coverFiles,totalCount):
    p1 = Process(target=StallSites,args=(f,cdsBedFile,gene_sp,total,up,down,by))
    p1.start()
    proc.append(p1)
for p in proc:
    p.join()
"""
#===============================================================================
#                         2. get genes encoding multiple proteins
#===============================================================================
"""
# 1. read genes and protein accession numbers
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt'
cds_pr_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False,usecols=[3,4])
cds_pr_df = cds_pr_df.drop_duplicates()
# 2. read target interested genes
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
df = pd.read_csv(gene_class_file,sep='\t',header=0)
gene_sp = df['gene_sp'].dropna().astype(str).tolist()
gene_no_sp = df['gene_no_sp'].dropna().astype(str).tolist()
genes = gene_sp + gene_no_sp
# 3. filter genes
criteria = cds_pr_df['GeneID'].map(lambda x: str(x) in genes)
cds_pr_df = cds_pr_df[criteria]
# 4. output to a file
gene_multiPr = cdsFile[:-10] + 'gene_multiPrs.txt'
handle = open(gene_multiPr,'w')
gene_pr_dict = {k:list(v) for k,v in cds_pr_df.groupby('GeneID')['Pr_Access']}
for key in gene_pr_dict:
    if len(gene_pr_dict[key])>1:
        handle.write(key+'\t' + '\t'.join(gene_pr_dict[key]) + '\n')
handle.close()
gene_multiPr = changeFileName(gene_multiPr,1)   # 11_gene_multiPrs.txt
"""
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
#                         4. get coverage for all position of codons in codon level or nt level
#===============================================================================
# 1. get coverage for each codon of a gene
def codonCovPr(coverFile,cdsFile,chr_len_file,gene_class_file):
    """
    This function calculates codon coverage for each gene in the gene_class_file
    
    * coverFile: str. Chromosome position coverage file. In the folder 01_cov5end
    * cdsFile: str. ['Chr','cds_start','cds_end','GeneID','Pr_Access','Strand']
    * chr_len_file: str. A file with length for each chromosome
    * gene_class_file: str. ['gene_sp','gene_no_sp']
    """
    # 1. read gene with all proteins file
    cds_pr_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False,names=['Chr','cds_start','cds_end','GeneID','Pr_Access','Strand'])
    cds_pr_df['cds_start'] = cds_pr_df['cds_start']-1
    gene_pr_df = cds_pr_df[['GeneID','Pr_Access']].drop_duplicates()
    gene_pr_dict = {k:list(v) for k,v in gene_pr_df.groupby('GeneID')['Pr_Access']}  # {gene:[pr1,pr2]}
    # 2. read chromosome length file
    df = pd.read_csv(chr_len_file,sep='\t',header=0)
    chr_len_dict = df.set_index('Chr')['Length'].to_dict()
    # 3. read target genes
    df = pd.read_csv(gene_class_file,sep='\t',header=0)
    gene_sp = df['gene_sp'].dropna().tolist()
    gene_no_sp = df['gene_no_sp'].dropna().tolist()
    genes = gene_sp + gene_no_sp
    # 4. get position
    cov_df = pd.read_csv(coverFile,header=0,names=['coverage','Chr','pos'],delim_whitespace=True)
    handle = open(coverFile.split('/')[-1][:3]+'.txt','w')
    for gene in genes:
        # (1). get position of longest protein
        pos,chrome,pr = posLongestPr(gene,cds_pr_df,gene_pr_dict)
        # (2). get coverage of the longest protein
        pos = [p-15 for p in pos]
        gene_cov_df = cov_df[cov_df['Chr'].values==chrome[0]]
        if gene_cov_df.empty:
            pr_cov = len(pos)*[0]
        else:
            pr_cov = getCDSCov(gene_cov_df,pos,chr_len_dict[chrome[0]])
        # (3). get coverage for each codon of protein
        for i in range(len(pr_cov)):
            if pr_cov[i]!='-':
                pr_cov[i]=str(pr_cov[i])
        """The next command get coverage at the codon level
        """
        #pr_codon_cov = [str(sum(p)) for p in chunks(pr_cov,3)]  
        handle.write(gene+'\t'+ pr + '\t' + '\t'.join(pr_cov)+'\n')
    handle.close()
"""
# 1. read coverage files
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_cov5end'
os.chdir(path)
coverFiles = [os.path.join(path,f) for f in os.listdir(path) if f.endswith('5endCov.txt')]
# 2. define other files
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt'
chr_len_file = '/data/shangzhong/RibosomeProfiling/Database/combined_ChrLen.txt'
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
# 3. define target pathway
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/10_pr_nt_cov'
if not os.path.exists(path):
    os.mkdir(path)
os.chdir(path)

# 4. parallel
proc = []
for f in coverFiles:
    p1 = Process(target=codonCovPr,args=(f,cdsFile,chr_len_file,gene_class_file))    
    p1.start()
    proc.append(p1)
for p in proc:
    p.join()
"""

#===============================================================================
#                         5. list all the start sites and stop sites in the chromosome for every protein
#===============================================================================
"""
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt'
cds_pr_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
proteins = list(set(cds_pr_df['Pr_Access'].tolist()))
chr_pr_start_end_sites = pd.DataFrame()
Chr = [];GeneID=[];Start=[];Stop=[];Pr=[];Strand=[]
for pr in proteins:
    Pr.append(pr)
    pr_df = cds_pr_df[cds_pr_df['Pr_Access'].values==pr]
    Chr.append(pr_df['Chr'].tolist()[0])
    GeneID.append(pr_df['GeneID'].tolist()[0])
    strand = pr_df['Strand'].tolist()
    if strand[0] == '+':
        Start.append(pr_df['cds_start'].tolist()[0])
        Stop.append(pr_df['cds_end'].tolist()[-1])
        
    else:
        Stop.append(pr_df['cds_start'].tolist()[0])
        Start.append(pr_df['cds_end'].tolist()[-1])
    Strand.append(strand[0])
chr_pr_start_end_sites['Chr']=pd.Series(Chr)
chr_pr_start_end_sites['GeneID']=pd.Series(GeneID)
chr_pr_start_end_sites['Start'] = pd.Series(Start)
chr_pr_start_end_sites['Stop'] = pd.Series(Stop)
chr_pr_start_end_sites['Pr_Access'] = pd.Series(Pr)
chr_pr_start_end_sites['Strand'] = pd.Series(Strand)
chr_pr_start_end_sites.to_csv('/data/shangzhong/RibosomeProfiling/cho_pr/14_chr_pr_start_end.txt',index=False,sep='\t')
"""
#===============================================================================
#                         6. get stats for the stalling sites
#===============================================================================
"""
1. protein cds coverage at nt level.
2. protein cds coverage at codon level
3. AA sequence
4. nt sequence
5. check how many proteins have different sequence with chromosome extracted sequence
6. for stalling sites, remember to remove the start and stop sites of other proteins which are covered by the analyzed transcript
"""

"""
# ======== find the protein sequence inconsistancy between refseq and ref genome =======================
# 1. check the inconsistancy between the sequence and the protein sequence
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt'
cds_pr_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
cds_pr_df['cds_start'] = cds_pr_df['cds_start'] -1
refFaFile = '/data/shangzhong/RibosomeProfiling/Database/combined.fa'
refDNA_dic = SeqIO.index(refFaFile,'fasta')
# 2. read correspond refseq protein sequences
prFaFile = '/data/shangzhong/RibosomeProfiling/Database/combined_pr.fa'
prRef_dic = SeqIO.index(prFaFile,'fasta')
# 3. get genes
gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
df = pd.read_csv(gene_class_file,sep='\t',header=0)
gene_sp = df['gene_sp'].dropna().astype(str).tolist()
gene_no_sp = df['gene_no_sp'].dropna().astype(str).tolist()
genes = gene_sp + gene_no_sp
# 3. get protein positions
criteria = cds_pr_df['GeneID'].map(lambda x: str(x) in genes)
cds_pr_df = cds_pr_df[criteria]
proteins = list(set(cds_pr_df['Pr_Access'].tolist()))
outFile = '/data/shangzhong/RibosomeProfiling/cho_pr/16_prInconsistant.txt'
outFile1 = '/data/shangzhong/RibosomeProfiling/cho_pr/16_prInconsistantID.txt'
handle = open(outFile,'w')
handle1 = open(outFile1,'w')
for pr in proteins:
    pr_df = cds_pr_df[cds_pr_df['Pr_Access'].values==pr]
    # get chrome sequence
    chrome = pr_df['Chr'].tolist()[0]
    chr_seq = str(refDNA_dic[chrome].seq)
    # get pr nt sequence
    pos = getGenepos(pr_df,feature_type='cds')
    pr_nt = Seq(''.join([chr_seq[n-1].upper() for n in pos]),generic_dna)
    # get pr AA sequence
    if pos[0]<pos[1]:
        AA = pr_nt.translate()
    else:
        AA = pr_nt.complement().translate()
    if AA.endswith('*'):
        AA = AA[:-1]
    try:
        if AA != str(prRef_dic[pr].seq):
            handle1.write(pr+'\n')
            handle.write('\n'.join(['>'+pr,str(AA)])+'\n')
    except:
        print pr,'not in the annotation file'
handle.close()
handle1.close()    
"""

#=====================  codon frequencies ======================
"""
# 1. read gene with all proteins file
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt'
cds_pr_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False,names=['Chr','cds_start','cds_end','GeneID','Pr_Access','Strand'])
cds_pr_df['cds_start'] = cds_pr_df['cds_start']-1
gene_pr_df = cds_pr_df[['GeneID','Pr_Access']].drop_duplicates()
gene_pr_dict = {k:list(v) for k,v in gene_pr_df.groupby('GeneID')['Pr_Access']}  # {gene:[pr1,pr2]}
# 2. read genome file
refFaFile = '/data/shangzhong/RibosomeProfiling/Database/combined.fa'
refDNA_dic = SeqIO.index(refFaFile,'fasta')
# 3). pr in consistant id
prInconID = '/data/shangzhong/RibosomeProfiling/cho_pr/16_prInconsistantID.txt'
prInconID_df = pd.read_csv(prInconID,header=None,names=['ID'])
prInconIDs = prInconID_df['ID'].tolist()
# 4) calculate frequencies 
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/04_allGenePos_cov'
os.chdir(path)
stallFiles = [os.path.join(path,f) for f in os.listdir(path) if f.endswith('stallsites.txt')]
stallFiles = natsorted(stallFiles)
target_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/11_codon_AA_freq'
if not os.path.exists(target_path):
    os.mkdir(target_path)
os.chdir(target_path)
for f in stallFiles:
    outFile = f.split('/')[-1][:3] + '_codon_freq.txt'
    outHandle = open(outFile,'w')
    stall_df = pd.read_csv(f,sep='\t',header=0,usecols=[1,2]) # Chr    GeneID    Chr_pos    ratio2median
    stall_df = stall_df.drop_duplicates()
    gene_stall_dic = {k:list(v) for k,v in stall_df.groupby('GeneID')['Chr_pos']}  # {gene: sites}
    genes = list(set(stall_df['GeneID'].tolist()))
    for gene in genes:
        cds_pr_df['len'] = cds_pr_df['cds_end'] - cds_pr_df['cds_start']
        proteins = gene_pr_dict[gene]
        length = []; position = []
        # 1. build the start and stop sites
        start_stop_sites = []
        for pr in proteins:
            if pr in prInconIDs:
                print pr,'sequence is inconsistant with refseq'
                continue
            # get cds positions
            gene_pr_df = cds_pr_df[cds_pr_df['Pr_Access'].values==pr]
            length.append(gene_pr_df['len'].sum())
            pos = getGenepos(gene_pr_df,feature_type='cds')
            if len(pos) < 75: continue
            position.append(pos)
            # get chromosome
            chrome = list(set(gene_pr_df['Chr'].tolist()))[0]
            # buld the dictionary
            start_stop_sites.extend(pos[:45]+pos[-30:])
        start_stop_sites = list(set(start_stop_sites))
        # 2. get the longest cds
        if position == []: continue
        max_len = max(length)
        index = length.index(max_len)
        pos = position[index]
        chr_seq = refDNA_dic[chrome].seq
        pr_nt = Seq(''.join([chr_seq[n-1].upper() for n in pos]),generic_dna)  # protein sequence
        if pos[0]>pos[1]:
            pr_nt = pr_nt.complement()
        AA = pr_nt.translate()
        if AA.endswith('*'):
            AA = AA[:-1]
        # 3. loop for each site
        for stall in gene_stall_dic[gene]:
            if stall in start_stop_sites:
                continue
            try:
                index = pos.index(stall)
                n = int((index+1)/3)
                m = (index+1)%3
                if m !=0: m = 1
                num_codon = n + m # how many codons are there from start to the stalling sites
                pr_nt_seq = pr_nt[:num_codon*3]
                AA = AA[:num_codon]
                outHandle.write(str(pr_nt_seq[-15:])+'\n')
            except:
                print stall, 'not in longest transcript of',gene
    outHandle.close()
"""
#=================== change nt to AA and calculate the t-test===============
"""
from plot_fig import getcodon_AAFreq
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/11_codon_AA_freq'
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
print aa_df[aa_df['p_value'].values>0]
"""






 



