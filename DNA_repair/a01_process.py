import sys,os
import pandas as pd
import subprocess
sys.path.append('/home/shangzhong/Codes/Pipeline')

from RibosomeProfilePipeline.f02_RiboDataModule import trpr
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Modules.p05_ParseGff import get_diff_pr_from_refseq
from Bio import Entrez
Entrez.email = 'shl198@eng.ucsd.edu'


old_gffFile = '/data/shangzhong/DNArepair/hamster.gff'
new_gffFile = '/data/shangzhong/DNArepair/super_scaffold/hamster.super.gff'
gene_file = '/data/shangzhong/DNArepair/Database/01_DNA_repair_genes.hamster.txt'
yeast_gene_file = '/data/shangzhong/DNArepair/Database/02_newGenesFromYeastGCR.txt'
all_id_file = '/data/shangzhong/DNArepair/Database/super_scaffold/hamster_super_AllIDS.txt'
old_all_id_file = '/data/shangzhong/DNArepair/Database/hamster_AllIDS.txt'
old_cds = '/data/shangzhong/DNArepair/Database/hamster.cds.txt'
old_ref = '/data/shangzhong/DNArepair/Database/hamster.fa'
new_cds = '/data/shangzhong/DNArepair/Database/super_scaffold/hamster.super.cds.txt'
new_ref = '/data/shangzhong/DNArepair/Database/super_scaffold/hamster.super.fa'

inconsistent_folder = '/data/shangzhong/DNArepair/f03_inconsistent_seq'
#===============================================================================
#                         1. extract the chromosomes of repair genes
#===============================================================================
# gene_df = pd.read_csv(gene_file,header=0,names=['geneid'])
# genes = gene_df['geneid'].tolist()
# 
# id_file = '/data/shangzhong/Database/cho/gff_chok1_all_ID.txt'
# id_df = pd.read_csv(id_file,sep='\t',header=0,low_memory=False)
# chr_df = id_df[id_df['GeneSymbol'].isin(genes)]
# chr_df['Chr'].to_csv('/data/shangzhong/DNArepair/Database/DNA_repair_chrom_hamster.txt',index=False)
#===============================================================================
#                         2. Test the consistancy between new gff and old gff on AA sequence (this is test for super gff and ncbi hamster fa)
#===============================================================================
def AA_sequence(refDNA_dic,cds_df,gene,seq_type='AA'):
    pr_seqs = []
    tr_seqs = []
    # 1. get all proteins
    gene_df = cds_df[cds_df['geneid'].values==gene]
    prs = list(set(gene_df['access'].tolist()))
    prs = sorted(prs)
    obj = trpr(gene_df)
    # 2. loop for each pr
    for pr in prs:
        # 1) get chromosome
        chrom = obj.get_chrom(pr,id_type='access')
        pos = obj.get_trpr_pos(pr)
        ref_seq = refDNA_dic[chrom].seq
        sequence = ''.join([ref_seq[p-1] for p in pos])
        nt_seq = Seq(sequence,generic_dna)
        if pos[0]>pos[1]:
            nt_seq = nt_seq.complement()
        AA = str(nt_seq.translate())
        tr_seqs.append(str(nt_seq))
        pr_seqs.append(AA)
    if seq_type=='AA':
        return pr_seqs,prs
    else:
        return tr_seqs,prs
# # 1. build {gene symbol:gene id dictionary}
# all_id_file = '/data/shangzhong/DNArepair/Database/super_scaffold/hamster_super_AllIDS.txt'
# id_df = pd.read_csv(all_id_file,sep='\t',header=0)
# symbol_id_dic = id_df.set_index('GeneSymbol')['GeneID'].to_dict()
# # 2. get all repair proteins
# gene_df = pd.read_csv(gene_file,header=0,names=['name'])
# gene_symbols = gene_df['name'].tolist()
# gene_ids = [symbol_id_dic[name] for name in gene_symbols]
# # 3. get protein sequence
# old_refDNA_dic = SeqIO.index(old_ref, "fasta")
# new_refDNA_dic = SeqIO.index(new_ref, "fasta")
#  
# columns = ['chr','start','end','geneid','access','strand']
# old_cds_df = pd.read_csv(old_cds,sep='\t',header=0,names=columns)
# new_cds_df = pd.read_csv(new_cds,sep='\t',header=0,names=columns)
# for g in gene_ids:
#     old_seq = AA_sequence(old_refDNA_dic,old_cds_df,g)
#     new_seq = AA_sequence(new_refDNA_dic,new_cds_df,g)
#     if old_seq != new_seq:
#         print g, 'has different sequence between old and new'
#===============================================================================
#                         3. Get all AA sequence of repair proteins
#===============================================================================
def get_repair_AA(all_id_file,gene_file,new_ref,new_cds,outPath,seq_type):
    # 1. build {gene symbol:gene id dictionary}, {protein access: transcript id}
    id_df = pd.read_csv(all_id_file,sep='\t',header=0)
    g_symbol_id_dic = id_df.set_index('GeneSymbol')['GeneID'].to_dict()
    g_id_symbol_dic = id_df.set_index('GeneID')['GeneSymbol'].to_dict()
    prac_trid_dic = id_df.set_index('PrAccess')['TrID'].to_dict()
    # 2. get all repair gene ids
    gene_df = pd.read_csv(gene_file,header=0,names=['name'])
    gene_symbols = gene_df['name'].tolist()
    gene_ids = [g_symbol_id_dic[name] for name in gene_symbols]
    # 3. get refDNA_dic
    columns = ['chr','start','end','geneid','access','strand']
    refDNA_dic = SeqIO.index(new_ref, "fasta")
    cds_df = pd.read_csv(new_cds,sep='\t',header=0,names=columns)
    # 3. get all protein sequences
    if not os.path.exists(outPath): os.mkdir(outPath)
    for g in gene_ids:
        AAs,prs = AA_sequence(refDNA_dic,cds_df,g,seq_type)
        for pr,aa in zip(prs,AAs):
            fn = ''.join([g_id_symbol_dic[g],'_',prac_trid_dic[pr],'.protein.fa'])
            handle = open(outPath+'/'+fn,'w')
            handle.write('>'+fn[:-11]+'_'+pr+'\n'+aa+'\n')
            handle.close()
# outPath = '/data/shangzhong/DNArepair/f01_Repair_AA'
# get_repair_AA(all_id_file,gene_file,old_ref,old_cds,outPath,'AA')
# 
# outPath = '/data/shangzhong/DNArepair/f02_Repair_nt'
# get_repair_AA(all_id_file,gene_file,old_ref,old_cds,outPath,'nt')

#===============================================================================
#                         4. find inconsistent AA with refseq sequence
#===============================================================================
# outPath = '/data/shangzhong/DNArepair/f03_inconsistent_seq'
# if not os.path.exists(outPath): os.mkdir(outPath)
# outputFile = outPath + '/hamster_pr_different_from_refseq.fa'
# get_diff_pr_from_refseq(outputFile,old_ref,old_gffFile,gene_file) # get inconsistant for old file
#===============================================================================
#                         5. get refseq cds seq, gff cds seq of inconsistent pr nt sequence
#===============================================================================
def get_cds_of_tr(tr):
    handle = Entrez.efetch(db = 'nucleotide',id=tr,rettype='gb',retmode='text')
    record = SeqIO.read(handle, "gb")
    record.features = [f for f in record.features if f.type == 'CDS']
    if len(record.features) > 1:
        print tr,'has fragmented cds'
    else:
        record.seq = record.features[0].extract(record.seq)
    return record
# inco_id_file = inconsistent_folder + '/hamster_pr_different_from_refseq.id.fa'
# #-------------- 1. build {prac:trac} dictionary ---------------------------
# id_df = pd.read_csv(old_all_id_file,sep='\t',header=0)
# pr_tr_dict = id_df.set_index('PrAccess')['TrAccess'].to_dict()
# pr_trid_dict = id_df.set_index('PrAccess')['TrID'].to_dict()
# #-------------- 2. output CDS of refseq tr---------------------------
# outFile = inconsistent_folder + '/hamster_pr_nt_diff_refseq.fa'
# df = pd.read_csv(inco_id_file,sep='\t',header=None,names=['g_sym','trid','prac'])
# pr_symbol_dict = df.set_index('prac')['g_sym'].to_dict()
# prs = df['prac'].tolist()
# out_handle = open(outFile,'w')
# for pr in prs:
#     tr = pr_tr_dict[pr]
#     record = get_cds_of_tr(tr)
#     record.id = pr_symbol_dict[pr]+'_'+pr_trid_dict[pr]+'_'+pr 
#     SeqIO.write(record,out_handle,'fasta')
# out_handle.close()
# #-------------- 3. output CDS of gff tr---------------------------
# outPath = '/data/shangzhong/DNArepair/f02_Repair_nt'
# os.chdir(outPath)
# all_pr_nt_fa_files = [f[:-11] for f in os.listdir(outPath) if f.endswith('protein.fa')]
# pr_files = []
# for pr in prs:
#     trid = pr_trid_dict[pr]
#     g_sym = pr_symbol_dict[pr]
#     if g_sym+'_'+trid in all_pr_nt_fa_files:
#         pr_files.append(g_sym+'_'+trid+'.protein.fa')
# outFile = inconsistent_folder + '/hamster_pr_nt_diff_gff.fa'
# cmd = 'cat {f} > {out}'.format(f=' '.join(pr_files),out=outFile)
# subprocess.call(cmd,shell=True)
    
#===============================================================================
#                     6. get DNA repair genes mRNA sequences
#===============================================================================
# # 1. get genes
# fn = gene_file
# df = pd.read_csv(fn,header=0,names=['gene'])
# hamster_genes = list(set(df['gene'].tolist()))
# 
# yeast_df = pd.read_csv(yeast_gene_file,header=0,names=['gene'])
# yeast_genes = list(set(yeast_df['gene'].tolist()))
# genes = hamster_genes + yeast_genes
# # 2. get sequence
# rna_fa = '/data/shangzhong/DNArepair/Database/hamster_rna.fa'
# handle = open(rna_fa,'r')
# outHandle = open('/data/shangzhong/DNArepair/Database/repair_rna.fa','w')
# for record in SeqIO.parse(handle,'fasta'):
#     for g in genes:
#         if g in record.description:
#             SeqIO.write(record,outHandle,'fasta')
#             break
# outHandle.close()
#===============================================================================
#                     7. compare gene count between modified genome and non modified genome 
#===============================================================================
# id_file = '/data/shangzhong/DNArepair/Database/hamster_AllIDS.txt'
# f1 = '/data/shangzhong/DE/test/before_htseq/1.txt'
# f2 = '/data/shangzhong/DE/test/before_htseq/4.txt'
# gene_file = '/data/shangzhong/DNArepair/Database/01_DNA_repair_genes.hamster.txt'
# gene_df = pd.read_csv(gene_file,header=0,names=['gene'])
# genes = gene_df['gene'].tolist()
#     
# df1 = pd.read_csv(f1,sep='\t',header=None,names=['id','1'],index_col=0)
# df1 = df1[:-5]
# df2 = pd.read_csv(f2,sep='\t',header=None,names=['id','2'],index_col=0)
# df2 = df2[:-5]
# df = pd.concat([df1,df2],axis=1)
# res_df = df[df['1']!= df['2']]
#     
# id_df = pd.read_csv(id_file,sep='\t',usecols=[0,1],names=['id','symbol'])
# id_df = id_df.drop_duplicates()
# dic = id_df.set_index('id')['symbol'].to_dict()
# res_df['symbol'] = res_df.index.map(lambda x: dic[str(x)] if x in dic.keys() else x)
#     
# res_df = res_df[res_df['symbol'].isin(genes)]
# res_df.to_csv('/data/shangzhong/DE/test/cmp.txt',sep='\t')
# print(res_df)
#===============================================================================
#                     8. get new gff removing DNArepair genes and add new annotations
# (the results are same for removing and not removing previous DNArepair genes annotations)
#===============================================================================
# gff = '/data/genome/hamster/ncbi_refseq_masked_repair_genes/hamster.masked_discrepancies.repair_refseq.gff'
# gene_file = '/data/shangzhong/DNArepair/Database/01_DNA_repair_genes.hamster.txt'
# gene_df = pd.read_csv(gene_file,header=0,names=['gene'])
# genes = gene_df['gene'].tolist()
# 
# new_gff = '/data/shangzhong/DE/test/new.gff'
# out = open(new_gff,'w')
# 
# handle = open(gff)
# for line in handle:
#     remove=False
#     for g in genes:
#         if 'gene='+g+';' in line:
#             remove = True
#             break
#     if remove==False:
#         out.write(line)
# handle.close()
# out.close()
#===============================================================================
#                     9. get cds start and end position given rna accession
#===============================================================================
def get_rna_cds_start_end(rna_access):
    handle = Entrez.efetch(db = 'nucleotide',id=rna_access,rettype='gb',retmode='text')
    record = handle.read()
    lines = record.split('\n')
    for line in lines:
        item = line.split()
        if 'CDS' in item:
            index = item.index('CDS')
            if '..' in item[index+1]:
                pos = item[index+1].split('..')
                start = pos[0]
                end = pos[1]
    try:
        return [start,end]
    except:
        return 'none'

# # get rna accessions
# refseq_rna = '/data/shangzhong/DNArepair/Database/hamster_rna.fa'
# rna = []
# for record in SeqIO.parse(open(refseq_rna,'r'),'fasta'):
#     rna.append(record.id)
# rna = list(set(rna))
# # get start and end postions
# outFile = '/data/shangzhong/DNArepair/f04_new_annotation/rna_cds_pos.txt'
# outHandle = open(outFile,'w')
# for r in rna:
#     pos = get_rna_cds_start_end(r)
#     if pos != 'none':
#         outHandle.write('\t'.join([r,pos[0],pos[1]])+'\n')
# outHandle.close()
#===============================================================================
#                     10. Merged all repair mRNA into one scaffold, test whether
# the CDS sequences are the same with refseq AA sequences.
#===============================================================================
# #-------------- 1. build {trac:prac} dictionary ---------------------------
# id_df = pd.read_csv(old_all_id_file,sep='\t',header=0)
# tr_pr_dict = id_df.set_index('TrAccess')['PrAccess'].to_dict()
# #-------------- 2. build {trac:gene} dictionary ---------------------------
# tr_gene_file = '/data/hooman/DNARepair/03_DNA_repair_genes.hamster.customscaffold.accession.txt'
# tr_gene_df = pd.read_csv(tr_gene_file,sep='\t',header=None,names=['gene','trac'])
# tr_gene_dict = tr_gene_df.set_index('trac')['gene'].to_dict()
# trs = tr_gene_df['trac'].tolist()
# #--------------- 3. get AA sequence ----------------------------------------
# new_ref_for_rna = '/data/genome/hamster/ncbi_refseq_masked_repair_genes/hamster.masked_discrepancies.repair_refseq.fa'
# new_gff_for_rna = '/data/genome/hamster/ncbi_refseq_masked_repair_genes/hamster.masked_discrepancies.repair_refseq.gff'
# new_gff_for_rna_df = pd.read_csv(new_gff_for_rna,sep='\t',header=None,usecols=[0,2,3,4,8],names=['chr','feature','start','end','anno'],comment='#')
# new_gff_for_rna_df = new_gff_for_rna_df[(new_gff_for_rna_df['chr'].values=='repair_gene_scaff') & (new_gff_for_rna_df['feature'].values=='CDS')]
# new_ref_dict = SeqIO.index(new_ref_for_rna,'fasta')
# repair_seq = str(new_ref_dict['repair_gene_scaff'].seq)

# for tr in trs:
#     # get ref_seq AA sequence
#     tr = 'XM_007625949.1'
#     pr = tr_pr_dict[tr]
#     handle = Entrez.efetch(db = 'protein',id=pr,rettype='fasta',retmode='text')
#     record = handle.read()
#     sequence = ''.join(record.split('\n')[1:-2])
#     # get gff AA sequence
#     g = tr_gene_dict[tr]
#     cri = new_gff_for_rna_df['anno'].map(lambda x: 'gene=' + g in x)
#     g_df = new_gff_for_rna_df[cri]
#     seq = repair_seq[int(g_df['start'].values[0])-1:int(g_df['end'].values[0])]
#     AA = str(Seq(seq,generic_dna).translate())
#     AA = ''.join([f  if f !='*' else 'X' for f in AA])
#     if sequence != AA[:-1]:
#         print g,'sequence is not consistent'
