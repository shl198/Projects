import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rosetta.parallel.pandas_easy import groupby_to_series_to_frame
from joblib import Parallel,delayed
import multiprocessing
#===============================================================================
#                     1. Merge scaffolds
#===============================================================================
fa = '/media/lewislab/Dropbox (UCSD SBRG)/LewisPub/Data/RefGenomes/hamster_2013/2015_PBJ_temp_assembly/ch_illumina_pbj.fasta'
outFile = '/data/shangzhong/LargeDeletion/hamster_pb_illu_super_scaffold.fa'
chr_junc_file = '/data/shangzhong/LargeDeletion/hamster_pb_super_scaffold_junc.txt'
"""
in_handle = open(fa,'r')
out_handle = open(outFile,'w')
chr_junc_handle = open(chr_junc_file,'w')
sequence = ''
n = 1
pos = 0
for record in SeqIO.parse(in_handle,'fasta'):
    sequence += str(record.seq)
    pos += len(record.seq)
    if len(sequence) >= 200000000:
        item = SeqRecord(Seq(sequence), id = 'chr'+str(n),description="")
        SeqIO.write(item,out_handle,'fasta')
        sequence = ''
        n += 1
        pos = 0
    else:
        sequence += 'N'*500
        chr_junc_handle.write('\t'.join(['chr'+str(n),str(pos),str(pos+500)]) + '\n')
        pos += 500
chr_junc_handle.close()
if sequence != '':
    item = SeqRecord(Seq(sequence[:-500]), id = 'chr'+str(n),description="")
    SeqIO.write(item,out_handle,'fasta')
    # remove last line of junctions
    df = pd.read_csv(chr_junc_file,sep='\t',header=None,names=['chr','start','end'])
    df = df[:-1]
    df.to_csv(chr_junc_file,sep='\t',index=False)
out_handle.close()
# check the length of each chromosome
handle = open(outFile)
for record in SeqIO.parse(handle,'fasta'):
    print len(record.seq)
# check the inserted 500 Ns
df = pd.read_csv(chr_junc_file,sep='\t',header=0)
for record in SeqIO.parse(open(outFile),'fasta'):
    chr_name = record.id
    chr_seq = str(record.seq)
    chr_df = df[df['chr'].values==chr_name]
    for start,end in zip(chr_df['start'],chr_df['end']):
        insert = chr_seq[int(start):int(end)]
        if insert != 'N'*500:
            raise ValueError('insert seqeunce pos extraction wrong')
"""
#===============================================================================
#                 2. filter largedel results (make sure the deletion doesn't span over scaffolds)
#===============================================================================
def del_overlap(chr_name,start,end,junc_df):
    """start,end are deletion start end"""
    junc_df = junc_df[junc_df['chr'].values=='chr' + str(chr_name)]
    filter_junc = junc_df[(junc_df['start']>=int(start)) & (junc_df['end']<=int(end))]
    if filter_junc.empty:
        return True
    else:
        return False
# large_del = '/data/shangzhong/LargeDeletion/01_largedel_res.txt'
# filter_large_del = '/data/shangzhong/LargeDeletion/02_filter_large_del.txt'
# 
# junc_df = pd.read_csv(chr_junc_file,sep='\t',header=0,names=['chr','start','end'])
# lar_del_df = pd.read_csv(large_del,sep='\t',header=0,names=['chr','start','end','cov'])
# lar_del_df['chr'] = lar_del_df['chr'].astype('str')
# print(lar_del_df.shape)
# #del_overlap('1','783881','785184',junc_df)
# cri = lar_del_df.apply(lambda x: del_overlap(x['chr'],x['start'],x['end'],junc_df),axis=1)
# res_df = lar_del_df[cri]
# print(res_df.shape)
# res_df.to_csv(filter_large_del,sep='\t',index=False)
#===============================================================================
#                 3. get the overlap between the deletion and gene 
#===============================================================================
def keep_map(single_df):
    """decide to keep the blast hits or not"""
    index = single_df.index
    if len(index) == 1:
        return single_df
    for i in range(len(index))[:-1]:
        if single_df.loc[index[i],'keep'] == 'F':
            continue
        max_index = max(single_df.loc[index[i],'sstart'],single_df.loc[index[i],'send'])
        min_index = min(single_df.loc[index[i],'sstart'],single_df.loc[index[i],'send'])
        fst_range = range(min_index,max_index+1)
        for j in range(i+1,len(index)):
            max_in = max(single_df.loc[index[j],'sstart'],single_df.loc[index[j],'send'])
            min_in = min(single_df.loc[index[j],'sstart'],single_df.loc[index[j],'send'])
            snd_range = range(min_in,max_in+1)
            if set(fst_range).intersection(snd_range) != set():
                single_df.loc[index[j],'keep'] = 'F'
    return single_df
#----------------- 1) filter blast results -------------------------------------
def applyParallel(dfGrouped,func,thread):
    retLst = Parallel(thread)(delayed(func)(group) for name,group in dfGrouped)
    return pd.concat(retLst)


# del_threshold = 500
# blast = '/data/shangzhong/LargeDeletion/03_cho_blast_low_sensitivity28_filter.txt'
# blst_df = pd.read_csv(blast,sep='\t',header=None,names=['qid','chr','iden','len',
#                     'mis','gap','qstart','qend','sstart','send','evalue','bit'])
# blst_df['keep'] = pd.Series(['T']*blst_df.shape[0])
# # remove overlapped regions
# res = applyParallel(blst_df.groupby(['qid','chr']),keep_map,10)
# result = res.sort_index()
# result = result[result['keep'].values=='T']
# result = result.drop('keep',1)
# result.to_csv('/data/shangzhong/del_blast2_hamster.txt',sep='\t',header=0,index=False)

#===============================================================================
#                     4. get overlap between deletion and genes
#===============================================================================








