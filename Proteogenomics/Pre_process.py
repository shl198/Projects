import pandas as pd
import os
import sarge
from Bio import SeqIO,Seq
from natsort import natsorted
import shutil
"""This file prepare the bam file for spliceDB"""

def chunk(l,n):
    n = max(1,n)
    result = [l[i:i+n] for i in range(0,len(l),n)]
    return result
#===============================================================================
#                     1. find the longest scaffold and the one that has largest number of genes for old hamster genome
#===============================================================================
#-------- 1. get gene number ------------------
# allID = '/data/shangzhong/Database/hamster/hamster_AllIDS.txt'
# id_df = pd.read_csv(allID,sep='\t',header=0)
# chr_gene_dic = {k:list(v) for k,v in id_df.groupby('Chrom')['GeneID']}
# for key in chr_gene_dic:
#     chr_gene_dic[key] = len(chr_gene_dic[key])
# 
# ref_fa = '/data/genome/hamster/ncbi_refseq/hamster.fa'
# len_chr = {}
# for record in SeqIO.parse(open(ref_fa,'r'),'fasta'):
#     len_chr[record.id] = len(record.seq)
# 
# df = pd.DataFrame({'chr':len_chr.keys()})
# df['chr_len'] = df['chr'].map(lambda x: len_chr[x])
# df['gene_num'] = df['chr'].map(lambda x: chr_gene_dic[x] if x in chr_gene_dic.keys() else 0)
# 
# max_chr_len = max(len_chr.values())
# print df[df['chr_len'].values==max_chr_len]
# max_gene_num = max(df['gene_num'].tolist())
# print df[df['gene_num'].values==max_gene_num]
# df.to_csv('/data/shangzhong/Proteogenomics/Results/01_chr_len_gene.txt',sep='\t',index=False)
#===============================================================================
#                 2. extract alignment based on chromosome
#===============================================================================
def copy_files(path,target_path,folders):
    for folder in folders:
        fq_path = path + '/' + folder
        os.chdir(fq_path)
        fqFiles = [f for f in os.listdir(fq_path) if f.endswith('.fastq.gz')]
        fqFiles = natsorted(fqFiles)
        fst = [f for f in fqFiles if 'R1' in f]
        snd = [f for f in fqFiles if 'R2' in f]
        
        cmd = ('cat {input} > {tar}/{folder}_1.fq.gz').format(
                    input=' '.join(fst),folder=folder,tar=target_path)
        print(cmd)
        sarge.run(cmd)
        
        cmd = ('cat {input} > {tar}/{folder}_2.fq.gz').format(
                    input=' '.join(snd),folder=folder,tar=target_path)
        print(cmd)
        sarge.run(cmd)

def rnaseq_map_and_extract_by_chr(path,target_path,batch,rna_pipeline_file,rna_pipeline_param,chrom=''):
    os.chdir(path)
    folders = [f for f in os.listdir(path) if os.path.isdir(f)]
    folders = natsorted(folders)
    sub_folders = chunk(folders,batch)
    
    for sub_dir in sub_folders:
        # 1. copy files
        copy_files(path,target_path,sub_dir)
        # 2. map using STAR
        os.chdir(target_path)
        cmd = ('python {pipe} {pipe_param}').format(pipe=rna_pipeline_file,pipe_param=rna_pipeline_param)
        sarge.run(cmd)
        # 3. extract chromosome
        if chrom != '':
            bam_path = target_path + '/sortBam'
            os.chdir(bam_path)
            if not os.path.exists(bam_path+'/chr'): os.mkdir(bam_path+'/chr')
            bams = [f for f in os.listdir(bam_path) if f.endswith('.sort.bam')]
            for bam in bams:
                out = bam[5:]
                cmd = ('samtools view {input} {chr} > chr/{out}').format(input=bam,chr=chrom,out=out)
                print(cmd)
                sarge.run(cmd)
                os.remove(bam)
                os.remove(bam+'.bai')
            shutil.rmtree(target_path+'/bam')
        # 4. clear files
#         fqs = [os.path.join(target_path,f) for f in os.listdir(target_path) if f.endswith('.fq.gz')]
#         for fq in fqs: os.remove(fq)
            
#path = '/media/lewislab/Dropbox (UCSD SBRG)/LewisPub/Data/DataForProteogenomics/DBH_CHO_mass_spec_data/YS_transfered_fastq_files'
path = '/media/lewislab/Dropbox (UCSD SBRG)/LewisPub/Data/RNAseq/2016_05_HH_L/CHO_natelewis-30316287'
target_path = '/data/shangzhong/DE/ercc'
chrom = ''#'NW_006870216.1'
batch= 100
rna_pipeline_file = '/home/shangzhong/Codes/NewPipeline/RNAseq_count.py'
rna_pipeline_param = '/data/shangzhong/DE/ercc/RNAseq_count.yaml'
rnaseq_map_and_extract_by_chr(path,target_path,batch,rna_pipeline_file,rna_pipeline_param,chrom='')
    
