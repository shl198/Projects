import pandas as pd
import os
import ipdb
import numpy
from Bio import SeqIO
import subprocess
from Modules.f00_Message import Message
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib as mpl
from multiprocessing import Process
mpl.style.use('ggplot')
#===============================================================================
#  read and merge fastq files: usually applied to CHO data from Denmark
#===============================================================================
# import os,subprocess
# path = '/media/lewislab/Dropbox (UCSD SBRG)/LewisPub/Data/RNAseq/2014_11_CHOS_baseline/B.Voldborg_14_05'
# os.chdir(path)
# folders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
# folders = natsorted(folders)
# target_path = '/data/shangzhong/DE/chos_base'
# folders = folders[20:40]
# for f in folders:
# 	subfolder = os.listdir(f)
# 	files = os.listdir(f+'/'+subfolder[0])
# 	fir = []
# 	snd = []
# 	for filename in files:
# 		if filename.endswith('_1.fastq.gz'):
# 			fir.append(f + '/' + subfolder[0] + '/' + filename)
# 		if filename.endswith('_2.fastq.gz'):
# 			snd.append(f + '/' + subfolder[0] + '/' + filename)
# 	cmd = 'cat ' + ' '.join(fir) + ' > ' + target_path + '/' + f + '_1.fq.gz'
# 	subprocess.call(cmd,shell=True)
# 	cmd = 'cat ' + ' '.join(snd) + ' > ' + target_path + '/' + f + '_2.fq.gz'
# 	subprocess.call(cmd,shell=True)


#===============================================================================
#                 merge cufflinks
#===============================================================================

# filepath = '/data/shangzhong/DE/pgsa/cufflinks'
# 
# os.chdir(filepath)
# folders = [f for f in os.listdir(filepath) if os.path.isdir(os.path.join(filepath,f))]
# folders = natsorted(folders)
# ConvertFile = '/data/shangzhong/Database/gff_chok1_ID_symbol_accession.txt'
# # start merge
# data = pd.read_csv(folders[0] +'/genes.fpkm_tracking',sep='\t',header=0,usecols=[0,4,6,9],names=['tracking_id','gene_short_name','locus',folders[0]])
# data['index'] = data['tracking_id'].map(str) + data['gene_short_name'] + data['locus']
# data = data.set_index(['index'])
# for i in range(1,len(folders)):
#     df = pd.read_csv(folders[i] +'/genes.fpkm_tracking',sep='\t',header=0,usecols=[0,4,6,9],names=['tracking_id','gene_short_name','locus',folders[i]])
#     df['index'] = df['tracking_id'].map(str) + df['gene_short_name'] + df['locus']
#     df = df.set_index(['index'])
#     df = df[folders[i]]
#     data = pd.concat([data,df],axis=1)
# 
# data.to_csv(filepath + '/sample_FPKM.csv',sep='\t',index=False)

#===============================================================================
#               Read winzeler fastq files
#===============================================================================
# folders = [f for f in os.listdir(filepath) if os.path.isdir(os.path.join(filepath,f))]
# folders = natsorted(folders)
# n = 1
# for folder in folders:
#     os.chdir(os.path.join(filepath,folder))
#     files = [f for f in os.listdir(filepath+'/'+folder) if f.endswith('fastq.gz')]
#     
#     cmd = ('cat {f} > {out}').format(f=' '.join(files),out='/data/shangzhong/DE/Winzeler/sample0'+str(n)+'.fq.gz')
#     subprocess.call(cmd,shell=True)
#     print cmd
#     n = n +1


#===============================================================================
#      merge cufflink and htseq for pgsa
#===============================================================================
# from Modules.f05_IDConvert import addProduct2CufflinkResultWithNCBIAnnotation
# htseq = '/data/shangzhong/DE/pgsa/htseq/htseq.name.csv'
# cufflink = '/data/shangzhong/DE/pgsa/cufflinks/sample_FPKM.csv.txt.csv'
# 
# # remove duplicate gene id + gene symbol in cufflinks results
# cuff_df = pd.read_csv(cufflink,sep='\t',header=0)
# cuff_sort = cuff_df.sort(['gene_short_name','Entrez_GeneID'])
# cuff_rmdup = cuff_sort.drop_duplicates(['gene_short_name','Entrez_GeneID'])
# cuff_res = cuff_sort.groupby(['gene_short_name','Entrez_GeneID'],as_index=False).sum()
# cuff_res.to_csv('/data/shangzhong/inter.csv',index=False,sep='\t')
# cuff_rmdupFile = addProduct2CufflinkResultWithNCBIAnnotation('/data/shangzhong/Database/141028chok1gene_info.txt','/data/shangzhong/inter.csv')
# cuffres_df = pd.read_csv(cuff_rmdupFile,sep='\t',header=0)
# # merge cufflinks and htseq
# htseq = pd.read_csv(htseq,header=0,sep='\t')
# htseq['GeneID'] = htseq['GeneID'].astype(str)
# htseq = htseq.rename(columns={'GeneID':'Entrez_GeneID','GeneSymbol':'gene_short_name'})
# res = pd.merge(cuffres_df,htseq,how='outer',on=['gene_short_name','Entrez_GeneID'])
# res.to_csv('/data/shangzhong/pgsa.csv',index=False,sep='\t')


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
my_dna = Seq("AGTACACTGGT", generic_dna)
print my_dna.complement()
print my_dna.translate()



