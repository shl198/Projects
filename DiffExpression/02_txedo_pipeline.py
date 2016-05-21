"""
This file run the pipeline of doing differential expression analysis of 
Contaminated CHO samples.
"""
import os
import sys
import subprocess
import pandas as pd
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.f02_aligner_command import gsnap,STAR,STAR_Db,bowtie
from Modules.f03_samtools import sam2bam_sort
from Modules.f05_IDConvert import addEntrezID2CufflinkResultWithNCBIAnnotation,addEntrezGeneID2CufflinkResultWithEnsemblAnnotation
from Modules.p01_FileProcess import get_parameters
from Modules.f06_txedo import cufflinks
#=========== define some parameters ===================
#parFile = '/data/shangzhong/DE/fpkm/02_txedo_Parameters.txt'
parFile = sys.argv[1]
param = get_parameters(parFile)

thread = param['thread']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']
ref_fa = param['refSequence']
file_path = param['filePath']
db_path = param['alignerDb']
trim = param['trim']
phred = param['phred']
picard = param['picard']
trimmomatic = param['trimmomatic']
trimmoAdapter = param['trimmoAdapter']

aligner = param['aligner']
annotation = param['annotation']
db_name = param['gsnapDbName']
gsnap_annotation = param['gsnapAnnotation']

Dict = param['symbolIDFile']
source = param['IDresource']
inputpath = file_path

#=========== 0. enter the path =================
os.chdir(file_path)
Message(startMessage,email)
"""
#========  (1) read files and trim ================================
fastqFiles = list_files(file_path)
print 'list file succeeded'
if trim == 'True':
    try:
        fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred,trimmoAdapter,batch=6)
        print 'trim succeeded'
        print 'fastqFiles is: ',fastqFiles
    except:
        print 'trim failed'
        raise
#========  (2) align reads ================================
if aligner == 'STAR':
    if not os.path.exists(db_path):
        os.mkdir(db_path)
    if os.listdir(db_path) == []:
        STAR_Db(db_path,ref_fa,thread)
    
    try:
        bam_files =  STAR(fastqFiles,db_path,thread,annotation,['--outSAMstrandField intronMotif',
                                                    '--outFilterType BySJout',
                                                    '--outFilterIntronMotifs RemoveNoncanonicalUnannotated',
                                                    '--outSAMtype BAM SortedByCoordinate',
                                                    '--limitBAMsortRAM 1417527247'])
                                                    #'--quantMode GeneCounts'])
        print 'align succeeded'
        print 'sam_files is: ',bam_files
    except:
        print 'align failed'
        raise
elif aligner =='bowtie':
    try:
        sam_files = bowtie(fastqFiles,db_path,thread=1,otherParameters=[''])
        print 'align succeed'
        print 'sam files are:',sam_files
    except:
        print 'align failed'
        raise"""
# #======== (3) sort sam files ==============================
# STAR already has the option to output sorted bam files, this part is not necessary
# try:
#     sorted_bams = sam2bam_sort(sam_files,thread)
#     print 'sort bam succeeded'
#     print 'sorted_bams is: ',sorted_bams
# except:
#     print 'sort sam failed'
#     raise
#======== (4) cufflinks ===================================
bam_files = ['trim_CHO1-35695882_1.bam', 'trim_CHO2-35697883_1.bam', 'trim_CHO3-35723733_1.bam', 'trim_CHO4-35723735_1.bam', 'trim_CHO5-35716755_1.bam', 'trim_CHO6-35724722_1.bam', 'trim_CHO7-35717764_1.bam', 'trim_CHO8-35712778_1.bam', 'trim_CHO9-35698856_1.bam', 'trim_CHO10-35723738_1.bam', 'trim_CHO11-35713776_1.bam', 'trim_CHO12-35711756_1.bam', 'trim_CHO13-35710779_1.bam', 'trim_CHO14-35722751_1.bam', 'trim_CHO15-35710781_1.bam', 'trim_CHO16-35703833_1.bam', 'trim_CHO18-35714763_1.bam', 'trim_CHO20-35727704_1.bam', 'trim_CHO21-35715768_1.bam', 'trim_CHO22-35721748_1.bam', 'trim_CHO24-35721745_1.bam', 'trim_CHO25-35726711_1.bam', 'trim_CHO26-35728707_1.bam', 'trim_CHO28-35714778_1.bam', 'trim_CHO29-35703840_1.bam', 'trim_CHO30-35729705_1.bam', 'trim_CHO31-35712790_1.bam', 'trim_CHO32-35726718_1.bam', 'trim_CHO33-35727712_1.bam', 'trim_CHO34-35730696_1.bam', 'trim_CHO35-35698875_1.bam', 'trim_CHO36-35724745_1.bam', 'trim_CHO37-35741706_1.bam', 'trim_CHO38-35742709_1.bam', 'trim_CHO39-35740708_1.bam', 'trim_CHO40-35737706_1.bam', 'trim_CHO42-35737702_1.bam', 'trim_CHO43-35734705_1.bam', 'trim_CHO44-35736705_1.bam', 'trim_CHO45-35733703_1.bam', 'trim_CHO46-35740715_1.bam', 'trim_CHO47-35731702_1.bam', 'trim_CHO48-35745711_1.bam']
try:
    folders = cufflinks(bam_files,annotation,thread)
    print 'cufflink succeeded'
    print 'folers is: ',folders
except:
    print 'cufflink failed'
#======== (5) id convertion ===============================
if source == 'ncbi':
    for f in folders:
        geneFpkmFile = f + '/genes.fpkm_tracking'
        new_res_file = addEntrezID2CufflinkResultWithNCBIAnnotation(Dict,geneFpkmFile)
        df = pd.read_csv(new_res_file,sep='\t',header=0)
        df = df[['Entrez_ID','FPKM']].groupby(['Entrez_ID']).sum()
        df.to_csv(new_res_file,sep='\t')
        print 'new file are: ', new_res_file
    print 'id conversion succeeded'
    
elif source == 'ensembl':
    for f in folders:
        geneFpkmFile = f + '/genes.fpkm_tracking'
        new_res_file = addEntrezGeneID2CufflinkResultWithEnsemblAnnotation(Dict,geneFpkmFile)
        df = pd.read_csv(new_res_file,sep='\t',header=0)
        df = df[['EntrezID','FPKM']].groupby(['EntrezID']).sum()
        df.to_csv(new_res_file,sep='\t')
        print 'new file are: ', new_res_file
    print 'id conversion succeeded'
else:
    pass
    
Message(endMessage,email)
