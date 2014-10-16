"""
This file run the pipeline of doing differential expression analysis of 
Contaminated CHO samples.
"""
import os
import sys 
import subprocess
sys.path.append('/home/shangzhong/Codes/Pipeline')
from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.f02_aligner_command import gsnap,STAR
from Modules.f03_samtools import sam2bam_sort
from Modules.f04_htseq import htseq_count
from Modules.f05_GeneIDConvert import ID_Convert
from Modules.FileProcess import remove,get_parameters
#=========== define some parameters ===================
#parFile = '/home/shangzhong/Codes/Pipeline/DiffExpression/DE_Parameters.txt'
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
trimmomatic = param['trimmomatic']
trimmoAdapter = param['trimmoAdapter']

aligner = param['aligner']
annotation = param['annotation']
output_path = param['htseqOutPath']
db_name = param['gsnapDbName']
gsnap_annotation = param['gsnapAnnotation']

Dict = param['symbolIDFile']
inputpath = file_path

#=========== (0) enter the directory ================
Message(startMessage,email)
os.chdir(file_path)
#=========== (1) reads files and trim ===============
fastqFiles = list_files(file_path)
if trim == 'True':
    fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred,trimmoAdapter)
print 'list file succeed'
#=========== (2) run gsnap to do the mapping ========
if aligner == 'gsnap':
    map_files = gsnap(fastqFiles,db_path, db_name,gsnap_annotation,thread)
else:
    map_files = STAR(fastqFiles,db_path,thread)
print 'align succeed'
#=========== (3) samtools to sort the file ==========
sorted_bam = sam2bam_sort(map_files,thread)
print 'sorted succeed'
#=========== (4) htseq_count ========================
htseq_count(sorted_bam,annotation,file_path)
print 'htseq count succeed'
#=========== (5) htseq symbol to id =================
ID_Convert(Dict,output_path,inputpath)
print 'id convert succeed'
Message(endMessage,email)