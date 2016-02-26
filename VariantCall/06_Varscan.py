"""
this file does variant calling using Varscan
"""
#=============  import required packages  =================
import os
import sys
import subprocess
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f00_Message import Message
from Modules.p01_FileProcess import remove,get_parameters,rg_bams
#=============  define some parameters  ===================
"""these parameters and read group names are different for 
   different samples, should only change this part for 
   running pipeline
"""
#parFile = '/data/shangzhong/DNArepair/RNA/GATK_parameters4DNAandRNA.txt'
parFile = sys.argv[1]
param = get_parameters(parFile)
thread = param['thread']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']

ref_fa = param['refSequence']
file_path = param['filePath']
starDb = param['alignerDb']
trim = param['trim']
phred = param['phred']

picard = param['picard']
trimmomatic = param['trimmomatic']
trimmoAdapter = param['trimmoAdapter']
gatk = param['gatk']
read_group = param['readGroup']
organism = param['organism']

##*****************  Part 0. Build index file for bwa and GATK ******
##*****************  Part I. Preprocess  ============================
#========  1. map and dedupping =====================================
#========  (0) enter the directory ========================
os.chdir(file_path)
Message(startMessage,email)
#========  (1) read files  ================================
fastqFiles = list_files(file_path)
if trim == 'True':
    trim_fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred,trimmoAdapter,batch=6)
    remove(fastqFiles)
else:
    trim_fastqFiles = fastqFiles
sys.stdout.write('list file succeed\n')
sys.stdout.write('fastqFiles is: {fq}\n'.format(fq=trim_fastqFiles))