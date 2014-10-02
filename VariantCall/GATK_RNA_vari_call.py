"""
this file does variant calling for RNAseq
"""
#=============  import required packages  =================
import os
import sys
sys.path.append('/home/shangzhong/Codes/Pipeline')
from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.f02_aligner_command import STAR2Pass
from Modules.f03_samtools import sam2bam_sort
from Modules.f07_picard import markduplicates,addReadGroup
from Modules.f08_GATK import *
from Modules.FileProcess import remove,get_parameters,rg_bams
#=============  define some parameters  ===================
"""these parameters and read group names are different for 
   different samples, should only change this part for 
   running pipeline
"""

parFile = '/home/shangzhong/Codes/Pipeline/VariantCall/RNA_Parameters.txt'
# parFile = sys.argv[1]
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
gatk = param['gatk']
read_group = param['readGroup']
organism = param['organism']

##*****************  Part 0. Build index file for bwa and GATK ******
##*****************  Part I. Preprocess  ============================
#========  1. map and dedupping =====================================
#========  (0) enter the directory ========================
os.chdir(file_path)
#Message(startMessage,email)
#========  (1) read files  ================================
fastqFiles = list_files(file_path)
if trim == 'True':
    fastqFiles = Trimmomatic(fastqFiles,phred)
print 'list file succeed'
#========  (2) align using 2 pass STAR ====================
map_sams= STAR2Pass(fastqFiles,starDb,ref_fa,thread)

#========  2. Add read groups, sort,mark duplicates, and create index
#========  (1) sort and add group =========================
sort_bams = addReadGroup(picard,map_sams,read_group)
#========  (2) mark duplicates ============================
dedup_bams = markduplicates(picard,sort_bams)

#========  3. Split 'N' Trim and reassign mapping qualiteies
