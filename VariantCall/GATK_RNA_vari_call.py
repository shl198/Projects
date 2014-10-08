"""
this file does variant calling for RNAseq
"""
#=============  import required packages  =================
import os
import sys
import subprocess
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

#parFile = '/home/shangzhong/Codes/Pipeline/VariantCall/RNA_Parameters.txt'
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
    fastqFiles = Trimmomatic(fastqFiles,phred)
subprocess.call('echo \"list file succeed\"',shell=True)
#========  (2) align using 2 pass STAR ====================
map_sams= STAR2Pass(fastqFiles,starDb,ref_fa,thread)
subprocess.call('echo \"align succeed\"',shell=True)

#========  2. Add read groups, sort,mark duplicates, and create index
#========  (1) sort and add group =========================
sort_bams = sam2bam_sort(map_sams,thread)
subprocess.call('echo \"sort bam succeed\"',shell=True)

group_bams = addReadGroup(picard,sort_bams,read_group)
subprocess.call('echo \"add group succeed\"',shell=True)
remove(map_sams)
#========  (2) mark duplicates ============================
dedup_bams = markduplicates(picard,group_bams)
subprocess.call('echo \"mark duplicate succeed\"',shell=True)
remove(sort_bams)

#========  3. Split 'N' Trim and reassign mapping qualiteies
split_bams = splitN(gatk,dedup_bams,ref_fa)
subprocess.call('echo \"split N succeed\"',shell=True)

#========  4. Indel realignment ===========================
#========  (1) generate intervals =========================
interval = RealignerTargetCreator(gatk,split_bams,ref_fa,thread)
subprocess.call('echo \"RealignerTarget Creator succeed\"',shell=True)
#========  (2) realignment of target intervals ============
realign_bams = IndelRealigner(gatk,split_bams,ref_fa,interval)
subprocess.call('echo \"IndelRealigner succeed\"',shell=True)
remove(split_bams)

#========  5. Base quality recalibration  =================

# since we don't have dbsnp for CHO, we need to:
# 1. find snp without recalibration, got vcf file
# 2. extract the snps we think are real snps, into a real_vcf file.
# 3. use the file in 2 to do the recalibration.

roundNum = 1
#========  6. Variant Calling =============================
vcf_files = HaplotypeCaller_RNA_VCF(gatk,realign_bams,ref_fa)
subprocess.call('echo \"1 round call succeed\"',shell=True)

#========  7. Variant filtering ===========================
gold_varis = RNA_Vari_Filter(gatk,vcf_files,ref_fa)
subprocess.call('echo \"1 round gold variant succeed\"',shell=True)

#========  5. Base quality recalibration  =================
recal_bams = RNA_BaseRecalibrator(gatk,realign_bams,ref_fa,
            gold_varis,roundNum,thread)
subprocess.call('echo \"recalibration succeed\"',shell=True)

#========  !!! merge lanes for the same sample ============
roundNum = '1'
if len(recal_bams) !=1:
    merged_bams = rg_bams(read_group,recal_bams)
    remove(recal_bams)
    dedup_files = markduplicates(picard,merged_bams)
    remove(merged_bams)
    interval = RealignerTargetCreator(gatk,dedup_files,ref_fa,thread)
    realign_bams = IndelRealigner(gatk,dedup_files,ref_fa,interval)
    remove(dedup_files)
subprocess.call('echo \"merge lanes succeed\"',shell=True)

roundNum = 2
#========  6. Variant Calling =============================
vcf_files = HaplotypeCaller_RNA_VCF(gatk,realign_bams,ref_fa)
subprocess.call('echo \"2 round call succeed\"',shell=True)
#========  7. Variant filtering ===========================
gold_varis = RNA_Vari_Filter(gatk,vcf_files,ref_fa)
subprocess.call('echo \"2 round gold variant succeed\"',shell=True)
subprocess.call('echo \" job finished succeessfully\"',shell=True)