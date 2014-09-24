"""
this file does variant calling for DNAseq
"""
#=============  import required packages  =================
import os
import sys
sys.path.append('/home/shangzhong/Codes/Pipeline')
from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.f02_aligner_command import bwa_vari
from Modules.f03_samtools import sam2bam_sort
from Modules.f07_picard import markduplicates
from Modules.f08_GATK import *
from Modules.FileProcess import remove,get_parameters,rg_bams
#=============  define some parameters  ===================
"""these parameters and read group names are different for 
   different samples, should only change this part for 
   running pipeline
"""

parFile = '/home/shangzhong/Codes/Pipeline/VariantCall/Parameters.txt'
#parFile = sys.argv[2]
param = get_parameters(parFile)
thread = param['thread']
email = param['email']

ref_fa = param['refSequence']
file_path = param['filePath']
bwaDb = param['alignerDb']
trim = param['trim']
phred = param['phred']

read_group = param['readGroup']
organism = param['organism']
##*****************  Part 0. Build index file for bwa and GATK ******
##=================  Part I. Preprocess  ============================
#========  1. map and dedupping =====================================
#========  (0) enter the directory ========================
os.chdir(file_path)
Message(param['startMessage'],email)
#========  (1) read files  ================================
fastqFiles = list_files(file_path)
if trim == 'True':
    fastqFiles = Trimmomatic(fastqFiles,phred)
print 'list file succeed'
#========  (2) define group ===============================
#defined above
#========  (3) align using bwa ============================
map_sam = bwa_vari(read_group,fastqFiles,bwaDb,thread)
#========  (4) Convert sam to sorted bam ==================
sort_bams = sam2bam_sort(map_sam,thread)
#========  (5) Markduplicates using picard ================
dedup_files = markduplicates(sort_bams)
remove(sort_bams)
#========  2. Indel realignment  ====================================
#========  (6) Create a target list of intervals===========
interval = RealignerTargetCreator(dedup_files,ref_fa)
#========  (7) realignment of target intervals ============
realign_bams = IndelRealigner(dedup_files,ref_fa,interval)
remove(dedup_files)
#========  3. Base quality recalibration  =================

# since we don't have dbsnp for CHO, we need to:
# 1. find snp without recalibration, got vcf file
# 2. extract the snps we think are real snps, into a real_vcf file.
# 3. use the file in 2 to do the recalibration.

##=================  Part II. Variant Calling  ======================
#========  1. call raw variant using HaplotypeCaller  =====
#========  (1) determine parameters  ======================
#========  (2) call variant  ==============================
#BQSR_round = 1
#for roundNum in range(BQSR_round):
roundNum = 1
raw_gvcf_files = HaplotypeCaller_DNA_gVCF(realign_bams,ref_fa,thread)
#========  (3) Joint Genotyping ===========================
joint_gvcf_file = JointGenotype(raw_gvcf_files,ref_fa,organism)
#*********** since we don't have the dbsnp for CHO, we need to repeat 
#*********** base reaclibration until it converge.
#========  (4) Variant hard filter  =======================
gold_files = HardFilter(joint_gvcf_file,ref_fa)
#========  (5) Base Recalibration  ========================
recal_bam_files = BaseRecalibrator(realign_bams,ref_fa,gold_files[0],
                                       gold_files[1],roundNum)
#======== second round ====================================
roundNum = 2
raw_gvcf_files = HaplotypeCaller_DNA_gVCF(recal_bam_files,ref_fa,thread)
joint_gvcf_file = JointGenotype(raw_gvcf_files,ref_fa,organism)

gold_files = HardFilter(joint_gvcf_file,ref_fa)
recal_bam_files = BaseRecalibrator(realign_bams,ref_fa,gold_files[0],
                                       gold_files[1],roundNum)
#========  !!! merge lanes for the same sample ============
if len(recal_bam_files!=1):
    merged_bams = rg_bams(read_group,recal_bam_files)
    remove(recal_bam_files)
    dedup_files = markduplicates(merged_bams)
    remove(merged_bams)
    realign_bams = IndelRealigner(dedup_files,ref_fa,interval)
    remove(dedup_files)
    
#========  (6) call variant  ==============================
raw_gvcf_files = HaplotypeCaller_DNA_gVCF(realign_bams,ref_fa,thread)
#========  (7) Joint Genotyping  ==========================
joint_gvcf_file = JointGenotype(raw_gvcf_files,ref_fa,organism)

#========  (8) VQSR  ======================================
#gold_files = HardFilter(joint_gvcf_file,ref_fa)
recal_variant = VQSR(joint_gvcf_file,gold_files[0],gold_files[1],ref_fa)

Message(param['endMessage'],email)
##=================  Part III. Analyze Variant  =====================
