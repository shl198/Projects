"""
this file does variant calling for DNAseq
"""
#=============  import required packages  =================
import os
import sys,subprocess
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.f02_aligner_command import bwa_vari
from Modules.f03_samtools import sam2bam_sort
from Modules.f07_picard import markduplicates
from Modules.f08_GATK import *
from Modules.p01_FileProcess import remove,get_parameters,rg_bams
#=============  define some parameters  ===================
"""these parameters and read group names are different for 
   different samples, should only change this part for 
   running pipeline
"""

#parFile = '/data/shangzhong/VariantCall/CHOS/DNA_Parameters.txt'
parFile = sys.argv[1]
param = get_parameters(parFile)
thread = param['thread']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']

ref_fa = param['refSequence']
file_path = param['filePath']
bwaDb = param['alignerDb']
trim = param['trim']
phred = param['phred']
picard = param['picard']
trimmomatic = param['trimmomatic']

gatk = param['gatk']
read_group = param['readGroup']
organism = param['organism']
##*****************  Part 0. Build index file for bwa and GATK ******
##=================  Part I. Preprocess  ============================
#========  1. map and dedupping =====================================
#========  (0) enter the directory ========================
os.chdir(file_path)
Message(startMessage,email)
#========  (1) read files  ================================

fastqFiles = list_files(file_path)
if trim == 'True':
    fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred)
print 'list file succeed'
print 'fastqFiles is: ',fastqFiles
#========  (2) define group ===============================
#defined above
#========  (3) align using bwa ============================
try:
    map_sam = bwa_vari(read_group,fastqFiles,bwaDb,thread)
    print 'align succeed'
    print 'map_sam is: ',map_sam
except:
    print 'align failed'
    sys.exit(1)
#========  (4) Convert sam to sorted bam ==================
try:
    sort_bams = sam2bam_sort(map_sam,thread)
    print 'sort bam files succeed'
    print 'sort_bams is: ',sort_bams
except:
    print 'sort bam files failed'
    sys.exit(1)
#========  (5) Markduplicates using picard ================
try:
    dedup_files = markduplicates(picard,sort_bams)
    print 'mark duplicates succeed'
    print 'dedup_files is: ',dedup_files
    remove(sort_bams)
except:
    print 'mark duplicates failed'
    sys.exit(1)
#========  2. Indel realignment  ====================================
#========  (6) Create a target list of intervals===========
try:
    interval = RealignerTargetCreator(gatk,dedup_files,ref_fa,thread)
    print 'RealignerTarget Creator succeed'
    print 'interval is: ',interval
except:
    print 'RealignerTarget Creator failed'
    sys.exit(1)
#========  (7) realignment of target intervals ============
try:
    realign_bams = IndelRealigner(gatk,dedup_files,ref_fa,interval)
    print 'IndexRealigner succeed'
    print 'realign_bams is: ',realign_bams
    remove(dedup_files)
except:
    print 'IndelRealigner failed'
    sys.exit(1)
#========  3. Base quality recalibration  =================

# since we don't have dbsnp for CHO, we need to:
# 1. find snp without recalibration, got vcf file
# 2. extract the snps we think are real snps, into a real_vcf file.
# 3. use the file in 2 to do the recalibration.

##=================  Part II. Variant Calling  ======================
#========  1. call raw variant using HaplotypeCaller  =====
#========  (1) determine parameters  ======================
#========  (2) call variant  ==============================
roundNum = 1
try:
    raw_gvcf_files = HaplotypeCaller_DNA_gVCF(gatk,realign_bams,ref_fa,thread)
    print 'round 1 call succeed'
    print 'raw_gvcf_files is: ',raw_gvcf_files
except:
    print 'round 1 call failed'
    sys.exit(1)
#========  (3) Joint Genotyping ===========================
try:
    joint_gvcf_file = JointGenotype(gatk,raw_gvcf_files,ref_fa,organism,thread)
    print 'round 1 join vcf succeed'
    print 'joint_gvcf_file is: ',joint_gvcf_file
except:
    print 'round 1 join vcf failed'
    sys.exit(1)
#*********** since we don't have the dbsnp for CHO, we need to repeat 
#*********** base reaclibration until it converge.
#========  (4) Variant hard filter  =======================
try:
    gold_files = HardFilter(gatk,joint_gvcf_file,ref_fa,thread)
    print 'round 1 gold files succeed'
    print 'gold_files is: ',gold_files
except:
    print 'round 1 gold files failed'
    sys.exit(1)
#========  (5) Base Recalibration  ========================
try:
    recal_bam_files = BaseRecalibrator(gatk,realign_bams,ref_fa,gold_files[0],
                                           gold_files[1],roundNum,thread)
    print 'round 1 recalibration succeed'
    print 'recal_bam_files is: ',recal_bam_files
except:
    print 'round 1 recalibration failed'
    sys.exit(1)
#======== second round ====================================
roundNum = 2
try:
    raw_gvcf_files = HaplotypeCaller_DNA_gVCF(gatk,recal_bam_files,ref_fa,thread)
    print 'round 2 call succeed'
    print 'raw_gvcf_files is:',raw_gvcf_files
except:
    print 'round 2 call failed'
    sys.exit(1)
#------- Joint Genotyping --------
try:
    joint_gvcf_file = JointGenotype(gatk,raw_gvcf_files,ref_fa,organism,thread)
    print 'round 2 join vcf succeed'
    print 'joint_gvcf_file is: ',joint_gvcf_file
except:
    print 'round 2 join vcf failed'
    sys.exit(1)
#------- Hard filter -------------
try:
    gold_files = HardFilter(gatk,joint_gvcf_file,ref_fa,thread)
    print 'round 2 gold files succeed'
    print 'gold_files is: ',gold_files
except:
    print 'round 2 gold files failed'
    sys.exit(1)
#------- Recalibration -----------
try:
    recal_bam_files = BaseRecalibrator(gatk,realign_bams,ref_fa,gold_files[0],
                                           gold_files[1],roundNum,thread)
    print 'round 2 recalibration succeed'
    print 'recal_bam_files is: ',recal_bam_files
except:
    print 'round 2 recalibration failed'
    sys.exit(1)
#========  !!! merge lanes for the same sample ============
if len(recal_bam_files) !=1:
    try:
        merged_bams = rg_bams(read_group,recal_bam_files)
        print 'merged succeed'
        print 'merged_bams is: ',merged_bams
        remove(recal_bam_files)
    except:
        print 'merged filed'
        sys.exit(1)
    try:
        dedup_files = markduplicates(picard,merged_bams)
        print 'dedup succeed'
        print 'merged dedup_files is: ',dedup_files
        remove(merged_bams)
    except:
        print 'merged dedup failed'
        sys.exit(1)
    try:
        interval = RealignerTargetCreator(gatk,dedup_files,ref_fa,thread)
        realign_bams = IndelRealigner(gatk,dedup_files,ref_fa,interval)
        print 'merged indelrealigner succeed'
        print 'merged realign_bams is: ',realign_bams
        remove(dedup_files)
    except:
        print 'merged realign failed'
        sys.exit(1)
    #========  (6) call variant  ==============================
    try:
        raw_gvcf_files = HaplotypeCaller_DNA_gVCF(gatk,realign_bams,ref_fa,thread)
        print 'merged final call succeed'
        print 'raw_gvcf_files is:',raw_gvcf_files
    except:
        print 'final call failed'
        sys.exit(1)
    #========  (7) Joint Genotyping  ==========================
    try:
        joint_gvcf_file = JointGenotype(gatk,raw_gvcf_files,ref_fa,organism,thread)
        print 'final joint succeed'
        print 'joint_gvcf_file is: ',joint_gvcf_file
    except:
        print 'final joint succeed'
        sys.exit(1)
else:
    # for only one file, just run calling with recalibration bam file
    try:
        raw_gvcf_files = HaplotypeCaller_DNA_gVCF(gatk,realign_bams,ref_fa,thread)
        print 'final call succeed'
        print 'raw_gvcf_files is:',raw_gvcf_files
    except:
        print 'fianl call succeed'
#========  (8) VQSR  ======================================
#gold_files = HardFilter(joint_gvcf_file,ref_fa)
try:
    recal_variant = VQSR(gatk,joint_gvcf_file,gold_files[0],gold_files[1],ref_fa,thread)
    print 'vcf recalibration succeed'
    print 'recal_variant is: ',recal_variant
except:
    print 'final vcf recalibration failed'
    sys.exit(1)
Message(endMessage,email)
##=================  Part III. Analyze Variant  =====================
