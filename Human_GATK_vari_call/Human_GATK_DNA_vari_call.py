"""
this file does variant calling for DNAseq
"""
#=============  import required packages  =================
import os
import sys,subprocess
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer so that in the log file all information is printed in order.
from f00_Message import Message
from f01_list_trim_fq import list_files_human,Trimmomatic
from f02_aligner_command import bwa_vari
from f03_samtools import sam2bam_sort
from f07_picard import markduplicates
from f08_GATK import *
from p01_FileProcess import remove,get_parameters,rg_bams
#=============  define some parameters  ===================
"""these parameters and read group names are different for 
   different samples, should only change this part for 
   running pipeline
"""

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
trimmoAdapter = param['trimmoAdapter']
gold_snp = param['dbSNP']
phaseINDEL= param['phase1INDEL']
gold_indel= param['MillINDEL']
omni = param['omni']
hapmap = param['hapMap']

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
fastqFiles = list_files_human(file_path)
if trim == 'True':
    fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred,trimmoAdapter)
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
    Message('align failed',email)
    sys.exit(1)
#========  (4) Convert sam to sorted bam ==================
try:
    sort_bams = sam2bam_sort(map_sam,thread)
    print 'sort bam files succeed'
    print 'sort_bams is: ',sort_bams
except:
    print 'sort bam files failed'
    Message('sort bam files failed',email)
    sys.exit(1)
#========  (5) Markduplicates using picard ================
try:
    dedup_files = markduplicates(picard,sort_bams)
    print 'mark duplicates succeed'
    print 'dedup_files is: ',dedup_files
    remove(sort_bams)
except:
    print 'mark duplicates failed'
    Message('mark duplicates failed',email)
    sys.exit(1)
#========  2. Indel realignment  ====================================
#========  (6) Create a target list of intervals===========
try:
    interval = RealignerTargetCreator(gatk,dedup_files,ref_fa,thread,phaseINDEL,gold_indel)
    print 'RealignerTarget Creator succeed'
    print 'interval is: ',interval
except:
    print 'RealignerTarget Creator failed'
    Message('RealignerTarget Creator failed',email)
    sys.exit(1)
#========  (7) realignment of target intervals ============
try:
    realign_bams = IndelRealigner(gatk,dedup_files,ref_fa,interval,phaseINDEL,gold_indel)
    print 'IndexRealigner succeed'
    print 'realign_bams is: ',realign_bams
    remove(dedup_files)
except:
    print 'IndelRealigner failed'
    Message('IndelRealigner failed',email)
    sys.exit(1)
#========  3. Base quality recalibration  =================
roundNum = '1'
try:
    recal_bam_files = BaseRecalibrator(gatk,realign_bams,ref_fa,gold_snp,
                                           gold_indel,roundNum,thread)
    print 'round 1 recalibration succeed'
    print 'recal_bam_files is: ',recal_bam_files
except:
    print 'round 1 recalibration failed'
    Message('round 1 recalibration failed',email)
    sys.exit(1)

##=================  Part II. Variant Calling  ======================
#========  1. call raw variant using HaplotypeCaller  =====
#========  (1) determine parameters  ======================
#========  (2) call variant  ==============================
#========  !!! merge lanes for the same sample ============
if len(recal_bam_files) !=1:
    #========= (3) merge samples  =========================
    try:
        merged_bams = rg_bams(read_group,recal_bam_files)
        print 'merged succeed'
        print 'merged_bams is: ',merged_bams
        remove(recal_bam_files)
    except:
        print 'merged failed'
        Message('merged failed',email)
        sys.exit(1)
    #========= (4) mark duplicates ========================
    try:
        dedup_files = markduplicates(picard,merged_bams)
        print 'dedup succeed'
        print 'merged dedup_files is: ',dedup_files
        remove(merged_bams)
    except:
        print 'merged dedup failed'
        Message('merged dedup failed',email)
        sys.exit(1)
    #========= (5) Realignment ============================
    try:
        interval = RealignerTargetCreator(gatk,dedup_files,ref_fa,
                                          thread,phaseINDEL,gold_indel)
        realign_bams = IndelRealigner(gatk,dedup_files,ref_fa,
                                      interval,phaseINDEL,gold_indel)
        print 'merged indelrealigner succeed'
        print 'merged realign_bams is: ',realign_bams
        remove(dedup_files)
    except:
        print 'merged realign failed'
        Message('merged realign failed',email)
        sys.exit(1)
    #=========  (6) call variant  ==============================
    try:
        raw_gvcf_files = HaplotypeCaller_DNA_gVCF(gatk,realign_bams,ref_fa,thread)
        print 'merged final call succeed'
        print 'raw_gvcf_files is:',raw_gvcf_files
    except:
        print 'final call failed'
        Message('final call failed',email)
        sys.exit(1)
    #========  (7) Joint Genotyping  ==========================
    try:
        joint_gvcf_file = JointGenotype(gatk,raw_gvcf_files,ref_fa,organism,thread)
        print 'final joint succeed'
        print 'joint_gvcf_file is: ',joint_gvcf_file
    except:
        print 'final joint failed'
        Message('final joint failed',email)
        sys.exit(1)
    #========  (8) VQSR  ======================================
    try:
        recal_variant = VQSR_human(gatk,joint_gvcf_file,ref_fa,thread,hapmap,omni,phaseINDEL,gold_snp,gold_indel)
        print 'vcf recalibration succeed'
        print 'recal_variant is: ',recal_variant
    except:
        print 'final vcf recalibration failed'
        Message('final vcf recalibration failed',email)
        sys.exit(1)
else:
    
    # for only one file, just run calling with recalibration bam file
    #======== Calling variant =================================
    try:
        raw_vcf_file = HaplotypeCaller_DNA_VCF(gatk,recal_bam_files[0],ref_fa,thread)  
        print 'final call succeed'
        print 'raw_gvcf_files is:',raw_vcf_file
    except:
        print 'final call failed'
        Message('final call failed',email)
        sys.exit(1)
#======== Hard filtering ==================================
    try:
        final_filtered_files = HardFilter(gatk,raw_vcf_file,ref_fa,thread)
        print 'final filter succeed'
        print 'final_filtered_files is: ',final_filtered_files
    except:
        print 'final filter failed'
        Message('final filter failed',email)
        sys.exit(1)
    #======== Combine snp and indel ===========================
    try:
        combinedVcf = CombineSNPandINDEL(gatk,ref_fa,final_filtered_files,'--assumeIdenticalSamples --genotypemergeoption UNSORTED')
        print 'combine snp and indel succeed'
        print 'combineVcf file is: ',combinedVcf
    except:
        print 'combine snp and indel failed'
        sys.exit(1)

Message(endMessage,email)