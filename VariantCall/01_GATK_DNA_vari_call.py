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
from Modules.f02_aligner_command import bwa_vari,bwa_Db
from Modules.f03_samtools import sam2bam_sort
from Modules.f07_picard import markduplicates,sortVCF
from Modules.f08_GATK import *
from Modules.p01_FileProcess import remove,get_parameters,rg_bams
from Modules.p02_ParseFasta import divide_scaffold_by_len
#=============  define some parameters  ===================
"""these parameters and read group names are different for
   different samples, should only change this part for
   running pipeline
"""

#parFile = '/data/shangzhong/DNArepair/GATK_parameters4DNAandRNA.txt'
parFile = sys.argv[1]
param = get_parameters(parFile)
thread = param['thread']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']

ref_fa = param['refSequence']
file_path = param['filePath']
bwaIndex = param['alignerDb']
trim = param['trim']
phred = param['phred']
picard = param['picard']
trimmomatic = param['trimmomatic']
trimmoAdapter = param['trimmoAdapter']

gatk = param['gatk']
read_group = param['readGroup']
organism = param['organism']
##*****************  Part 0. Build index file for bwa and GATK ******
##=================  Part I. Preprocess  ============================
#========  1. map and dedupping =====================================
Message(startMessage,email)
#========  (0) enter the directory ========================
bwa_path = bwaIndex[:bwaIndex.rfind('/')]
if not os.path.exists(bwa_path): os.mkdir(bwa_path)
if os.listdir(bwa_path) == []:
    bwa_Db(bwa_path,ref_fa)
os.chdir(file_path)
#========  (1) read files  ================================
fastqFiles = list_files(file_path)
if trim == 'True':
    trim_fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred,trimmoAdapter,batch=6) 
    remove(fastqFiles)
else:
    trim_fastqFiles = fastqFiles
print 'list file succeed'
print 'fastqFiles is: ',trim_fastqFiles
#========  (2) define group ===============================
#defined above
#========  (3) align using bwa ============================
try:
    map_sam = bwa_vari(read_group,trim_fastqFiles,bwaIndex,thread)
    print 'align succeed'
    print 'map_sam is: ',map_sam
except:
    print 'align failed'
    Message('align failed',email)
    raise
#========  (4) Convert sam to sorted bam ==================
try:
    sort_bams = sam2bam_sort(map_sam,thread)
    print 'sort bam files succeed'
    print 'sort_bams is: ',sort_bams
except:
    print 'sort bam files failed'
    Message('sort bam files failed',email)
    raise
#========  (5) Markduplicates using picard ================
try:
    dedup_files = markduplicates(picard,sort_bams)
    print 'mark duplicates succeed'
    print 'dedup_files is: ',dedup_files
    remove(sort_bams)
except:
    print 'mark duplicates failed'
    Message('mark duplicates failed',email)
    raise
#========  2. Indel realignment  ====================================
#========  (1) Create a target list of intervals===========
try:
    interval = RealignerTargetCreator(gatk,dedup_files,ref_fa,thread)
    print 'RealignerTarget Creator succeed'
    print 'interval is: ',interval
except:
    print 'RealignerTarget Creator failed'
    Message('RealignerTarget Creator failed',email)
    raise
#========  (2) realignment of target intervals ============
try:
    realign_bams = IndelRealigner(gatk,dedup_files,ref_fa,interval,9)  # (gatk,dedupbams,reference,intervals,batch=1,*gold_indels)
    print 'IndexRealigner succeed'
    print 'realign_bams is: ',realign_bams
    remove(dedup_files)
except:
    print 'IndelRealigner failed'
    Message('IndelRealigner failed',email)
    raise
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
    raw_gvcf_files = HaplotypeCaller_DNA_gVCF(gatk,realign_bams,ref_fa,thread,batch=int(thread))
    print 'round 1 call succeed'
    print 'raw_gvcf_files is: ',raw_gvcf_files
except:
    print 'round 1 call failed'
    Message('round 1 call failed',email)
    raise
#========  (3) Joint Genotyping ===========================
try:
    joint_gvcf_file = JointGenotype(gatk,raw_gvcf_files,ref_fa,organism,thread)
    print 'round 1 join vcf succeed'
    print 'joint_gvcf_file is: ',joint_gvcf_file
    remove(raw_gvcf_files)
except:
    print 'round 1 join vcf failed'
    Message('round 1 join vcf failed',email)
    raise
#*********** since we don't have the dbsnp for CHO, we need to repeat 
#*********** base reaclibration until it converge.
#========  (4) Variant hard filter  =======================
try:
    gold_files = HardFilter(gatk,joint_gvcf_file,ref_fa,thread)
    print 'round 1 gold files succeed'
    print 'gold_files is: ',gold_files
    remove(joint_gvcf_file)
except:
    print 'round 1 gold files failed'
    Message('round 1 gold files failed',email)
    raise
#========  (5) Base Recalibration  ========================
try:
    recal_bam_files = BaseRecalibrator(gatk,realign_bams,ref_fa,gold_files[0],
                                           gold_files[1],roundNum,thread,6)
    print 'round 1 recalibration succeed'
    print 'recal_bam_files is: ',recal_bam_files
    remove(realign_bams)
except:
    print 'round 1 recalibration failed'
    Message('round 1 recalibration failed',email)
    raise
# #======== second round ====================================
# roundNum = 2
# try:
#     raw_gvcf_files = HaplotypeCaller_DNA_gVCF(gatk,recal_bam_files,ref_fa,thread)
#     print 'round 2 call succeed'
#     print 'raw_gvcf_files is:',raw_gvcf_files
# except:
#     print 'round 2 call failed'
#     Message('round 2 call failed',email)
#     raise
# #------- Joint Genotyping --------
# try:
#     joint_gvcf_file = JointGenotype(gatk,raw_gvcf_files,ref_fa,organism,thread)
#     print 'round 2 join vcf succeed'
#     print 'joint_gvcf_file is: ',joint_gvcf_file
#     remove(raw_gvcf_files)
# except:
#     print 'round 2 join vcf failed'
#     Message('round 2 join vcf failed',email)
#     raise
# #------- Hard filter -------------
# try:
#     gold_files = HardFilter(gatk,joint_gvcf_file,ref_fa,thread)
#     print 'round 2 gold files succeed'
#     print 'gold_files is: ',gold_files
#     remove(joint_gvcf_file)
# except:
#     print 'round 2 gold files failed'
#     Message('round 2 gold files failed',email)
#     raise
# #------- Recalibration -----------
# try:
#     recal_bam_files = BaseRecalibrator(gatk,realign_bams,ref_fa,gold_files[0],
#                                            gold_files[1],roundNum,thread)
#     print 'round 2 recalibration succeed'
#     print 'recal_bam_files is: ',recal_bam_files
#     remove(realign_bams)
# except:
#     print 'round 2 recalibration failed'
#     Message('round 2 recalibration failed',email)
#     raise
#========  !!! merge lanes for the same sample ============
if len(recal_bam_files) !=1:
    #========= merge samples  =========================
    try:
        merged_bams = rg_bams(read_group,recal_bam_files)
        print 'merged succeed'
        print 'merged_bams is: ',merged_bams
        remove(recal_bam_files)
    except:
        print 'merged failed'
        Message('merged failed',email)
        raise
    #========= mark duplicates ========================
    try:
        dedup_files = markduplicates(picard,merged_bams)
        print 'dedup succeed'
        print 'merged dedup_files is: ',dedup_files
        remove(merged_bams)
    except:
        print 'merged dedup failed'
        Message('merged dedup failed',email)
        raise
    #========= Realignment ============================
    try:
        interval = RealignerTargetCreator(gatk,dedup_files,ref_fa,thread)
        realign_bams = IndelRealigner(gatk,dedup_files,ref_fa,interval)
        print 'merged indelrealigner succeed'
        print 'merged realign_bams is: ',realign_bams
        remove(dedup_files)
    except:
        print 'merged realign failed'
        Message('merged realign failed',email)
        raise
    #========  (6) call variant  ==============================
    try:
        raw_gvcf_files = HaplotypeCaller_DNA_gVCF(gatk,realign_bams,ref_fa,thread,batch=int(thread))
#         raw_gvcf_files,par_L_files = par_HaplotypeCaller_DNA_gVCF(gatk,realign_bams,ref_fa,L_path)
#         remove(par_L_files)
        print 'merged final call succeed'
        print 'raw_gvcf_files is:',raw_gvcf_files
    except:
        print 'final call failed'
        Message('final call failed',email)
        raise
    #========  (7) Joint Genotyping ===========================
    try:
        joint_gvcf_file = JointGenotype(gatk,raw_gvcf_files,ref_fa,organism,thread)
        print 'final joint succeed'
        print 'joint_gvcf_file is: ',joint_gvcf_file
        remove(raw_gvcf_files)
    except:
        print 'final joint failed'
        Message('final joint failed',email)
        raise
else:
    # for only one file, just run calling with recalibration bam file
    try:
        joint_gvcf_file = HaplotypeCaller_DNA_VCF(gatk,recal_bam_files[0],ref_fa,thread)
        print 'final call succeed'
        print 'raw_gvcf_files is:',joint_gvcf_file
    except:
        print 'final call failed'
        Message('final call failed',email)
        raise
#========  (8) VQSR or Hard filter  ======================================
# since for CHO samples we don't have enough samples and snp resources, the VQSR step cannot give a very good prediction.
# we choose to use hardFilter.
try:
    final_filtered_files = HardFilter(gatk,joint_gvcf_file,ref_fa,thread)
    print 'final filter succeed'
    print 'final_filtered_files is: ',final_filtered_files
except:
    print 'final filter failed'
    Message('final filter failed',email)
    raise

# try:
#     recal_variant = VQSR(gatk,joint_gvcf_file,gold_files[0],gold_files[1],ref_fa,thread)
#     print 'vcf recalibration succeed'
#     print 'recal_variant is: ',recal_variant
# except:
#     print 'final vcf recalibration failed'
#     Message('final vcf recalibration failed',email)
#     raise
#========  (9) combine snp and indel ======================================
try:
    combinedVcf = CombineSNPandINDEL(gatk,ref_fa,final_filtered_files,'--assumeIdenticalSamples --genotypemergeoption UNSORTED')
    print 'combine snp and indel succeed'
    print 'combineVcf file is: ',combinedVcf
    remove(final_filtered_files)
except:
    print 'combine snp and indel failed'
    raise
Message(endMessage,email)
##=================  Part III. Analyze Variant  =====================
