"""
this file does variant calling for RNAseq
"""
#=============  import required packages  =================
import os
import sys
import subprocess
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.f02_aligner_command import STAR2Pass
from Modules.f03_samtools import sam2bam_sort
from Modules.f07_picard import markduplicates,addReadGroup
from Modules.f08_GATK import *
from Modules.p01_FileProcess import remove,get_parameters,rg_bams
#=============  define some parameters  ===================
"""these parameters and read group names are different for 
   different samples, should only change this part for 
   running pipeline
"""
#parFile = '/data/shangzhong/VariantCall/GlycoKO/rest/GATK_parameters4DNAandRNA.txt'
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
    fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred,trimmoAdapter,batch=6)
sys.stdout.write('list file succeed\n')
sys.stdout.write('fastqFiles is: {fq}\n'.format(fq=fastqFiles))

#========  (2) align using 2 pass STAR ====================
try:
    map_sams= STAR2Pass(fastqFiles,starDb,ref_fa,thread)
    sys.stdout.write('align succeed\n')
    sys.stdout.write('map_sams is: {map}\n'.format(map=map_sams))
except:
    sys.stdout.write('align failed\n')
    Message('align failed',email)
    raise
#========  2. Add read groups, sort,mark duplicates, and create index
#========  (1) sort and add group =========================
try:
    sort_bams = sam2bam_sort(map_sams,thread)
    sys.stdout.write('sort bam succeed\n')
    sys.stdout.write('sort_bams is: {bam}\n'.format(bam=sort_bams))
except:
    sys.stdout.write('sort bam failed\n')
    Message('sort bam failed',email)
    raise

try:
    group_bams = addReadGroup(picard,sort_bams,read_group,batch=9) 
    sys.stdout.write('add group succeed\n')
    sys.stdout.write('group_bams is: {group}\n'.format(group=group_bams))
except:
    sys.stdout.write('add group failed\n')
    Message('add group failed',email)
    raise

#========  (2) mark duplicates ============================
try:
    dedup_bams = markduplicates(picard,group_bams)
    sys.stdout.write('mark duplicate succeed\n')
    sys.stdout.write('dedup_bams is: {dedup}\n'.format(dedup=dedup_bams))
    remove(group_bams)
except:
    sys.stdout.write('mark duplicate failed\n')
    Message('mark duplicate failed',email)
    raise
#========  3. Split 'N' Trim and reassign mapping qualiteies
try:
    split_bams = splitN(gatk,dedup_bams,ref_fa,batch=5)
    sys.stdout.write('split N succeed\n')
    sys.stdout.write('split N is: {N}\n'.format(N=split_bams))
    remove(dedup_bams)
except:
    sys.stdout.write('split N failed\n')
    Message('split N failed',email)
    raise
#========  4. Indel realignment ===========================
#========  (1) generate intervals =========================
try:
    interval = RealignerTargetCreator(gatk,split_bams,ref_fa,thread,batch=6)
    sys.stdout.write('RealignerTarget Creator succeed\n')
    sys.stdout.write('interval is: {int}\n'.format(int=interval))
except:
    sys.stdout.write('RealignerTarget Creator failed\n')
    Message('RealignerTarget Creator failed',email)
    raise
#========  (2) realignment of target intervals ============
try:
    realign_bams = IndelRealigner(gatk,split_bams,ref_fa,interval,batch=6)
    sys.stdout.write('IndelRealigner succeed\n')
    sys.stdout.write('realign bams is: {reali}\n'.format(reali=realign_bams))
    remove(split_bams)
except:
    sys.stdout.write('IndelRealigner failed\n')
    Message('IndelRealigner failed',email)
    raise

#========  5. Base quality recalibration  =================

# since we don't have dbsnp for CHO, we need to:
# 1. find snp without recalibration, got vcf file
# 2. extract the snps we think are real snps, into a real_vcf file.
# 3. use the file in 2 to do the recalibration.

roundNum = 1
#========  6. Variant Calling =============================
try:
    vcf_files = HaplotypeCaller_RNA_VCF(gatk,realign_bams,ref_fa,thread,batch=9)
    sys.stdout.write('1 round call succeed\n')
    sys.stdout.write('vcf_files is: {vcf}\n'.format(vcf=vcf_files))
except:
    sys.stdout.write('1 round call failed\n')
    Message('1 round call failed',email)
    raise
#========  7. Variant filtering ===========================
try:
    gold_varis = RNA_Vari_Filter(gatk,vcf_files,ref_fa,batch=9)
    sys.stdout.write('1 round gold variant succeed\n')
    sys.stdout.write('gold_varis is: {gold}\n'.format(gold=gold_varis))
except:
    sys.stdout.write('1 round gold variant failed\n')
    Message('1 round gold variant failed',email)
    raise
#========  5. Base quality recalibration  =================
try:
    recal_bams = RNA_BaseRecalibrator(gatk,realign_bams,ref_fa,
                gold_varis,roundNum,thread,batch=3)
    sys.stdout.write('recalibration succeed\n')
    sys.stdout.write('recal_bams is: {recal}\n'.format(recal=recal_bams))
except:
    sys.stdout.write('recalibration failed\n')
    Message('recalibration failed',email)
    raise

#========  !!! merge lanes for the same sample ============
if len(recal_bams) !=1:
    try:
        merged_bams = rg_bams(read_group,recal_bams)
        sys.stdout.write('merge_bams is: {mer}\n'.format(mer=merged_bams))
        remove(recal_bams)
    except:
        sys.stdout.write('merge failed\n')
        Message('merge failed',email)
        raise
    try:
        dedup_files = markduplicates(picard,merged_bams)
        sys.stdout.write('dedup_files is: {dedup}\n'.format(dedup=dedup_files))
        remove(merged_bams)
    except:
        sys.stdout.write('mark duplicate merged failed\n')
        Message('mark duplicate merged failed',email)
        raise
    try:
        interval = RealignerTargetCreator(gatk,dedup_files,ref_fa,thread,batch=6)
        realign_bams = IndelRealigner(gatk,dedup_files,ref_fa,interval,batch=6)
        remove(dedup_files)
        sys.stdout.write('realign_bams is: {reali}\n'.format(reali=realign_bams))
        sys.stdout.write('merge lanes succeed\n')
    except:
        sys.stdout.write('realign merged failed\n')
        Message('realign merged failed',email)
        raise
else:
    realign_bams = recal_bams

roundNum = 2
#========  6. Variant Calling =============================
try:
    vcf_files = HaplotypeCaller_RNA_VCF(gatk,realign_bams,ref_fa,thread,batch=9)
    sys.stdout.write('2 round call succeed\n')
    sys.stdout.write('vcf_files is: {vcf}\n'.format(vcf=vcf_files))
    remove(realign_bams)
except:
    sys.stdout.write('2 round call failed\n')
    Message('2 round call failed',email)
    raise
#========  7. Variant filtering ===========================
try:
    gold_varis = RNA_Vari_Filter(gatk,vcf_files,ref_fa,batch=9)
    sys.stdout.write('2 round gold variant succeed\n')
    sys.stdout.write('gold_varis is: {gold}\n'.format(gold=gold_varis))
    sys.stdout.write('job finished succeessfully\n')
except:
    sys.stdout.write('2 round gold vairant failed')
    Message('2 round gold variant failed',email)
    raise

Message(endMessage,email)
