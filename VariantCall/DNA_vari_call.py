"""
this file does variant calling for DNAseq
"""
#=============  import required packages  =================
import os
from Modules.f01_list_files import list_fastq
from Modules.f02_aligner_command import bwa_vari
from Modules.f03_samtools import sam2bam_sort
from Modules.f07_picard import markduplicates
from Modules.f08_GATK import *
#=============  define some parameters  ===================
thread = '10'
file_path = '/data/shangzhong/VariantCall/DNAseq'
Trim = 'False'
bwa_chok1Db = '/opt/genome/cho/bwa_chok1Db/bwachok1'
chok1_fa = '/opt/genome/cho/cgr_ref_CriGri_1.0_chrUn.fa'
samplename = 'CHOS'
#samplename = sys.argv[2]
##=================  Part I. Preprocess  ============================
#========  1. map and dedupping =====================================
#========  (0) enter the directory ========================
os.chdir(file_path)
#========  (1) read files  ================================
fastqFiles = list_fastq(file_path,Trim)
print 'list file succeed'
#========  (2) define group ===============================
read_group = """'@RG\\tID:chosgroup1\\tSM:sample1\\tPL:illumina\\tLB:lib1\\tPU:unit1'"""
#========  (3) align using bwa ============================
#map_sam = bwa_vari(read_group,fastqFiles,bwa_chok1Db,thread)
#========  (4) Convert sam to sorted bam ==================
#sort_bams = sam2bam_sort(map_sam)
#========  (5) Markduplicates using picard ================
#dedup_files = markduplicates(sort_bams)
#========  2. Indel realignment  ====================================
#========  (6) Create a target list of intervals===========
#interval = RealignerTargetCreator(dedup_files,chok1_fa)
#========  (7) realignment of target intervals ============
#realign_bams = IndelRealigner(dedup_files,chok1_fa,interval)
#========  3. Base quality recalibration  =================
"""
since we don't have dbsnp for CHO, we need to:
1. find snp without recalibration, got vcf file
2. extract the snps we think are real snps, into a real_vcf file.
3. use the file in 2 to do the recalibration.
"""
##=================  Part II. Variant Calling  ======================
#========  1. call raw variant using HaplotypeCaller  =====
#========  (1) determine parameters  ======================
#========  (2) call variant  ==============================
#raw_gvcf_files = HaplotypeCaller_DNA_gVCF(realigned_files,chok1_fa,thread)
#========  (3) Joint Genotyping ===========================
#joint_gvcf_files = JointGenotype(raw_gvcf_files,chok1_fa,samplename)
raw_gvcf = 'CHOS.raw.g.vcf'
realignbams = ['trim_1_140605_BC42FFACXX_P1135_201_1.reali.bam']
#*********** since we don't have the dbsnp for CHO, we need to repeat 
#*********** base reaclibration until it converge.
BQSR_round = 2
for i in range(BQSR_round):
#========  (4) Variant hard filter  =======================
    gold_files = HardFilter(raw_gvcf,chok1_fa)
#========  (5) Base Recalibration  ========================
    recal_bam_files = BaseRecalibrator(realignbams,chok1_fa,gold_files[0],
                                       gold_files[1])
#========  (6) call variant  ==============================
    raw_gvcf_files = HaplotypeCaller_DNA_gVCF(recal_bam_files,chok1_fa,thread)
#========  (7) Joint Genotyping  ==========================
    joint_gvcf_files = JointGenotype(raw_gvcf_files,chok1_fa,samplename)
#========  (8) VQSR  ======================================
recal_variant = VQSR(raw_gvcf,gold_files[0],gold_files[1],chok1_fa)


##=================  Part III. Analyze Variant  =====================
