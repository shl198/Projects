"""
this pipeline quantify expression level of both host and virus
"""
#=============  import required packages  =================
import sys,os
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f00_Message import Message
from Modules.p01_FileProcess import get_parameters
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.f02_aligner_command import gsnap,STAR
from Modules.f03_samtools import sam2bam_sort,extract_bam
from Modules.f04_htseq import htseq_count
from Modules.f07_picard import sam2fastq
from Modules.p01_FileProcess import remove

#=============  define some parameters  ===================
parFile = sys.argv[1]
#parFile = '/data/shangzhong/DE/Winzler/02_QuantifyVirus_Parameters.txt'
param = get_parameters(parFile)
thread = param['thread']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']
seqType = param['seqType']
#-------------  sequence parameters  ----------------------
host_fa = param['hostRefSequence']
host_annotation = param['hostAnnotation']
host_AnnotationSource = param['hostAnnotationSource']
virus_ref_fa = param['virusRefSequence']
virus_annotation = param['virusAnnotation']
virus_AnnotationSource = param['virusAnnotationSource']
#-------------  software parameters  ----------------------
file_path = param['filePath']
aligner = param['aligner']
host_alignerDb = param['host_alignerDb']
virus_alignerDb = param['virus_alignerDb']
trim = param['trim']
phred = param['phred']

picard = param['picard']
trimmomatic = param['trimmomatic']
trimmoAdapter = param['trimmoAdapter']
organism = param['organism']
#-------------  specific parameters for gsnap  ------------
host_gsnapDbName = param['host_gsnapDbName']
host_gsnapAnnotation = param['host_gsnapAnnotation']
virus_gsnapDbName = param['virus_gsnapDbName']
virus_gsnapAnnotation = param['virus_gsnapAnnotation']
#-------------  htseqCount parameters  --------------------
host_htseqFolder = param['host_htseqOutputFolder']
virus_htseqFolder = param['virus_htseqOutputFolder']

#===============================================================================
#          1. map to host reference genome
#===============================================================================
#========  (0) enter the directory ========================
os.chdir(file_path)
Message(startMessage,email)
#========  (1) read files  ================================
fastqFiles = list_files(file_path)
if trim == 'True':
    fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred)     # [[fq.gz]]
print 'list file succeed'
print 'fastqFiles is: ',fastqFiles
#========  (2) align fastq files to host ================================
try:
    if aligner == 'gsnap':
        map_files = gsnap(fastqFiles,host_alignerDb, host_gsnapDbName,host_gsnapAnnotation,thread)
    else:
        map_files = STAR(fastqFiles,host_alignerDb,thread)
    print 'host align succeed'
    print 'map_files is: ',map_files                          # [file.sam]
except:
    print 'host align failed'
    Message('host align failed',email)
    raise
#========  (3) sam to bam and sort  ================================
try:
    sorted_bams = sam2bam_sort(map_files,thread)            # [file.sort.bam]
    print 'host sorted succeed'
    print 'sorted_bam is: ',sorted_bams
except:
    print 'host sorted failed'
    Message('host sorted failed',email)
    raise
#========  (4) get htseq Count to host  ============================
try:
    htseq_count(sorted_bams,host_annotation,host_htseqFolder,host_AnnotationSource)
    print 'host htseqCount succeed'
except:
    print 'host htseq count failed'
    Message('host htseq count failed',email)
    raise
#========  (5) extract unmapped reads  =============================
try:
    unmap2host_bams = extract_bam(sorted_bams,'unmap',seqType,thread)  # [file.sort.unmap.bam]
    print 'extract unmap2host_bams succeed'
    print 'unmap2host_bams is: ',unmap2host_bams
    remove(sorted_bams)
    # rename files
    for f in unmap2host_bams: os.rename(f,f[:-4]+'2host.bam')
    unmap2host_bams = [f[:-4]+'2host.bam' for f in unmap2host_bams]    # [file.sort.unmap2host.bam]
except:
    print 'extract unmap2host_bams failed'
    Message('extract unmap2host_bams failed',email)
    raise
#========  (6) unmap2host_bams to fastq ============================
try:
    unmap2host_fqs = sam2fastq(picard,unmap2host_bams,seqType)    # [[file.sort.unmap2host.fq.gz]]
    print 'unmap2host_fq succeed'
    print 'unmap2host_fqs is: ',unmap2host_fqs
    remove(unmap2host_bams)
except:
    print 'unmap2host_fq failed'
    Message('unmap2host_fq failed',email)
    raise
#===========================================================================
#                 2. Map to virus reference genome
#===========================================================================
#========  (1) align unmap2host_fq to virus  =========================
try:
    if aligner == 'gsnap':
        map_files = gsnap(unmap2host_fqs,virus_alignerDb, virus_gsnapDbName,virus_gsnapAnnotation,thread)
    else:
        map_files = STAR(unmap2host_fqs,virus_alignerDb,thread)  # [file.sort.unmap2host.sam]
    print 'virus align succeed'
    print 'map_files is: ',map_files
except:
    print 'virus align failed'
    Message('virus align failed',email)
    sys.exit(1)
#========  (2) sam to bam and sort  ================================
try:
    sorted_bams = sam2bam_sort(map_files,thread)   # [file.sort.unmap2host.sort.bam]
    print 'virus sorted succeed'
    print 'sorted_bam is: ',sorted_bams
except:
    print 'virus sorted failed'
    Message('virus sorted failed',email)
    raise
#========  (3) get htseq Count to virus  ============================
try:
    htseq_count(sorted_bams,virus_annotation,virus_htseqFolder,virus_AnnotationSource)
    print 'virus htseqCount succeed'
except:
    print 'host htseq count failed'
    Message('host htseq count failed',email)
    raise

Message(endMessage,email)