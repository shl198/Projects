"""
this pipeline quantify expression level of both host and virus
"""
#=============  import required packages  =================
import sys,os
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f00_Message import Message
from Modules.p01_FileProcess import get_parameters,remove
from Modules.p02_ParseFasta import fq2fa
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.f02_aligner_command import gsnap,STAR,gsnap_Db,blastn
from Modules.f03_samtools import sam2bam_sort,extract_bam
from Modules.f07_picard import sam2fastq
from Modules.f12_trinity import Trinity
from a01_blast2fasta import blast2fa,merge_fa
from Modules.p04_ParseBlast import annotateBlast
#=============  define some parameters  ===================
parFile = sys.argv[1]
#parFile = '/data/shangzhong/DetectVirus/Brain/01_DetectVirus_Parameters.txt'
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
#-------------  specific parameters for gsnap  ------------
host_gsnapDbName = param['host_gsnapDbName']
host_gsnapAnnotation = param['host_gsnapAnnotation']
virus_gsnapDbName = param['virus_gsnapDbName']
virus_gsnapAnnotation = param['virus_gsnapAnnotation']
#------------- blast parameters ---------------------------
blast_Db = param['blast_Db']
runBlastAfterAssemble = param['runBlastAfterAssemble']
blast_nt_DB = param['blast_nt_DB']
#===============================================================================
#          1. map to host reference genome
#===============================================================================
#========  (0) enter the directory ========================
os.chdir(file_path)
Message(startMessage,email)
#========  (1) read files  ================================
fastqFiles = list_files(file_path)
if trim == 'True':
    fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred)
print 'list file succeed'
print 'fastqFiles is: ',fastqFiles
#========  (2) align fastq files to host ================================
try:
    if aligner == 'gsnap':
        # check index
        if os.listdir(host_alignerDb) == []:
            gsnap_Db(host_fa,host_alignerDb,host_gsnapDbName,host_gsnapAnnotation)
        map_files = gsnap(fastqFiles,host_alignerDb, host_gsnapDbName,host_gsnapAnnotation,thread) # [file.sam]
    else:
        map_files = STAR(fastqFiles,host_alignerDb,thread)
    print 'host align succeed'
    print 'map_files is: ',map_files
except:
    print 'host align failed'
    Message('host align failed',email)
    raise
#========  (3) sam to bam and sort  ================================
try:
    sorted_bams = sam2bam_sort(map_files,thread)  # [file.sort.bam]
    print 'host sorted succeed'
    print 'sorted_bam is: ',sorted_bams
except:
    print 'host sorted failed'
    Message('host sorted failed',email)
    raise
#========  (4) extract reads that unmapped to host  =====
try:
    unmap2host_bams = extract_bam(sorted_bams,'unmap',seqType,thread) # [file.sort.unmap.bam]
    print 'extract unmap2host_bams succeed'
    print 'unmap2host_bams is: ',unmap2host_bams
    remove(sorted_bams)
    # rename files
    for f in unmap2host_bams: os.rename(f,f[:-4]+'2host.bam')
    unmap2host_bams = [f[:-4]+'2host.bam' for f in unmap2host_bams] # [file.sort.unmap2host.bam]
except:
    print 'extract unmap2host_bams failed'
    Message('extract unmap2host_bams failed',email)
    raise
#========  (6) unmap2host_bams to fastq.gz =========================
try:
    unmap2host_fq_gzs = sam2fastq(picard,unmap2host_bams,seqType) # [[file.sort.unmap2host.fq.gz]]
    print 'unmap2host_fq_gzs succeed'
    print 'unmap2host_fq_gzs is: ',unmap2host_fq_gzs
    remove(unmap2host_bams)
except:
    print 'unmap2host_fq_gzs failed'
    Message('unmap2host_fq_gzs failed',email)
    raise
#===========================================================================
#                 2. Map to virus reference genome
#===========================================================================
#========  (1) align unmap2host_fq to virus  =========================
try:
    if aligner == 'gsnap':
        # check index
        if os.listdir(virus_alignerDb) == []:
            gsnap_Db(virus_ref_fa,virus_alignerDb,virus_gsnapDbName,virus_gsnapAnnotation)
        # align the reads
        map_files = gsnap(unmap2host_fq_gzs,virus_alignerDb, virus_gsnapDbName,virus_gsnapAnnotation,thread)
    else:
        map_files = STAR(unmap2host_fq_gzs,virus_alignerDb,thread)  # [file.sort.unmap2host.sam]
    print 'virus align succeed'
    print 'map_files is: ',map_files
    remove(unmap2host_fq_gzs)
except:
    print 'virus align failed'
    Message('virus align failed',email)
    raise
#========  (2) sam to bam and sort  ================================
try:
    sorted_bams = sam2bam_sort(map_files,thread)    # [file.sort.unmap2host.sort.bam]
    print 'virus sorted succeed'
    print 'sorted_bams is: ',sorted_bams
except:
    print 'virus sorted failed'
    Message('virus sorted failed',email)
    raise
#========  (3) extract reads that mapped and unmapped to virus  ============================
try:
    map2virus_bams = extract_bam(sorted_bams,'map',seqType,thread)     # [file.sort.unmap2host.sort.map.bam]
    unmap2virus_bams = extract_bam(sorted_bams,'unmap',seqType,thread) # [file.sort.unmap2host.sort.unmap.bam]
    # rename files
    for f in map2virus_bams: os.rename(f,f[:-29]+'.only2virus.bam')
    map2virus_bams = [f[:-29]+'.only2virus.bam' for f in map2virus_bams]     # [file.only2virus.bam]
    for f in unmap2virus_bams: os.rename(f,f[:-31]+'.map2neither.bam')        # [file.map2neither.bam]
    unmap2virus_bams = [f[:-31]+'.map2neither.bam' for f in unmap2virus_bams]
    print 'extract map and unmap2virus_bams succeed'
    print 'map2virus_bams is: ',map2virus_bams,'unmap2virus_bams is: ',unmap2virus_bams
    remove(sorted_bams)
except:
    print 'extract map and unmap2virus_bams failed'
    Message('extract map and unmap2virus_bams failed',email)
    raise
#========  (4) transfer the mapped and unmapped to virus bam to fastq  =======
try:
    map2virus_fq_gzs = sam2fastq(picard,map2virus_bams,seqType)     # [[file.only2virus.fq.gz]]
    unmap2virus_fq_gzs = sam2fastq(picard,unmap2virus_bams,seqType) # [[file.map2neither.fq.gz]]
    print 'transfer from bam to fq succeed'
    print 'map2virus_fq_gzs is: ',map2virus_fq_gzs,'unmap2virus_fq_gzs',unmap2virus_fq_gzs
    remove(map2virus_bams);remove(unmap2virus_bams)
except:
    print 'transfer from bam to fq failed'
    Message('transfer from bam to fq failed',email)
    raise
#========  (5) transfer the unmapped to virus fastq to fasta =================
try:
    only2virus_faFiles = fq2fa(map2virus_fq_gzs)      # [file.only2virus.fa.gz]
    map2neither_faFiles = fq2fa(unmap2virus_fq_gzs)       # [file.map2neither.fa.gz]
    print 'fq files to fa files succeed'
    print 'only2virus_faFiles is: ',only2virus_faFiles,'map2neither_faFiles is: ',map2neither_faFiles
    remove(map2virus_fq_gzs); remove(unmap2virus_fq_gzs)
except:
    print 'fq files to fa files failed'
    Message('fq files to fa files failed',email)
    raise
#========  (6) map the unmapped to virus fa files with blastn ======================================
try:
    virus_blast = blastn(map2neither_faFiles,blast_Db,thread)  # [file.map2niether.txt]
    print 'blast to virus succeed'
    print 'virus_blast is: ',virus_blast
except:
    print 'blast to virus failed'
    Message('blast to virus failed',email)
    raise
#========  (7) extract reads from blast results ================================
try:
    blastFas = blast2fa(virus_blast,map2neither_faFiles,seqType)   # file.map2neither.blast.fa.gz
    print 'extract reads from blast tabular txt succeed'
    print 'blastFas is: ', blastFas
    remove(virus_blast)
except:
    print 'extract reads from blast tabular txt failed'
    Message('extract reads from blast tabular txt failed',email)
    raise
#========  (8) merge blast mapped with gsnap/STAR mapped ==========================
try:
    mergedFa = merge_fa(seqType,only2virus_faFiles,blastFas)
    print 'merge fa files for Trinity succeed'
    print 'mergedFa is: ',mergedFa
except:
    print 'merge fa files for Trinity failed'
    Message('merge fa files for Trinity',email)
    raise
#========  (9) assemble the reads using trinity ==============================
try:
    Trinity(mergedFa,thread)
    print 'assemble using Trinity succeed'
    remove(mergedFa)
except:
    print 'assemble using Trinity failed'
    Message('assemble using Trinity failed',email)
    raise
#========  (10) run blast to get homology ==============================
if runBlastAfterAssemble == 'Yes':
    ### 1. blast assembled fa file
    faFiles = ['trinity_out_dir/Trinity.fasta']
    try:
        blast_assemble_result = blastn(faFiles,blast_nt_DB,thread)
        print 'blast assembled fa succeed'
        print 'blast_assemble_result is: ',blast_assemble_result
    except:
        print 'blast assembled fa failed'
        Message('blast assembled fa failed',email)
        raise
    ### 2. annotate the blast result
    try:
        annotateBlast(blast_assemble_result[0],'nucleotide')
        print 'annotate blast succeed'
    except:
        print 'annotate blast failed'
        Message('annotate blast failed',email)
        raise
Message(endMessage,email)