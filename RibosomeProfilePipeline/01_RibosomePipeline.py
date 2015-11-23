"""
this file process ribosome profiling data
"""
#=============  import required packages  =================
import os
import sys,subprocess
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.p01_FileProcess import remove,get_parameters,rg_bams
from Modules.f02_aligner_command import bowtie2,tophat
from Modules.f06_txedo import cufflinks
from Modules.f03_samtools import sortBam
#=============  define some parameters  ===================
"""these parameters and read group names are different for 
   different samples, should only change this part for 
   running pipeline
"""

#parFile = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_RiboPro_Parameters.txt'
parFile = sys.argv[1]
param = get_parameters(parFile)
thread = param['thread']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']

file_path = param['filePath']
alignerrRNADb = param['alignerrRNADb']
alignerDb = param['alignerDb']

trim = param['trim']
phred = param['phred']
trimmomatic = param['trimmomatic']
trimmoAdapter = param['trimmoAdapter']

annotation = param['annotation']

##*****************  Part 0. Build index file for bowtie ******
##=================  Part I. Preprocess  ============================
#========  1. map reads ==================================
#========  (0) enter the directory =======================
os.chdir(file_path)
Message(startMessage,email)
#========  (1) read and trim files  ======================
fastqFiles = list_files(file_path)
if trim == 'True':
    fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred,batch=6)  # [[filename.fq.gz]]
print 'list file succeed'
print 'fastqFiles is: ',fastqFiles
#========  (2) align to rRNA =============================
try:
    noRNA_fqs = bowtie2(fastqFiles,alignerrRNADb,thread,otherParameters=['--un-gz']) # [[filename.norna.fq.gz]]
    print 'extract norRNA succeed'
    print 'noRNA_fqs is: ',noRNA_fqs
except:
    print 'extract norRNA failed'
    Message('extract norRNA failed',email)
    raise
#========  (3) align to reference using tophat ===========
try:
    map_folders = tophat(noRNA_fqs,alignerDb,annotation,thread,['--no-novel-juncs','-g 1']) # [filename.norna.tophat]
    print 'tophat alignment succeed'
    print 'map_folders is: ',map_folders
except:
    print 'tophat alignment failed'
    Message('tophat alignment failed',email)
    raise
#========  (3) extract the perfect match =================
try:
    bam4cufflinks = []
    for f in map_folders:
        bam = f + '/accepted_hits.bam'
        better_bam = f[:-3] + '.bam'                         
        bam4cufflinks.append(better_bam)                    # [filename.norna.top.bam]
        cmd = ('samtools view -h {fst_map} | grep -E {pattern} | '
               'samtools view -bS - > {out}').format(fst_map=bam,
                pattern='\'(NM:i:[012])|(^@)\'',out=better_bam)
        subprocess.call(cmd,shell=True)
    print 'extract perfect match succeed'
    print 'perfect matched bam is:',bam4cufflinks
except:
    print 'extract perfect match failed'
    Message('extract perfect match failed',email)
    raise
#========  (4) sort the bam file ==================
try:
    sortBams = sortBam(bam4cufflinks,thread,sortType='')
    print 'sort bam succeed'
    print 'sort bams is:',sortBams
except:
    print 'sort bam failed'
    Message('sort bam failed',email)
"""
#========  (4) quantify the transcripts ==================
try:
    result = cufflinks(bam4cufflinks,annotation,thread)  # [filename.nor_cufflinks]
    print 'cufflinks succeed'
    print 'result is: ',result
except:
    print 'cufflinks failed'
    Message('cufflinks failed',email)
"""
Message(endMessage,email)