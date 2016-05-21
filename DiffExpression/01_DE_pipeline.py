"""
This file run the pipeline of doing differential expression analysis of 
Contaminated CHO samples.
"""
import os
import sys 
import subprocess
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.f02_aligner_command import gsnap,STAR,STAR_Db,bowtie
from Modules.f03_samtools import sam2bam_sort,flagstat
from Modules.f04_htseq import htseq_count
from Modules.f05_IDConvert import geneSymbol2EntrezID
from Modules.p01_FileProcess import get_parameters
from Modules.p01_FileProcess import remove
#=========== define some parameters ===================
#parFile = '/data/shangzhong/DE/FDA/01_DE_Parameters.txt'
parFile = sys.argv[1]
param = get_parameters(parFile)

thread = param['thread']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']

dataSource = param['dataSource']
ref_fa = param['refSequence']
file_path = param['filePath']
db_path = param['alignerDb']
trim = param['trim']
phred = param['phred']
trimmomatic = param['trimmomatic']
trimmoAdapter = param['trimmoAdapter']

aligner = param['aligner']
annotation = param['annotation']
output_path = param['htseqOutPath']
htseqBatch = param['htseqBatch']
db_name = param['gsnapDbName']
gsnap_annotation = param['gsnapAnnotation']

Dict = param['symbolIDFile']
inputpath = file_path

#=========== (0) enter the directory ================
os.chdir(file_path)
Message(startMessage,email)
#=========== (1) reads files and trim ===============

fastqFiles = list_files(file_path)
print 'list file succeed'
if trim == 'True':
    try:
        trim_fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred,trimmoAdapter,batch=6)
        print 'trim succeed'
        print 'fastqFiles is: ',fastqFiles
        remove(fastqFiles)
    except:
        print 'trim failed'
        Message('trim failed',email)
        raise
else:
    trim_fastqFiles = fastqFiles
#=========== (2) run STAR to do the mapping ========
try:
    if aligner == 'gsnap':
        map_files = gsnap(trim_fastqFiles,db_path, db_name,gsnap_annotation,thread)
    elif aligner == 'STAR':
        if not os.path.exists(db_path): os.mkdir(db_path)
        if os.listdir(db_path) == []:
            STAR_Db(db_path,ref_fa,thread)
        map_files = STAR(trim_fastqFiles,db_path,thread,annotation,['--outSAMtype BAM SortedByCoordinate','--quantMode GeneCounts'])
    elif aligner == 'bowtie':
        map_files = bowtie(trim_fastqFiles,db_path,thread=1,otherParameters=[''])
    print 'align succeed'
    print 'map_files is: ',map_files
except:
    print 'align failed'
    Message('align failed',email)
    raise
#=========== (3) samtools to sort the file ==========
try:
    sorted_bams = sam2bam_sort(map_files,thread,'name')
    print 'sorted succeed'
    print 'sorted_bam is: ',sorted_bams
except:
    print 'sorted failed'
    Message('sorted failed',email)
    raise
#=========== (4) get mapping stats ==================
try:
    flagstat(sorted_bams)
    print 'flagstat succeed'
except:
    print 'flagstat failed'
    Message('flagstat failed',email)
    raise
#=========== (4) htseq_count ========================
try:
    htseq_count(sorted_bams,annotation,output_path,dataSource,htseqBatch)
    print 'htseq count succeed'
except:
    print 'htseq count failed'
    Message('htseq count failed',email)
    raise
#=========== (5) htseq symbol to id =================
if dataSource == 'ncbi':
    try:
        geneSymbol2EntrezID(Dict,output_path,output_path)
        print 'id convert succeed'
    except:
        print 'id convert failed'
        Message('id convert failed',email)
        raise
Message(endMessage,email)