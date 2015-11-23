"""
This file builds pipeline for detecting any contaminations in cho RNAseq samples.
"""
#=============  import required packages  =============
import os
from os import listdir
import sys
sys.path.append('/home/shangzhong/Codes/Pipeline')

from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.f02_aligner_command import gsnap,gsnap_Db,STAR
from Modules.p01_FileProcess import get_parameters
from Modules.f03_samtools import sam2bam_sort
#=============  define some parameters  ===============
parFile = sys.argv[1]
#parFile = '/data/shangzhong/DE/Autism/01_DetectContamination_Parameters.txt'
param = get_parameters(parFile)
thread = param['thread']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']
seqType = param['seqType']
#-------------  sequence parameters  ----------------------
ref_fa = param['ref_fa']
ref_gff = param['ref_gff']
#-------------  software parameters  ----------------------
file_path = param['filePath']
aligner = param['aligner']
alignerDb = param['alignerDb']
trim = param['trim']
phred = param['phred']

picard = param['picard']
trimmomatic = param['trimmomatic']
trimmoAdapter = param['trimmoAdapter']
#-------------  specific parameters for gsnap  ------------
gsnapDbName = param['gsnapDbName']
gsnapAnnotation = param['gsnapAnnotation']
#========  (0) enter the directory ================
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
        if os.listdir(alignerDb) == []:
            gsnap_Db(ref_fa,alignerDb,gsnapDbName,gsnapAnnotation)
        map_files = gsnap(fastqFiles,alignerDb,gsnapDbName,gsnapAnnotation,thread) # [file.sam]
    else:
        map_files = STAR(fastqFiles,alignerDb,thread)
    print 'align succeed'
    print 'map_files is: ',map_files
except:
    print 'align failed'
    Message('host align failed',email)
    raise
#========  (3) sam to bam and sort  ================================
try:
    sorted_bams = sam2bam_sort(map_files,thread)  # [file.sort.bam]
    print 'bam sorted succeed'
    print 'sorted_bam is: ',sorted_bams
except:
    print 'bam sorted failed'
    Message('bam sorted failed',email)
    raise
"""
#========  (2) use gsnap map to bacteria ========
map_files = gsnap(fastqFiles,microDb_path,microDb_name,'',thread)
print 'mapping succeed'
#========  (3) sam2Bam and sort bam =============
sorted_bam = sam2bam_sort(map_files)
print 'sortting succeed'
#========  (4) extract mapped reads =============
mapped_files = extract_mapped(sorted_bam)
print 'extracting succeed'
#========  (5) bam file to sam ==================
samfiles = bam2sam(mapped_files)
print 'bam2sam succeed'

#========  (6) get stats for sam file ===========
outputSam(samfiles) 

#========  (6.1) get stats for sam file =========
#The following steps are for read in sam file directly and do
#the same thing as (6)

allFiles = [f for f in listdir(file_path) if f.endswith(".sam")]
allFiles.sort()
outputSam(allFiles)
print 'output sam result succeed'
#========  (7) view the mapping result in IGV tools =======
#========  (8) blast unique sequence if result from (7)====
#==================  don't make sense  ====================
faFiles = [f for f in listdir(file_path) if f.endswith(".fa")]
blastn(faFiles,ntdb,thread=4)

#========  (9) integrate annotation to blast result =======
blast_files = [f for f in listdir(file_path) if f.endswith('.blast.txt')]
blast_files.sort()
annotateBlast(blast_files,'nucleotide')
print 'annotate blast result succeed'
print 'folder ' + file_path + ' done'
"""