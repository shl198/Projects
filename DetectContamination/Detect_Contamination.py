"""
This file builds pipeline for detecting any contaminations in cho RNAseq samples.
"""
#=============  import required packages  =============
from Modules.f01_list_files import list_fastq
from Modules.f02_aligner_command import gsnap,blastn
from Modules.f03_samtools import *
from Modules.f06_reads_stat import *
import os
from os import listdir
#file_path = sys.argv[1]
file_path = '/data/RNAseq/seq/108' 
chok1_annotation = '/opt/genome/cho/chok1.gff3'
ntdb = '/data/140822ntdatabase/nt'
#=============  define some parameters  ===============
thread = '6'
Trim = 'False'
microDb_path = '/home/shangzhong/Database/bacteria_fungi/gsnap_bacteriaDb'
microDb_name ='microbial'
#========  (0) enter the directory ================
os.chdir(file_path)
"""
#========  (1) read files  ======================
fastqFiles = list_fastq(file_path,Trim)
print 'list file succeed'
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
"""
#========  (9) integrate annotation to blast result =======
blast_files = [f for f in listdir(file_path) if f.endswith('.blast.txt')]
blast_files.sort()
annotateBlast(blast_files,'nucleotide')
print 'annotate blast result succeed'
print 'folder ' + file_path + ' done'