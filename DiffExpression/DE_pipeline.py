"""
This file run the pipeline of doing differential expression analysis of 
Contaminated CHO samples.
"""
import os 
import subprocess
from Modules.f01_list_files import list_fastq
from Modules.f02_aligner_command import gsnap,STAR
from Modules.f03_samtools import sam2bam_sort
from Modules.f04_htseq import htseq_count
from Modules.f05_GeneIDConvert import ID_Convert
#=========== define some pathways ===================
#file_path = sys.argv[1]
file_path = '/data/shangzhong/Diff_express/star_test'
db_path = '/opt/genome/cho/gsnap_chok1Db'
db_name = 'chok1'
gsnap_annotation = '/opt/genome/cho/gsnap_chok1Db/chok1.maps/chok1.splicesites.iit'
thread = '1'
annotation = '/opt/genome/cho/chok1.gff3'

Dict = '/home/shangzhong/Database/chok1_geneID_symbol.txt'
output_path = '/data/RNAseq/htseq_count' 
inputpath = file_path

subprocess.call('mail shl198@eng.ucsd.edu < /home/shangzhong/start.txt',shell=True)
#=========== (0) enter the directory ================
os.chdir(file_path)
#=========== (1) reads files and trim ===============
fastqFiles = list_fastq(file_path,'False')
#=========== (2) run gsnap to do the mapping ========
#map_files = gsnap(fastqFiles,db_path, db_name,gsnap_annotation,thread)

starDb_path = '/opt/genome/cho/STAR_chok1Db'
map_files = STAR(fastqFiles,starDb_path,thread)
#=========== (3) samtools to sort the file ==========
sorted_bam = sam2bam_sort(map_files)
#=========== (4) htseq_count ========================
htseq_count(sorted_bam,annotation,file_path)
#=========== (5) htseq symbol to id =================
ID_Convert(Dict,output_path,inputpath)

subprocess.call('mail shl198@eng.ucsd.edu < /home/shangzhong/end.txt',shell=True)