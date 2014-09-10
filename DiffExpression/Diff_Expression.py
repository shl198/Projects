"""
This file does differential expression.
"""
#============ import required packages ================
import os
import sys
sys.path.append('/home/shangzhong/Scripts/Modules')
from f01_list_files import *
from f02_aligner_command import tophat
#from reads_stat import combine_results,  sam_result, write_file
database = '/home/shangzhong/Database/mouse/mouse'
#file_path = "/data/shangzhong/Diff_express/Esko_project/Processed_data/" + \
#             "Sample_mbmec_Ext1_Ext2_M09"
annotation = '/home/shangzhong/Database/mouse/mouse.gff3'
htseq_dir = '/data/shangzhong/Diff_express/Esko_project/htseqcount'
#=============  define some parameters  ===============

#========== (1) read files =================
#========== (2) concatenate files ===============
def combinePEfastq(fastqFiles):
    """
    This function combine the technical paired end replicates together
    And will return the names of combined files
    """
    first = []  ;   second = []
    for fastq in fastqFiles:
        first.append(fastq[0])
        second.append(fastq[1])
        output1 = file_path[-3:] + '_1.fq.gz'
        output2 = file_path[-3:] + '_2.fq.gz'
    combineCmd = 'cat {files} > {output}'.format(files = (' ').join(first), 
                                                 output = output1)
    os.system(combineCmd)
    combineCmd = 'cat {files} > {output}'.format(files = (' ').join(second), 
                                                 output = output2)
    os.system(combineCmd)
    return [[output1,output2]]

#========== (3) map using tophat ================
#========== (4) sort bam files and run htseq ====
def sort_htseq(map_result):
    for bam_path in map_result:
        
        os.chdir(file_path + '/' + bam_path)
        sortCmd = 'samtools sort accepted_hits.bam hits.sort'
        indexCmd = 'samtools index hits.sort.bam'
        os.system(sortCmd)
        os.system(indexCmd)
    #---------- run htseq ----------------
        output = htseq_dir + '/' + 'htseq_' + bam_path[:-4] + '.txt'
        htseqCmd = ('htseq-count -f bam -r pos -s no -i Dbxref hits.sort.bam '
                    '{annotation} > {outputfile}').format(annotation=annotation,
                                outputfile = output)
        os.system(htseqCmd)
#========== run the pipeline ====================================
def diff_exp(file_path):
    fastqFiles = list_fastq(file_path)
    print 'list fastq files done'
    run_files = combinePEfastq(fastqFiles)
    print 'combine fastq files done'
    map_result = tophat(run_files)
    print 'mapping done'
    sort_htseq(map_result)
    print 'sort and htseq done'
    print map_result[0][-4:] + 'done'


file_path = sys.argv[1]
diff_exp(file_path)
print 'done'