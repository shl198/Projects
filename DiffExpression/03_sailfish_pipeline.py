"""
This function runs sailfish
"""
import subprocess,os,sys
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.p01_FileProcess import get_parameters
##=============== some functions ===============
def sailfish_index(transcript_fa,kmer,outdir,thread):
    """
    This function build index for sailfish
    
    * transcript_fa: fasta file that stores all mRNA of target organism, should not be zipped
    * kmer: integer to set kmer
    * outdir: pathway to store the index files
    * thread: integer to set # of threads to use
    """
    cmd = ('sailfish index -t {transcript} -k {kmer} -o {outdir} '
           '-p {thread} -f').format(transcript=transcript_fa,kmer=str(kmer),
                                outdir=outdir,thread=str(thread))
    subprocess.call(cmd.split(' '))
    
    
def sailfish_run(fqFiles,index_dir,libtype,thread):
    """
    This function runs quantification step of sailfish
    
    * fqFiles: a list of fastq files. should be [[1.fq.gz,2.fq.gz]] for pair end, 
                or [[1.fq.gz]] for single end
    * index_dir: the pathway of the index files
    * libtype: library type
    * out_dir: output directory
    * thread: int. number of thread to use
    """
    # check the output directory
    cmd = ''
    for fq in fqFiles:
        if fq[0].endswith('fq.gz'):
            out_dir = fq[0][:-6]
        else:
            out_dir = fq[0][:-9]
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        # run cmd
        if len(fq) == 2:
            cmd = cmd + ('sailfish quant -i {index} -l {libtype} '
                   '-1 <(gunzip -c {fq1}) -2 <(gunzip -c {fq2}) -o {out} -p {thread} && ').format(index=index_dir,
                libtype=libtype,fq1=fq[0],fq2=fq[1],out=out_dir,thread=thread)
        else:
            cmd = cmd + ('sailfish quant -i {index} -l {libtype} '
                    '-r <(gunzip -c {fq}) -o {out} -p {thread} && ').format(index=index_dir,
                libtype=libtype,fq1=fq[0],out=out_dir,thread=thread)
    subprocess.call(cmd[:-3],shell=True,executable='/bin/bash')


#=========== define some parameters ===================
#parFile = '/data/shangzhong/DE/Autism/03_sailfish_Parameters.txt'
parFile = sys.argv[1]
param = get_parameters(parFile)

thread = param['thread']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']

file_path = param['filePath']
sailfish_index = param['sailfish_index']
libtype = '\"' + param['libtype'] + '\"'

trim = param['trim']
phred = param['phred']
trimmomatic = param['trimmomatic']
trimmoAdapter = param['trimmoAdapter']


#===============================================================================
#                 Pipeline
#===============================================================================
#=========== (0) enter the directory ================
Message(startMessage,email)
os.chdir(file_path)
#=========== (1) reads files and trim ====================
fastqFiles = list_files(file_path)
print 'list file succeeded'
if trim == 'True':
    try:
        fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred,trimmoAdapter)
        print 'trim succeed'
        print 'fastqFiles is: ',fastqFiles
    except:
        print 'trim failed'
        Message('trim failed',email)
        raise
#=========== (2) run sailfish mapping ====================
try:
    sailfish_run(fastqFiles,sailfish_index,libtype,thread)
    print 'sailfish quant succeed'
except:
    print 'sailfish quant failed'
    raise

Message(endMessage,email)    