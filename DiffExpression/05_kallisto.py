import subprocess,os,sys
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f00_Message import Message
from Modules.f01_list_trim_fq import list_files,Trimmomatic
from Modules.p01_FileProcess import get_parameters

def kallisto_index(transcript_fa,index_name,kmer,index_dir):
    """
    This function build index for sailfish
    
    * transcript_fa: fasta file that stores all mRNA of target organism, should not be zipped
    * index_name: index name
    * kmer: integer to set kmer
    * outdir: pathway to store the index files
    """
    if not os.path.exists(index_dir):
        os.mkdir(index_dir)
    os.chdir(index_dir)
    cmd = ('kallisto index -i {index} -k {kmer} {fa}'
           ).format(index=index_name,fa=transcript_fa,kmer=str(kmer))
    subprocess.call(cmd.split(' '))

# transcript_fa = '/data/shangzhong/RibosomeProfiling/Database/kallisto_index/combined_mRNA.fa'
# index_name = 'cho'
# kmer = 25
# index_dir = '/data/shangzhong/RibosomeProfiling/Database/kallisto_index'
# kallisto_index(transcript_fa,index_name,kmer,index_dir)
#============  kallisto quant ===============================
def kallisto_quant(fqFiles,index,thread=1):
    """
    This function runs kallisto
    
    * fqFiles: list. a list of fq files can be either paired end or sing end
    * index: str. Full path plust index name
    * thread: int. Number of thread
    """
    cmd = ''
    # define output path
    for fq in fqFiles:
        if fq[0].endswith('fq.gz'):
            out_dir = fq[0][:-6]
        else:
            out_dir = fq[0][:-9]
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        # run quant
        if len(fq) == 1:
            cmd = cmd + ('kallisto quant -i {index} -o {outPath} -t {thread} --plaintext '
                         '--single -l 50 {fq} && ').format(
                        index=index,outPath=out_dir,thread=thread,fq=fq[0])
    subprocess.call(cmd[:-3],shell=True)
    
    
#=========== define some parameters ===================
parFile = '/home/shangzhong/Codes/Pipeline/DiffExpression/05_kallisto_Parameters.txt'
#parFile = sys.argv[1]
param = get_parameters(parFile)

thread = param['thread']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']

file_path = param['filePath']
kallisto_index = param['kallisto_index']

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
        fastqFiles = Trimmomatic(trimmomatic,fastqFiles,phred,trimmoAdapter,batch=6)
        print 'trim succeed'
        print 'fastqFiles is: ',fastqFiles
    except:
        print 'trim failed'
        Message('trim failed',email)
        raise
#=========== (2) run sailfish mapping ====================
try:
    kallisto_quant(fastqFiles,kallisto_index,thread)
    print 'kallisto quant succeed'
except:
    print 'kallisto quant failed'
    raise

Message(endMessage,email)   
