"""
this pipeline annotation variant calling results in vcf file and then 
use provean to predict the effect
"""
import sys,subprocess,os
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f11_snpEff_provean import *
from Modules.p01_FileProcess import get_parameters
from Modules.f00_Message import Message
from Modules.p05_ParseGff import *
from multiprocessing import Process
parFile = sys.argv[1]
#parFile = '/data/shangzhong/DNArepair/Annotation/Annotation_Parameters.txt'
param = get_parameters(parFile)
# parameters
thread = param['thread']
pathway = param['pathway']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']
# database reference
fastaFile = param['reference']
record_dict = SeqIO.index(fastaFile,'fasta')
gffFile = param['annotation']
genome = param['genome']
CodonFile = param['CodonFile']
# software parameters
snpSift = param['snpSift']
snpEff = param['snpEff']
provean = param['provean']
support_set_path = param['support_set']
provean_res_path = param['provean_results']
# other parameters
gene_file = param['gene_file']

#===============================================================================
#        Variant analysis pipeline
#===============================================================================
def chunk(l,n):
    n = max(1,n)
    res = [l[i:i+n] for i in range(0,len(l),n)]
    return res

def get_genes_from_file(gene_file):
    """read gene list from the file and return a list of gene symbols"""
    if gene_file == '':
        genes = ['']
    else:
        genes = []
        gene_df = pd.read_csv(gene_file,header=0,names=['GeneID'])
        genes = gene_df['GeneID'].tolist()
    return genes

def get_all_folders(pathway):
    """put each pair of vcf,vcf.idx files into separate folder, return folders"""
    folders = []
    files = [f for f in os.listdir(pathway) if f.endswith('.merged.filter.vcf')]
    files = natsorted(files)
    for f in files:
        fp = f[:-18]
        folders.append(fp)
        if not os.path.exists(fp): os.mkdir(fp)
        os.rename(f,fp+'/'+f)
        os.rename(f+'.idx',fp+'/'+f+'.idx')
    print 'list directories succeeds'
    print 'folders are:',folders
    return folders

def prepare_fa_vari(workdir,snpEff,snpSift,email,genome,genes,record_dict,gffFile,CodonFile):
    """
    Prepare files for running provean, each folder should only have vcf and vcf.idx file
    * workdir: the folder that has vcf files
    * snpEff: path to snpEff
    * snpSift: path to snpSift
    * email: email or phone number (number@txt.att.net)
    * genome: genome name defined in snpEff
    * genes: A list of gene symbols
    * record_dict: 
    """
    os.chdir(workdir) # set work directory
    vcfFiles = [f for f in os.listdir(workdir) if f.endswith('filter.vcf')]
    vcfFile = vcfFiles[0]
    proteinFiles = [];variantFiles = []
    #============= 1. Annotate vcf results using snpEff ================
    annotatedVCF = snpEff_annotateVCF(vcfFile,snpEff,genome)  # annotated: filename.eff.vcf
    #============= 2. Loop for every genes ================================
    for gene in genes:
        print gene,'start to get input files for provean'
        if gene == '':
            try:
                filteredVCF = snpSift_filterVCF(annotatedVCF,snpSift,
                            ['((ANN[*].IMPACT=\'HIGH\') | (ANN[*].IMPACT=\'MODERATE\'))'])
            except:
                print gene,'snpSift filter failed'
                Message('snpSift filter failed',email)
        else:
            gene_if = ('(ANN[*].GENE=\'{gene}\')').format(gene=gene)
            #============= (1). Filter the annotated file ========================
            try:
                filteredVCF = snpSift_filterVCF(annotatedVCF,snpSift,[gene_if,'&'
                            '((ANN[*].IMPACT=\'HIGH\') | (ANN[*].IMPACT=\'MODERATE\'))'])
                print 'filteredVCF is: ',filteredVCF
            except:
                print gene,'snpSift filter failed'
                Message('snpSift filter failed',email)
        #============= (2). Get input files for provean ======================
        try:
            [protein_files,variant_files] = vcf2input4provean(filteredVCF,record_dict,gffFile,CodonFile)
        except:
            print gene,'fail to get provean inputs'
            Message('fail to get provean inputs',email)
            raise
        if protein_files != '':
            proteinFiles.extend(protein_files)
            variantFiles.extend(variant_files)
            print gene,'prepare for provean input succeeds'
        else:
            print gene,'does not have interested variants'
            raise

Message(startMessage,email)
genes = get_genes_from_file(gene_file)
#================= 0. list directories =========================================
os.chdir(pathway) # set work directory
folders = get_all_folders(pathway)
#============= 2. prepare input files for provean ======================================
batch_folders = chunk(folders,int(thread))
for batch in batch_folders:
    proc = [Process(target=prepare_fa_vari,args=(pathway+'/'+f,snpEff,snpSift,email,genome,genes,record_dict,gffFile,CodonFile,)) for f in batch]
    for p in proc:
        p.start()
    for p in proc:
        p.join()
#============= 3. Run provean ======================================    
# # support set for provean, it can help proven skip the time consuming blast step
support_set = [f for f in os.listdir(support_set_path) if f.endswith('.sss')]
for folder in folders:
    workdir = pathway+'/'+folder
    os.chdir(workdir)
    proteinFiles = sorted([f for f in os.listdir(workdir) if f.endswith('protein.fa')])
    variantFiles = sorted([f for f in os.listdir(workdir) if f.endswith('variant.txt')])
    if not os.path.exists(provean_res_path): os.mkdir(provean_res_path)
    provean_result = provean_res_path +'/'+folder+'_proveanScore.txt'
    try:
        capture_provean_scores(provean_result,provean,proteinFiles,variantFiles,support_set_path,support_set,thread)
        print folder,'folder analysis succeeds'
    except:
        print 'capture provean scores failed'
        Message('capture provean scores failed',email)
        raise
    #============= 4. move the sss support to the standard pathway ======================================
    new_support_set = [f for f in os.listdir(pathway+'/'+folder) if f.endswith('.sss')]
    for f in new_support_set:
        os.rename(f,support_set_path+'/'+f)
        os.rename(f+'.fasta',support_set_path+'/'+f+'.fasta')
    for p in proteinFiles: os.remove(p)
    for v in variantFiles: os.remove(v)
#============= 4. Merge provean results ======================================
outFile = pathway+'provean_final_result.txt'
try:
    merge_provean_results(provean_res_path,outFile)
    print 'merge succeed'
except:
    print 'merge failed'
Message(endMessage,email)


