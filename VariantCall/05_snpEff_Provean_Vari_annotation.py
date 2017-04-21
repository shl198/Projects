"""
this pipeline annotation variant calling results in vcf file and then 
use provean to predict the effect
"""
import sys,subprocess,os
sys.path.append('/home/shangzhong/Codes/Projects')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f11_snpEff_provean import *
from Modules.p01_FileProcess import get_parameters
from Modules.f00_Message import Message
from Modules.p05_ParseGff import *
from multiprocessing import Pool,Process
import glob
parFile = sys.argv[1]
#parFile = '/data/hooman/LdhaCombiClones/Annotation_Parameters.txt'
param = get_parameters(parFile)
# parameters
thread = param['thread']
pathway = param['pathway']
pr_path = param['protein_path']
email = param['email']
startMessage = param['startMessage']
endMessage = param['endMessage']
# database reference
fastaFile = param['reference']
record_dict = SeqIO.index(fastaFile,'fasta')
gffFile = param['annotation']
genome = param['genome']
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
        gene_df = pd.read_csv(gene_file,header=None,names=['GeneID'])
        genes = gene_df['GeneID'].tolist()
    return genes

def get_all_folders(pathway):
    """put each pair of vcf,vcf.idx files into separate folder, return folders"""
    folders = []
    files = [f for f in os.listdir(pathway) if f.endswith('.merged.filter.vcf')]
    if files != []:
        files = natsorted(files)
        for f in files:
            fp = f[:-18]
            folders.append(fp)
            if not os.path.exists(fp): os.mkdir(fp)
            os.rename(f,fp+'/'+f)
            os.rename(f+'.idx',fp+'/'+f+'.idx')
    else:
        all_folders = [fo for fo in os.listdir(pathway) if os.path.isdir(fo)]
        for folder in all_folders:
            fns = [f for f in os.listdir(folder) if f.endswith('merged.filter.vcf')]
            if fns != []:
                folders.append(folder)
    print 'list directories succeeds'
    print 'folders are:',folders
    return folders

def get_all_pr_seq(gff,fa,genes,folder):
    '''
    this function get all protein seqeunces of all the genes provided, each protein stores in 
    a separate file and the file anme is gene:rnaid
    * folder is the full path
    '''
    if not os.path.exists(folder): os.mkdir(folder)
    os.chdir(folder)
    # index fa 
    ref_dic = SeqIO.index(fa,'fasta')
    # read gff
    df = pd.read_csv(gff,sep='\t',header=None,comment='#')
    gff_obj = ncbi_gff(df)
    all_id_df = gff_obj.get_all_id()
    gene_id_df = all_id_df[all_id_df['GeneSymbol'].isin(genes)]
    gene_id_df = gene_id_df[gene_id_df['PrAccess'].values != '-'].drop_duplicates()
    for g,r,p in zip(gene_id_df['GeneSymbol'],gene_id_df['TrID'],gene_id_df['PrAccess']):
        AA = gff_obj.get_gene_seq(ref_dic, p, id_type='pr')
        AA = str(AA[1])
        fn = g+':'+r
        if not os.path.exists(os.path.join(folder,fn+'.protein.fa')):
            with open(fn+'.protein.fa','w') as f:
                f.write('>'+fn+'\n'+AA[:-1])
    

def prepare_vari(workdir,snpEff,snpSift,email,genome,genes):
    """
    Prepare files for running provean, each folder should only have vcf and vcf.idx file
    * workdir: the folder that has vcf files
    * snpEff: path to snpEff
    * snpSift: path to snpSift
    * email: email or phone number (number@txt.att.net)
    * genome: genome name defined in snpEff
    * genes: A list of gene symbols
    """
    
    vcfFiles = glob.glob(workdir+'/*.filter.vcf')
    vcfFile = vcfFiles[0]
    #============= 1. Annotate vcf results using snpEff ================
    annotatedVCF = vcfFile[:-3] + 'eff.vcf'
    if not os.path.exists(annotatedVCF):
        annotatedVCF = snpEff_annotateVCF(vcfFile,snpEff,genome)  # annotated: filename.eff.vcf
    #============= 2. Loop for every gene ================================
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
            vari_files = vari_input4provean(filteredVCF)
        except:
            print gene,'fail to get provean inputs'
            Message('fail to get provean inputs',email)
            raise
        if vari_files == '':
            print gene,'does not have any interested variants'
    print workdir,'provean input succeed'

# Message(startMessage,email)
genes = get_genes_from_file(gene_file)
#================= 0. list directories =========================================
os.chdir(pathway) # set work directory
folders = get_all_folders(pathway)
folders = natsorted(folders)
#============= 1. get all protein sequences ============================================
get_all_pr_seq(gffFile,fastaFile,genes,pr_path)
#============= 2. prepare input files for provean ======================================
# prepare_vari(pathway+'/'+folders[0],snpEff,snpSift,email,genome,genes)
batch_folders = chunk(folders,int(thread))
for batch in batch_folders:
    proc = [Process(target=prepare_vari,args=(pathway+'/'+f,snpEff,snpSift,email,genome,genes,)) for f in batch]
    for p in proc:
        p.start()
    for p in proc:
        p.join()
#============= 3. Run provean ======================================    
for folder in folders:
    # support set for provean, it can help provean skip the time consuming blast step
    support_set = glob.glob(support_set_path+'/*.sss') #[f for f in os.listdir(support_set_path) if f.endswith('.sss')]
    workdir = pathway+'/'+folder
    variantFiles = sorted(glob.glob(workdir+'/vari/*variant.txt'))
    gene_rna = [os.path.split(f)[1][:-12] for f in variantFiles]
    proteinFiles = sorted([f for f in glob.glob(pr_path+'/*protein.fa') if os.path.split(f)[1][:-11] in gene_rna])
    if not os.path.exists(provean_res_path): os.mkdir(provean_res_path)
    provean_result = provean_res_path +'/'+folder+'_proveanScore.txt'
    try:
        capture_provean_scores(provean_result,provean,proteinFiles,variantFiles,support_set_path,support_set,thread)
        print folder,'folder analysis succeeds'
    except:
        print 'capture provean scores failed'
        Message('capture provean scores failed',email)
        raise
#     #============= 4. move the sss support to the standard pathway ======================================
#     new_support_set = glob.glob(pathway+'/'+folder+'/*.sss') #[f for f in os.listdir(pathway+'/'+folder) if f.endswith('.sss')]
#     for f in new_support_set:
#         if os.path.exists(f+'.fasta'):
#             os.rename(f,support_set_path+'/'+f)
#             os.rename(f+'.fasta',support_set_path+'/'+f+'.fasta')
#============= 4. Merge provean results ======================================
outFile = pathway+'/provean_final_result.txt'
try:
    merge_provean_results(provean_res_path,outFile)
    print 'merge succeed'
except:
    print 'merge failed'
Message(endMessage,email)

