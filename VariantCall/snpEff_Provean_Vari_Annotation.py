import sys,subprocess,os
from os import listdir
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f11_snpEff_provean import *
from Modules.p01_FileProcess import get_parameters
from Modules.f00_Message import Message

parFile = sys.argv[1]
#parFile = '/data/shangzhong/VariantCall/CHOS_ToT_VS_chok1/filteredVcfResults/CHOS_TOT_A8/test/Annotation_Parameters.txt'
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
# other parameters
gene_file = param['gene_file']

Message(startMessage,email)

#===============================================================================
#        Variant analysis pipeline
#===============================================================================
## read gene list
genes = []
geneFile = open(gene_file,'r')
for line in geneFile:
    if line[0].isalpha():
        genes.append(line[:-1])
geneFile.close()
genes = list(set(genes))
genes.sort()
#================= 0. list directories =========================================
os.chdir(pathway) # set work directory
folders = [f for f in listdir(pathway) if os.path.isdir(os.path.join(pathway, f))]
print 'list directories succeeds'
print 'folders are:',folders

support_set = [f for f in os.listdir(support_set_path) if f.endswith('.sss')]
for folder in folders:
    workdir = pathway + '/' + folder
    os.chdir(workdir) # set work directory
    vcfFiles = [f for f in listdir(workdir) if f.endswith('.vcf')]
    vcfFile = vcfFiles[0]
    proteinFiles = [];variantFiles = []
    #============= 1. Annotate vcf results using snpEff ================
    annotatedVCF = snpEff_annotateVCF(vcfFile,snpEff,genome)  # annotated: filename.eff.vcf
    #cmd = ('sed -i \'s/{original}/{new}/g\' {file}').format(original='5_prime_UTR_truncation+',
    #                                                        new='5_prime_UTR_truncation\&',file=annotatedVCF)
    #subprocess.call(cmd,shell=True)
    #============= Loop for every genes ================================
    for gene in genes:
        print gene,'start to get input files for provean'
        gene_if = ('(ANN[*].GENE=\'{gene}\')').format(gene=gene)
        #============= 2. Filter the annotated file ========================
        try:
            filteredVCF = snpSift_filterVCF(annotatedVCF,snpSift,[gene_if,'&'
                        '((ANN[*].IMPACT=\'HIGH\') | (ANN[*].IMPACT=\'MODERATE\'))'])
        except:
            print gene,'snpSift filter failed'
            Message('snpSift filter failed',email)
        #============= 3. Get input files for provean ======================
        try:
            [protein_files,variant_files] = vcf2input4provean(filteredVCF,record_dict,gffFile,CodonFile)
        except:
            print gene,'fail to get provean inputs'
            sys.exit(1)
        if protein_files != '':
            proteinFiles.extend(protein_files)
            variantFiles.extend(variant_files)
            print gene,'prepare for provean input succeeds'
        else:
            print gene,'does not have interested variants'
            sys.exit(1)
    #============= 4. Run provean ======================================
    provean_result = 'proveanScore.txt'
    #support_set = [f for f in os.listdir(support_set_path) if f.endswith('.sss')]
    try:
        capture_provean_scores(provean_result,provean,proteinFiles,variantFiles,support_set_path,support_set,thread)
        print folder,'folder analysis succeeds'
    except:
        print 'capture provean scores failed'
        Message('capture provean scores failed',email)
        sys.exit(1)
Message(endMessage,email)
 


