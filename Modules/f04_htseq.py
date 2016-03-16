import subprocess,os
import pysam
from natsort import natsorted
import pandas as pd
import HTSeq
from RibosomeProfilePipeline.f02_RiboDataModule import trpr

def chunk(l,n):
    n = max(1,n)
    result = [l[i:i+n] for i in range(0,len(l),n)]
    return result


def htseq_count(sortedBamFiles,annotation,outputpath,annotationSource,batch=1):
    """
    This function use htseq_count to count reads
    * sortedBamFiles: a list of bam files
    * annotation: annotation gff3 files
    * outputpath: pathway to store the hteseq count results
    * annotationSource:     'ncbi'(default)  'ensembl'(alternative) 'genedb'(alternative)
    """
    # 1. check whether outputpath exist
    if not os.path.exists(outputpath):
        os.mkdir(outputpath)
    # 2. check the annotation source
    if annotationSource == 'ncbi':
        seqType = 'exon'
        id_attr = 'gene'
    elif annotationSource == 'ensembl':
        seqType = 'exon'
        id_attr = 'gene_id'
    elif annotationSource == 'genedb':
        seqType = 'CDS'
        id_attr = 'Parent'
    # 3. run htseq-count
    bamFiles = chunk(sortedBamFiles,int(batch))
    for bams in bamFiles:
        cmd = ''
        for bamfile in bams: 
            outputfile = outputpath + '/' + bamfile[:-3] + 'txt'
            htseqCmd = ('htseq-count -f bam -s no -t {type} -i {gene} {bam} {annotation} '
            '> {output} & ').format(type=seqType,gene=id_attr,bam=bamfile, 
                                    annotation=annotation,output=outputfile)
            cmd = cmd + htseqCmd
        print cmd
        subprocess.check_output(cmd + 'wait',shell=True) # %1; echo $?
#         if call[0]=='1':
#             raise
# bam_path = '/data/shangzhong/RibosomeProfiling/TotalRNA_align'
# os.chdir(bam_path)
# sortedBamFiles = [f for f in os.listdir(bam_path) if f.endswith('.bam')]
# annotation = '/data/shangzhong/RibosomeProfiling/Database/new.gff'
# outputpath = '/data/shangzhong/RibosomeProfiling/TotalRNA_align/htseq'
# annotationSource = 'ncbi'
# htseq_count(sortedBamFiles,annotation,outputpath,annotationSource)

def merge_htseqCount(path):
    """
    This function merge all the htseqCount into one matrix
    """
    os.chdir(path)
    files = [f for f in os.listdir(path) if f.endswith('.txt')]
    files = natsorted(files)
    dfs = []
    for f in files:
        df = pd.read_csv(f,sep='\t',header=None,names=['id',f[:-4]],index_col=0)
        dfs.append(df)
    res = pd.concat(dfs,axis=1)
    res.to_csv(path+'/Count.csv',sep=',')
    
def add_rep_count(files,outFile):
    """This function adds up technical replicate for each gene
    * files: list. A list of htseq count files
    """
    dfs = []
    for f in files:
        df = pd.read_csv(f,sep='\t',header=None,index_col=0,names=['geneid',f])
        dfs.append(df)
    merge_df = pd.concat(dfs,axis=1)
    merge_df['count'] = merge_df.sum(axis=1)
    res_df = merge_df['count']
    res_df.to_csv(outFile,sep='\t',header=0)
    

def htseq_count_py(gffFile,bamFile):
    gff_handle = HTSeq.GFF_Reader(gffFile)
    exons = HTSeq.GenomicArrayOfSets('auto',stranded=False)
    counts = {}
    for feature in gff_handle:
        if (feature.type == 'exon') and ('gene' in feature.attr):
            exons[feature.iv] +=feature.attr['gene']
            counts[feature.attr['gene']]=0        
    
    bam_handle = HTSeq.BAM_Reader(bamFile)
    for aln in bam_handle:
        iset = None
        for iv2,step_set in exons[aln.iv].steps():
            if iset is None:
                iset = step_set.copy()
            else:
                iset.intersection_update(step_set)
        if len(iset)==1:
            counts[list(iset)[0]] += 1
    df = pd.DataFrame(counts.items(),columns=['gene_name','count'])
    df = df.sort()
    outFile = bamFile.split('.')[0]+'_count.txt'
    df.to_csv(outFile,sep='\t',index=False)
    
    
def fpkm_from_htseq(bam_path,ruv_path,exn_file):
    """
    This function calculates fpkm from the htseq-count results.
    * bam_path: pathway that has bam files. Used to get total mapped reads.
    * ruv_path: pathway that has ruvseq corrected count data.
    * exn_file: 6 columns. including ['chr','start','end','geneid','traccess','strand'].
    output file that ends with .fpkm.
    """
    os.chdir(bam_path)
    bams = [f for f in os.listdir(bam_path) if f.endswith('.bam')]
    bams = natsorted(bams)
    # 1. get total count
    totalCount = []
    for b in bams:
        bamHandle = pysam.AlignmentFile(b,'rb')
        totalCount.append(bamHandle.mapped)
    # 2. get rna_obj
    rna_df = pd.read_csv(exn_file,sep='\t',header=0,low_memory=False)
    rna_obj = trpr(rna_df)
    # 3. get length for each gene
    os.chdir(ruv_path)
    norm_count_files = [f for f in os.listdir(ruv_path) if f.endswith('.txt')]
    norm_count_files = natsorted(norm_count_files)
    for fn,total in zip(norm_count_files,totalCount):
        df = pd.read_csv(fn,sep=' ',header=None,names=['geneid','count'],index_col=0,low_memory=False)
        df['len'] = df.index.map(lambda x: rna_obj.get_gene_trpr_len(x,multi_chrom='Y'))
        df['fpkm'] = df['count']/float(total)/df['len']*10**9
        df['fpkm'].ix[:-20].to_csv(fn[:-3]+'fpkm.txt',sep='\t')