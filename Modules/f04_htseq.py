import subprocess,os
from natsort import natsorted
import pandas as pd
import HTSeq

def htseq_count(sortedBamFiles,annotation,outputpath,annotationSource):
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
    cmd = ''
    for bamfile in sortedBamFiles: 
        outputfile = outputpath + '/' + bamfile[:-3] + 'txt'
        htseqCmd = ('htseq-count -f bam -s no -t {type} -i {gene} {bam} {annotation} '
        '> {output} & ').format(type=seqType,gene=id_attr,bam=bamfile, 
                                annotation=annotation,output=outputfile)
        cmd = cmd + htseqCmd
    subprocess.call(cmd + 'wait',shell=True)
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
    res.to_csv(path+'/Count.txt',sep='\t')

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