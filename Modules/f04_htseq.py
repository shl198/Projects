import subprocess,os
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
