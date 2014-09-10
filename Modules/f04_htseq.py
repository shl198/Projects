def htseq_count(sortedBamFiles,annotation,outputpath):
    """
    This function use htseq_count to count reads
    """
    import subprocess
    for bamfile in sortedBamFiles: 
        outputfile = outputpath + '/' + bamfile[:-3] + 'txt'
        htseqCmd = ('htseq-count -f bam -r pos -s no -t exon -i gene {bam} {annotation} '
        '> {output}').format(bam=bamfile, annotation = annotation,output=outputfile)
        subprocess.call(htseqCmd,shell=True)
