import subprocess
"""
this file includes all functions about txedo pipeline. (cufflinks,cufmerge,cuffdiff) 
"""
def cufflinks(bamFiles,annotation,thread,otherParameters=['']):
    """
    this function run cufflinks.
    bamFiles:           a list of bam files. [1.bam,2.bam,...]
    annotation:         gtf or gff files
    thread:             umber of cores
    otherParameters:    additional parameters in a list [param1,param2,...]
    """
    cmd = ''
    result = []
    for bam in bamFiles:
        out_dir = bam[:-9]
        result.append(out_dir)
        cufflinkCmd = ('cufflinks -o {out_dir} -p {thread} -G {annotation} {bam}').format(
                        out_dir=out_dir,thread=thread,annotation=annotation,bam=bam)
        cmd = cmd + cufflinkCmd + ' ' + ' '.join(otherParameters) + ' && '
    subprocess.call(cmd[:-3],shell=True)
    
    return result