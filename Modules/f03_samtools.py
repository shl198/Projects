import subprocess
def sam2bam_sort(samfiles,thread=1):
    """
    This function will change sam file to bam file and sort bam files
    """
    sorted_files = []
    sam2bamCmd=''
    rmSamCmd=''
    sortBamCmd=''
    indexCmd=''
    rmBamCmd=''
    for sam in samfiles:
        bam = sam[:-3] + 'bam'
        sort = bam[:-3] + 'sort'
        sort_bam = sort + '.bam'
        sorted_files.append(sort_bam)
        # sam file to bam file
        sam2bamCmd = sam2bamCmd + ('samtools view -@ {thread} -bS {sam} '
                                   '-o {bam} && ').format(sam=sam, bam=bam,thread=thread)
        # remove sam file
        rmSamCmd = rmSamCmd + 'rm {sam} & '.format(sam=sam)
        # sort bam file
        sortBamCmd = sortBamCmd + ('samtools sort -m 4G -@ {thread} -T {sort} {bam} '
                                   '-o {sortbam} && ').format(thread=thread,sort=sort,
                                    bam=bam,sortbam=sort_bam)
        # index bam file
        indexCmd = indexCmd + 'samtools index {sortbam} & '.format(sortbam=sort_bam)
        # remove unsorted bam file
        rmBamCmd = rmBamCmd + 'rm {bam} & '.format(bam=bam)
        
    subprocess.call(sam2bamCmd[:-3],shell=True)    
    subprocess.call(rmSamCmd[:-3],shell=True)
    subprocess.call(sortBamCmd[:-3],shell=True)
    subprocess.call(indexCmd[:-3],shell=True)
    subprocess.call(rmBamCmd[:-3],shell=True)
    return sorted_files
 
def bam2sam(bamfiles):
    """
    This function will transfer bam file to sam files
    """
    samFiles = []
    cmd = ''
    for bam in bamfiles:
        samfile = bam[:-3] + 'sam'
        samFiles.append(samfile)
        # command
        cmd = cmd + 'samtools view {bam} > {sam} & '.format(bam=bam,sam=samfile)
    subprocess.call(cmd[:-3],shell=True)
    
    return samFiles
#===========================================================================
#===========================================================================
def extract_mapped(map_result):
    """
    This function extract the mapped reads from the sam file which 
    is got from mapping using aligners. input map_result is a list
    of mapped sam file
    """
    # define returned files
    returnFile = []
    cmd = ''
    for mapfile in map_result:
        filename = mapfile[:-9] + '.mapped.sort.bam'
        returnFile.append(filename)
        cmd = cmd + 'samtools view -F 4 -bh {input} > {output} & '.format(input=mapfile,output=filename)
    subprocess.call(cmd[:-3],shell=True)
    
    return returnFile

def merge_bam(bamfiles,outputbam):
    """
    this function merges bam files into one
    """
    bam = ' '.join(bamfiles)
    cmd = ('samtools merge -f {output} {input}').format(output=outputbam,input=bam)
    subprocess.call(cmd,shell=True)
    
    
def index_bam(bamFiles):
    """
    This function indexes bam files for easy of access
    """
    cmd = ''
    for bam in bamFiles:
        cmd = cmd + 'samtools index {bam} & '.format(bam=bam)
    subprocess.call(cmd[:-3],shell=True)
    print 'done'
    
    
