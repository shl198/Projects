import subprocess
def sam2bam_sort(samfiles,thread=1,sortType=''):
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
        if sortType == 'name':
            tag = ' -n'
        else:
            tag = ''
        sortBamCmd = sortBamCmd + ('samtools sort{tag} -m 4G -@ {thread} -T {sort} {bam} '
                                   '-o {sortbam} && ').format(tag=tag,thread=thread,sort=sort,
                                    bam=bam,sortbam=sort_bam)
        # index bam file
        indexCmd = indexCmd + 'samtools index {sortbam} & '.format(sortbam=sort_bam)
        # remove unsorted bam file
        rmBamCmd = rmBamCmd + 'rm {bam} & '.format(bam=bam)
        
    subprocess.check_call(sam2bamCmd[:-3],shell=True)    
    subprocess.check_call(rmSamCmd[:-3],shell=True)
    subprocess.check_call(sortBamCmd[:-3],shell=True)
    subprocess.check_call(indexCmd + 'wait',shell=True)
    subprocess.check_call(rmBamCmd + 'wait',shell=True)
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
    subprocess.check_call(cmd + 'wait',shell=True)
    
    return samFiles
#===========================================================================
#===========================================================================
def extract_bam(sortedBamFiles,extractType,seqType='pair',thread=1):
    """
    This function extract the mapped/unmapped reads from the bam file which 
    is got from mapping using aligners. 
    
    * sortedBamFiles: a list of bam files
    * extractType: default(map), other(unmap)
    * thread: number of threads
    """
    if seqType == 'single':
        sam_tag = '4'
    elif seqType == 'pair':
        sam_tag = '12'
    # define returned files
    returnFile = []
    cmd = ''
    if extractType == 'map':
        for bam in sortedBamFiles:
            filename = bam[:-3] + 'map.bam'
            returnFile.append(filename)
            cmd = cmd + 'samtools view -@ {thread} -F {tag} -bh {input} > {output} && '.format(
                        thread=str(thread),tag=sam_tag,input=bam,output=filename)
    else:
        for bam in sortedBamFiles:
            filename = bam[:-3] + 'unmap.bam'
            returnFile.append(filename)
            cmd = cmd + 'samtools view -@ {thread} -f {tag} -bh {input} > {output} && '.format(
                        thread=str(thread),tag=sam_tag,input=bam,output=filename)
    subprocess.call(cmd[:-3],shell=True)
    
    return returnFile

def merge_bam(bamfiles,outputbam):
    """
    this function merges bam files into one
    """
    bam = ' '.join(bamfiles)
    cmd = ('samtools merge -f {output} {input}').format(output=outputbam,input=bam)
    subprocess.check_call(cmd,shell=True)
    
    
def index_bam(bamFiles):
    """
    This function indexes bam files for easy of access
    """
    cmd = ''
    for bam in bamFiles:
        cmd = cmd + 'samtools index {bam} & '.format(bam=bam)
    subprocess.check_call(cmd + 'wait',shell=True)
    print 'done'
    
    
