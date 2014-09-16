import subprocess
def sam2bam_sort(samfiles):
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
        sam2bamCmd = sam2bamCmd + 'samtools view -bS {sam} > {bam} & '.format(sam=sam, bam=bam)
        # remove sam file
        rmSamCmd = rmSamCmd + 'rm {sam} & '.format(sam=sam)
        # sort bam file
        sortBamCmd = sortBamCmd + 'samtools sort {bam} {bamsort} & '.format(bam=bam,bamsort=sort)
        # index bam file
        indexCmd = indexCmd + 'samtools index {sortbam} & '.format(sortbam=sort_bam)
        # remove unsorted bam file
        rmBamCmd = rmBamCmd + 'rm {bam} & '.format(bam=bam)
        
    subprocess.call(sam2bamCmd[:-2],shell=True)    
    subprocess.call(rmSamCmd[:-2],shell=True)
    subprocess.call(sortBamCmd[:-2],shell=True)
    subprocess.call(indexCmd[:-2],shell=True)
    subprocess.call(rmBamCmd[:-2],shell=True)
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
    subprocess.call(cmd[:-2],shell=True)
    
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
    subprocess.call(cmd[:-2],shell=True)
    
    return returnFile