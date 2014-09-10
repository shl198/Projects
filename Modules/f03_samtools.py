import subprocess
def sam2bam_sort(samfiles):
    """
    This function will change sam file to bam file and sort bam files
    """
    sorted_files = []
    for sam in samfiles:
        bam = sam[:-3] + 'bam'
        sort = bam[:-3] + 'sort'
        sort_bam = sort + '.bam'
        sorted_files.append(sort_bam)
        # sam file to bam file
        sam2bamCmd = 'samtools view -bS {sam} > {bam}'.format(sam=sam, bam=bam)
        subprocess.call(sam2bamCmd,shell=True)
        # remove sam file
        removeCmd = 'rm {sam}'.format(sam=sam)
        subprocess.call(removeCmd,shell=True)
        # sort bam file
        sortBamCmd = 'samtools sort {bam} {bamsort}'.format(bam=bam,bamsort=sort)
        subprocess.call(sortBamCmd,shell=True)
        # index bam file
        indexCmd = 'samtools index {sortbam}'.format(sortbam=sort_bam)
        subprocess.call(indexCmd,shell=True)
        
    return sorted_files
 
def bam2sam(bamfiles):
    """
    This function will transfer bam file to sam files
    """
    samFiles = []
    for bam in bamfiles:
        samfile = bam[:-3] + 'sam'
        samFiles.append(samfile)
        # command
        cmd = 'samtools view {bam} > {sam}'.format(bam=bam,sam=samfile)
        subprocess.call(cmd,shell=True)
    
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
    for mapfile in map_result:
        filename = mapfile[:-9] + '.mapped.sort.bam'
        returnFile.append(filename)
        cmd = 'samtools view -F 4 -bh {input} > {output}'.format(input=mapfile,output=filename)
        subprocess.call(cmd,shell=True)
    
    return returnFile