import subprocess
import re
import os
def extractBamWithLength(sam,length):
    handle = open(sam,'r')
    out = sam[:-3]+str(length)+'.sam'
    outHandle = open(out,'w')
    for line in handle:
        sum_nt = 0
        item = line[:-1].split('\t')
        cigar = item[5]
        num_map = re.findall('\d+|\D+', cigar)
        for i in range(len(num_map)):
            if num_map[i] == 'M':
                sum_nt = sum_nt + int(num_map[i-1])
        if sum_nt == length:
            outHandle.write(line)
    outHandle.close()
    
    
def sam2bam_sort(samfiles,thread=1,sortType=''):
    """
    This function will change sam/bam file to sorted bam file
    
    * samfiles:list. A list of sam/bam files
    * thread: int. Number of thread
    * sortType: str. Default: '' which means sort by position. Alter: name
    
    return sorted bam file
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
        if sam.endswith('sam'):
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

    if samfiles[0].endswith('sam'):
        print sam2bamCmd
        subprocess.check_call(sam2bamCmd[:-3],shell=True)
    print rmSamCmd
    subprocess.check_call(rmSamCmd[:-3],shell=True)
    print sortBamCmd
    subprocess.check_call(sortBamCmd[:-3],shell=True)
    if tag == '':
        subprocess.check_call(indexCmd + 'wait',shell=True)
    print rmBamCmd
    subprocess.check_call(rmBamCmd + 'wait',shell=True)
    return sorted_files


def sortBam(bamFiles,thread=1,sortType=''):
    """
    This function sort bam files
    """
    if sortType == 'name':
            tag = ' -n'
    else:
        tag = ''
    cmd = ''
    sortBams = []
    for bam in bamFiles:
        sort = bam[:-3] + 'sort'
        sort_bam = sort + '.bam'
        sortBams.append(sort_bam)
        sortCmd = ('samtools sort{tag} -m 4G -@ {thread} -T {sort} {bam} '
                   '-o {sortbam} && ').format(tag=tag,thread=thread,sort=sort,
                                              bam=bam,sortbam=sort_bam)
        cmd = cmd + sortCmd
    subprocess.call(cmd[:-3],shell=True)
    return sortBams



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
def extract_bam(sortedBamFiles,extractType,seqType='',thread=1):
    """
    This function extract the mapped/unmapped reads from the bam file which 
    is got from mapping using aligners. 
    
    * sortedBamFiles: list. A list of bam files
    * extractType: str. default(map), other(unmap)
    * seqType: str. 'pair' or 'single'
    * thread: int. number of threads
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
    subprocess.call(cmd + 'wait',shell=True)
    

def flagstat(bamFiles):
    """
    This function calculates basic mapping results of bam file(% of reads mapping,etc)
    
    * bamFiles: list. A list of bam files
    """
    cmd = ''
    for bam in bamFiles:
        cmd = cmd + 'samtools flagstat {bam} && '.format(bam=bam)
    subprocess.call(cmd[:-3],shell=True)
    