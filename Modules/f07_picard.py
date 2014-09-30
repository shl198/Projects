import subprocess

def read_group(ID,sample,platform,library,platunit):
    return ('@RG\\tID:'+ID+'\\tSM:'+sample+'\\tPL:'+platform+'\\tLB:'+library
            +'\\tPU:'+platunit)
    #@RG\\tID:chosgroup1\\tSM:sample1\\tPL:illumina\\tLB:lib1\\tPU:unit1
def sam2sortbam(picard,samfiles):
    """
    this function change samfile to sorted bam file
    """
    SortSam = picard + '/' + 'SortSam.jar'
    sorted_files = []
    for sam in samfiles:
        sort_bam = sam[:-3] + 'sort.bam'
        sorted_files.append(sort_bam)
        cmd = ('java -jar {SortSam} INPUT={input} OUTPUT={output} '
               'SORT_ORDER=coordinate').format(SortSam=SortSam,
                input=sam,output=sort_bam)
        subprocess.call(cmd,shell=True)
    return sorted_files

def markduplicates(picard,sortbam):
    """
    this function mark duplicates of the sorted bam file
    """
    mark = picard + '/' + 'MarkDuplicates.jar'
    dedup_files = []
    cmd = ''
    for bam in sortbam:
        dedup = bam[:-8] + 'dedup.bam'
        dedup_files.append(dedup)
        cmd = cmd + ('java -jar {mark} I={input} O={output} CREATE_INDEX=true '
        'METRICS_FILE=metrics.txt MAX_RECORDS_IN_RAM=2000000 '
        'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 & ').format(mark=mark,input=bam,
        output=dedup)
    subprocess.call(cmd[:-2],shell=True)
    return dedup_files


        