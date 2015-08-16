import subprocess,os

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
        subprocess.check_call(cmd,shell=True)
    return sorted_files

def markduplicates(picard,sortBams):
    """
    this function mark duplicates of the sorted bam file
    """
    mark = picard + '/' + 'MarkDuplicates.jar'
    dedup_files = []
    cmd = ''
    for bam in sortBams:
        dedup = bam[:-8] + 'dedup.bam'
        dedup_files.append(dedup)
        cmd = cmd + ('java -jar {mark} I={input} O={output} CREATE_INDEX=true '
        'METRICS_FILE=metrics.txt MAX_RECORDS_IN_RAM=8000000 '
        'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 '
        'VALIDATION_STRINGENCY=LENIENT && ').format(mark=mark,input=bam,
        output=dedup)
    subprocess.check_call(cmd[:-3],shell=True)
    return dedup_files

def addReadGroup(picard,sortBamFiles,readgroups):
    """
    This function add readgroup to a list of samfiles
    """
    add = picard + '/' + 'AddOrReplaceReadGroups.jar'
    sortBams = []
    cmd = ''
    for sam,rg in zip(sortBamFiles,readgroups):
        sortbam = sam[:-3] + 'adrg.bam'
        sortBams.append(sortbam)
        readgroup = rg.split('\\t')
        ID = readgroup[1][3:-1]
        SM = readgroup[2][3:-1]
        PL = readgroup[3][3:-1]
        LB = readgroup[4][3:-1]
        PU = readgroup[5][3:]
        cmd = cmd + ('java -jar {addGp} I={input} O={sortbam} '
                     'RGID={ID} RGSM={SM} RGPL={PL} RGLB={LB} RGPU={PU} & ').format(
                    addGp=add,input=sam,sortbam=sortbam,ID=ID,SM=SM,PL=PL,LB=LB,
                    PU=PU)
    subprocess.check_call(cmd + 'wait',shell=True)
#     # the file name in sortBams is filename.sort.sort.bam, need to change to filename.sort.bam
#     final_sort_bams = []
#     for bam in sortBams:
#         finalSortBam = bam[:-13] + 'sort.bam'
#         final_sort_bams.append(finalSortBam)
#         os.remove(finalSortBam)
#         renameCmd = ('mv {sortBam} {finalSortBam}').format(sortBam=bam,finalSortBam=finalSortBam)
#         subprocess.check_call(renameCmd,shell=True)
#     
    return sortBams

def sam2fastq(picard,samFiles,endType):
    """
    * samFiles is a list of sam/bam files
    * Type: 'single' or 'pair'
    """
    fqs = []
    cmd = ''
    sam2fq = picard + '/' + 'SamToFastq.jar'
    if endType == 'pair':
        for sam in samFiles:
            fq1 = sam[:-4] + '_1.fq.gz'
            fq2 = sam[:-4] + '_2.fq.gz'
            fqs.append([fq1,fq2])
            sam2fqCmd = ('java -jar {sam2fq} I={input} F={fq1} F2={fq2} '
                         'VALIDATION_STRINGENCY=LENIENT').format(
                        sam2fq=sam2fq,input=sam,fq1=fq1,fq2=fq2)
            cmd = cmd + sam2fqCmd + ' & '
    else:
        for sam in samFiles:
            fq = sam[:-4] + '.fq.gz'
            fqs.append(fq)
            sam2fqCmd = ('java -jar {sam2fq} I={intput} F={fq} '
                         'VALIDATION_STRINGENCY=LENIENT').format(
                        sam2fq=sam2fq,input=sam,fq=fq)
            cmd = cmd + sam2fqCmd + ' & '
    subprocess.check_call(cmd + 'wait',shell=True)
    return fqs

