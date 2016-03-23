import subprocess,os

def chunk(l,n):
    n = max(n,1)
    res = [l[i:i+n] for i in range(0,len(l),n)]
    return res


def read_group(ID,sample,platform,library,platunit):
    return ('@RG\\tID:'+ID+'\\tSM:'+sample+'\\tPL:'+platform+'\\tLB:'+library
            +'\\tPU:'+platunit)
    #@RG\\tID:chosgroup1\\tSM:sample1\\tPL:illumina\\tLB:lib1\\tPU:unit1
def sam2sortbam(picard,samfiles):
    """
    this function change samfile to sorted bam file
    """
    sorted_files = []
    for sam in samfiles:
        sort_bam = sam[:-3] + 'sort.bam'
        sorted_files.append(sort_bam)
        cmd = ('java -jar {picard} SortSam INPUT={input} OUTPUT={output} '
               'SORT_ORDER=coordinate').format(picard=picard,
                input=sam,output=sort_bam)
        subprocess.check_call(cmd,shell=True)
    return sorted_files

def markduplicates(picard,sortBams):
    """
    this function mark duplicates of the sorted bam file
    """
    if not os.path.exists('tmp'): 
        os.makedirs('tmp')
    dedup_files = []
    cmd = ''
    for bam in sortBams:
        dedup = bam[:-8] + 'dedup.bam'
        dedup_files.append(dedup)
        cmd = cmd + ('java -Djava.io.tmpdir=tmp -jar {picard} MarkDuplicates I={input} O={output} CREATE_INDEX=true '
        'METRICS_FILE=metrics.txt MAX_RECORDS_IN_RAM=8000000 '
        'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 '
        'VALIDATION_STRINGENCY=LENIENT && ').format(picard=picard,input=bam,
        output=dedup)
    subprocess.check_call(cmd[:-3],shell=True)
    return dedup_files

def addReadGroup(picard,sortBamFiles,readgroups,batch=1):
    """
    This function adds readgroup to a list of samfiles
    """
    batch = min(batch,len(sortBamFiles))
    subBams = chunk(sortBamFiles,batch)
    subGroups = chunk(readgroups,batch)
    sortBams = []
    for Bams,Groups in zip(subBams,subGroups):
        cmd = ''
        for sam,rg in zip(Bams,Groups):
            sortbam = sam[:-3] + 'adrg.bam'
            sortBams.append(sortbam)
            readgroup = rg.split('\\t')
            ID = readgroup[1][3:-1]
            SM = readgroup[2][3:-1]
            PL = 'illumina' #readgroup[3][3:-1]
            LB = 'lib20000' #readgroup[4][3:-1]
            PU = 'unit1' #readgroup[5][3:]
            cmd = cmd + ('java -jar {picard} AddOrReplaceReadGroups I={input} O={sortbam} SO=coordinate '
                         'RGID={ID} RGSM={SM} RGPL={PL} RGLB={LB} RGPU={PU} & ').format(
                        picard=picard,input=sam,sortbam=sortbam,ID=ID,SM=SM,PL=PL,LB=LB,
                        PU=PU)
        print cmd
        subprocess.call(cmd + 'wait',shell=True)
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
    This function transfer sam/bam to fastq files
    For paired end data, return [['f1.fq.gz','f2.fq.gz'],...]
    for single end data, return [['f1.fq.gz'],...]
    
    * samFiles is a list of sam/bam files
    * Type: 'single' or 'pair'
    """
    fqs = []
    cmd = ''
    if endType == 'pair':
        for sam in samFiles:
            fq1 = sam[:-4] + '_1.fq.gz'
            fq2 = sam[:-4] + '_2.fq.gz'
            fqs.append([fq1,fq2])
            sam2fqCmd = ('java -jar {picard} SamToFastq I={input} F={fq1} F2={fq2} '
                         'VALIDATION_STRINGENCY=LENIENT').format(
                        picard=picard,input=sam,fq1=fq1,fq2=fq2)
            cmd = cmd + sam2fqCmd + ' & '
    else:
        for sam in samFiles:
            fq = sam[:-4] + '.fq.gz'
            fqs.append([fq])
            sam2fqCmd = ('java -jar {picard} SamToFastq I={input} F={fq} '
                         'VALIDATION_STRINGENCY=LENIENT').format(
                        picard=picard,input=sam,fq=fq)
            cmd = cmd + sam2fqCmd + ' & '
    subprocess.call(cmd + 'wait',shell=True)
    return fqs


def sortVCF(picard,vcfFiles,fa_dict,batch=1):
    """This function reorders chromosome in vcf, making it the same as reference
    * picard
    * vcfFile: str.vcf file name
    * fa_dict: reference.dict that GATK used
    """
    VCFs = chunk(vcfFiles,batch)
    outFiles = []
    for vcfs in VCFs:
        cmd = ''
        outFiles = []
        for vcf in vcfs:
            outVCF = vcf[:-3] + 'sort.vcf'
            outFiles.append(outVCF)
            cmd = cmd + ('java -jar {picard} SortVcf I={input} O={output} SEQUENCE_DICTIONARY={fa_dict} & ').format(
                            picard=picard,input=vcf,output=outVCF,fa_dict=fa_dict)
        print cmd
        subprocess.call(cmd + 'wait',shell=True)
        # remove the original files
        for f_in,f_out in zip(vcfs,outFiles):
            os.remove(f_in)
            os.rename(f_out,f_in)
