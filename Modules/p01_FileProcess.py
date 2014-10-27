import subprocess
from f03_samtools import merge_bam
def remove(files):
    """
    this function can remove files provided
    Arguments:  1. files: a list of files to be removed
    
    files: a list of files to be removed. [f1,f2,f3,...]
    """
    if isinstance(files,str):
        subprocess.call('rm {file}'.format(file=files),shell=True)
    if isinstance(files,list):
        cmd = ''
        for f in files:
            cmd = cmd + ('rm {file} & ').format(file=f)
        subprocess.call(cmd[:-3],shell=True)
    

def get_parameters(parFile):
    """
    This function list all parameters for all pipelines.
    And return a dictionary
    
    parFiles: filename of the parmeter file.
    """
    res = open(parFile)
    dic = {}
    for line in res:
        if line[0].isalpha():
            item = line.split('\t')
            value = item[1] # combine all values into a single string
            if ',' in value:
                rg = value.split(',')
                rg[-1] = rg[-1][:-1]
                dic[item[0]] = rg
            else:
                dic[item[0]] = item[1][:-1]
        else:
            continue
    if isinstance(dic['readGroup'],str):
                dic['readGroup'] = [dic['readGroup'][:-1]]
    return dic

def rg_bams(rgs,bamfiles):
    """
    this function merge bam files by read group names
    bam files belong to same sample will be merged into one bam file
    
    rgs: a list of read group names. [rg1,rg2,rg3,...]
    
    bamfiles: a list of bam files. [f1.bam,f2.bam,f3.bam,...]
    """
    readic = {}
    # build a dic, key is sample name, value is bam file
    for rg,bam in zip(rgs,bamfiles):
        start = rg.index('SM:')
        end = rg.index('PL:')
        sample = rg[start+3:end-3]
        if sample in readic:
            readic[sample].append(bam)
        else:
            readic[sample] = [bam]
    merged = []
    for sample in readic:
        output = sample + '.merged.sort.bam'
        merged.append(output)
        if len(readic[sample]) == 1:
            renameCmd = ('mv {before} {after}').format(
                    before=readic[sample][0],after=output)
            subprocess.call(renameCmd,shell=True)
        else:
            merge_bam(readic[sample],output)
    return merged

def changeFastqReadName(fastqFiles):
    """
    this function changes fastq file names to numbers in order to make sure fastq
    files be aligned without error
    
    fastqFiles: a list of fastq files. [f1.fq.gz,f2.fq.gz,f3.fq.gz,...]
    """
    cmd = ''
    for fq in fastqFiles:
        output = fq[0][5:]
        renameCmd = ('gunzip -c {fq} | fastx_renamer -n COUNT '
                     '-z -o {output}').format(fq=fq[0],output=output)
        cmd = cmd + renameCmd + ' & '
    subprocess.call(cmd[:-3],shell=True)


def extract_gene2ref(filename,taxID,columnNums = []):
    """
    This file extracts rows by taxonomy id, and columns. 
    
    filename: 'gene2ref'
    
    taxID: a list of IDs. eg: [id1,id2,...]
    
    columnNums: a list of columns to extract
    """
    # define output file
    taxID = map(str,taxID)
    outputfile = filename + '_sub.txt'
    result = open(filename,'r')
    output = open(outputfile,'w')
    # process each line
    for line in result:
        line = line[:-1]
        if line[0] == '#':
            continue 
        item = line.split('\t')
        # extract by tax ID
        if item[0] in taxID:
            if columnNums == []:
                output.write(line)
            else:
                row = []
                for index in columnNums:
                    row.append(item[index-1])
                output.write('\t'.join(row) + '\n')
    result.close()
    output.close()

def extract_from_gene2ref(filename,taxID,columnNum=[]):
    """
    This file extracts rows by taxonomy id, and columns. 
    
    filename: 'gene2ref'
    
    taxID: integer indicates taxonomy id of interested organism.
    
    columnNums: a list of columns to extract
    """
    taxID = str(taxID)
    columnNum = map(str,columnNum)
    outputfile = filename + '_sub.txt'
    if columnNum ==[]:
        cmd = ('awk -f {sep} {printout} {inputfile} | uniq > {output}').format(sep='\"\\t\"',
                printout='\'{print $0}\'',inputfile=filename,output=outputfile)
    else:
        printout = 'print '
        for index in columnNum:
            printout = printout + '$' + index + ' ' + '\"\\t\"' + ' '
        printout = printout[:-6]
        cmd = ('awk -F {sep} {parse} {printcmd} {inputfile} | uniq > {output}').format(sep='\"\\t\"',
                    printcmd='{' + printout + '}\'',inputfile=filename,output=outputfile,
                    parse =   '\'$1 == ' + taxID)
    subprocess.call(cmd,shell=True)

