import subprocess
import sys
import os
from f03_samtools import merge_bam
def remove(files):
    """
    this function can remove files provided
    Arguments:  1. files: a list of files to be removed
    
    files: a list of files to be removed. [f1,f2,f3,...] or [[f1,f2],[f3],...] or with any depth of list layers
    """
    if isinstance(files,str):
        os.remove(files)
        try:
            os.remove(files + '.bai')
        except:
            pass
    if isinstance(files,list):
        for f in files:
            remove(f)
            try:
                remove(f+'.bai')
            except:
                continue


def get_parameters(parFile):
    """
    This function list all parameters for all pipelines.
    And return a dictionary
    
    parFiles: filename of the parmeter file.
    """
    res = open(parFile,'r')
    dic = {}
    for line in res:
        if line[0].isalpha(): # the first character is a letter
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
    if 'readGroup' in dic:
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
            renameCmd = ('cp {before} {after}').format(
                    before=readic[sample][0],after=output)
            subprocess.check_call(renameCmd,shell=True)
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
    subprocess.check_call(cmd[:-3],shell=True)



