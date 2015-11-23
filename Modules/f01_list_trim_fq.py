from os import listdir
import subprocess
import copy
from natsort import natsorted

def chunk(l,n):
    n = max(1,n)
    result = [l[i:i+n] for i in range(0,len(l),n)]
    return result


def list_files_human(file_path):
    """
    This function lists all fastq files into a list for human samples,
    it is only for paired and strand specific pairs. R1 is for forward fastq,
    R2 for forward fastq.
    """
    allFiles = [f for f in listdir(file_path) if f.endswith(".fastq.gz") or f.endswith(".fq.gz")]
    allFiles = natsorted(allFiles)
    R1 = []; R2 = []
    for filename in allFiles:
        if 'R1' in filename:
            R1.append(filename)
        if 'R2' in filename:
            R2.append(filename)
    
    fastqFiles = []
    for r1,r2 in zip(R1,R2):
        fastqFiles.append([r1,r2])
    
    return fastqFiles


def list_files(file_path):
    """
    This function list all fastq files into a list
    """
    allFiles = [f for f in listdir(file_path) if f.endswith(".fastq.gz") or f.endswith(".fq.gz")]
    allFiles = natsorted(allFiles)
    fastqFiles = []  # this list is going to stroe the paired or single file for running aligner
    while len(allFiles) > 1:           # this is to append the single end or pair end files into a list.
        if allFiles[0].endswith(".fastq.gz"):
            index = allFiles[0].index(".fastq.gz")
            if allFiles[1][index-2:index] == '_2':
                fastqFiles.append(allFiles[:2])
                del allFiles[:2]
            else:
                fastqFiles.append(allFiles[:1])
                del allFiles[:1]
        
        if len(allFiles) != 0:        
            if allFiles[0].endswith(".fq.gz"):
                index = allFiles[0].index(".fq.gz")
                if allFiles[1][index-2:index] == '_2':
                    fastqFiles.append(allFiles[:2])
                    del allFiles[:2]
                else:
                    fastqFiles.append(allFiles[:1])
                    del allFiles[:1]
    if len(allFiles) == 1:
        fastqFiles.append(allFiles)

    return fastqFiles

def Trimmomatic(Trimmomatic,fastqFiles,phred ='33',adapter_file='',batch=1):  # batch means number of analysis run together in parallel in each batch
    """
    this function trims fastq files using Trimmomatic
    """
    trimmedFiles = copy.deepcopy(fastqFiles)
    for i in range(len(fastqFiles)):
            for j in range(len(fastqFiles[i])):
                trimmedFiles[i][j] = 'trim_' + fastqFiles[i][j]
    batch=min(batch,len(fastqFiles))
    subFqs = chunk(fastqFiles,batch)
    subTrims = chunk(trimmedFiles,batch)
    for Fqs,Trims in zip(subFqs,subTrims):
        cmd = ''
        for fastq, trimFastq in zip(Fqs,Trims):
            if len(fastq) ==2:
                trimCmd1st = ('java -jar {Trim} PE -threads {thread} -phred{type} {fastq1} {fastq2} '
                '{Trimmed1} unpair1.fq.gz {Trimmed2} unpair2.fq.gz ').format(Trim=Trimmomatic,
                    thread='3',type=phred,fastq1 = fastq[0], fastq2=fastq[1], 
                    Trimmed1 = trimFastq[0], Trimmed2 = trimFastq[1])
                
                if adapter_file != '':
                    adaptCmd = 'ILLUMINACLIP:{adapter}:2:30:10 '.format(adapter=adapter_file)
                else:
                    adaptCmd = ''
                trimCmd2nd = 'SLIDINGWINDOW:5:10 LEADING:15 TRAILING:10 MINLEN:36 TOPHRED33 & '
                cmd = cmd + trimCmd1st + adaptCmd + trimCmd2nd
            else:
                trimCmd1st = ('java -jar {Trim} SE -threads {thread} -phred{type} '
                              '{input} {output} ').format(Trim = Trimmomatic,thread = '3',
                            input = fastq[0],output=trimFastq[0],type=phred)
                if adapter_file != '':
                    adaptCmd = 'ILLUMINACLIP:{adapter}:2:30:10 '.format(adapter=adapter_file)
                else:
                    adaptCmd = ''
                trimCmd2nd = 'SLIDINGWINDOW:5:10 LEADING:5 TRAILING:3 MINLEN:22 TOPHRED33 & '
                cmd = cmd + trimCmd1st + adaptCmd + trimCmd2nd
        print cmd
        subprocess.call(cmd + 'wait',shell=True)
    return trimmedFiles