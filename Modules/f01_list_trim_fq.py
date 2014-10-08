from os import listdir
import subprocess
import copy
def list_fastq(file_path,Trim):
    """
    This function can put all the fastq files in a folder into a list
    """

    Trimmomatic = '/home/shangzhong/Installation/Trimmomatic-0.32/trimmomatic-0.32.jar'
    adapter_file = '/home/shangzhong/Installation/Trimmomatic-0.32/adapters/TruSeq3-SE.fa'
    allFiles = [f for f in listdir(file_path) if f.endswith(".fastq.gz") or f.endswith(".fq.gz")]
    allFiles.sort()
    fastqFiles = []  # this list is going to stroe the paired or single file for running aligner
    while len(allFiles) > 1:           # this is to append the single end or pair end files into a list.
        if allFiles[0].endswith(".fastq.gz"):
            index = allFiles[0].index(".fastq.gz")
            if allFiles[1][index - 1] == '2':
                fastqFiles.append(allFiles[:2])
                del allFiles[:2]
            else:
                fastqFiles.append(allFiles[:1])
                del allFiles[:1]
        
        if len(allFiles) != 0:        
            if allFiles[0].endswith(".fq.gz"):
                index = allFiles[0].index(".fq.gz")
                if allFiles[1][index - 1] == '2':
                    fastqFiles.append(allFiles[:2])
                    del allFiles[:2]
                else:
                    fastqFiles.append(allFiles[:1])
                    del allFiles[:1]
    if len(allFiles) == 1:
        fastqFiles.append(allFiles)
        
    if Trim == 'False':
        trimmedFiles = fastqFiles
    if Trim == 'True':
        #-----   begin trimming -------
        trimmedFiles = copy.deepcopy(fastqFiles)
        for i in range(len(fastqFiles)):
            for j in range(len(fastqFiles[i])):
                trimmedFiles[i][j] = 'trim_' + fastqFiles[i][j]
        for fastq, trimFastq in zip(fastqFiles, trimmedFiles):
            if len(fastq) ==2:
                TrimCommand = ('java -jar {Trim} PE -threads {thread} -phred33 {fastq1} {fastq2} '
                '{Trimmed1} unpair1.fq.gz {Trimmed2} unpair2.fq.gz SLIDINGWINDOW:5:10 LEADING:15 '
                'TRAILING:15 MINLEN:36').format(Trim = Trimmomatic,thread = '5', 
                fastq1 = fastq[0], fastq2=fastq[1], Trimmed1 = trimFastq[0], Trimmed2 = trimFastq[1])
                subprocess.call(TrimCommand,shell=True)
            else:
                TrimCommand = ('java -jar {Trim} SE -threads {thread} -phred33 {input} {output} '
                'ILLUMINACLIP:{adapter}:2:30:10 SLIDINGWINDOW:5:10 LEADING:2 TRAILING:2 MINLEN:36').format(
                Trim = Trimmomatic,thread = '5',input = fastq[0],output =trimFastq[0],
                adapter = adapter_file)
                subprocess.call(TrimCommand,shell=True)
    return trimmedFiles

def list_files(file_path):
    """
    This function list all fastq files into a list
    """
    allFiles = [f for f in listdir(file_path) if f.endswith(".fastq.gz") or f.endswith(".fq.gz")]
    allFiles.sort()
    fastqFiles = []  # this list is going to stroe the paired or single file for running aligner
    while len(allFiles) > 1:           # this is to append the single end or pair end files into a list.
        if allFiles[0].endswith(".fastq.gz"):
            index = allFiles[0].index(".fastq.gz")
            if allFiles[1][index - 1] == '2':
                fastqFiles.append(allFiles[:2])
                del allFiles[:2]
            else:
                fastqFiles.append(allFiles[:1])
                del allFiles[:1]
        
        if len(allFiles) != 0:        
            if allFiles[0].endswith(".fq.gz"):
                index = allFiles[0].index(".fq.gz")
                if allFiles[1][index - 1] == '2':
                    fastqFiles.append(allFiles[:2])
                    del allFiles[:2]
                else:
                    fastqFiles.append(allFiles[:1])
                    del allFiles[:1]
    if len(allFiles) == 1:
        fastqFiles.append(allFiles)
    
    return fastqFiles

def Trimmomatic(fastqFiles,phred = '33'):
    """
    this function trims fastq files using Trimmomatic
    """
    Trimmomatic = '/home/shangzhong/Installation/Trimmomatic-0.32/trimmomatic-0.32.jar'
    adapter_file = '/home/shangzhong/Installation/Trimmomatic-0.32/adapters/TruSeq3-SE.fa'
    trimmedFiles = copy.deepcopy(fastqFiles)
    for i in range(len(fastqFiles)):
            for j in range(len(fastqFiles[i])):
                trimmedFiles[i][j] = 'trim_' + fastqFiles[i][j]
    
    cmd = ''
    for fastq, trimFastq in zip(fastqFiles, trimmedFiles):
        if len(fastq) ==2:
            cmd = cmd + ('java -jar {Trim} PE -threads {thread} -phred{type} {fastq1} {fastq2} '
            '{Trimmed1} unpair1.fq.gz {Trimmed2} unpair2.fq.gz SLIDINGWINDOW:5:10 LEADING:15 '
            'TRAILING:10 MINLEN:36 TOPHRED33 & ').format(Trim = Trimmomatic,thread = '5', type=phred,
            fastq1 = fastq[0], fastq2=fastq[1], Trimmed1 = trimFastq[0], Trimmed2 = trimFastq[1])
        else:
            cmd = cmd + ('java -jar {Trim} SE -threads {thread} -phred{type} {input} {output} '
            'ILLUMINACLIP:{adapter}:2:30:10 SLIDINGWINDOW:5:10 LEADING:2 TRAILING:2 '
            'MINLEN:36 TOPHRED33 & ').format(
            Trim = Trimmomatic,thread = '5',input = fastq[0],output =trimFastq[0],type=phred,
            adapter = adapter_file)
    subprocess.call(cmd[:-2],shell=True)
    return trimmedFiles