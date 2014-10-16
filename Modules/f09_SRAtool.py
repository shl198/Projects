import subprocess
from os import listdir
def list_sra(path):
    """
    this function list all sra files in a folder
    """
    allFiles = [f for f in listdir(path) if f.endswith('.sra')]
    allFiles.sort()
    return allFiles

def sra2fastq(sraFiles):
    """
    This file transfer sra files to fastq files
    """
    cmd = ''
    for sra in sraFiles:
        sraCmd = ('fastq-dump --split-files --gzip {sra}').format(sra=sra)
        cmd = cmd + sraCmd + ' & '
    subprocess.call(cmd[:-3],shell=True)