import os,subprocess,sys

def Trinity(files,thread):
    """
    This function runs Trinity function
    """
    if len(files) == 1:
        if files[0].endswith('fa.gz'):
            seqtype = 'fa'
        else:
            seqtype = 'fq'
        cmd = ('Trinity --seqType {seqtype} --max_memory 50G --single {f} --CPU {t}').format(seqtype=seqtype,f=files[0],t=thread)
    else:
        if files[0].endswith('fa.gz') or files[0].endswith('.fa'):
            seqtype = 'fa'
        else:
            seqtype = 'fq'
        cmd = ('Trinity --seqType {seqtype} --max_memory 50G --left {f1} --right {f2} --CPU {t}').format(
                                                        seqtype=seqtype,f1=files[0],f2=files[1],t=thread)
    subprocess.call(cmd.split(' '))