import subprocess
def remove(files):
    """
    this function can remove files provided
    """
    if isinstance(files,str):
        subprocess.call('rm {file}'.format(file=files),shell=True)
    if isinstance(files,list):
        cmd = ''
        for f in files:
            cmd = cmd + 'rm {file} & '.format(file=f)
        subprocess.call(cmd[:-2],shell=True)
    
        