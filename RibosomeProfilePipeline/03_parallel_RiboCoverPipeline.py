import os,subprocess
from natsort import natsorted

pythonFile = '/data/shangzhong/Codes/Pipeline/RibosomeProfilePipeline/03_RiboCoverPipeline.py'
path = '/data/shangzhong/RibosomeProfiling/Ribo_align'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
coverFiles = natsorted(coverFiles)
 
cmd = ''
for f in coverFiles:
    logFile = f[:-3]+'log.txt'
    inter = ('python {py} {file} > {log} & ').format(py=pythonFile,file=os.path.join(path,f),log=logFile)
    cmd = cmd + inter
subprocess.call(cmd[:-3],shell=True)