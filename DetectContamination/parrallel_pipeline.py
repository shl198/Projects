import os
import sys
import subprocess
path = '/data/RNAseq/seq'
sub_folders = os.walk(path).next()
folders = [f for f in sub_folders[1]]
folders.sort()

python_file = '/home/shangzhong/Scripts/DetectContamination/Detect_Contamination.py'
command =''
for syntax in folders:
    command = command + 'python {python} {path} & '.format(python=python_file,path=path + 
                                          '/' + syntax)
command = command[:-2]
subprocess.call(command,shell=True)
try:
    sys.stdout.close()
except:
    pass
try:
    sys.stderr.close()
except:
    pass
print 'done'