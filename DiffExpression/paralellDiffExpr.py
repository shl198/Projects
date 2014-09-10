import os
import subprocess
htseq_dir = '/data/shangzhong/Diff_express/Esko_project/htseqcount'
directory = '/data/shangzhong/Diff_express/Esko_project/Processed_data'
annotation = '/home/shangzhong/Database/mouse/mouse.gff3'
#========  get all folders ================
sub_folders = os.walk(directory).next()
folders = [f for f in sub_folders[1]]
folders.sort()
"""
folders.remove('Sample_mbmec_Ext1_Ext2_M09')
python_file = '/home/shangzhong/Scripts/Diff_Expression/Diff_Expression.py'
command =''
for syntax in folders:
    command = command + 'python {python} {path} & '.format(python=python_file,path=directory + 
                                          '/' + syntax)
command = command[:-2]
os.system(command)
"""


