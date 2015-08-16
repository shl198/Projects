#===============================================================================
#  read and merge file for Autism fastq files
#===============================================================================
import pandas as pd
path = '/media/lewislab/Dropbox (UCSD SBRG)/LewisPub/Data/RNAseq/2015_03_Autism/58_samples'
import os
os.chdir(path)
filenames = [f for f in os.listdir(path) if f.endswith('fastq.gz')]
df = pd.DataFrame({'58_samples':filenames})

filename = '/data/shangzhong/DE/Autism/readin.txt'

names = pd.read_csv(filename,header=None,names=["name"])
import subprocess
for name in names['name']:
    name_fileobject = df[df['58_samples'].map(lambda x: x.endswith(name+'.fastq.gz'))]
    first_obj = name_fileobject[name_fileobject['58_samples'].map(lambda x: '_1_' in x)]
    snd_obj = name_fileobject[name_fileobject['58_samples'].map(lambda x: '_2_' in x)]
    
    fir_list = first_obj['58_samples'].tolist()
    snd_list = snd_obj['58_samples'].tolist()
    cmd = 'cat ' + ' '.join(fir_list) + ' > ' + '/data/shangzhong/DE/Autism/' + name + '_1.fq.gz'
    subprocess.call(cmd,shell=True)
    print cmd
    cmd = 'cat ' + ' '.join(snd_list) + ' > ' + '/data/shangzhong/DE/Autism/' + name + '_2.fq.gz'
    subprocess.call(cmd,shell=True)
    print cmd