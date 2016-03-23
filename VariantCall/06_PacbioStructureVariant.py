from natsort import natsorted
import os
from Modules.f02_aligner_command import blasr
path = '/data/shangzhong/Pacbio/Analysis_Results'
ref_fa = '/data/moji/real_data/cho-k1-control/hamster_pb_illu.fa'
thread = 4

#========  1. list files ==================================
os.chdir(path)
try:
    fa_files = [f for f in os.listdir(path) if f.endswith('.fa') or f.endswith('.fasta')]
    fa_files = natsorted(fa_files)
    print 'list files succeed'
    print 'fa_files is:',fa_files
except:
    print 'list files failed'
    raise
#======== 2. align using blasr ============================
try:
    samFiles = blasr(fa_files,ref_fa,thread,otherParameters=['-clipping soft'])
    print 'align succeed'
    print 'samFiles',samFiles
except:
    print 'align failed'
    raise




