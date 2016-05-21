# import pbcore.io as pbi 
# 
# pbFn = '/data/shangzhong/Pacbio/Analysis_Results/m151203_130739_42149_c100885822550000001823202004021691_s1_p0.1.bax.h5'
# handle = pbi.BasH5Reader(pbFn)
import sarge

def bash5tool(hdf5File,seqtype):
    cmd = ('bash5tools.py --outFilePrefix {h5} --readType subreads --outType {type} '
           '--minReadScore 0.8 {input}').format(h5=hdf5File[:-7],type=seqtype,input=hdf5File)
    sarge.run(cmd)

