import os,sarge
from natsort import natsorted
import shutil
import pandas as pd

def chunk(l,n):
    n = max(1,n)
    result = [l[i:i+n] for i in range(0,len(l),n)]
    return result

def copy_fq_merge_bio_replicate(lst_fn,target_path,index,datafolder='/media/lewislab/Dropbox (UCSD SBRG)/LewisPub'):
    '''
    This function automate copy of samples in lst_fn and merge the all fq files belong to the same biological replicate sample.
    * index: needs to be a list
    '''
    os.chdir(datafolder)
    smp_info_dic = {}
    with open(lst_fn) as f:
        for line in f:
            item = line.strip().split('\t')
            smp = item[0]
            if smp in smp_info_dic:
                smp_info_dic[smp].append(item[1:])
            else:
                smp_info_dic[smp] = [item[1:]]
    
    for idx in index:  # idx is sample id
        files = []
        for info in smp_info_dic[str(idx)]:
            cmd = 'find ./{path} -name {name}*.gz'.format(path=info[1],name=info[0])
            fns = sarge.get_stdout(cmd)
            files.extend(fns.strip().split('\n'))
        files = natsorted(files)
        fst = [f for f in files if ('_1.' in f) or ('_R1_' in f)]
        snd = [f for f in files if ('_2.' in f) or ('_R2_' in f)]
        cmd = ('cat {in_f} > {target}/{out_f}_1.fq.gz').format(in_f=' '.join(fst),target=target_path,out_f='sp_'+str(idx))
        print cmd
        sarge.run(cmd)
        cmd = ('cat {in_f} > {target}/{out_f}_2.fq.gz').format(in_f=' '.join(snd),target=target_path,out_f='sp_'+str(idx))
        sarge.run(cmd)

# lst_fn = '/home/shangzhong/Codes/Projects/cp_lst.txt'
# target_path = '/data/shangzhong/DE/pgsa/fq'
# copy_fq_merge_bio_replicate(lst_fn,target_path,[1001])

def run_cufflink_pipeline(batch,lst_fn,run_path,store_path,ppl_fn,param_fn,indexes='all'):
    '''
    * run_path: pathway to run the pipeline
    * store_path: pathway store all the results
    '''
    # find indexes samples
    if indexes == 'all':
        lst_df = pd.read_csv(lst_fn,sep='\t',header=None)
        indexes = list(set(lst_df['0'].tolist()))
    batch_index = chunk(indexes,batch)
    if not os.path.exists(store_path):
        os.mkdir(store_path)
    for sub_index in batch_index:
        if os.path.exists(run_path):
            shutil.rmtree(run_path)
        os.mkdir(run_path)
        copy_fq_merge_bio_replicate(lst_fn,run_path,sub_index)
        cmd = ('python {py} {param}').format(py=ppl_fn,param=param_fn)
        sarge.run(cmd)
        cuff = [os.path.join(run_path,fo) for fo in os.listdir(run_path) if fo.endswith('cufflinks')]
        for c in cuff:
            shutil.move(c,store_path)
        shutil.rmtree(run_path)

# indexes = range(218,240)
# batch = 6
# lst_fn = '/data/hooman/DNARepair/fullRNASeqToProcess_FinalBatch.txt'
# run_path = '/data/shangzhong/DE/DNA_repair_de/fq'
# store_path = '/data/shangzhong/DE/DNA_repair_de/count'
# ppl_fn = '/home/shangzhong/Codes/Projects/DiffExpression/02_txedo_pipeline.py'
# param_fn = '/data/shangzhong/DE/DNA_repair_de/02_txedo_Parameters.txt'
# run_cufflink_pipeline(batch,lst_fn,run_path,store_path,ppl_fn,param_fn,indexes)

def run_vcf_pipeline(indexes,batch,lst_fn,run_path,store_path,ppl_fn,param_fn):
    batch_index = chunk(indexes,batch)
    if not os.path.exists(store_path):
        os.mkdir(store_path)
    for sub_index in batch_index:
        if os.path.exists(run_path):
            shutil.rmtree(run_path)
            os.mkdir(run_path)
        