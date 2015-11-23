import pysam
import os,subprocess,sys
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
import matplotlib
matplotlib.use('Agg')
from multiprocessing import Process,Pool
import scipy.stats as sp_stats
from natsort import natsorted
from f02_RiboDataModule import *
from Modules.p03_ParseSam import bam_parse
import shutil
import matplotlib.pyplot as plt
from Modules.p05_ParseGff import extractAllPr
import HTSeq as ht
#---------- parameters ---------------------------------------
# paths and files
bam_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/bam'
gffFile = '/data/shangzhong/RibosomeProfiling/Database/combined.gff'
db_path = '/data/shangzhong/RibosomeProfiling/Database'

fwd_rev_path = bam_path + '/01_cov'
tss_tse_cov_path = bam_path + '/02_TSS_TSE_cov'
pr_pos_cov_path = bam_path + '/03_pr_pos_cov_path'
# other parameters
up = 300; down = 300


os.chdir(bam_path)
bamFiles = [f for f in os.listdir(bam_path) if f.endswith('.sort.bam')]
bamFiles = natsorted(bamFiles)
#==============================================================================
#                     1. get the 01_pr_cds.txt
#==============================================================================
# if not os.path.exists(db_path): os.mkdir(db_path)
# outFile = db_path + '/01_pr_rna.txt'
# extractAllPr(gffFile,feature='exon',out=outFile)                                # 01_rna_exon.txt
# df = pd.read_csv(outFile,sep='\t',header=0)
# df = df.sort(['Chr','ex_start'])
# df.to_csv(outFile,sep='\t',index=False)
# outFile = db_path + '/01_pr_cds.txt'
# extractAllPr(gffFile,feature='CDS',out=outFile)
# df = pd.read_csv(outFile,sep='\t',header=0)
# df = df.sort(['Chr','cds_start'])
# df.to_csv(outFile,sep='\t',index=False)
#===============================================================================
#                     2. distribution of alignment length (around 5 minutes)
#===============================================================================
def align_len_dist(bamFile):
    bamHandle = pysam.AlignmentFile(bamFile,'rb')
    Handle = bam_parse(bamHandle)
    len_dist = Handle.align_len_distribution()
    return len_dist

# # align_len_dist(bamFiles[0])
# #------------- 1) put length distribution into a dataframe -------------
# p = Pool(processes=4)
# len_dist_list = [p.apply(align_len_dist,args=(bam,)) for bam in bamFiles]
# dfs = []
# for dic,bam in zip(len_dist_list,bamFiles):
#     name = bam[:-4]
#     df = pd.DataFrame.from_dict(dic,orient='index')
#     df.columns = [name]
#     dfs.append(df)
# len_dist_df = pd.concat(dfs,axis=1)
# #------------- 2) plot the distribution ------------------------
# ax = len_dist_df.plot(kind = 'line')
# ax.set_xlabel('alignment read length')
# ax.set_ylabel('percentage of total mapped reads %')
# ax.set_title('distribution of alignment read length')
#  
# figure_path = bam_path + '/figures'
# if not os.path.exists(figure_path): os.mkdir(figure_path)
# figure = '01_alignment_length_distribution'
# plt.savefig(figure_path + '/' + figure +'.png')
# plt.savefig(figure_path + '/' + figure +'.svg')
# # decide the length to include
# len_dist_df['mean'] = len_dist_df.mean(axis=1)
# len_dist_df['cumu'] = len_dist_df['mean'].cumsum()
# max_len = len_dist_df[len_dist_df['cumu']<0.996].index[-1]

#===============================================================================
#                     3. get forward and reverse count
#===============================================================================

#fwd_rev_cov(bamFiles[0])

# proc = [Process(target=fwd_rev_cov,args=(bam,)) for bam in bamFiles]              # in folder 01_cov
# for p in proc:
#     p.start()
# for p in proc:
#     p.join()
# # move file into 01_cov
# fwd_rev_path = bam_path + '/01_cov'
# if not os.path.exists(fwd_rev_path): os.mkdir(fwd_rev_path)
# for f in os.listdir(bam_path):
#     if f.endswith('.txt'):
#         if os.path.exists(fwd_rev_path + '/' + f):
#             os.remove(fwd_rev_path + '/' + f)
#         shutil.move(f,fwd_rev_path)

#===============================================================================
#                     3. coverage around TSS TSE  
#===============================================================================
def TSS_TSE_coverage(covFile,cdsFile,outpath,up=50,down=50,align_len=0):
    """
    This function calculates coverage around TSS and TSE.
    * covFile: filename of cover df. columns are: ['num','chr','end5','end3','strand','len']
    * cdsFile: 01_pr_cds.txt. columns are:['chr','start','end','geneid','praccess','strand']
    """
    # 1). read coverage data
    cov_df = pd.read_csv(covFile,sep='\t',header=None,names=['num','chr','end5','end3','strand','len'])
    if align_len !=0:
        cov_df = cov_df[cov_df['len'].values==align_len]
    # 2). read protein position file
    cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)   # Chr    cds_start    cds_end    GeneID    Pr_Access    Strand
    trpr_obj = trpr(cds_df)  # read the class
    # 3). create array set by htseq. This is for retriving positions that are covered by many genes
    gene_array_set = trpr_obj.htseq_genome_array_set(False)
    # 3). get all genes
    genes = trpr_obj.all_gene_ids()
    genes = natsorted(genes)
    # 4). get coverage for each gene
    tss_cov = pd.DataFrame()
    tse_cov = pd.DataFrame()
    for g in genes:
        # get coverage
        try:
            gene_df = cds_df[cds_df['geneid'].values==g]
            g_obj = trpr(gene_df)
            long_pr = g_obj.get_longest_trpr(g)  # longest protein of the gene
            pos = g_obj.get_trpr_pos(long_pr,up,down) # position of the protein
            cov_pos = posCoverage(gene_df,cov_df,pos,gene_array_set)
        except:
            print g,'maps to multiple scaffold'
            continue
        tss_cov[g] = pd.Series(cov_pos[:(up+down)])
        tse_cov[g] = pd.Series(cov_pos[-(up+down):])
    if not os.path.exists(outpath):os.mkdir(outpath)
    tss_cov.T.to_csv(outpath+'/'+covFile.split('.')[0]+'_tss.txt',sep='\t')
    tse_cov.T.to_csv(outpath+'/'+covFile.split('.')[0]+'_tse.txt',sep='\t')
    
    
def groupCovForPosRep(files,totalCounts,up,down,calType='total'):
    """
    This function groups the replicates of samples' coverage around TSS and TSE sites.
    return two columns: ['mean'/'total'/'median'/'geometricmean','std']
    
    * files: list. A list of replicate files.
    * totalCounts: list. A list of int with same length of files. Stores the total mapped reads in files.
    * calType: str. Calculation type.
    return a dataframe with two columns. ['mean_value','standard deviation']
    """
    res_df = pd.DataFrame()
    mean = pd.DataFrame()
    std = pd.DataFrame()
    for f,total in zip(files,totalCounts):
        df = pd.read_csv(f,sep='\t',header=0,index_col=0,low_memory=False)
        # remove the antibody, signal peptide and non signal peptide are separate
        try:
            df = df.drop('heavychain');df=df.drop('lightchain');df=df.drop('NeoRKanR')
            #df = df.loc[['heavychain','lightchain']]
        except:
            pass
        # filter by median of gene
        df['median'] = df.median(axis=1)
        #df = df[df['median']>=0.5]
        del df['median']
        up_df = df.iloc[:,0:up]
        down_df = df.iloc[:,up:up+down]
        try:
            up_df = up_df.replace('-',np.nan).dropna().astype('int').T
            down_df = down_df.replace('-',np.nan).dropna().astype('int').T
        except:
            pass
        up_df = up_df/float(total)*(10**6)   # row: position. col: gene
        down_df = down_df/float(total)*(10**6)
        if calType == 'total':
            up_df[calType] = up_df.sum(axis=1);down_df[calType] = down_df.sum(axis=1)
        if calType == 'mean':
            up_df[calType] = up_df.mean(axis=1);down_df[calType] = down_df.mean(axis=1)
        if calType == 'median':
            up_df[calType] = up_df.median(axis=1);down_df[calType] = down_df.median(axis=1)
        if calType == 'geoMean':
            up_df = up_df.replace([0],[1]);down_df=down_df.replace([0],[1])
            up_df[calType] = sp_stats.gmean(up_df.values,axis=1)
            down_df[calType] = sp_stats.gmean(down_df.values,axis=1)
        df = df.T
        df[calType]= pd.concat([up_df[calType],down_df[calType]])  # add caltype column
        up_df['std']=up_df.std(axis=1);down_df['std']=down_df.std(axis=1)
        df['std'] = pd.concat([up_df['std'],down_df['std']])
        mean[f+calType] = df[calType]  # stores the results
        std[f+'std'] = df['std']
    res_df['mean'] = mean.mean(axis=1)
    res_df['std'] = ((std**2).sum(axis=1)/std.shape[1]).apply(np.sqrt)
    return res_df

def align_lengths(covFiles):
    """This function get all the align ment lengths
    * covFiles: list. A list of covFiles.
    """
    lens=[]
    for f in covFiles:
        df = pd.read_csv(f,sep='\t',header=None,names=['num','chr','start','end','strand','len'])
        length = list(set(df['len'].tolist()))
        lens = lens + length
    return natsorted(list(set(lens)))

def plot_tsse(df,xlabel,ylabel,title,file_len,up,down,savfig):
    """ This function plot the pandas dataframe around TSS TSE. 
    
    * df: dataframe. rows are each position around TSS TSE, columns are sample names.
    * file_len: int. indicates how many subplot the figure contains.
    * up,down: int. indicates how many nts of upstream and downstream of TSS TSE to include.
    * savfig: str. Filename to sotre the output figure.
    """
    f,ax=plt.subplots(file_len,sharex=True)
    f.set_size_inches(14.5,8)
    x = np.array(range(-up,down))
    for i in range(file_len):
        ax[i].bar(x,df.iloc[:,i],align='center',alpha=0.89)
        ax[i].set_xlim([-up,down])
        ax[i].set_title(df.columns.values[i])
    ax[file_len-1].set_xlabel(xlabel,color='black',fontsize=18)
    f.text(0.07, 0.5, ylabel, ha='center', va='center', rotation='vertical',fontsize=18)
    plt.xticks(x,range(-up,down),rotation=90)
    plt.suptitle(title,fontsize=26,color='black')
    plt.savefig(savfig)
# #====================== 1. get coverage (results in folder 02_TSS_TSE_cov)===================
os.chdir(fwd_rev_path)
covFiles = [f for f in os.listdir(fwd_rev_path) if f.endswith('cov.txt')]
covFiles = natsorted(covFiles)
cdsFile = db_path + '/01_pr_cds.txt'
#TSS_TSE_coverage(covFiles[0],cdsFile,tss_tse_cov_path,up,down,20)
lens = align_lengths(covFiles) # get all lengths
lens.append(0)
for l in lens:
    tsse_path = tss_tse_cov_path + '/' + str(l)
    proc = [Process(target=TSS_TSE_coverage,args=(covFile,cdsFile,tsse_path,up,down,l,)) for covFile in covFiles]
    for p in proc:
        p.start()
    for p in proc:
        p.join()
#======================== 2. plot the coverage (results in folder figures)===============    

#-------------- 1) get total count --------------------
os.chdir(bam_path)
calType = 'total'
ylabel = 'total count (RPM)'
totalCount = []
for bam in bamFiles:
    handle = pysam.AlignmentFile(bam,'rb')
    totalCount.append(handle.mapped)

for l in lens:
    #-------------- 2) merge replicates ------------------- 
    tsse_path = tss_tse_cov_path + '/' + str(l)
    if not os.path.exists(tsse_path):os.mkdir(tsse_path)
    os.chdir(tsse_path)
    tss_files = [f for f in os.listdir(tsse_path) if f.endswith('tss.txt')]
    tss_files = natsorted(tss_files)
    tse_files = [f for f in os.listdir(tsse_path) if f.endswith('tse.txt')]
    tse_files = natsorted(tse_files)
    day3_TSS = groupCovForPosRep(tss_files[0:3],totalCount[0:3],up,down,calType)
    day6_TSS = groupCovForPosRep(tss_files[3:6],totalCount[3:6],up,down,calType)
    day3_TSE = groupCovForPosRep(tse_files[0:3],totalCount[0:3],up,down,calType)
    day6_TSE = groupCovForPosRep(tse_files[3:6],totalCount[3:6],up,down,calType)
    # mean
    tsse = pd.DataFrame();tsse_err = pd.DataFrame()
    tsse['day3_start'] = day3_TSS['mean'];tsse['day6_start']=day6_TSS['mean']
    tsse['day3_end'] = day3_TSE['mean'];tsse['day6_end']=day6_TSE['mean']
    # error
    tsse_err['day3_start']=day3_TSS['std'];tsse_err['day6_start']=day6_TSS['std']
    tsse_err['day3_end']=day3_TSS['std'];tsse_err['day6_end']=day6_TSS['std']
    #-------------- 3) plot --------------------------------
    x = np.array(range(-up,down))
    f, ax = plt.subplots(4, sharex=True)
    f.set_size_inches(14.5,8)
    colors = ['#FF4500','#008B8B','#6495ED','#808080']
    for i in range(4):
        ax[i].bar(x,tsse.iloc[:,i],color=colors[i],align='center',alpha=0.89)#,yerr=[tuple([0]*(TSS.shape[0])),tuple(TSS_err.iloc[:,i].values)])
        #ax[i].set_ylim([0,200])
        ax[i].set_xlim([-up,down])
        ax[i].set_title(tsse.columns.values[i])
    ax[3].set_xlabel('distance from TSS or TSE',color='black',fontsize=18)
    f.text(0.07, 0.5, ylabel, ha='center', va='center', rotation='vertical',fontsize=18)
    plt.xticks(x,range(-up,down),rotation=90)
    if l == 0:
        title = 'TSS TSE coverage'
    else:
        title = 'TSS TSE coverage at length {l}'.format(l=l)
    plt.suptitle(title,fontsize=26,color='black')
    plt.savefig(bam_path+'/figures/'+str(l)+'.png')
    plt.savefig(bam_path+'/figures/'+str(l)+'.svg')

#======================== 3. plot the coverage for each bam file (results in folder ) ===============

#-------------- 1) get total count --------------------
os.chdir(bam_path)
calType = 'total'
ylabel = 'total count (RPM)'
totalCount = []
for bam in bamFiles:
    handle = pysam.AlignmentFile(bam,'rb')
    totalCount.append(handle.mapped)
  
for l in lens:
    tsse_path = tss_tse_cov_path + '/' + str(l)
    os.chdir(tsse_path)
    tss_files = [f for f in os.listdir(tsse_path) if f.endswith('tss.txt')]
    tss_files = natsorted(tss_files)
    tse_files = [f for f in os.listdir(tsse_path) if f.endswith('tse.txt')]
    file_len = len(tss_files)
    tse_files = natsorted(tse_files)
    tss_df = pd.DataFrame()
    tse_df = pd.DataFrame()
    for f1,f2,total in zip(tss_files,tse_files,totalCount):
        df1 = pd.read_csv(f1,sep='\t',header=0,index_col=0,low_memory=False)
        df1 = df1.drop('heavychain');df1=df1.drop('lightchain');df1=df1.drop('NeoRKanR')
        df2 = pd.read_csv(f2,sep='\t',header=0,index_col=0,low_memory=False)
        df2 = df2.drop('heavychain');df2=df2.drop('lightchain');df2=df2.drop('NeoRKanR')
        tss_df[f1.split('.')[0]] = df1.T.sum(axis=1)/float(total)*(10**6)
        tse_df[f2.split('.')[0]] = df2.T.sum(axis=1)/float(total)*(10**6)
    plot_tsse(tss_df,'distance from TSS',ylabel,'TSS coverage at lenght {len}'.format(len=str(l)),file_len,up,down,tss_tse_cov_path+'/'+str(l)+'_tss.png')
    plot_tsse(tse_df,'distance from TSE',ylabel,'TSE coverage at lenght {len}'.format(len=str(l)),file_len,up,down,tss_tse_cov_path+'/'+str(l)+'_tse.png')
    print 'done'










