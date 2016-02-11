from __future__ import division
import subprocess,os,sys
import pandas as pd
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('ggplot')
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import pysam
import HTSeq as ht
from Modules.p03_ParseSam import bam_parse

def chunk(l,n):
    n = max(1,n)
    res = [l[i:i+n] for i in range(0,len(l),n)]
    return res

def signalP(inputFa):
    """
    This function runs signalp for files
    
    * inputFa: str. Input fasta file that stores protein sequence
    """
    outFile = inputFa[:-3] + '_sp.txt'
    cmd = ('signalp {input} > {out}').format(input=inputFa,out=outFile)
    subprocess.call(cmd,shell=True)
    return outFile
    
def removeIDVersion(ID):
    """
    This function remove version information in all kinds of id versions
    
    * ID: str. should be genome mRNA or protein accession number with version information. eg: NW0123.1
    """
    if '.' in ID:
        return ID[:ID.index('.')]
    else:
        return ID

def changeFileName(filename,n=1):
    """
    This function changes the index of file names
    
    * filename: str. File name
    * n: int. index increase
    """
    index = filename.rindex('/')
    name = filename[index+1:]
    newName = filename[:index+1] + str(int(name[:2]) + n).zfill(2) + name[2:]
    cmd = ('mv {old} {new}').format(old=filename,new=newName)
    subprocess.call(cmd,shell=True)
    return newName
    
    
def addGeneIDSymbolChr2signalP(signalP,mapFile,organism='',mapSource='gff'):
    """
    This function adds gene id in the first column in the signalP result
    
    * signalP: str. Filename of signalP result
    * mapFile: str. All information of cho extract from gene2refseq
    * organism: str. 'others' or 'chok1, or 'hamster'
    * mapSource: str. Default is 'gff'. Alter is 'gene2refseq' 
    """
    # read map data
    
    mapdf = pd.read_csv(mapFile,sep='\t',header=0,names=['GeneID','GeneSymbol','Chr','TrAccess','PrAccess'],low_memory=False)
    # extract either cho or hamster
    if mapSource=='gene2refseq':
        mapdf = pd.read_csv(mapFile,sep='\t',skiprows=[0],header=None,usecols=[1,3,5,7],names=['GeneID','TrAccess','PrAccess','Chr'],low_memory=False)
        if organism == 'chok1':
            criteria = mapdf['Chr'].map(lambda x: 'NW_0036' in x)
            mapdf = mapdf[criteria]
        if organism == 'hamster':
            criteria = mapdf['Chr'].map(lambda x: 'NW_0068' in x)
            mapdf = mapdf[criteria]
        
    mapdf = mapdf.drop_duplicates()
    # read signalP result
    signalP_df = pd.read_csv(signalP,header=None,skiprows=[0,1],delim_whitespace=True,names=['name','Cmax','posC','Ymax','posY','Smax','posS','Smean','D','?','Dmaxcut','Networks-used'],low_memory=False)
    # build dictionary
    gene_pr_dic = mapdf.set_index('PrAccess')['GeneID'].to_dict()
    mrna_pr_dic = mapdf.set_index('PrAccess')['TrAccess'].to_dict()
    try:
        del gene_pr_dic['-'],mrna_pr_dic['-']
    except:
        pass
    # add rna accessions
    rnas = []
    # add geneIDs
    geneIDs = []
    for pr in signalP_df['name']:
        if pr in gene_pr_dic:
            geneIDs.append(gene_pr_dic[pr])
        else:
            #geneIDs.append(mRNA_prID2geneIDRemote(pr))
            geneIDs.append('-')
        if pr in mrna_pr_dic:
            rnas.append(mrna_pr_dic[pr])
        else:
            #rnas.append(mRNA_prID2geneIDRemote(pr,target='rna'))
            rnas.append('-')
    signalP_df = signalP_df.rename(columns={'name':'Pr_Access'})
    signalP_df.insert(0,'RNA_Access',pd.Series(rnas))    
    signalP_df.insert(0,'GeneID',pd.Series(geneIDs))
    outFile = signalP[:-3] + 'gene.txt'
    signalP_df = signalP_df.sort(columns='GeneID')
    signalP_df.to_csv(outFile,index=False,sep='\t',float_format='%.3f')
    return outFile
# addGeneIDSymbolChr2signalP('/data/shangzhong/RibosomeProfiling/cho_pr/02_choPrRefseq_sp.txt',
#                         '/data/shangzhong/Database/gff_chok1_all_ID.txt')

#===============================================================================
#      Calculate how many genes map to different scaffolds
#===============================================================================
def genesMap2diffChrom(inputFile,gffMapFile):
    """
    This function finds how many genes in inputFile map to multiple chromosome/scaffolds
    
    * inputFile: str. Filename of input file. Should have a column named GeneID that 
                 corresponds to gene ids.
    * gffMapFile: str. Filename of id mappings extracted from annotation file. First 
                    3 columns should be ['GeneID','GeneSymbol','Chr']
    
    return a file with GeneID as 1st column,chromosomes at the rest columns
    """
    # get gene ids
    df = pd.read_csv(inputFile,header=0,sep='\t',low_memory=False)
    geneIDs = df['GeneID'].astype(str).tolist()
    geneIDs = list(set(geneIDs))
    # read gff file
    gff = pd.read_csv(gffMapFile,header=0,sep='\t',usecols=[0,2],names=['gene','chrom'])
    gff = gff.drop_duplicates()
    gff['gene'] = gff['gene'].astype(str)
    gff_dic = {k:list(v) for k,v in gff.groupby('gene')['chrom']}
    # get genelist
    output = inputFile[:-3] + 'multichr.txt'
    f = open(output,'w')
    for gene in geneIDs:
        if gene in gff_dic:
            if len(gff_dic[gene]) > 1:
                f.write(gene+'\t'+'\t'.join(gff_dic[gene])+'\n')
        else:
            print gene,'not in dict'
        
    f.close()
    return output
# inputFile = '/data/shangzhong/RibosomeProfiling/cho_pr/03_choPrRefseq_sp.gene.txt'
# gffMapFile = '/data/shangzhong/Database/gff_chok1_all_ID.txt'
# genesMap2diffChrom(inputFile,gffMapFile)

#===============================================================================
#             classify gene ids regarding signal peptide
#===============================================================================
def gene_sp_classify(sp_gene_res):
    """
    This function classfies which gene encodes protein only with signal peptide,
    which gene encodes protein only without signal peptide and which gene encods both.
    
    * sp_gene_res: str. Filename of signalP prediction results with gene id inserted with column name 'GeneID'.
    return a file with 3 columns: both (genes encodes sp protein and non sp protien), no_sp, with_sp
    """
    df = pd.read_csv(sp_gene_res,header=0,sep='\t')
    dic = {k:list(v) for k,v in df.groupby('GeneID')['?']}
    # categorize gene ids
    both = []; with_sp=[];no_sp=[]
    for key in dic:
        if ('Y' in dic[key]) and ('N' in dic[key]):
            both.append(key)
            with_sp.append('-')
            no_sp.append('-')
        elif ('Y' in dic[key]) and ('N' not in dic[key]):
            both.append('-')
            with_sp.append(key)
            no_sp.append('-')
        else:
            both.append('-')
            with_sp.append('-')
            no_sp.append(key)
    res = list(set(both)); res.remove('-')
    print 'both: ',len(res)
    res = list(set(with_sp)); res.remove('-')
    print 'with_sp: ',len(res)
    res = list(set(no_sp)); res.remove('-')
    print 'no_sp: ',len(res)
    df = pd.DataFrame({'both':both,'with_sp':with_sp,'no_sp':no_sp})
    outFile = sp_gene_res[:-3]+'classify.txt'
    df.to_csv(outFile,sep='\t',index=False)
    
    return outFile
# filename = '/data/shangzhong/RibosomeProfiling/cho_pr/03_choPrRefseq_sp.gene.txt'
# gene_sp_classify(filename)

#===============================================================================
#     generate bed file
#===============================================================================

def gff2Bed4Genes(gffFile,genes,feature='exon',outBed=''):
    """
    This function converts gff file to bed file for a list of genes in the input.
    Only extract exons or CDS of genes.
    
    * gffFile: str. Filename of gff file
    * genes: list. A list of gene IDs.
    * feature: str. exon or CDS
    """
    gff_handle = open(gffFile,'r')
    if outBed =='':
        outBed = gffFile[:-4] + '.bed'
    out_handle = open(outBed,'w')
    for line in gff_handle:
        if line.startswith('#'):
            continue
        item = line[:-1].split('\t')
        if item[2] != feature:
            continue
        # get gene id
        try:
            geneid = ''
            index = line.index('GeneID:')
            index = index + 7
            while line[index] !=';' and line[index] != ',':
                    geneid = geneid + line[index]
                    index = index + 1
        except:
            geneid = '-'
            print line, 'dont have geneid'
        # output to file
        if geneid in genes:
            outline = '\t'.join([item[0],str(int(item[3])),item[4],geneid,'0',item[6]]) + '\n'
            out_handle.write(outline)
    gff_handle.close()
    out_handle.close()
    df = pd.read_csv(outBed,sep='\t',header=None,names=['chr','start','stop','geneid','None','Strand'],low_memory=False)
    df = df.drop_duplicates()
    df.to_csv(outBed,sep='\t',header=None,index=False)
    return outBed


def mergeBed(bed,feature='CDS'):
    """
    This function merge the bed file, if the interested feature in a gene overlaps.
    
    bed: str. Filename input bed file.
    feature: str. Feature of interest.
    """
    finalBed = bed[:-3] + feature + '.bed'
    cmd = ('sort -k1,1 -k2,2n {input} | bedtools merge -i stdin -s -c 4,5,6 '
           '-o distinct,distinct,distinct > {out}').format(input=bed,out=finalBed)
    subprocess.call(cmd,shell=True)
    df = pd.read_csv(finalBed,sep='\t',header=None,names=['chr','start','stop','geneid','None','Strand'],low_memory=False)
    df['start']=df['start']-1  # this is because bedtools would ignore the first position when calculating coverage
    df = df.sort(['chr','geneid','start'])
    df.to_csv(finalBed,sep='\t',header=None,index=False)
    return finalBed

def overlapCDS(bed,feature='CDS'):
    finalBed = bed[:bed.rindex('/')+3]+'_overlap_'+feature+'.txt'
    cmd = ('sort -k1,1 -k2,2n {input} | bedtools merge -i stdin -c 4,6 '
           '-o distinct,distinct > {out}').format(input=bed,out='inter.txt')
    subprocess.call(cmd,shell=True)
    handle = open('inter.txt')
    outHandle = open(finalBed,'w')
    outHandle.write('\t'.join(['Gene1','Gene2','Strand','\n']))
    for line in handle:
        item = line.split('\t')
        if ',' in item[3]:
            outHandle.write('\t'.join(item[3].split(','))+'\t'+item[4])
    outHandle.close()
    handle.close()
    os.remove('inter.txt')
    return finalBed
        
# df = pd.read_table('/data/shangzhong/RibosomeProfiling/cho_pr/03_choPrRefseq_sp.gene.txt',header=0)
# genes = df['GeneID'].astype(str).tolist()
# genes = list(set(genes))
# gff2Bed4Genes('/data/shangzhong/RibosomeProfiling/Database/combined.gff',genes,'exon')


    
#===============================================================================
#              calculate coverage
#===============================================================================
def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]

#===============================================================================
#             gene length statistics
#===============================================================================
def geneCDSStats(cds_File,geneIDs):
    """
    This function calculates the frequencies of length of the sp genes and no_sp_genes
    return a file with two columns, sp_genes and non_sp_genes. each row represents a position.
    values are how many genes have reached the length of row index.
    
    * cds_File: str. Filename of the file that stores 
    * geneIDs: list. A list of gene ids
    """
    gene_lens =[]
    cds_df = pd.read_csv(cds_File,header=None,sep='\t',names=['Chr','cds_start','cds_end','GeneID','None','Strand'],low_memory=False)
    for gene in geneIDs:
        gene_cds_df = cds_df[cds_df['GeneID']==gene]
        gene_obj = trpr(gene_cds_df)
        pos = gene_obj.get_trpr_pos(gene,level='gene')
        gene_lens.append(len(pos))
    len_stats = [0]*max(gene_lens)
    for g in gene_lens:
        for i in range(g):
            len_stats[i] = len_stats[i]+1
    percent = [x/len(geneIDs) for x in len_stats]
    return len_stats


def getCDSCov(gene_cov_df,pos,chr_len):
    """
    This function gets the coverage at each position 
    
    * gene_cov_df: df. Pandas df with 'coverage','Chr','pos'
    * pos: list. A list of positions
    * chr_len: int. Length of chromosome
    """
    cov_target = []
    dic = gene_cov_df.set_index('pos')['coverage'].to_dict()
    for p in pos:
        if (p > chr_len) or (p < 0):
            cov_target.append('-')
        else:
            if p in dic:
                cov_target.append(dic[p])
            else:
                cov_target.append(0)
    return cov_target

#===============================================================================
#        divide gene into 100 bin and calculate coverage in each bin
#===============================================================================
def ribo_percent_cov(exon_cov_df,cds_df,geneIDs):
    """
    This function calculates coverage at each integer percentile length of each gene in geneIDs.
    
    * exon_cov_df: dataframe. Stores coverage at each position of exons, has columns 
                        ['Chr','ex_start','ex_end','GeneID','None','Strand','pos','coverage']
    * gene_cds_df: dataframe. Stores 6 columns. ['Chr','cds_start','cds_end','GeneID','None','Strand']
    return a tabole with gene as row, each percentile as column
    """
    
    res_df = pd.DataFrame()  # this stores the final table, column is each geen, row is each percentage
    for gene in geneIDs:
        # get coverage of each position at all exon of gene
        gene_exon_cov_df = exon_cov_df[exon_cov_df['GeneID'].values==gene].copy() # gene exon coverage
        if gene_exon_cov_df.empty:
            print gene,'exon overlap with other genes'
            continue
        gene_exon_cov_df.loc[:,'real_pos'] = gene_exon_cov_df['pos'].add(gene_exon_cov_df['ex_start'])
        # get gene longest cds length
        gene_cds_df = cds_df[cds_df['GeneID']==gene]
        
        pos = [] # store the cds positions
        for start,end,stra in zip(gene_cds_df['cds_start'],gene_cds_df['cds_end'],gene_cds_df['Strand']):
            inter = range(start+1,end+1)
            if stra == '-':
                inter.reverse()
                pos = inter + pos
            else:
                pos.extend(inter)
        
        criteria = gene_exon_cov_df['real_pos'].map(lambda x: x in pos)
        target = gene_exon_cov_df[criteria]  # has coverage at each postion for target gene
        cov_target = target['coverage'].tolist()
        strand = target['Strand'].tolist()
        if '+' not in strand:
            cov_target.reverse()
        if ('+' in strand) & ('-' in strand):
            print gene,'has both strands'
            continue
        # divide into 100 bin
        #gene_len = gene_cds_df['length'].sum(0)
        gene_array = np.array(cov_target)
        gene_bin_array = np.array_split(gene_array,100)
        total_cov = np.sum(gene_array)
        # calculate average at each bin
        percent_cov = []
        for i in gene_bin_array:
            percent_cov.append(np.sum(i)/total_cov)
        res_df[gene]= pd.Series(percent_cov)
        print gene,'analysis done'
#     # 8. calculate means
#     res_df = res_df.fillna(0)
#     res_df['mean'] = res_df.mean(axis=1)
    res_df = res_df.transpose()
    return res_df


def ribo_percent5endCov(exon_cov_df,cds_df,geneIDs):
    """
    This function calculates coverage at each integer percentile length of each gene in geneIDs.
    It only count 5'end of a mapping read
    
    * exon_cov_df: dataframe. Stores coverage at each position of exons, has columns 
                        ['Chr','ex_start','ex_end','GeneID','None','Strand','pos','coverage']
    * gene_cds_df: dataframe. Stores 6 columns. ['Chr','cds_start','cds_end','GeneID','None','Strand']
    return a tabole with gene as row, each percentile as column
    """
    res_df = pd.DataFrame()
    for gene in geneIDs:
        # 1). get the chromosome for gene
        gene_cds_df = cds_df[cds_df['GeneID'].values==gene]
        chrome = list(set(gene_cds_df['Chr'].tolist()))
        if len(chrome) > 1:
            print gene,'maps to many chromosome'
            continue
        # 2) get the coverage for the specific chromosome
        gene_cov_df = exon_cov_df[exon_cov_df['Chr'].values==chrome[0]]  # columns= [count,chr,pos]
        # 3) list all CDS position for the gene
        pos = [] # store the cds positions
        for start,end,stra in zip(gene_cds_df['cds_start'],gene_cds_df['cds_end'],gene_cds_df['Strand']):
            inter = range(start+1,end+1)
            if stra == '-':
                inter.reverse()
                pos = inter + pos
            else:
                pos.extend(inter)
        strand = gene_cds_df['Strand'].tolist()
        if ('-' in strand) and ('+' in strand):
            assert False,gene+' have both strand'
        elif '-' in strand:
            pos = range(pos[0]+13,pos[0],-1) + pos[:-13]
            if 0 in pos:
                index = pos.index(0)
                pos = pos[:index]
        else:
            pos = range(pos[0]-13,pos[0]) + pos[:-13]
            if 0 in pos:
                index = pos.index(0)
                pos = pos[index:]
        # 4) get coverage for all CDS positions
        cov_target = []
        dic = gene_cov_df.set_index('pos')['coverage'].to_dict()
        for p in pos:
            if p in dic:
                cov_target.append(dic[p])
            else:
                cov_target.append(0)
        # 5) split into 100 bins
        gene_array = np.array(cov_target)
        gene_bin_array = np.array_split(gene_array,100)
        #total_cov = np.sum(gene_array)
        # calculate average at each bin
        percent_cov = []
        for i in gene_bin_array:
            percent_cov.append(np.sum(i))
        res_df[gene]= pd.Series(percent_cov)
        print gene,'analysis done'
#     # 8. calculate means
#     res_df = res_df.fillna(0)
#     res_df['mean'] = res_df.mean(axis=1)
    res_df = res_df.transpose()
    return res_df

def mRNApercent5endCov(exon_cov_df,exon_df,geneIDs):
    """
    This function calculates coverage at each integer percentile length of each gene in geneIDs.
    It considers all the exons in each gene and only count 5'end of a mapping read
    
    * exon_cov_df: dataframe. Stores coverage at each position of exons, has columns 
                        ['Chr','ex_start','ex_end','GeneID','None','Strand','pos','coverage']
    * gene_exon_df: dataframe. Stores 6 columns. ['Chr','cds_start','cds_end','GeneID','None','Strand']
    return a tabole with gene as row, each percentile as column
    """
    res_df = pd.DataFrame()
    for gene in geneIDs:
        # 1). get the chromosome for gene
        gene_cds_df = exon_df[exon_df['GeneID'].values==gene]
        chrome = list(set(gene_cds_df['Chr'].tolist()))
        if len(chrome) > 1:
            print gene,'maps to many chromosome'
            continue
        # 2) get the coverage for the specific chromosome
        gene_cov_df = exon_cov_df[exon_cov_df['Chr'].values==chrome[0]]  # columns= [count,chr,pos]
        # 3) list all exon positions for the gene
        pos = [] # store the cds positions
        for start,end,stra in zip(gene_cds_df['ex_start'],gene_cds_df['ex_end'],gene_cds_df['Strand']):
            inter = range(start+1,end+1)
            if stra == '-':
                inter.reverse()
                pos = inter + pos
            else:
                pos.extend(inter)
        # 4) get coverage for all CDS positions
        cov_target = []
        dic = gene_cov_df.set_index('pos')['coverage'].to_dict()
        for p in pos:
            if p in dic:
                cov_target.append(dic[p])
            else:
                cov_target.append(0)
        # 5) split into 100 bins
        gene_array = np.array(cov_target)
        gene_bin_array = np.array_split(gene_array,100)
        #total_cov = np.sum(gene_array)
        # calculate average at each bin
        percent_cov = []
        for i in gene_bin_array:
            percent_cov.append(np.sum(i))
        res_df[gene]= pd.Series(percent_cov)
        print gene,'analysis done'
    res_df = res_df.transpose()
    return res_df

#===============================================================================
#             coverage calculation at each position
#===============================================================================
def normGeneCov(line,up,down,by='mean'):
    """
    This function normalize gene coverage by divie count at each position by total count and mean.
    The mean value is calculated after removing the 1st 15 codon and last 10 codon.
    
    * line: str. One line in a gene coverage file. In the file, each line starts with gene id and 
    * up: int. How many positions are there before the TSS sites.
    * down: int. How many positions are there after the TSE sites.
    * total: int. Total count in the sample.
    * by: str. Normalize by mean or median.
    return normalized coverage for each gene
    """
    item = line[:-1].split('\t')
    gene = item[0];cov=item[1:]
    if gene == 'heavychain':
        pass
    df = pd.DataFrame({gene:cov})
    df.index = df.index-up       # change index to those relative to start sites
    try:
        df = df.replace('-',np.nan)
    except:
        pass
    df = df.astype(float)
    # filter by median = 1, no trimming of base pairs
    endIndex = df.index[-1]
    median = (df.loc[-15:endIndex-down-14]).median().values[0]
    if median < 1: return ''
    # filter out first 15 and last 10 codons for calculating mean
    if by == 'mean':
        mean = (df.loc[-15+45:endIndex-down-14-30]).mean().values[0]
        df = df/float(mean)
    if by == 'median':
        median = (df.loc[-15+45:endIndex-down-14-30]).median().values[0]
        df = df/float(median)
    return df


def normGeneCovWindow(geneCovFile,gene_sp,total,up,down,center):
    """
    This function normalize gene coverage by divie count at each position by total count and mean.
    The mean value is calculated after removing the 1st 15 codon and last 10 codon. 
    
    * geneCovList: str. Filename that stores coverage of all positon of all genes.
    * gene_sp: list. A list of genes that has signal peptide.
    * total: int. Number of total count in that sample.
    * up: int. Number of nts upstream of center position
    * down: int. Number of nts downstream of center position
    """
    handle = open(geneCovFile,'r')
    geneCov_df = pd.DataFrame()
    for line in handle:
        df = normGeneCov(line,up,down,total)
        if type(df) == str: continue
        gene = df.columns[0]
        if gene in gene_sp:
            df = df.rename(columns={gene:gene+'_sp'})
        # get sequence in the window
        if center == 0:
            geneCov_df = pd.concat([geneCov_df,df.loc[center-up:center+down]],axis=1)
        else:
            data_df = df.loc[df.shape[0]-2*up-down-1:df.shape[0]-down]
            data_df.index = range(-50,51)
            geneCov_df = pd.concat([geneCov_df,data_df],axis=1)
    geneCov_df = geneCov_df.T
    return geneCov_df


def normCodonCov(line,rm_start=15,rm_end=10,by='median'):
    """
    This function normalize the coding coverage.
    
    * line: str. First item is gene id, second column is longest protein. The rest are coverage at each codon.
    * by: str. 'median','mean'.
    * rm_start: int. remove the first rm_start nts of the amino acids.
    * rm_end: int. remove the rm_end nts of the amino acids.
    """
    item = line[:-1].split('\t')
    cov = item[2:]
    cov = [int(p) for p in cov]
    filter_cov = cov[rm_start:-rm_end]
    if filter_cov==[]: return ''   # length is less than 25 codons, skip
    median = np.median(np.array(cov[rm_start:-rm_end]))
    if median < 3: return ''        # median of trimmed cov is less than 1
    if by == 'median':
        norm = cov/np.median(np.array(cov[rm_start:-rm_end]))
    if by == 'mean':
        norm = cov/np.mean(np.array(cov[rm_start:-rm_end]))
    return norm
    

def CodonStallSites(f,cds_df,target_path,rm_start,rm_end):
    """
    This function detects stall sites in codon level for each gene
    
    * f: str. Filename, each line is a gene, protein with coverage of each codon.
    * cdsFile: str. Filename of cds information for each protein
    * target_path: str. Pathway of the output.
    """
    if not os.path.exists(target_path): os.mkdir(target_path)
    handle = open(f,'r')
    res = pd.DataFrame()
    chr_id=[]; gene=[]; chr_pos=[]; pr=[]; pr_relative_pos=[]; ratio=[]
    for line in handle:
        item = line[:-1].split('\t')  # 1st is gene id, second is protein access
        g = item[0] # gene
        p = item[1] # protein access
        norm = normCodonCov(line,rm_start,rm_end)
        if norm == '':
            continue
        pr_df = cds_df[cds_df['access'].values==p]
        pr_obj = trpr(pr_df)
        pr_pos = pr_obj.get_trpr_pos(p)
        for i in range(rm_start,len(norm)-rm_end):
            if norm[i] > 25:
                pr_relative_pos.append(i+1)
                chr_pos.append(pr_pos[i*3])
                chr_id.append(pr_obj.get_chrom(p,'protein'))
                gene.append(g)
                pr.append(p)
                ratio.append(norm[i])
    res['Chr'] = pd.Series(chr_id)
    res['GeneID'] = pd.Series(gene)
    res['Chr_pos'] = pd.Series(chr_pos)
    res['Pr_Access'] = pd.Series(pr)
    res['Pr_pos'] = pd.Series(pr_relative_pos)
    res['ratio'] = pd.Series(ratio)
    handle.close()
    res.to_csv(os.path.join(target_path,f.split('.')[0] + '.stallsites.txt'),sep='\t',index=False)
    

def StallSeq(stallFiles,cdsFile,refDNA_dic,prInconID,chr_pr_start_end_file,outputFile):
    """
    This function extracts nt/AA sequence before the stall sites.
    
    * stallFiles: list. A list of stall site files.
    * outputFile: str. Target file 
    """
    # 1. read gene with all proteins file
    cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False,names=['Chr','cds_start','cds_end','GeneID','Pr_Access','Strand'])
    # 2. pr inconsistant id
    prInconID_df = pd.read_csv(prInconID,header=None,names=['ID'])
    prInconIDs = prInconID_df['ID'].tolist()
    # 3 list chromosome protein start and end sites
    start_end_df = pd.read_csv(chr_pr_start_end_file,sep='\t',header=0)  # 'Chr' 'GeneID' 'Start' 'Stop' 'Pr_Access' 'Strand'
    # 4 define output file names 
    codonFile = outputFile.split('/')[-1][:4] + 'codon_freq.txt'
    codonHandle = open(codonFile,'w')
    aaFile = outputFile.split('/')[-1][:4] + 'aa_freq.txt'
    aaHandle = open(aaFile,'w')
    # 5. merge all sites into one df
    all_dfs = []
    for f in stallFiles:
        stall_df = pd.read_csv(f,sep='\t',header=0,usecols=[1,2,3]) # Chr    GeneID    Chr_pos    Pr_Access Pr_pos    ratio
        stall_df['GeneID'] = stall_df['GeneID'].astype(str)
        all_dfs.append(stall_df)
    merged_df = pd.concat(all_dfs)
    merged_df = merged_df.drop_duplicates()
    # 6. get proteins, and necessary dictionaries
    pr_gene_dic = merged_df.set_index('Pr_Access')['GeneID'].to_dict() # {protein: gene}
    pr_stall_dic = {k:list(v) for k,v in merged_df.groupby('Pr_Access')['Chr_pos']}  # {protein: sites}
    proteins = list(set(merged_df['Pr_Access'].tolist()))
    # 7 loop for each protein
    # (1). build the start and stop sites of proteins in the same gene into a list
    start_stop_sites = []
    for pr in proteins:
        if pr in prInconIDs:
            print pr,'sequence is inconsistant with refseq'
            continue
        # remove the start and stop sites in the same gene
        gene = pr_gene_dic[pr]
        gene_start_end_df = start_end_df[start_end_df['GeneID'].values==gene]
        start_stop_sites.extend(gene_start_end_df['Start'].tolist())
        start_stop_sites.extend(gene_start_end_df['Stop'].tolist())
        pr_stall_sites = [s for s in pr_stall_dic[pr] if s not in start_stop_sites]
        # get cds positions
        pr_df = cds_df[cds_df['Pr_Access'].values==pr]
        pr_obj = trpr(pr_df)
        pos = pr_obj.get_trpr_pos(pr)
        # get protein nt and AA sequence
        chrome = pr_obj.get_chrom(pr, id_type='protein') # get chromosome
        chr_seq = refDNA_dic[chrome].seq
        pr_nt = Seq(''.join([chr_seq[n-1].upper() for n in pos]),generic_dna)  # protein sequence
        if pos[0]>pos[1]:
            pr_nt = pr_nt.complement()
        AA = pr_nt.translate()
        if AA.endswith('*'):
            AA = AA[:-1]
        # get 15 nt and 5 AA before the stalling sites
        for stall in pr_stall_sites:
            try:
                index = pos.index(stall)
                pr_nt_seq = pr_nt[:index+3]
                AA_seq = AA[:int(index/3+1)]
                codonHandle.write(str(pr_nt_seq[-15:])+'\n')
                aaHandle.write(str(AA_seq[-5:])+'\n')
            except:
                print stall, 'not in protein of',pr
    codonHandle.close()
    aaHandle.close()
#===============================================================================
#                calculate at each position
#===============================================================================
def merge_percent_cov(files,sample,calculate_type='mean'):
    """
    This function extracts mean values of each file and calculate average of these mean values
    at each percentile. In each file, column is percentile, row are gene ids, the last one is mean value
     
    * files: list. A list of files.
    * sample: str. sample name for the files.
    * calculate_type: str. Default: mean.
    """
    res_df = pd.DataFrame()
    for f in files:
        df = pd.read_csv(f,sep='\t',header=0,index_col=0)
        df = df.transpose()
        if calculate_type=='mean':
            df[calculate_type] = df.mean(axis=1)
        if calculate_type=='median':
            df[calculate_type] = df.median(axis=1)
        if calculate_type=='std':
            df[calculate_type] = df.std(axis=1)
        res_df[f] = df[calculate_type]
#     if calculate_type=='mean':
#         res_df[sample] = res_df.mean(axis=1)
#     if calculate_type=='median':
#         res_df[sample] = res_df.median(axis=1)
    return res_df.loc[:,[sample]]

#===============================================================================
#         coverage at gene level
#===============================================================================
def mergeGeneExpress(outFile,coverFiles,chrPosCovFiles,covType='rawCount'):
    """
    This function calculates the rpkm for all input files and merge into one file
    
    * coverFiles: list. Filename of raw count data. ['GeneID','rawCount']
    * chrPosCovFiles: list. Filename of chromosome count data.
    * covType: str. type of calculation to represent the expression level.
    """
    gene_count_df = pd.DataFrame()
    for f in coverFiles:
        # get gene count data
        count_df = pd.read_csv(f,header=0,sep='\t',names=['GeneID','rawCount','length'])
        index = f.rfind('/')
        gene_count_df[f[index+1:index+4]] = count_df['rawCount']
    if covType!='rawCount':   # rpkm or rpm
        totalCounts = []
        length = count_df['length'].tolist()
        # get total count
        for t in chrPosCovFiles:
            df = pd.read_csv(t,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
            total = df['coverage'].sum()
            totalCounts.append(total)
        gene_count_df = gene_count_df.div(totalCounts)*(10**6)
        if covType=='rpkm': # calculate rpkm
            gene_count_df = gene_count_df.div(length,axis=0)*(1000)
    gene_count_df.insert(0,'GeneID',count_df['GeneID'])
    gene_count_df.to_csv(outFile,sep='\t',index=False)
    
    
#===============================================================================
#             motif stalling sites
#===============================================================================
def getcodon_AAFreq(freqFile):
    """
    This functions gets the frequencies of AA and codon. return them as dictionaries: {codon/AA:percentage}
    
    * freqFile: str. Has only one column of nt sequence of prtein.
    """
    codon_num_dic = {}
    AA_num_dic = {}
    handle = open(freqFile,'r')
    total = 0.0
    for line in handle:
        total = total + 1
        codon = line[-4:-1] # the last one is '\n' so it is -4:-1
        aa = Seq(codon,generic_dna).translate()
        if codon in codon_num_dic:
            codon_num_dic[codon] = codon_num_dic[codon] + 1
        else:
            codon_num_dic[codon] = 1
        if aa in AA_num_dic:
            AA_num_dic[aa] = AA_num_dic[aa] + 1
        else:
            AA_num_dic[aa] = 1
    for key in codon_num_dic:
        codon_num_dic[key] =  codon_num_dic[key]/total
    for key in AA_num_dic:
        AA_num_dic[key] = AA_num_dic[key]/total
    return codon_num_dic,AA_num_dic


#===============================================================================
#                 class defination
#===============================================================================
class trpr(object):
    """
    Input is a pandas dataframe from 01_pr_cds.txt. Chr    cds_start    cds_end    GeneID    PrAccess    Strand
    """
    def __init__(self,df):
        self.df = df
        self.df.columns = ['chr','start','end','geneid','access','strand']
    # id part
    def get_chrom(self,g,id_type='gene'):
        """
        This function get chromosome from gene name
        """
        if id_type == 'gene':
            q_id = 'geneid'
        else:
            q_id = 'access'
        g_df = self.df[self.df[q_id].values==g]
        chrom = list(set(g_df['chr'].tolist()))
        if len(chrom)==1:
            return chrom[0]
        else:
            assert False,chrom + 'is wrong,morethan one mapping'
    
    def all_gene_ids(self):
        """get all the gene ids in the dataframe"""
        genes = list(set(self.df['geneid'].tolist()))
        return genes
    
    def fwd_gene_ids(self):
        """get all gene ids in the forward strand"""
        genes = list(set(self.df[self.df['strand'].values=='+'].tolist()))
        return genes
    
    def rev_gene_ids(self):
        """get all gene ids in the reverse strand"""
        genes = list(set(self.df[self.df['strand'].values=='-'].tolist()))
        return genes
    
    def genes2multi_chr(self):
        """This function finds genes that map to multiple chromosomes"""
        gene_chr_df = self.df[['chr','geneid']].drop_duplicates()
        gene_chr_dic = gene_chr_df.groupby('geneid').groups
        genes = []
        for key in gene_chr_dic:
            if len(gene_chr_dic[key])>1:
                genes.append(key)
        return genes
        
    def trpr_gene_dic(self):
        """This function generates a function of {praccess:geneid}"""
        pr_g_dic = self.df.set_index('access')['geneid'].to_dict()
        return pr_g_dic
    # position part
    def get_longest_trpr(self,gene):
        """get the longest protein access or transcript access"""
        gene_df = self.df[self.df['geneid'].values==gene].copy()
        gene_df.loc[:,'len'] = gene_df['end'] - gene_df['start']
        trprs = list(set(gene_df['access'].tolist()))
        if len(trprs)!=1:
            lens = {}
            for trpr in trprs:
                trpr_len = gene_df[gene_df['access'].values==trpr]['len'].sum()
                lens[trpr_len] = trpr
            return lens[max(lens.keys())]
        else:
            return trprs[0]
        
    def get_trpr_pos(self,trpr,up=0,down=0,level='trpr'):
        """
        This function gets all position of transcript or protein or gene
        """
        # decide whether to proceed
        trpr = str(trpr)
        if level == 'trpr':
            trpr_df = self.df[self.df['access'].values==trpr]
        elif level == 'gene':
            trpr_df = self.df[self.df['geneid'].values==trpr]
        chrom = list(set(trpr_df['chr'].tolist()))
        if len(chrom) !=1:
            assert False,'Protein {trpr} map to multiple chromosome'.format(trpr=trpr)
        strand = list(set(trpr_df['strand'].tolist()))
        if len(strand)!=1: 
            assert False,'Protein {trpr} map to both strand'.format(trpr=trpr)
        # get position
        pos = [] # store the cds positions
        for start,end,stra in zip(trpr_df['start'],trpr_df['end'],trpr_df['strand']):
            inter = range(int(start)+1,int(end)+1)
            if stra == '-':
                inter.reverse()
                pos = inter + pos
            else:
                pos.extend(inter)
        if pos[0]>pos[1]: # - strand
            pos = natsorted(list(set(pos)))[::-1]
            new_pos = range(pos[0]+up,pos[0],-1) + pos + range(pos[-1]-1,pos[-1]-down,-1)
        else:
            pos = natsorted(list(set(pos)))
            new_pos = range(pos[0]-up,pos[0]) + pos + range(pos[-1]+1,pos[-1]+down)
        return new_pos
    
    def get_ribo_trpr_pos(self,trpr,offset,level='trpr'):
        """get position of A site for ribosome profiling.
        includes the offset of the distance between 5' end and the A site. 
        """
        pos = self.get_trpr_pos(trpr,level=level)
        off = offset + 1
        if pos[0]>pos[1]: # - strand
            new_pos = range(pos[0]+off,pos[0],-1) + pos[:-off]
        else:
            new_pos = range(pos[0]-off,pos[0]) + pos[:-off]
        return new_pos
    
    def htseq_genome_array_set(self,stranded=False):
        """
        This function generates a genomic array of set for finding the position that are ovelapped by genes
        * stranded: logical. 
        """
        gene_array_set = ht.GenomicArrayOfSets('auto',stranded=stranded)
        for row in self.df.itertuples(index=False):
            try:
                gene_array_set[ht.GenomicInterval(row[0],row[1],row[2],row[5])] += row[3]
            except:
                print row,'has problem'
                continue
        return gene_array_set
       
    def genes_cov_target_pos(self,chrom,pos):
        """
        This function returns the genes that cover the position.
        * chrom: str. chromosome.
        * pos: int. position. 
        return list of genes containing that position
        """
        chr_df = self.df[self.df['chr'].values==chrom]
        chr_df = chr_df[(chr_df['start']< pos) & (chr_df['end']>=pos)]
        genes = list(set(chr_df['geneid'].tolist()))
        return genes
    
    def all_pr_start_end_pos(self):
        """This function gets all the start and end position of all proteins.
        For end position it is the start of the stop codon.
        """
        proteins = list(set(self.df['access'].tolist()))
        chr_pr_start_end_sites = pd.DataFrame()
        Chr = [];GeneID=[];Start=[];Stop=[];Pr=[];Strand=[]
        for pr in proteins:
            pr_df = self.df[self.df['access'].values==pr]
            pr_obj = trpr(pr_df)
            try:
                pos = pr_obj.get_trpr_pos(pr)
            except:
                print pr,'map to multiple chromosome'
                continue
            strand = pr_df['strand'].tolist()
            # assign
            Chr.append(pr_df['chr'].tolist()[0])
            GeneID.append(pr_df['geneid'].tolist()[0])
            Start.append(pos[0])
            Stop.append(pos[-3])
            Pr.append(pr)
            Strand.append(strand[0])
        chr_pr_start_end_sites['Chr']=pd.Series(Chr)
        chr_pr_start_end_sites['GeneID']=pd.Series(GeneID)
        chr_pr_start_end_sites['Start'] = pd.Series(Start)
        chr_pr_start_end_sites['Stop'] = pd.Series(Stop)
        chr_pr_start_end_sites['Pr_Access'] = pd.Series(Pr)
        chr_pr_start_end_sites['Strand'] = pd.Series(Strand)
        
        return chr_pr_start_end_sites
    
    def get_exn_cds_start_end_pos(self,gene):
        """
        This function gets the start and end position of exons or cds of a gene provided
        returns 2 list of start position and end position
        """
        gene_df = self.df[self.df['geneid'].values==gene]
        strand = list(set(gene_df['strand'].tolist()))[0]
        if strand == '+':
            exn_start = gene_df['start'].tolist()
            exn_end = gene_df['end'].tolist()
#             ex_start = [abs(n-exn_start[0]) for n in exn_start]
#             ex_end = [abs(n-exn_start[0]) for n in exn_end]
        else:
            exn_start = gene_df['end'].tolist()
            exn_end = gene_df['start'].tolist()
            exn_start.reverse()
            exn_end.reverse()
#             ex_start = [abs(n-exn_end[-1]) for n in exn_start]
#             ex_end = [abs(n-exn_end[-1]) for n in exn_end]
        return exn_start,exn_end
    
def write_dic(dic,File):
    """
    This function write the dictionary to a file with 'key tab value'
    * dic: dict. Dictionary that need to be write.
    * File: filename. File that stores the output.
    """
    Handle = open(File,'w')
    for key in natsorted(dic):
        Handle.write(str(dic[key])+'\t'+key+'\n')
    Handle.close()

def fwd_rev_cov(bamFile,max_len=37,seq_len=50):
    """
    This function calculates mapping 5' and 3' end position of each read and count the number for those with the same positions
    * bamFile: str. Bam file
    output a file with cover dataframe. The columns are: ['num','chr','end5','end3','strand','len']    
    """
    bamHandle = pysam.AlignmentFile(bamFile,'rb')
    Handle = bam_parse(bamHandle)
    count_dic = Handle.bam_fwd_rev_count(max_len,seq_len)
    covFile = bamFile[:3] + '_cov.txt'
    write_dic(count_dic,covFile)
    
    
class cover(object):
    """
    Input is a pandas dataframe with 6 columns. # num chr start end strand length
    """
    def __init__(self,df):
        self.df = df
        self.df.columns = ['num','chr','end5','end3','strand','len']
        
    def get_pos_coverage(self,chrom,pos,strand_specific='N'):
        """
        This function calculates coverage at each position of a gene, and returns a list with same lenght of pos.
        
        * chrom: str. chromosome.
        * pos: list. A list of integer positaion.
        * strand_specific: str. The RNA seq or DNAseq is strand specific or not
        """
        chr_df = self.df[self.df['chr'].values==chrom]
        if pos[0] > pos[1]: # - strand
            if strand_specific=='Y':
                neg_df = chr_df[chr_df['strand'].values=='-']
            else: neg_df = chr_df
            #neg_df = neg_df[(neg_df['end3']>=pos[-1]) & (neg_df['end3']<=pos[0])]
            dic = {k:list(v) for k,v in neg_df.groupby('end3')['num']}
        else:
            if strand_specific=='Y':
                neg_df = chr_df[chr_df['strand'].values=='+']
            else:
                neg_df = chr_df
            #neg_df = neg_df[(neg_df['end5']>=pos[0]) & (neg_df['end5']<=pos[-1])]
            dic = {k:list(v) for k,v in neg_df.groupby('end5')['num']}
        if neg_df.empty:
            cov = [0]* len(pos)
        else:
            cov = []
            for p in pos:
                if p in dic:
                    cov.append(sum(dic[p]))
                else:
                    cov.append(0)
        return cov
    
    def get_align_lengths(self):
        """This function get how many different alignment length a file has
        """
        length = list(set(self.df['len'].tolist()))
        return length
            
            
def read_2_multi_genes(gene_array,chrom,start,end):
    """
    This function test whether a read map to multiple genes
    * return True if read map to only one gene and should be kept
    """
    set5end = gene_array[ht.GenomicPosition(chrom,start)]
    set3end = gene_array[ht.GenomicPosition(chrom,end)]
    if (len(set5end)==1) and (len(set3end)==1) and (set5end!=set3end):
        return False
    elif(set5end==set()) and (set3end==set()):
        return False
    else:
        return True


def posCoverage(gene_df,cov_df,pos,gene_array_set):
    """
    This function calculates the coverage at each position of gene
    
    * gene_df: pd dataframe. with 6 columns # 'chr','start','end','geneid','praccess','strand'. Here it is used to get chromosoem.
    * cov_df: pd dataframe. with 6 columns # 'num','chr','end5','end3','strand','len'. end5 and end 3 are alignment end, not read end.
    * pos: list. List of position.
    * gene_array_set: htseq genearray set. Used for removing reads mapping to multiple genes.
    """
    gene = list(set(gene_df['geneid'].tolist()))
    g_obj = trpr(gene_df) 
    chrom = g_obj.get_chrom(gene[0])  # chromosome
    gCov_df = cov_df[cov_df['chr'].values==chrom]  # coverage dataframe of the chromosome
    if gCov_df.empty:
        return [0]*len(pos)
    max_len = max(gCov_df['len'].tolist())
    if pos[0]>pos[1]:  # - strand
        large = pos[0] + max_len
        small = pos[-1] - max_len
        new_pos = range(large,pos[0],-1) + pos
        gCov_df = gCov_df[(gCov_df['end3'].isin(new_pos))]
    else:
        small = pos[0] - max_len
        large = pos[-1] + max_len
        new_pos = range(small,pos[0]) + pos
        gCov_df = gCov_df[(gCov_df['end5'].isin(new_pos))]
    if gCov_df.empty:
        return [0]*len(pos)
    cri = gCov_df.apply(lambda row: read_2_multi_genes(gene_array_set,row['chr'],row['end5'],row['end3']),axis=1)
    fil_cov_df = gCov_df[cri]
    cov_obj = cover(fil_cov_df)            
    cov_pos = cov_obj.get_pos_coverage(chrom, pos)
    return cov_pos      


def gene_pr_cov(covFile,cdsFile,ribo_offset_file,outpath,covType='pos',genes=[],pr_cov_level='codon'):
    """
    This function gets all protein position coverage in a coverfile, longest protein for a gene is chosen
    
    * covFile: str. Filename with 6 columns.
    * cdsFile: str. Filename with 6 columns.
    * ribo_offset_file: str. Filename with 2 columns: length,offset
    * covType: str. get coverage at position level or at the gene level.
    * genes: list. If setted it will only get coverage for the genes provided, if not it will calculate all genes.
    If covType is pos: each line of output would be "geneid    protein access    coverage at each position"
    If covType is gene: each line of output would be "geneid    total count    length of total CDS positions"
    """
    if not os.path.exists(outpath): os.mkdir(outpath)
    # 1) read coverage files
    cov_df = pd.read_csv(covFile,sep='\t',header=None,names=['num','chr','end5','end3','strand','len'],low_memory=False)
    # 2) read protein position file
    cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)   # Chr    cds_start    cds_end    GeneID    Pr_Access    Strand
    trpr_obj = trpr(cds_df)
    # 3). create array set by htseq. This is for retriving positions that are covered by many genes
    gene_array_set = trpr_obj.htseq_genome_array_set(False)
    # 4) get the offset length
    ribo_offset_df = pd.read_csv(ribo_offset_file,sep='\t',header=0,names=['len','off'],low_memory=False)
    off_len_dic = {k:list(v) for k,v in ribo_offset_df.groupby('off')['len']}
    # 5) get all genes
    if genes==[]:
        genes = trpr_obj.all_gene_ids()
    genes = natsorted(genes)
    # 6) get coverage for each gene
    if covType == 'pos':
        handle = open(outpath+'/'+covFile.split('.')[0]+'_prpos.txt','w')
        for g in genes:
            gene_cov_df = pd.DataFrame()
            try:
                gene_df = cds_df[cds_df['geneid'].values==g]
                g_obj = trpr(gene_df)
                long_pr = g_obj.get_longest_trpr(g)  # longest protein of the gene
                pos = g_obj.get_trpr_pos(long_pr) # position of the protein
                for offset in natsorted(off_len_dic):  # alignment length
                    off_cov_df = cov_df[cov_df['len'].isin(off_len_dic[offset])]
                    off = offset + 1
                    # if gene is NeoR, offset change to 16 for all alignment length
                    if g == 'NeoRKanR':
                        off = 17
                    # end
                    if pos[0]>pos[1]: # - strand
                        new_pos = range(pos[0]+off,pos[0],-1) + pos[:-off]
                    else:
                        new_pos = range(pos[0]-off,pos[0]) + pos[:-off]
                    cov_pos = posCoverage(gene_df,off_cov_df,new_pos,gene_array_set)
                    gene_cov_df[str(offset)] = pd.Series(cov_pos)
            except:
                print g,'maps to multiple scaffold'
                continue
            gene_cov = gene_cov_df.sum(axis=1).tolist()
            
            if pr_cov_level == 'codon':
                pr_cov = [str(sum(p)) for p in chunks(gene_cov,3)]
            elif pr_cov_level == 'nt':
                pr_cov = [str(p) for p in gene_cov]
            handle.write(g+'\t'+long_pr+'\t'+'\t'.join(pr_cov)+'\n')
    else:  # covType == 'gene'
        handle = open(outpath+'/'+covFile.split('.')[0]+'_geneCount.txt','w')
        handle.write('\t'.join(['GeneID','count','length'])+'\n')
        for g in genes:
            gene_cov_df = pd.DataFrame()
            try:
                gene_df = cds_df[cds_df['geneid'].values==g]
                g_obj = trpr(gene_df)
                pos = g_obj.get_trpr_pos(g,level='gene') # position of the protein
                for offset in natsorted(off_len_dic):  # alignment length
                    off_cov_df = cov_df[cov_df['len'].isin(off_len_dic[offset])]
                    off = offset + 1
                    if pos[0]>pos[1]: # - strand
                        new_pos = range(pos[0]+off,pos[0],-1) + pos[:-off]
                    else:
                        new_pos = range(pos[0]-off,pos[0]) + pos[:-off]                    
                    cov_pos = posCoverage(gene_df,off_cov_df,new_pos,gene_array_set)
                    gene_cov_df[str(off)] = pd.Series(cov_pos)
            except:
                print g,'maps to multiple scaffold'
                continue
            gene_cov = gene_cov_df.sum(axis=1).tolist()
            handle.write('\t'.join([g,str(sum(gene_cov)),str(len(gene_cov))])+'\n')
    handle.close()
            

def gene_tr_cov(exnFile,cdsFile,all_id_file,covFile,outpath,genes=[],offset='no',ribo_offset_file=''):
    """This function get coverage at position level for RNAseq. gene transcript coverage. The transcript with longest CDS length is chosen.
    
    * exnFile: str. Filename has 6 columns. ['Chr','Start','End,'GeneID','TrAccess','Strand']
    * cdsFile: str. Filename has 6 columns. ['Chr','Start','End,'GeneID','PrAccess','Strand']
    * all_id_file: str. Filename has all kinds ids
    * covFile: str. Filename
    * outpath: str. Path that store the results. Each line in each file, first column is gene id, second is tr id, third is pr id, the rest is coverage at each position of transcript.
    * genes: list. A list of genes.
    * offset: 'yes' or 'no'. For RNAseq should be no. 'yes' is for getting riboseq coverage at transcript level
    * ribo_offset_file: str. File that stores the offset for each alignment length.
    """
    #------------- build protein:transcript dictionary -----------------------------
    df = pd.read_csv(all_id_file,sep='\t',header=0,names=['id','symbol','chr','tr','pr'])
    dic_df = df[['tr','pr']].drop_duplicates()
    pr_tr_dic = dic_df.set_index('pr')['tr'].to_dict()
    
    if not os.path.exists(outpath): os.mkdir(outpath)
    # 1) read coverage files
    cov_df = pd.read_csv(covFile,sep='\t',header=None,names=['num','chr','end5','end3','strand','len'],low_memory=False)
    # 2) read transcript position file
    exn_df = pd.read_csv(exnFile,sep='\t',header=0,low_memory=False)   # Chr    cds_start    cds_end    GeneID    Pr_Access    Strand
    exn_obj = trpr(exn_df)
    # 3) read cds position file
    cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
    cds_obj = trpr(cds_df)
    # 4). create array set by htseq. This is for retriving positions that are covered by many genes
    gene_array_set = exn_obj.htseq_genome_array_set(False)
    # 5) get all protein coding genes 
    if genes == []:  # consider all genes
        genes = cds_obj.all_gene_ids()
    genes = natsorted(genes)
    # 6) get coverage for each gene
    handle = open(outpath+'/'+covFile.split('.')[0]+'_trpos.txt','w')
    if offset == 'no':
        for g in genes:
            try:
                gene_cds_df = cds_df[cds_df['geneid'].values==g]
                gene_exn_df = exn_df[exn_df['geneid'].values==g]
                g_cds_obj = trpr(gene_cds_df)
                g_exn_obj = trpr(gene_exn_df)
                long_pr = g_cds_obj.get_longest_trpr(g)  # longest protein of the gene
                tr = pr_tr_dic[long_pr]  # correspond transcript
                pos = g_exn_obj.get_trpr_pos(tr) # position of the protein
                cov_pos = posCoverage(gene_exn_df,cov_df,pos,gene_array_set)
                cov_pos = [str(p) for p in cov_pos]
            except:
                print g,'maps to multiple scaffold'
                continue
            handle.write(g+'\t'+tr+'\t'+long_pr+'\t'+'\t'.join(cov_pos)+'\n')
    else:  # consider the offset for riboseq coverage at transcript level
        # 1) get the offset length
        ribo_offset_df = pd.read_csv(ribo_offset_file,sep='\t',header=0,names=['len','off'],low_memory=False)
        off_len_dic = {k:list(v) for k,v in ribo_offset_df.groupby('off')['len']}
        # 2) get coverage for each gene
        for g in genes:
            gene_cov_df = pd.DataFrame()
            try:
                gene_cds_df = cds_df[cds_df['geneid'].values==g]
                gene_exn_df = exn_df[exn_df['geneid'].values==g]
                g_cds_obj = trpr(gene_cds_df)
                g_exn_obj = trpr(gene_exn_df)
                long_pr = g_cds_obj.get_longest_trpr(g)  # longest protein of the gene
                tr = pr_tr_dic[long_pr]  # correspond transcript
                pos = g_exn_obj.get_trpr_pos(tr) # position of the protein
                for offset in natsorted(off_len_dic):  # alignment length
                    off_cov_df = cov_df[cov_df['len'].isin(off_len_dic[offset])]
                    off = offset + 1
                    # if gene is NeoR, offset change to 16 for all alignment length
                    if g == 'NeoRKanR':
                        off = 17
                    # end
                    if pos[0]>pos[1]: # - strand
                        new_pos = range(pos[0]+off,pos[0],-1) + pos[:-off]
                    else:
                        new_pos = range(pos[0]-off,pos[0]) + pos[:-off]
                    cov_pos = posCoverage(gene_exn_df,off_cov_df,new_pos,gene_array_set)
                    gene_cov_df[str(offset)] = pd.Series(cov_pos)
            except:
                print g,'maps to multiple scaffold'
                continue
            gene_cov = gene_cov_df.sum(axis=1).tolist()  # sum all alignment lengths
            pr_cov = [str(p) for p in gene_cov]
            handle.write(g+'\t'+tr+'\t'+long_pr+'\t'+'\t'.join(pr_cov)+'\n')
    handle.close()        


def merge_rep_gene_pos_cov(totalCount,covFiles,gene,exn_df,covType='trpt'):
    """
    This function calculates the ribo seq coverage at each position and . and then merge the replicates.
    * totalCount: list. A list of total count for replicated bam files.
    * covFiles: list. A list of filenames with 6 columns. # 'num','chr','end5','end3','strand','len'
    * gene: str. Target gene
    * exn_df: dataframe. File name of 6 columns.
    * covType: str. 'gene' or 'trpt'.
    """
    gene_cov_df = pd.DataFrame()
    for f in covFiles:
        cov_df = pd.read_csv(f,sep='\t',header=None,names=['num','chr','end5','end3','strand','len'])
        gene_df = exn_df[exn_df['GeneID'].values==gene]
        gene_obj = trpr(gene_df)
        trpr_id = list(set(gene_df['access'].tolist()))[0]
        pos = gene_obj.get_trpr_pos(trpr_id)
        if covType == 'gene':
            tr_last = pos[-1]
            if pos[0]>pos[1]:  # - strand
                pos = range(pos[0],pos[-1]-1,-1)
            else:
                pos = range(pos[0],pos[-1]+1)
        gene_array_set = gene_obj.htseq_genome_array_set(False)
        cov_pos = posCoverage(gene_df,cov_df,pos,gene_array_set)
        gene_cov_df[f] = pd.Series(cov_pos)
    print 'position of start of 3 UTR relative to gene:',pos.index(tr_last)
    gene_cov_df = gene_cov_df.div(totalCount) * (10**6)
    gene_cov_df['mean'] = gene_cov_df.mean(axis=1)
    return gene_cov_df['mean']


def utr3_cov(exnFile,cdsFile,all_id_file,covFile,outpath,genes=[]):
    """This function calculates the coverage at the 3'UTR region
    """
    if not os.path.exists(outpath): os.mkdir(outpath)
    #------------- 1. build protein:transcript dictionary -----------------------------
    df = pd.read_csv(all_id_file,sep='\t',header=0,names=['id','symbol','chr','tr','pr'])
    dic_df = df[['tr','pr']].drop_duplicates()
    pr_tr_dic = dic_df.set_index('pr')['tr'].to_dict()
    # 1) read coverage files
    cov_df = pd.read_csv(covFile,sep='\t',header=None,names=['num','chr','end5','end3','strand','len'],low_memory=False)
    # 2) read transcript position file
    exn_df = pd.read_csv(exnFile,sep='\t',header=0,low_memory=False)   # Chr    cds_start    cds_end    GeneID    Pr_Access    Strand
    exn_obj = trpr(exn_df)
    # 3) read cds position file
    cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
    cds_obj = trpr(cds_df)
    # 4). create array set by htseq. This is for retriving positions that are covered by many genes
    gene_array_set = exn_obj.htseq_genome_array_set(False)
    # 5) get all protein coding genes 
    if genes == []:  # consider all genes
        genes = cds_obj.all_gene_ids()
    genes = natsorted(genes)
    # 6) get coverage for each gene
    handle = open(outpath+'/'+covFile.split('.')[0]+'_3utr.txt','w')
    for g in genes:
        try:
            gene_cds_df = cds_df[cds_df['geneid'].values==g]
            gene_exn_df = exn_df[exn_df['geneid'].values==g]
            g_cds_obj = trpr(gene_cds_df)
            g_exn_obj = trpr(gene_exn_df)
            long_pr = g_cds_obj.get_longest_trpr(g)  # longest protein of the gene
            tr = pr_tr_dic[long_pr]  # correspond transcript
            cds_pos = g_cds_obj.get_trpr_pos(long_pr)  # position of cds
            rna_pos = g_exn_obj.get_trpr_pos(tr) # position of the rna
            index = rna_pos.index(cds_pos[-1])
            pos = rna_pos[index+1:]  # utr position
            if pos == []:
                print g,'has no 3utr'
                continue
            cov_pos = posCoverage(gene_exn_df,cov_df,pos,gene_array_set)
            cov_pos = [str(p) for p in cov_pos]
        except:
            print g,'maps to multiple scaffold'
            continue
        handle.write(g+'\t'+tr+'\t'+long_pr+'\t'+'\t'.join(cov_pos)+'\n')
    handle.close()


def refseq_reffa_inconsist_pr(cdsFile,refFaFile,prFaFile,pr_id_file,pr_seq_file):
    """
    This function test whether the AA extracted from annotation file are the same with those in refseq.
    # 1. protein cds coverage at nt level.
    # 2. protein cds coverage at codon level
    # 3. AA sequence
    # 4. nt sequence
    # 5. check how many proteins have different sequence with chromosome extracted sequence
    # 6. for stalling sites, remember to remove the start and stop sites of other proteins which are covered by the analyzed transcript
    
    * cdsFile: str. File name of cds information. has 6 columns.
    * refFaFile: str. Reference fasta file.
    * pr_id_file: str. Filename that stores the inconsistant protein accessions.
    * pr_seq_file: str. fasta file that stores the inconsistant protein sequences.
    """
    # 1. check the inconsistancy between the genome translated sequence and the refseq protein sequence
    cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
    refDNA_dic = SeqIO.index(refFaFile,'fasta')
    # 2. read correspond refseq protein sequences
    prRef_dic = SeqIO.index(prFaFile,'fasta')
    # 3. get inconsistant genes
    proteins = list(set(cds_df['Pr_Access'].tolist()))
    handle = open(pr_id_file,'w')
    handle1 = open(pr_seq_file,'w')
    for pr in proteins:
        pr_df = cds_df[cds_df['Pr_Access'].values==pr]
        pr_obj = trpr(pr_df)
        # get chrome sequence
        try:
            chrome = pr_obj.get_chrom(pr, 'protein')
            chr_seq = str(refDNA_dic[chrome].seq)
        except:
            print pr,'maps to multiple chromosome'
            continue
        # get pr nt sequence
        pos = pr_obj.get_trpr_pos(pr)
        pr_nt = Seq(''.join([chr_seq[n-1].upper() for n in pos]),generic_dna)
        # get pr AA sequence
        if pos[0]<pos[1]:
            AA = pr_nt.translate()
        else:
            AA = pr_nt.complement().translate()
        if AA.endswith('*'):
            AA = AA[:-1]
        try:
            if AA != str(prRef_dic[pr].seq):
                handle.write(pr+'\n')
                handle1.write('\n'.join(['>'+pr,str(AA)])+'\n')
        except:
            print pr,'not in the annotation file'
    handle.close()
    handle1.close()

#===============================================================================
#             frame count
#===============================================================================
def frame_count(cov,calType='sum'):
    """This function calculcate and compare the frame coverage.
    * cov: list. A list of coverage at each position.
    * calType: str. If set as sum, the frame coverage would be the summation.
                    If set as median, the frame coverage would be the median of all in the same frame.
                    If set as scale, the framme coverage would be scale inside each codon, largest value is 1, the other is the fraction.
    """
    if len(cov)%3 != 0:
        print 'sequence length is not times of 3'
    frames = [0] * 3
    if calType == 'median':
        full_list = [[] for i in range(3)]
        for i in range(len(cov)):
            if i%3 == 0:
                full_list[2].append(int(cov[i]))
            elif i%3 == 1:
                full_list[0].append(int(cov[i]))
            else:
                full_list[1].append(int(cov[i]))
        for i in range(3):
            frames[i] = np.median(full_list[i])
    elif calType == 'sum':
        for i in range(len(cov)):
            if i%3 == 0:
                frames[2] += int(cov[i])
            elif i%3 == 1:
                frames[0] += int(cov[i])
            else:
                frames[1] += int(cov[i])
        sum_cov = sum(frames)
        frames = [float(p)/max(1,sum_cov) for p in frames]
    elif calType == 'scale':
        cov_list = chunk(cov,3)
        for codon in cov_list:
            if len(codon) != 3:
                continue
            max_cov = max(codon)
            if max_cov == 0:
                continue
            for i in range(len(codon)):
                try:
                    codon[i] = codon[i]/max_cov
                except:
                    raise False, 'divide by 0'
#             codon = [c/float(codon[1]) for c in codon]
            frames[1] += float(codon[2])
            frames[2] += float(codon[0])
            frames[0] += float(codon[1])
    return frames


def genes_frame_cov(pr_pos_cov_file,out_path,genes=[],calType='sum'):
    """This function calculates frame coverage for all the genes that are listed.
    * pr_pos_cov_file: filename of coverage file.
    * out_path: str. Pathway of the results.
    * genes: list. A list of genes.
    * calType: str. can be 'sum','median','scale'.
    """
    if not os.path.exists(out_path): os.mkdir(out_path)
    handle = open(pr_pos_cov_file)
    outFile = out_path + '/' + pr_pos_cov_file[:-3]+'frame.txt'
    out_handle = open(outFile,'w')
    out_handle.write('\t'.join(['geneid','frame1','frame2','frame3\n']))
    for line in handle:
        item = line[:-1].split('\t')
        geneid = item[0]
        if (genes!=[]) and (geneid not in genes): continue
        cov = item[2:]  # coverage
        cov = [int(p) for p in cov]
        if sum(cov) < 128:
            continue
        frames = frame_count(cov,calType)
        frames = [str(p) for p in frames]
        out_handle.write(geneid + '\t' + '\t'.join(frames)+'\n')
    handle.close()
    out_handle.close()
        