from __future__ import division
import subprocess,os,sys
import pandas as pd
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('ggplot')
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

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
    
    mapdf = pd.read_csv(mapFile,sep='\t',header=0,names=['GeneID','GeneSymbol','Chr','TrAccess','PrAccess'])
    # extract either cho or hamster
    if mapSource=='gene2refseq':
        mapdf = pd.read_csv(mapFile,sep='\t',skiprows=[0],header=None,usecols=[1,3,5,7],names=['GeneID','TrAccess','PrAccess','Chr'])
        if organism == 'chok1':
            criteria = mapdf['Chr'].map(lambda x: 'NW_0036' in x)
            mapdf = mapdf[criteria]
        if organism == 'hamster':
            criteria = mapdf['Chr'].map(lambda x: 'NW_0068' in x)
            mapdf = mapdf[criteria]
        
    mapdf = mapdf.drop_duplicates()
    mapdf['TrAccess'] = mapdf['TrAccess'].apply(lambda x: removeIDVersion(x))
    mapdf['PrAccess'] = mapdf['PrAccess'].apply(lambda x: removeIDVersion(x))
    # read signalP result
    signalP_df = pd.read_csv(signalP,header=0,skiprows=[0],delim_whitespace=True)
    signalP_df['name'] = signalP_df['name'].apply(lambda x: removeIDVersion(x))
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
    df = pd.read_csv(inputFile,header=0,sep='\t')
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
    Only extract exons of genes.
    
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
    df = pd.read_csv(outBed,sep='\t',header=None,names=['chr','start','stop','geneid','None','Strand'])
    df = df.drop_duplicates()
    df.to_csv(outBed,sep='\t',header=None,index=False)
    finalBed = outBed[:-3] + feature + '.bed'
    cmd = ('sort -k1,1 -k2,2n {input} | bedtools merge -i stdin -s -c 4,5,6 '
           '-o distinct,distinct,distinct > {out}').format(input=outBed,out=finalBed)
    subprocess.call(cmd,shell=True)
    os.remove(outBed)
    df = pd.read_csv(finalBed,sep='\t',header=None,names=['chr','start','stop','geneid','None','Strand'])
    df['start']=df['start']-1  # this is because bedtools would ignore the first position when calculating coverage
    df = df.sort(['chr','geneid','start'])
    df.to_csv(finalBed,sep='\t',header=None,index=False)
    return finalBed

# df = pd.read_table('/data/shangzhong/RibosomeProfiling/cho_pr/03_choPrRefseq_sp.gene.txt',header=0)
# genes = df['GeneID'].astype(str).tolist()
# genes = list(set(genes))
# gff2Bed4Genes('/data/shangzhong/RibosomeProfiling/Database/combined.gff',genes,'exon')

# def checkExonCdsDistance(exon,cds,genes):
#     """
#     This function checks the distance of start and stop positions between CDS and exon of genes
#     
#     * exons: str. Filename of a bed file stores exon position. Columns: ['Chr','Start','End','GeneID','None','Strand']
#     * cds: str. Filename of a bed file stores cds position. Columns: ['Chr','Start','End','GeneID','None','Strand']
#     * genes: list. A list of genes.
#     """
#     exon_df = pd.read_csv(exon,header=None,sep='\t',names=['Chr','ex_start','ex_end','GeneID','None','Strand'])
#     exon_df['GeneID'] = exon_df['GeneID'].astype(str)
#     cds_df = pd.read_csv(cds,header=None,sep='\t',names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
#     cds_df['GeneID'] = cds_df['GeneID'].astype(str)
#     for g in genes:
#         gene_exon_df = exon_df[exon_df['GeneID'].values==g]
#         gene_cds_df = cds_df[cds_df['GeneID'].values==g]
#         exon_maximum = max(gene_exon_df['ex_end'].tolist())
#         exon_minimum = min(gene_exon_df['ex_start'].tolist())
#         cds_maximum = max(gene_cds_df['cds_end'].tolist())
#         cds_minimum = min(gene_cds_df['cds_start'].tolist())
#         strand = gene_exon_df['Strand'].tolist()
#         if ('-' in strand) & ('+' in strand):
#             assert False,g + 'has both strands'
#         if '-' in strand:
#             exon_start = exon_maximum
#             exon_end = exon_minimum
#             cds_start = cds_maximum
#             cds_end = cds_minimum
#         else:
#             exon_start = exon_minimum
#             exon_end = exon_maximum
#             cds_start = cds_minimum
#             cds_end = cds_maximum
#         start_dis = abs(exon_start - cds_start)
#         end_dis = abs(exon_end - cds_end)
#         if start_dis < 16:
#             print g,'start site distance is less than 16 bases'
#         if end_dis < 16:
#             print g,'end site distance is less than 16 bases'

# exon = '/data/shangzhong/RibosomeProfiling/cho_pr/06_combined.exon.bed'
# cds = '/data/shangzhong/RibosomeProfiling/cho_pr/06_combined.CDS.bed'
# genes = ['100761614']
# checkExonCdsDistance(exon,cds,genes)   
    
#===============================================================================
#              calculate coverage
#===============================================================================
def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]
        
def posCoverage(refBed,bamFiles,otherParameters=['']):
    """
    This function calculates coverage in the bamFiles
    
    * refBed: str. Filename that has genes whose coverage would be calculated
    * bamFiles: list. A list of bam file.
    """
    a = refBed
    cmd = ''
    for bam in bamFiles:
        out = bam[:-3]+'txt'
        cmd = cmd + ('coverageBed -a {a} -b {b} -d > {out} && ').format(a=a,b=bam,out=out)
    subprocess.call(cmd[:-3],shell=True)
        
# filepath = '/data/shangzhong/RibosomeProfiling/MergeRibo/total_RNA'
# os.chdir(filepath)
# bamFiles = [f for f in os.listdir(filepath) if f.endswith('bam')]
# bamFiles = natsorted(bamFiles)
# chr_geneFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_combined.exon.bed'
# posCoverage(chr_geneFile,bamFiles)

###===================== position 5'end coverage  ===============================
def pos5Coverage(bams,batch=1):
    """
    This function calculates coverage of 5' reads at each position
    
    * bams: list. A list of bam files. Alternative: sam
    return a list of files that have 3 columns: ['count','Chr','pos']
    """
    cmd = ''
    res = []
    sub_bams = chunks(bams,batch)
    for bams in sub_bams:
        for bam in bams:
            outFile = bam[:-3] +'5endCov.txt'
            res.append(outFile)
            if bam[-3:] == 'bam':
                cmd = cmd + ('samtools view {bam} | cut -f 3,4 | uniq -c > {outFile} & ').format(
                                bam=bam,outFile=outFile)
            else:
                cmd = cmd + ('cut -f 3,4 {sam} | uniq -c > {outFile} & ').format(
                                sam=bam,outFile=outFile)
        subprocess.call(cmd + 'wait',shell=True)
    return res
# filepath = '/data/shangzhong/RibosomeProfiling/Ribo_align/23'
# os.chdir(filepath)
# bamFiles = [f for f in os.listdir(filepath) if f.endswith('sam')]
# bamFiles = natsorted(bamFiles)
# pos5Coverage(bamFiles,6)

#===============================================================================
#             gene length statistics
#===============================================================================
def getRiboCDSpos(gene_cds_df):
    """
    This function gets all cds positions in a gene. This is for coverage of A site. 
    
    * gene_cds_df: df. Pandas dataframe that have column 'cds_start','cds_end','Strand'. The df is for only one gene.
    """
    # store the cds positions
    pos = []
    for start,end,stra in zip(gene_cds_df['cds_start'],gene_cds_df['cds_end'],gene_cds_df['Strand']):
        inter = range(start+1,end+1)
        if stra == '-':
            inter.reverse()
            pos = inter + pos
        else:
            pos.extend(inter)
    strand = gene_cds_df['Strand'].tolist()
    if ('-' in strand) and ('+' in strand):
            assert False,'gene have both strand'
    elif '-' in strand:
        pos = range(pos[0]+15,pos[0],-1) + pos[:-15]
        if 0 in pos:
            index = pos.index(0)
            pos = pos[:index]
    else:
        pos = range(pos[0]-15,pos[0]) + pos[:-15]
        if 0 in pos:
            index = pos.index(0)
            pos = pos[index:]
    return pos


def getGenepos(gene_df,feature_type='cds'):
    """
    This function gets all positions for a gene feature(cds or exon).
    
    * gene_df: dataframe. In the format of Bed files, has 6 coumns ['Chr','cds/ex_start','cds/ex_end','GeneID','None','Strand']
    * feature_type: str. Either cds or exon. 
    """
    if feature_type == 'cds':
        start = 'cds_start'
        end = 'cds_end'
    else:
        start = 'ex_start'
        end = 'ex_end'
    pos = [] # store the cds positions
    for start,end,stra in zip(gene_df[start],gene_df[end],gene_df['Strand']):
        inter = range(start+1,end+1)
        if stra == '-':
            inter.reverse()
            pos = inter + pos
        else:
            pos.extend(inter)
    return pos


def geneCDSStats(cds_File,geneIDs):
    """
    This function calculates the frequencies of length of the sp genes and no_sp_genes
    return a file with two columns, sp_genes and non_sp_genes. each row represents a position.
    values are how many genes have reached the length of row index.
    
    * cds_File: str. Filename of the file that stores 
    * geneIDs: list. A list of gene ids
    """
    gene_lens =[]
    cds_df = pd.read_csv(cds_File,header=None,sep='\t',names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
    for gene in geneIDs:
        gene_cds_df = cds_df[cds_df['GeneID']==gene]
        pos = getGenepos(gene_cds_df)
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
def covAllPos(exonCovFile,cdsBedFile,geneIDs,chr_len_file,up,down):
    """
    This function calculates how many reads map to each position of genes, with
    up number of nts and down number of nts.
    
    * exon_cov_df: df. Dataframe read by pandas, with 3 columns:['coverage','Chr','pos']
    * cdsBedFile: str. Filename, with 6 columns: ['Chr','cds_start','cds_end','GeneID','None','Strand']
    * geneIDs: list. A list of gene IDs you want to include.
    * chr_len_file: dict. Stores chromosome length information.
    * up: number of nucleotides upstream of TSS sites
    * down: number of nucleotides downstream of TSE sites.
    """
    # 1) read choromosome length
    df = pd.read_csv(chr_len_file,sep='\t',header=0)
    chr_len_dict = df.set_index('Chr')['Length'].to_dict()
    # 2) read coverage file
    exon_cov_df = pd.read_csv(exonCovFile,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
    cds_df = pd.read_csv(cdsBedFile,sep='\t',header=None,names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
    #res_df = pd.DataFrame() # column: genes; rows: positions
    #3) loop for each gene
    output = exonCovFile[:-11]+'geneAllposCov.txt'
    outHandle = open(output,'w')
    for gene in geneIDs:
        # 1). get the chromosome for gene
        gene_cds_df = cds_df[cds_df['GeneID'].values==gene]
        if gene_cds_df.empty:
            outHandle.write(gene + '\t' + '0'+'\n')
            continue
        chrome = list(set(gene_cds_df['Chr'].tolist()))
        if len(chrome) > 1:
            print gene,'maps to many chromosome'
            continue
        # 2) get the coverage for the specific chromosome
        gene_cov_df = exon_cov_df[exon_cov_df['Chr'].values==chrome[0]]  # columns= [count,chr,pos]
        # 3) list all CDS position for the gene and upstream downstream bases.
        pos = getGenepos(gene_cds_df,feature_type='cds')
        strand = gene_cds_df['Strand'].tolist()
        if '-' in strand:
            pos = range(pos[0]+up,pos[0],-1) + pos + range(pos[-1]-1,pos[-1]-down,-1)
        else:
            pos = range(pos[0]-up,pos[0]) + pos + range(pos[-1]+1,pos[-1]+down)  # pos[up] is TSS, pos[-down] is TSE
        # 4) get coverage for all positions including upstream and downstream locations.
        chr_len = chr_len_dict[chrome[0]]
        cov_target = getCDSCov(gene_cov_df,pos,chr_len)
        #cov_df = pd.Series(cov_target,name=gene)
        outHandle.write('\t'.join([gene]+map(str,cov_target))+'\n')
        #res_df = res_df.join(cov_df,how='outer')
        #res_df
#     res_df.index = res_df.index-up
#     output = exonCovFile[:-11]+'geneAllposCov.txt'
#     res_df.to_csv(output,sep='\t')
    outHandle.close()
    return output


def covNearTSS_TSE(exonCovFile,cdsBedFile,geneIDs,chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down):
    """
    This function calculates how many reads map to each position around TSS and TSE sites, with
    up number of nts and down number of nts.
    
    * exon_cov_df: df. Dataframe read by pandas, with 3 columns:['coverage','Chr','pos']
    * cdsBedFile: str. Filename, with 6 columns: ['Chr','cds_start','cds_end','GeneID','None','Strand']
    * geneIDs: list. A list of gene IDs you want to include.
    * chr_len_file: str. Filename stores chromosome length information. ['Chr','Length']
    * TSS_up,TSE_up: number of nucleotides upstream of TSS and TSE sites.
    * TSS_down,TSE_up: number of nucleotides downstream of TSS and TSE sites.
    """
    # 1) read choromosome length
    df = pd.read_csv(chr_len_file,sep='\t',header=0)
    chr_len_dict = df.set_index('Chr')['Length'].to_dict()
    # upstream and downstream nts around TSS and TSE
    exon_cov_df = pd.read_csv(exonCovFile,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
    cds_df = pd.read_csv(cdsBedFile,sep='\t',header=None,names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
    df_start = pd.DataFrame() # column: genes; rows: positions
    df_end = pd.DataFrame()
    for gene in geneIDs:
        # 1). get the chromosome for gene
        gene_cds_df = cds_df[cds_df['GeneID'].values==gene]
        chrome = list(set(gene_cds_df['Chr'].tolist()))
        if len(chrome) > 1:
            print gene,'maps to many chromosome'
            continue
        # 2) get the coverage for the specific chromosome
        gene_cov_df = exon_cov_df[exon_cov_df['Chr'].values==chrome[0]]  # columns= [count,chr,pos]
        # 3) list all CDS position for the gene and upstream downstream bases.
        pos = getGenepos(gene_cds_df,feature_type='cds')
        strand = gene_cds_df['Strand'].tolist()
        if '-' in strand:
            pos = range(pos[0]+TSS_up,pos[0],-1) + pos + range(pos[-1]-1,pos[-1]-TSE_down,-1)
        else:
            pos = range(pos[0]-TSS_up,pos[0]) + pos + range(pos[-1]+1,pos[-1]+TSE_down)  # pos[up] is TSS, pos[-down] is TSE
        # 4) get coverage for all positions including upstream and downstream locations.
        chr_len = chr_len_dict[chrome[0]]
        cov_target = getCDSCov(gene_cov_df,pos,chr_len)
        # 5) get coverage around start and stop sites
        start_cov = cov_target[0:TSS_up+TSS_down]
        end_cov = cov_target[-(TSE_up+TSE_down):]
        #start_sum = [x+y for x,y in zip(start_sum,start_cov)]
        #stop_sum = [x+y for x,y in zip(stop_sum,stop_cov)]
        df_start[gene] = pd.Series(start_cov)
        df_end[gene] = pd.Series(end_cov)
    df_start = df_start.T
    df_end = df_end.T
    return df_start,df_end

def covNearSpEnd(exonCovFile,id_file,spFile,cdsFile,chr_len_file,gene_sp,up,down):
    """
    This function claculates coverage around the end of sp genes
    
    * exon_cov_df: df. Dataframe read by pandas, with 3 columns:['coverage','Chr','pos']
    * id_file: str. Filename stores all ids in the organism
    * spFile: str. Filename stores sp resutls predicted by signalP.
    * cdsBedFile: str. Filename, with 6 columns: ['Chr','cds_start','cds_end','GeneID','None','Strand']
    * geneIDs: list. A list of gene IDs you want to include.
    * chr_len_file: str. Filename stores chromosome length information. ['Chr','Length']
    * gene_sp: list. A list of gene ids.
    * up: number of nucleotides upstream of end of signal peptide.
    * down: number of nucleotides downstream of end of signal peptide.
    """
    # 1) get {gene:chr} from all ID files
    all_id = pd.read_csv(id_file,sep='\t',header=0) # GeneID    GeneSymbol    Chrom    TrAccess    PrAccess
    all_id['GeneID'] = all_id['GeneID'].astype(str)
    gene_chr_df = all_id[['GeneID','Chrom']].drop_duplicates()
    gene_chr_dic = gene_chr_df.set_index('GeneID')['Chrom']
    # 2) read choromosome length
    df = pd.read_csv(chr_len_file,sep='\t',header=0)
    chr_len_dict = df.set_index('Chr')['Length'].to_dict()
    # 3) get sp predicted results
    sp_df = pd.read_csv(spFile,header=0,sep='\t')
    sp_df['GeneID'] = sp_df['GeneID'].astype(str)
    # 4) filter sp results by gene_sp,remove gene without signal peptides
    cri = sp_df['GeneID'].map(lambda x: x in gene_sp)
    sp_df = sp_df[cri]
    # get {gene: protein} dictionary
    gene_pr_dic = {k:list(v) for k,v in sp_df.groupby('GeneID')['Pr_Access']}
    # get {pr:Ymax} dictionary
    pr_Ymax_dic = sp_df.set_index('Pr_Access')['pos.Y'] 
    # 5) add chr to sp_result
    sp_df['Chr'] = sp_df['GeneID'].apply(lambda x: gene_chr_dic[x])
    # 6) read 10_pr_cds.txt file
    cds_df = pd.read_csv(cdsFile,header=0,sep='\t',low_memory=False)
    cds_df['GeneID'] = cds_df['GeneID'].astype(str)
    cds_df['cds_start'] = cds_df['cds_start'] - 1
    cri = cds_df['GeneID'].map(lambda x: x in gene_sp)
    cds_df = cds_df[cri]
    # do analysis
    gene_res = pd.DataFrame()
    cov_df = pd.read_csv(exonCovFile,header=0,names=['coverage','Chr','pos'],delim_whitespace=True)
    for gene in gene_sp:
        pr = gene_pr_dic[gene]
        cov_near_sp = pd.DataFrame()
        for p in pr:
            pr_cds_df = cds_df[cds_df['Pr_Access'].values==p]
            pr_pos = getGenepos(pr_cds_df)
            sp_len = pr_Ymax_dic[p]*3
            sp_end_pos = pr_pos[sp_len-1]
            res_pos = range(sp_end_pos-up,sp_end_pos)+range(sp_end_pos,sp_end_pos+down)
            # get coverage of that protein
            pr_cov_df = cov_df[cov_df['Chr'].values == gene_chr_dic[gene]]
            chrome = gene_chr_dic[gene]
            chr_len = chr_len_dict[chrome]
            sp_end_cov = getCDSCov(pr_cov_df,res_pos,chr_len)
            cov_near_sp[p] = pd.Series(sp_end_cov)
        cov_near_sp = cov_near_sp.T.drop_duplicates().T
        shape = cov_near_sp.shape
        if shape[1] ==1:
            cov_near_sp.columns = [gene]
        gene_res = pd.concat([gene_res,cov_near_sp],axis=1)
    gene_res.index = range(-up,down)
    gene_res = gene_res.T
    output = exonCovFile[:-11]+'5endNearSPendCov.txt'
    gene_res.to_csv(output,sep='\t')
    return output


def normGeneCov(line,up,down,total,by='mean'):
    """
    This function normalize gene coverage by divie count at each position by total count and mean.
    The mean value is calculated after removing the 1st 15 codon and last 10 codon.
    
    * line: str. One line in a gene coverage file. In the file, each line starts with gene id and 
    * up: int. How many positions are there before the TSS sites.
    * down: int. How many positions are there after the TSE sites.
    * total: int. Total count in the sample.
    * by: str. Normalize by mean or median.
    """
    item = line[:-1].split('\t')
    gene = item[0];cov=item[1:]
    df = pd.DataFrame({gene:cov})
    df.index = df.index-up
    try:
        df = df.replace('-',np.nan)
    except:
        pass
    df = df.astype(float)
    # filter by median = 1
    endIndex = df.index[-1]
    median = (df.loc[-15:endIndex-down-14]).median().values[0]
    if median < 1: return ''
    df = df/total*(10**6)
    # filter out first 15 and last 10 codons for calculating mean
    if by == 'mean':
        mean = (df.loc[-15+45:endIndex-down-14-30]).mean().values[0]
        df = df/mean
    if by == 'median':
        median = (df.loc[-15+45:endIndex-down-14-30]).median().values[0]
        df = df/median
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

def StallSites(geneCDSposCovFile,cdsBedFile,gene_sp,total,up,down,by='median'):
    """
    This function detects all the stalling sites
    
    * geneCovFile: str. 1st column is gene id, other columns is count for each position of the gene
    * gene_sp: list. A list of genes with signal peptide
    """
    handle = open(geneCDSposCovFile,'r')
    #geneCov_df = pd.DataFrame()
    out = geneCDSposCovFile[:3]+'.stallsites.txt'
    outHandle = open(out,'w')
    outHandle.write('\t'.join(['Chr','GeneID','Chr_pos','ratio2median'])+'\n')
    for line in handle:
        # 1. normalize the coverage
        df = normGeneCov(line,up,down,total,by)
        if type(df) == str: continue  # gene cds median ribo count < 1
        # 2. get ratio to median and cut the lenth, only retain the 
        gene = df.columns[0]
        ratio = df[gene].tolist() # sequence with count/median
        ratio = ratio[up-15:-down-14]
        # 3. get Chromosome pos for the gene
        cds_df = pd.read_csv(cdsBedFile,header=None,sep='\t',names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
        gene_cds_df = cds_df[cds_df['GeneID'].values==gene]
        gene_cds_pos = getGenepos(gene_cds_df,feature_type='cds')
        # 4. output the stalling sites. ['Chr','GeneID','Chr_pos','ratio']
        chrome = gene_cds_df['Chr'].tolist()[0]
        indexes = [i for i in range(len(ratio)) if ratio[i] > 25]  # index of the stalling sites in the ratios
        if indexes == []: continue # no pausing sites in this gene
        for index in indexes:
            outHandle.write('\t'.join([chrome,gene,str(gene_cds_pos[index]),str(ratio[index])]) + '\n')
    # output
    outHandle.close()
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
#         coverage at gene levle
#===============================================================================
def GeneExpressLevel(outFile,genePosCovFiles,chrPosCovFiles,up,down,seqType='ribo',expressType='rpkm'):
    """
    This function calculates expression levels for each gene, can be used for ribo seq data or
    mRNA seq data.
    
    * genePosCovFiles: list. A list of replicate files which store coverage for each gene at each position.
                        In each line, the first column is gene id, the rest are coverage at each position.
    * chrPosCovFiles: list. A list of replicate files which store coverage for chromosome. Has 3 columns,
                            ['coverage','Chr','pos']
    * up,down: int. Number of nucletides in upstream of start position, downstream of stop position of gene.
    * seqType: str. Can be 'ribo' or 'rna'
    * expressType: str. Can be 'rpkm' or 'raw'
    """
    gene_count_df = pd.DataFrame()
    for f,t in zip(genePosCovFiles,chrPosCovFiles):
        # get total count
        df = pd.read_csv(t,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
        total = df['coverage'].sum()
        # calculate coverage for each line
        geneIDs = []
        handle = open(f,'r')
        counts = []
        for line in handle:
            item = line[:-1].split('\t')
            gene = item[0]
            geneIDs.append(gene)
            cov = item[1:]
            for i in range(len(cov)):
                if cov[i]=='-': cov[i] = 0
                else: cov[i] = int(cov[i])
            if seqType == 'ribo':
                count = sum(cov[up-15:-down-15])
                if expressType == 'rpkm':
                    count = count/total/len(cov[up:-down+1])*(10**9)
            elif seqType == 'rna':
                count = sum(cov[up:-down])
                if expressType == 'rpkm':
                    count = count/total/len(cov)*(10**9)
            counts.append(count)
        handle.close()
        column = f.split('/')
        gene_count_df[column[-1][:3]] = pd.Series(counts)
    gene_count_df.insert(0,'GeneID',geneIDs)
    gene_count_df.to_csv(outFile,sep='\t',index=False)

#===============================================================================
#             motif stalling sites
#===============================================================================
def posLongestPr(gene,cds_pr_df,gene_pr_dict):
    """
    This function gets the positions of the longest transcripts of a gene
    
    * gene: str. gene id
    * cds_pr_df: dataframe. ['Chr','cds_start','cds_end','GeneID','Pr_Access','Strand']
    * gene_pr_dict: dict. A dictionary in the format: {gene:[pr1,pr2]}
    """
    cds_pr_df['len'] = cds_pr_df['cds_end'] - cds_pr_df['cds_start']
    proteins = gene_pr_dict[gene]
    length = []; position = []
    for pr in proteins:
        # get cds positions
        gene_pr_df = cds_pr_df[cds_pr_df['Pr_Access'].values==pr]
        length.append(gene_pr_df['len'].sum())
        position.append(getGenepos(gene_pr_df,feature_type='cds'))
        # get chromosome
        chrome = list(set(gene_pr_df['Chr'].tolist()))
    # get the longest cds
    max_len = max(length)
    index = length.index(max_len)
    pos = position[index]
    return pos,chrome,proteins[index]

def codonStallSites(pos,resolution='codon',by='median'):
    """
    This function normalize the coverage for detecting the stall sites
    """
    start = 15; end = 10
    if resolution != 'codon':
        start=start*3;end=end*3
    if by == 'median':
        median = np.median(np.array(pos[start:-end]))
    else:
        median = np.mean(np.array(pos[start:-end]))
    if median == 0:
        return 0
    norm = [f/median for f in pos]
    return norm

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
        codon = line[-4:-1]
        aa = Seq(codon,generic_dna).translate()
        if codon in codon_num_dic:
            codon_num_dic[codon] = codon_num_dic[codon] + 1
        else:
            codon_num_dic[codon] = 0
        if aa in AA_num_dic:
            AA_num_dic[aa] = AA_num_dic[aa] + 1
        else:
            AA_num_dic[aa] = 0
    for key in codon_num_dic:
        codon_num_dic[key] =  codon_num_dic[key]/total
    for key in AA_num_dic:
        AA_num_dic[key] = AA_num_dic[key]/total
    return codon_num_dic,AA_num_dic