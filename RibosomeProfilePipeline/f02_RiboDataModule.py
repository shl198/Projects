from __future__ import division
import subprocess,os,sys
import pandas as pd
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('ggplot')


def signalP(inputFa):
    """
    This function runs signalp for files
    
    * inputFa: str. Input fasta file that stores protein sequence
    """
    outFile = inputFa[:-3] + '_sp.txt'
    cmd = ('signalp {input} > {out}').format(input=inputFa,out=outFile)
    subprocess.call(cmd.split(' '))
    
def removeIDVersion(ID):
    """
    This function remove version information in all kinds of id versions
    
    * ID: str. should be genome mRNA or protein accession number with version information. eg: NW0123.1
    """
    if '.' in ID:
        return ID[:ID.index('.')]
    else:
        return ID
    
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

# inputFile = '/data/shangzhong/RibosomeProfiling/cho_pr/03_choPrRefseq_sp.gene.txt'
# gffMapFile = '/data/shangzhong/Database/gff_chok1_all_ID.txt'
# genesMap2diffChrom(inputFile,gffMapFile)

#===============================================================================
#             classify gene ids regarding signal peptide
#===============================================================================
def gene_sp_classify(sp_gene):
    """
    This function classfies which gene encodes protein only with signal peptide,
    which gene encodes protein only without signal peptide and which gene encods both.
    
    * sp_gene: str. Filename of signalP prediction results with gene id inserted with column name 'GeneID'.
    """
    df = pd.read_csv(sp_gene,header=0,sep='\t')
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
    outFile = sp_gene[:-3]+'classify.txt'
    df.to_csv(outFile,sep='\t',index=False)
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
        try:
            # get gene id
            geneid = ''
            index = line.index('GeneID:')
            index = index + 7
            while line[index] !=';' and line[index] != ',':
                    geneid = geneid + line[index]
                    index = index + 1
        except:
            geneid = '-'
            print line, 'dont have geneid'
        if geneid in genes:
            outline = '\t'.join([item[0],str(int(item[3])),item[4],geneid,'0',item[6]]) + '\n'
            out_handle.write(outline)
            #genes.remove(geneid)
    gff_handle.close()
    out_handle.close()
    df = pd.read_csv(outBed,sep='\t',header=None,names=['chr','start','stop','geneid','None','Strand'])
    df = df.drop_duplicates()
    #df['chr_id'] = df['chr'].map(str) + '_' + df['geneid'].astype(str)
    #res_df = df[['chr','start','stop','chr_id']]
    df.to_csv(outBed,sep='\t',header=None,index=False)
    finalBed = outBed[:-3] + feature + '.bed'
    cmd = ('sort -k1,1 -k2,2n {input} | bedtools merge -i stdin -s -c 4,5,6 '
           '-o distinct,distinct,distinct > {out}').format(input=outBed,out=finalBed)
    subprocess.call(cmd,shell=True)
    #os.remove(outBed)
    df = pd.read_csv(finalBed,sep='\t',header=None,names=['chr','start','stop','geneid','None','Strand'])
    df['start']=df['start']-1
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
def pos5Coverage(bams):
    """
    This function calculates coverage of 5' reads at each position
    
    * bams: list. A list of bam files
    return a list of files that have 3 columns: ['count','Chr','pos']
    """
    cmd = ''
    res = []
    for bam in bams:
        outFile = bam[:-3] +'5endCov.txt'
        res.append(outFile)
        cmd = cmd + ('samtools view {bam} | cut -f 3,4 | uniq -c > {outFile} & ').format(
                        bam=bam,outFile=outFile)
    subprocess.call(cmd[:-3],shell=True)
    return res
# filepath = '/data/shangzhong/RibosomeProfiling/MergeRibo/total_RNA'
# os.chdir(filepath)
# bamFiles = [f for f in os.listdir(filepath) if f.endswith('bam')]
# bamFiles = natsorted(bamFiles)
# pos5Coverage(bamFiles)

#===============================================================================
#             gene length statistics
#===============================================================================
def getGeneCDSpos(gene_cds_df,gene):
    pos = [] # store the cds positions
    for start,end,stra in zip(gene_cds_df['cds_start'],gene_cds_df['cds_end'],gene_cds_df['Strand']):
        inter = range(start+1,end+1)
        if stra == '-':
            inter.reverse()
            pos = inter + pos
        else:
            pos.extend(inter)
    return pos

def getGeneExonpos(gene_exon_df,gene):
    pos = [] # store the cds positions
    for start,end,stra in zip(gene_exon_df['ex_start'],gene_exon_df['exon_end'],gene_exon_df['Strand']):
        inter = range(start+1,end+1)
        if stra == '-':
            inter.reverse()
            pos = inter + pos
        else:
            pos.extend(inter)
    return pos

def geneCDSStats(cds_File,geneIDs):
    """
    """
    gene_lens =[]
    cds_df = pd.read_csv(cds_File,header=0,sep='\t',names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
    for gene in geneIDs:
        gene_cds_df = cds_df[cds_df['GeneID']==gene]
        pos = getGeneCDSpos(gene_cds_df,gene)
        gene_lens.append(len(pos))
    len_stats = [0]*max(gene_lens)
    for g in gene_lens:
        for i in range(g):
            len_stats[i] = len_stats[i]+1
    percent = [x/len(geneIDs) for x in len_stats]
    return len_stats


#===============================================================================
#        divide gene into 100 bin and calculate coverage in each bin
#===============================================================================
def percent_cov(exon_cov_df,cds_df,geneIDs):
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


def percent5endCov(exon_cov_df,cds_df,geneIDs):
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
#             coverage near TSS and TSE sites
#===============================================================================
def covNearTSS_TSE(exon_cov_df,cds_df,geneIDs,TSS_up,TSS_down,TSE_up,TSE_down):
    """
    This function calculates how many reads map to each position around TSS and TSE sites, with
    up number of nts and down number of nts.
    
    * exon_cov_df: df. Dataframe read by pandas, with 3 columns:['coverage','Chr','pos']
    * cds_df: df. Dataframe read from bed file by pandas, with 6 columns: ['Chr','cds_start','cds_end','GeneID','None','Strand']
    * GeneIDs: list. A list of gene IDs you want to include.
    * up: number of nucleotides upstream of TSS and TSE sites.
    * down: number of nucleotides downstream of TSS and TSE sites.
    """
#     up = 50; down=50  # upstream and downstream nts around TSS and TSE
    start_sum = [0]*(TSS_up+TSE_down); stop_sum = [0]*(TSE_up+TSE_down) # start stores coverage around TSS site, stop for TSE site.
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
        pos = [] # store the cds positions
        for start,end,stra in zip(gene_cds_df['cds_start'],gene_cds_df['cds_end'],gene_cds_df['Strand']):
            inter = range(start+1,end+1)
            if stra == '-':
                inter.reverse()
                pos = inter + pos
            else:
                pos.extend(inter)
        #TSS = pos[0]; TSE = [-1]
        strand = gene_cds_df['Strand'].tolist()
        if '-' in strand:
            pos = range(pos[0]+TSS_up,pos[0],-1) + pos + range(pos[-1]-1,pos[-1]-TSE_down,-1)
        else:
            pos = range(pos[0]-TSS_up,pos[0]) + pos + range(pos[-1]+1,pos[-1]+TSE_down)  # pos[up] is TSS, pos[-down] is TSE
        # 4) get coverage for all positions including upstream and downstream locations.
        cov_target = []
        dic = gene_cov_df.set_index('pos')['coverage'].to_dict()
        for p in pos:
            if p in dic:
                cov_target.append(dic[p])
            else:
                cov_target.append(0)
        # 5) get coverage around start and stop sites
        start_cov = cov_target[0:TSS_up+TSS_down]
        stop_cov = cov_target[-(TSE_up+TSE_down):]
        start_sum = [x+y for x,y in zip(start_sum,start_cov)]
        stop_sum = [x+y for x,y in zip(stop_sum,stop_cov)]
    return start_sum,stop_sum




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



