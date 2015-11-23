import pandas as pd
from f02_RiboDataModule import *
#===============================================================================
#         calculate coverage of house keeping genes at each position
#===============================================================================
# ========== 1. analyze ribo data =================
"""
# 1) read cds file
cdsFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_combined.CDS.bed'
cds_df = pd.read_csv(cdsFile,sep='\t',header=None,names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
# 2) read cov file
path = '/data/shangzhong/RibosomeProfiling/Ribo_align'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('5endCov.txt')]
coverFiles = natsorted(coverFiles)
# 3) define genes
genes = ['heavychain','lightchain','100689477','100765489','100689409','100736557','100757362','100689395',
         '100689467','100754509','100756201','100769768','100757453','100750732']
genes_label = ['heavychain','lightchain','Actb','Actr5','B2m','Gapdh','Gusb','Tfrc','Pgk1','Hsp90ab1','Rplp0','Hprt1','Nono','Sdha']
# 4) loop for each file
gene_ribo_cov_df = pd.DataFrame()
for f in coverFiles:
    gene_cov = []
    cov_df = pd.read_csv(f,header=0,names=['coverage','Chr','pos'],delim_whitespace=True)
    total = cov_df['coverage'].sum()
    for gene in genes:
        gene_cds_df = cds_df[cds_df['GeneID'].values==gene]
        chrome = list(set(gene_cds_df['Chr'].tolist()))
        if len(chrome) > 1:
            assert False,'gene has multi chromosomes'
        pos = getRiboCDSpos(gene_cds_df)
        gene_cov_df = cov_df[cov_df['Chr'].values==chrome[0]]
        cov = getCDSCov(gene_cov_df,pos)
        gene_count = sum(cov)
        gene_rpkm = gene_count*(10**9)/total/len(pos)
        gene_cov.append(gene_rpkm)
    gene_ribo_cov_df[f[:3]] = pd.Series(gene_cov)
gene_ribo_cov_df.index = genes
gene_ribo_cov_df['day3'] = gene_ribo_cov_df.iloc[:,0:3].mean(axis = 1)
gene_ribo_cov_df['day6'] = gene_ribo_cov_df.iloc[:,3:6].mean(axis = 1)
gene_ribo_cov_df['day3_heavy'] = gene_ribo_cov_df.loc['heavychain','day3']/gene_ribo_cov_df['day3']
gene_ribo_cov_df['day3_light'] = gene_ribo_cov_df.loc['lightchain','day3']/gene_ribo_cov_df['day3']
gene_ribo_cov_df['day6_heavy'] = gene_ribo_cov_df.loc['heavychain','day6']/gene_ribo_cov_df['day6']
gene_ribo_cov_df['day6_light'] = gene_ribo_cov_df.loc['lightchain','day6']/gene_ribo_cov_df['day6']
ribo_df = gene_ribo_cov_df.iloc[:,-4:]
ribo_df.index = genes_label
ribo_df = ribo_df[2:]
ribo_df = np.log2(ribo_df)
ax = ribo_df.plot(kind='bar',title='ribosome log2(Ab/house keeping gene)')
ax.set_ylabel('log2 ratio')
#=========== 2. analyze total mRNA =========
# 1) read exon file
exonFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_combined.exon.bed'
exon_df = pd.read_csv(exonFile,sep='\t',header=None,names=['Chr','ex_start','ex_end','GeneID','None','Strand'])
# 2) read cov file
path = '/data/shangzhong/RibosomeProfiling/Ribo_align/total_RNA'
os.chdir(path)
coverFiles = [f for f in os.listdir(path) if f.endswith('5endCov.txt')]
coverFiles = natsorted(coverFiles)
# 4) loop for each file
gene_rna_cov_df = pd.DataFrame()
for f in coverFiles:
    gene_cov = []
    cov_df = pd.read_csv(f,header=0,names=['coverage','Chr','pos'],delim_whitespace=True)
    total = cov_df['coverage'].sum()
    for gene in genes:
        gene_exon_df = exon_df[exon_df['GeneID'].values==gene]
        chrome = list(set(gene_exon_df['Chr'].tolist()))
        if len(chrome) > 1:
            assert False,'gene has multi chromosomes'
        pos = getGenepos(gene_exon_df)
        gene_cov_df = cov_df[cov_df['Chr'].values==chrome[0]]
        cov = getCDSCov(gene_cov_df,pos)
        gene_count = sum(cov)
        gene_rpkm = gene_count*(10**9)/total/len(pos)
        gene_cov.append(gene_rpkm)
    gene_rna_cov_df[f[:3]] = pd.Series(gene_cov)
gene_rna_cov_df.index = genes
gene_rna_cov_df['day3'] = gene_rna_cov_df.iloc[:,0:3].mean(axis = 1)
gene_rna_cov_df['day6'] = gene_rna_cov_df.iloc[:,3:6].mean(axis = 1)
gene_rna_cov_df['day3_heavy'] = gene_rna_cov_df.loc['heavychain','day3']/gene_rna_cov_df['day3']
gene_rna_cov_df['day3_light'] = gene_rna_cov_df.loc['lightchain','day3']/gene_rna_cov_df['day3']
gene_rna_cov_df['day6_heavy'] = gene_rna_cov_df.loc['heavychain','day6']/gene_rna_cov_df['day6']
gene_rna_cov_df['day6_light'] = gene_rna_cov_df.loc['lightchain','day6']/gene_rna_cov_df['day6']
rna_df = gene_rna_cov_df.iloc[:,-4:]
rna_df.index = genes_label
rna_df = rna_df[2:]
rna_df = np.log2(rna_df)
ax = rna_df.plot(kind='bar',title='total rna log2(Ab/house keeping gene)')
ax.set_ylabel('log2 ratio')
# 3. ====== translation efficiency =======
trans_df = pd.DataFrame((gene_ribo_cov_df[['day3','day6']].values) / (gene_rna_cov_df[['day3','day6']].values),columns=['day3','day6'],index=gene_rna_cov_df.index)
trans_df['day3_heavy'] = trans_df.loc['heavychain','day3']/trans_df['day3']
trans_df['day3_light'] = trans_df.loc['lightchain','day3']/trans_df['day3']
trans_df['day6_heavy'] = trans_df.loc['heavychain','day6']/trans_df['day6']
trans_df['day6_light'] = trans_df.loc['lightchain','day6']/trans_df['day6']
trans_df = trans_df.iloc[:,-4:]
trans_df.index = genes_label
trans_df = trans_df[2:]
ax = trans_df.plot(kind='bar',title='translation efficiency (Ab/house keeping gene)')
ax.set_ylabel('ratio')
plt.show()
print 'done'
"""

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
        # calculate coverage for each gene(line)
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
                count = sum(cov[up-15:-down-14])  # offsite for A site
                if expressType == 'rpkm':
                    count = count/total/len(cov[up:-down+1])*(10**9)
            elif seqType == 'rna':
                count = sum(cov)
                if expressType == 'rpkm':
                    count = count/total/len(cov)*(10**9)
            counts.append(count)
        handle.close()
        column = f.split('/')
        gene_count_df[column[-1][:3]] = pd.Series(counts)
    gene_count_df.insert(0,'GeneID',geneIDs)
    gene_count_df.to_csv(outFile,sep='\t',index=False)


def StallSites(geneCDSposCovFile,cdsBedFile,gene_sp,up,down,by='median'):
    """
    This function detects all the stalling sites
    
    * geneCovFile: str. 1st column is gene id, other columns is count for each position of the gene
    * cdsBedFile: str. cds bed file name. columns are ['Chr','cds_start','cds_end','GeneID','None','Strand']
    * gene_sp: list. A list of genes with signal peptide.
    * total: int. total number of reads mapping to all genes.
    * up: int. Number of additional base pairs upstream of TSS site.
    * down: int. Number of additional base pairs downstream of TSE site.
    """
    handle = open(geneCDSposCovFile,'r')
    #geneCov_df = pd.DataFrame()
    out = geneCDSposCovFile[:3]+'.stallsites.txt'
    outHandle = open(out,'w')
    outHandle.write('\t'.join(['Chr','GeneID','Chr_pos','ratio2median','Gene_pos'])+'\n')
    for line in handle:
        # 1. normalize the coverage
        df = normGeneCov(line,up,down,by)
        if type(df) == str: continue  # gene cds median ribo count < 1
        # 2. get ratio to median and cut the lenth, only retain the cds reagion
        gene = df.columns[0]  # get gene id
        ratio = df[gene].tolist() # sequence with count/median
        ratio = ratio[up-15:-down-14]
        # 3. get Chromosome pos for the gene
        cds_df = pd.read_csv(cdsBedFile,header=None,sep='\t',names=['Chr','cds_start','cds_end','GeneID','None','Strand'])
        gene_cds_df = cds_df[cds_df['GeneID'].values==gene]
        gene_cds_pos = getGenepos(gene_cds_df,feature_type='cds')
        # 4. output the stalling sites. ['Chr','GeneID','Chr_pos','ratio']
        chrome = gene_cds_df['Chr'].tolist()[0]
        indexes = [i for i in range(46,len(ratio)-30) if ratio[i] > 25]  # index of the stalling sites relative to TSS site
        if indexes == []: continue # no pausing sites in this gene
        for index in indexes:
            outHandle.write('\t'.join([chrome,gene,str(gene_cds_pos[index]),str(ratio[index]),str(index-15)]) + '\n')
    # output
    outHandle.close()
    

# #===============================================================================
# #                         5. generate bed file for CDS
# #===============================================================================
# choPrRefseq_sp_gene = signalP_path + '/cho_pr/03_cho_ab_pr_sp.gene.txt'
# df = pd.read_csv(choPrRefseq_sp_gene,header=0,sep='\t',low_memory=False)
# genes = df['GeneID'].astype(str).tolist()
# genes = list(set(genes))
# # exonBed = choPrRefseq_sp_gene[:-4]+'.bed'
# # exonBedFile = gff2Bed4Genes(gffFile,genes,'exon',exonBed)
# # exonMergeBed = mergeBed(exonBedFile,feature='exon')
# # outBed = changeFileName(exonMergeBed,3)                                 # 06_cho_ab_pr_sp.gene.exon.bed
# # os.remove(exonBedFile)
# cdsBed = choPrRefseq_sp_gene[:-4] + '.bed'
# cdsBedFile = gff2Bed4Genes(gffFile,genes,'CDS',cdsBed)
# cdsMergeBed = mergeBed(cdsBedFile,feature='CDS')
# outBed = changeFileName(cdsMergeBed,3)                                  # 06_cho_ab_pr_sp.gene.CDS.bed
# # get 07_overlap_CDS.txt
# overlap_CDS = overlapCDS(cdsBedFile,feature='CDS')
# overlapFile = changeFileName(overlap_CDS,4)  
# os.remove(cdsBedFile)
# # two genes encode two proteins at different strand direction.

# # # # # # # # # # # # #===============================================================================
# # # # # # # # # # # # #                         8. get cover data for mapping bam files, should delete
# # # # # # # # # # # # #===============================================================================
# # # # # # # # # # # # path = '/data/shangzhong/RibosomeProfiling/Ribo_align'
# # # # # # # # # # # # os.chdir(path)
# # # # # # # # # # # # bamFiles = [f for f in os.listdir(path) if f.endswith('.sort.bam')]
# # # # # # # # # # # # cov5end = pos5Coverage(bamFiles,batch=6)
# # # # # # # # # # # # # move to a folder
# # # # # # # # # # # # folder = '01_cov5end'
# # # # # # # # # # # # if not os.path.exists(folder):
# # # # # # # # # # # #     os.mkdir(folder)
# # # # # # # # # # # # for f in cov5end:
# # # # # # # # # # # #     cmd = ('mv {input} {out}').format(input=f,out=folder)
# # # # # # # # # # # #     subprocess.call(cmd.split(' '))                                                    # folder 01_cov5end

# # # # # # # # # # # # #===============================================================================
# # # # # # # # # # # # #                         10. calculate coverage around TSS and TSE sites obsolete
# # # # # # # # # # # # #===============================================================================
# # # # # # # # # # # # def wrap_covNearTSS_TSE(target_path,exonCovFile,cdsBedFile,genes,chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down,gene_type='sp'):
# # # # # # # # # # # #     """
# # # # # # # # # # # #     This function is a wrapper of covNearTSS_TSE
# # # # # # # # # # # #     """
# # # # # # # # # # # #     #index = exonCovFile.rindex('/')
# # # # # # # # # # # #     f = exonCovFile#[index+1:]
# # # # # # # # # # # #     df_start,df_end = covNearTSS_TSE(exonCovFile,cdsBedFile,genes,chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down)
# # # # # # # # # # # #     if gene_type == 'sp':
# # # # # # # # # # # #         startFile = 'start_sp_' + f[:-16]+'txt';endFile='end_sp_'+f[:-16]+'txt'
# # # # # # # # # # # #     else:
# # # # # # # # # # # #         startFile = 'start_nosp_' + f[:-16]+'txt';endFile='end_nosp_'+f[:-16]+'txt'
# # # # # # # # # # # #     df_start.to_csv(startFile,sep='\t')
# # # # # # # # # # # #     df_end.to_csv(endFile,sep='\t')
# # # # # # # # # # # #     cmd = ('mv {old} {new}').format(old=startFile,new=target_path)
# # # # # # # # # # # #     subprocess.call(cmd.split(' '))
# # # # # # # # # # # #     cmd = ('mv {old} {new}').format(old=endFile,new=target_path)
# # # # # # # # # # # #     subprocess.call(cmd.split(' '))
# # # # # # # # # # # # """
# # # # # # # # # # # # # 1).read chromosome length
# # # # # # # # # # # # chr_len_file = '/data/shangzhong/RibosomeProfiling/Database/combined_ChrLen.txt'
# # # # # # # # # # # # # 2).sp genes and no sp genes
# # # # # # # # # # # # gene_class_file = '/data/shangzhong/RibosomeProfiling/cho_pr/08_sp_no_sp_genes.txt'
# # # # # # # # # # # # df = pd.read_csv(gene_class_file,sep='\t',header=0)
# # # # # # # # # # # # gene_sp = df['gene_sp'].dropna()
# # # # # # # # # # # # gene_no_sp = df['gene_no_sp'].dropna()
# # # # # # # # # # # # # 3).read cds position file
# # # # # # # # # # # # cdsBedFile = '/data/shangzhong/RibosomeProfiling/cho_pr/06_cho_ab_pr_sp.gene.CDS.bed'
# # # # # # # # # # # # # 4). define target folder
# # # # # # # # # # # # target_path = '/data/shangzhong/RibosomeProfiling/Ribo_align/02_TSS_TSE_cov'
# # # # # # # # # # # # if not os.path.exists(target_path):
# # # # # # # # # # # #     os.mkdir(target_path)
# # # # # # # # # # # # # 5).define results
# # # # # # # # # # # # res_sp = pd.DataFrame()
# # # # # # # # # # # # res_no_sp = pd.DataFrame()
# # # # # # # # # # # # TSS_up = 50;TSS_down= 50
# # # # # # # # # # # # TSE_up = 50;TSE_down= 50
# # # # # # # # # # # # # 6).parallel
# # # # # # # # # # # # path = '/data/shangzhong/RibosomeProfiling/Ribo_align/01_cov5end'
# # # # # # # # # # # # os.chdir(path)
# # # # # # # # # # # # coverFiles = [f for f in os.listdir(path) if f.endswith('sort.5endCov.txt')]
# # # # # # # # # # # # coverFiles = natsorted(coverFiles)
# # # # # # # # # # # # #covNearTSS_TSE(coverFiles[0],cdsBedFile,['100689328'],chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down)
# # # # # # # # # # # # proc = []
# # # # # # # # # # # # for f in coverFiles:
# # # # # # # # # # # #     p1 = Process(target=wrap_covNearTSS_TSE,args=(target_path,f,cdsBedFile,gene_sp,chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down,'sp'))    
# # # # # # # # # # # #     p1.start()
# # # # # # # # # # # #     proc.append(p1)
# # # # # # # # # # # #     p2 = Process(target=wrap_covNearTSS_TSE,args=(target_path,f,cdsBedFile,gene_no_sp,chr_len_file,TSS_up,TSS_down,TSE_up,TSE_down,'nosp'))    
# # # # # # # # # # # #     p2.start()
# # # # # # # # # # # #     proc.append(p2)
# # # # # # # # # # # # for p in proc:
# # # # # # # # # # # #     p.join()
# # # # # # # # # # # # """

def RiboGeneCovAllPos(chrCovFile,cdsFile,geneIDs,chr_len_file):
    """
    This function calculates how many reads map to each gene for RiboSeq. Using 15 base pairs to offset.
    
    * exon_cov_df: df. Dataframe read by pandas, with 3 columns:['coverage','Chr','pos']
    * cdsBedFile: str. Filename, with 6 columns: ['Chr','cds_start','cds_end','GeneID','None','Strand']
    * geneIDs: list. A list of gene IDs you want to include.
    * chr_len_file: dict. Stores chromosome length information. ['Chr','Length']
    """
    # 1) read choromosome length
    df = pd.read_csv(chr_len_file,sep='\t',header=0)
    chr_len_dict = df.set_index('Chr')['Length'].to_dict()
    # 2) read coverage file
    exon_cov_df = pd.read_csv(chrCovFile,header=None,names=['coverage','Chr','pos'],delim_whitespace=True)
    cds_df = pd.read_csv(cdsFile,sep='\t',header=None,names=['Chr','cds_start','cds_end','GeneID','Pr_Access','Strand'],low_memory=False)
    # 3) loop for each gene
    output = chrCovFile[:-11]+'geneRawCount.txt'
    outHandle = open(output,'w')
    outHandle.write('\t'.join(['GeneID','rawCount','length'])+'\n')
    for gene in geneIDs:
        # 1). get the chromosome for gene
        gene_cds_df = cds_df[cds_df['GeneID'].values==gene]
        if gene_cds_df.empty: # gene doesn't encode protein
            outHandle.write(gene + '\t' + '0'+'\n')  
            continue
        chrome = list(set(gene_cds_df['Chr'].tolist()))
        if len(chrome) > 1:
            print gene,'maps to many chromosome'
            continue
        # 2) get the coverage for the specific chromosome
        gene_cov_df = exon_cov_df[exon_cov_df['Chr'].values==chrome[0]]  # columns= [count,chr,pos]
        # 3) list all CDS position for the gene and upstream downstream bases.
        proteins = list(set(gene_cds_df['Pr_Access'].tolist()))
        gene_pos = []
        for pr in proteins:
            pr_cds_df = gene_cds_df[gene_cds_df['Pr_Access'].values==pr]  # get protein df
            pos = getGenepos(pr_cds_df,feature_type='cds')                # get gene postion
            strand = pr_cds_df['Strand'].tolist()
            if '-' in strand:
                pos = range(pos[0]+15,pos[0],-1) + pos[:-15]
            else:
                pos = range(pos[0]-15,pos[0]) + pos[:-15]
            gene_pos.extend(pos)
        gene_pos = list(set(gene_pos))    # remove '-' in gene position
        try:
            gene_pos.remove('-')
        except:
            pass
        # 4) get coverage for all positions including upstream and downstream locations.
        chr_len = chr_len_dict[chrome[0]]
        cov_target = getCDSCov(gene_cov_df,gene_pos,chr_len)
        cov_target = [p for p in cov_target if p != '-']
        count = str(sum(map(int,cov_target)))
        length = str(len(cov_target))
        outHandle.write('\t'.join([gene,count,length])+'\n')
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
        df_start[gene] = pd.Series(start_cov)
        df_end[gene] = pd.Series(end_cov)
    df_start = df_start.T
    df_end = df_end.T
    return df_start,df_end

def covNearSpEnd(exonCovFile,id_file,spFile,cdsFile,chr_len_file,gene_sp,up,down,level='nt'):
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
    gene_res = pd.DataFrame()  # row is geneid, column is position, values are read count
    cov_df = pd.read_csv(exonCovFile,header=0,names=['coverage','Chr','pos'],delim_whitespace=True)
    for gene in gene_sp:
        pr = gene_pr_dic[gene]
        cov_near_sp = pd.DataFrame()  # one dimension
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
            if level == 'pr':
                sp_end_cov_new = chunks(sp_end_cov,3)
                sp_end_cov = []
                for i in range(len(sp_end_cov_new)):
                    if '-' in sp_end_cov_new[i]:
                        sp_end_cov.append('-')
                    else:
                        sp_end_cov.append(sum(sp_end_cov_new[i]))
            cov_near_sp[p] = pd.Series(sp_end_cov)
        cov_near_sp = cov_near_sp.T.drop_duplicates().T
        shape = cov_near_sp.shape
        if shape[1] ==1:                    # if a gene has multiple proteins, it keeps the protein accession number
            cov_near_sp.columns = [gene]
        gene_res = pd.concat([gene_res,cov_near_sp],axis=1)
    #gene_res.index = range(-up,down)
    gene_res = gene_res.T
    output = exonCovFile[:-11]+'5endNearSPendCov.txt'
    gene_res.to_csv(output,sep='\t')
    return output