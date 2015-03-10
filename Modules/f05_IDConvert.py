"""
This file includes functions of all kinds of ID converts, such as gene symbols to gene ids, protein ids to gene ids...
"""
import pandas as pd
import os
import subprocess
from Bio import Entrez
Entrez.email = "shl198@eng.ucsd.edu"
#===============================================================================
#  This part is for DE analysis. Convert geneSymbols to Entrez genes in the given files
#===============================================================================

def geneSymbol2EntrezID(Dict,output_path,inputpath,sym2ID='yes'):
    """
    This fuction converts gene symbol in htseq-count file to entriz_gene_ID, and ID to symbol
    
    * Dict: Convert file. 1st column is ID, 2nd column is symbol
    * output_path: a folder stores all the files end with sort.Count.txt 
    * inputpath: a folder stores all the htseq count files.
    * sym2ID: if yes: symbol to ID
              if no : ID to symbol
    """
    allFiles = os.listdir(inputpath)
    filelist = [f for f in allFiles if f.endswith('.sort.txt')]
    #=========== build dictionary ============
    dic = {}
    result = open(Dict,'r')
    if sym2ID =='yes':
        for line in result:
            item = line.split('\t')
            dic[item[1][:-1]] = item[0]
    else:
        for line in result:
            item = line.split('\t')
            dic[item[0]] = item[1][:-1] 
    for filename in filelist:
        outputName = filename[:-3] + 'Count.txt'
        result = open(inputpath + '/' + filename,'r')
        output = open(output_path + '/' + outputName,'w')
        for line in result:
            item = line.split('\t')
            try:
                item[0] = dic[item[0]]
                output.write('\t'.join(item))
            except:
                output.write('\t'.join(item))
                print item[0]
        result.close()
        output.close()


# def addEntrezGeneID2CufflinkResultWithEnsemblAnnotation(ConvertFile,geneFpkmFile):
#     """
#     This function insert entrez gene ids as the 1st column in the convertfile.
#     The cufflink results are got based on ensemble annotation file
#     
#     * ConvertFile: filename. has 2 columsn, 1st is gene id, 2nd is ensembl id.
#     * geneFpkmFile: the result file genes.fpkm_tracking returned by cufflinks 
#     """
#     mappedIDs = []
#     mapData = pd.read_csv(ConvertFile,sep='\t',names =['GeneID','EnsemblID'],header=0)
#     cuffData = pd.read_csv(geneFpkmFile,sep='\t')
#     ensemblGenes = mapData['EnsemblID'].tolist()
#     geneIDs = mapData['GeneID'].tolist()
#     for gene in cuffData['gene_id']:
#         index = gene.index('.')
#         ensembl_gene = gene[:index]
#         #geneID = mapData.ix[mapData["ensemblID"]==ensembl_gene,'geneID']
#         try:
#             index = ensemblGenes.index(ensembl_gene)
#             geneID = geneIDs[index]
#             mappedIDs.append(geneID)
#         except:
#             mappedIDs.append('-')
#     cuffData.insert(0, 'EntrezGeneID', pd.Series(mappedIDs) )
#     
#     outputFile = geneFpkmFile + '.csv'
#     cuffData.to_csv(outputFile)
#     print 'done'


def geneSymbol2ID4Cufflink(ConvertFile,geneFpkmFile):
    """
    This function convert the gene symbols to gene id in cufflinks results
    and then return a file with two columns:1st Entrez geneID, 2nd FPKM
    """
    cuffData = pd.read_csv(geneFpkmFile,sep='\t')
    mapData = pd.read_csv(ConvertFile,sep='\t',names =['GeneID','GeneSymbol','Accession'],header=None)
    
    outputFile = geneFpkmFile + '.txt'
    with open(outputFile,'w') as f:
        for symbol,fpkm,accession_locus in zip(cuffData['gene_short_name'],cuffData['FPKM'],cuffData['locus']):
            # first, get the accession number
            accession_number = accession_locus.split(':')[0]
            try:
                geneID = mapData[(mapData['GeneSymbol']==symbol) & (mapData['Accession']== accession_number)].GeneID.values[0]
                f.write(str(geneID) + '\t' + str(fpkm) + '\n')
            except:
                f.write(str(symbol) + '\t' + str(fpkm) + '\n')
    print 'done'
# pathway = '/data/shangzhong/DE/contamination_base_line/31to40'
# os.chdir(pathway)
# folders = [f for f in os.listdir(pathway) if os.path.isdir(os.path.join(pathway, f))]
# folders.sort()
# for folder in folders:
#     geneFpkmFile = folder + '/genes.fpkm_tracking'
#     geneSymbol2ID4Cufflink('/data/shangzhong/Database/chok1_ID_symbol.txt',
#                        geneFpkmFile)
#===============================================================================
# This part is more useful for blast results
#===============================================================================
def extractIDFromBlastName(name):
    """
    This function deals with names in the form: gi|625204682|ref|XP_007643745.1|,
    will return only the id, 625204682
    
    * name: name described above.
    """
    index = name.replace('|','x',1).index('|')
    ids = name[3:index]
    return ids


def extract_from_gene2ref(filename,taxID,organism,columnNum=[]):
    """
    This file extracts rows by taxonomy id, and columns from the ncbi gene2refseq file. 
    
    * filename: 'gene2ref.gz' or 'gene2ref'
    
    * taxID: integer indicates taxonomy id of interested organism.
    
    * organism: name of organism to be extracted
    
    * columnNums: a list of columns to extract
    """
    taxID = str(taxID)
    columnNum = map(str,columnNum)
    # gzipped file
    if filename.endswith('.gz'):
        outputfile = filename[:-2] + organism  + '.txt'
        if columnNum == []:
            cmd = ('gunzip -c {input} | awk -F {sep} {parse} {printout} | '
                   'uniq > {output}').format(sep='\"\\t\"',parse='\'$1 == ' + taxID,
                    printout='{print $0}\'',input=filename,
                    output=outputfile)
        else:
            printout = 'print '
            for index in columnNum:
                printout = printout + '$' + index + ' ' + '\"\\t\"' + ' '
            printout = printout[:-6]
            cmd = ('gunzip -c {input} | awk -F {sep} {parse} {printcmd} '
                   '| uniq > {output}').format(sep='\"\\t\"',parse = '\'$1 == ' + taxID,
                    printcmd='{' + printout + '}\'',
                    input=filename,output=outputfile)
    # unzipped file
    else:
        outputfile = filename + organism  + '.txt'
        if columnNum == []:
            cmd = ('awk -F {sep} {parse} {printout} {input} | '
                   'uniq > {output}').format(sep='\"\\t\"',parse='\'$1 == ' + taxID,
                    printout='{print $0}\'',input=filename,
                    output=outputfile)
        else:
            printout = 'print '
            for index in columnNum:
                printout = printout + '$' + index + ' ' + '\"\\t\"' + ' '
            printout = printout[:-6]
            cmd = ('awk -F {sep} {parse} {printcmd} {inputfile} | uniq > {output}').format(sep='\"\\t\"',
                        printcmd='{' + printout + '}\'',inputfile=filename,output=outputfile,
                        parse = '\'$1 == ' + taxID)
    subprocess.call(cmd,shell=True)
    
    return outputfile


def mRNA_prID2geneIDRemote(prID,IDtype='protein'):
    """
    This function searches latest gene ID using given mRNA or protein ID
    through biopython. It connects to romote ncbi latest database. 
    Remote search is slow,so you should only use this if you cannot map gene IDs 
    using local gene2refseq file.
    
    return correspond geneID
    
    * prID: protein or mRNA entrez id.
    
    * IDtype: defalut is protein.
              alternative is nucleotide.
    """
    handle = Entrez.efetch(db=IDtype,id=str(prID),rettype='gb')
    record = handle.read()
    try:
        geneIndex = record.rindex('GeneID:')  # get the last index
        endIndex = geneIndex + 7
        while record[endIndex] != '\"':
            endIndex = endIndex + 1
        geneId = record[geneIndex + 7:endIndex]
    except:
        geneId = 'NA'
        print prID,'is not traceble online'
    return geneId


def mRNA_prIDMap2geneIDMap(prMapFile,gene2refseq1,gene2refseq2,switch='False',IDtype = 'protein'):
    """
    This function change protein id mapping between two species into 
    gene id mappings, after changing to gene ID, extract the unique
    mapping.
    
    * prMapFile:  file name. The first two columns of this file must be protein ID, each is 
                protein id for respective organism. The other columns can be any string
                
    * gene2refseq1: filename of organism corresponds to the 1st column in 
                prMapFile. This file should have three columns, 1st is gene
                id, 2nd is protein accession id, 3rd is protein id.
                
    * gene2refseq2: filename of organism corresponds to the 2nd column in
                prMapFile. This file should have three columns, 1st is gene
                id, 2nd is protein accession id, 3rd is protien id.
    
    * switch: if False, get the uniqueline id mapping sorted by 1st column.
              if True,  get the uniqueline id mapping sorted by 2nd column.
              
    * IDtype: default is IDtype, alternative is nucleotide.
    
    return filename.topn.gene.uniqline.txt
    """
    # build library for each gene2refseq
    protein_gene_lib = {}
    accession_gene_lib = {}
    genePr1 = open(gene2refseq1,'r')
    genePr2 = open(gene2refseq2,'r')
    
    # * library for the 1st organism
    for line1 in genePr1:
        line1 = line1[:-1]
        item1 = line1.split('\t')
        if item1[1] == '-' or item1[2] == '-':
            continue
        protein_gene_lib[item1[2]] = item1[0]
        accession_gene_lib[item1[1][:-2]] = item1[0]
    # * library for the 2nd organism
    for line2 in genePr2:
        line2 = line2[:-1]
        item2 = line2.split('\t')
        if item2[1] == '-' or item2[2] == '-':
            continue
        protein_gene_lib[item2[2]] = item2[0]
        accession_gene_lib[item2[1][:-2]] = item2[0]
    
    # define output file
    outputfile = prMapFile[:-3] + 'gene.txt' # stores the files with gene id mapping
    # change protein id to gene id
    prMap = open(prMapFile,'r')
    geneMap = open(outputfile,'w')
    for line in prMap:
        line = line[:-1]
        item = line.split('\t')
        # test whether there is outdated ids
        try:
            item[0] = protein_gene_lib[item[0]]
        except:
            # get updated ids.
            item[0] = mRNA_prID2geneIDRemote(item[0],IDtype)

        try:
            item[1] = protein_gene_lib[item[1]]
        except:
            # get updated ids.
            item[1] = mRNA_prID2geneIDRemote(item[1],IDtype)

        geneMap.write('\t'.join(item) + '\n')
    geneMap.close()
    prMap.close()

    return outputfile

def uniq1stGene(prMappingFile):
    """
    This file list unique gene ids in the first column, the second column would list 
    all genes mapping to the gene id in the 1st column and can be repeated.
    
    * prMappingFile: a file with two columns. each column is gene id.
    """
    # # # get unique gene ids in the first column
    uniq1stGene = prMappingFile[:-3] + 'uniq1stgene.txt'
    inputFile = open(prMappingFile,'r')
    # create library for uniqline.txt
    library = {}
    for line in inputFile:
        item = line[:-1].split('\t')
        if item[0] in library:
            library[item[0]].append(item[1])
        else:
            library[item[0]] = [item[1]]
    inputFile.close()
    # write to file
    output = open(uniq1stGene,'w')
    for line in library:
        output.write(line + '\t' + ','.join(library[line]) + '\n')
    output.close()
    
    return uniq1stGene

def indexUniqline(uniqFile,switch='False'):
    """
    This function index the id mapping file which has uniqlines. For one id in org1 that has many ids in org2 mapping to, the
    ids in org2 are indexed by blast order.
    
    * uniqFile: string. filename that have uniquelines of id mapping.
    
    * switch: if False, it means index the 2nd column,
              if True, it means index the 1st column.
    """
    res = open(uniqFile,'r')
    outputfile = uniqFile[:-3] + 'index.txt'
    output = open(outputfile,'w')
    i = 0
    before = ''
    for line in res:
        item = line[:-1].split('\t')
        if switch == 'False':   # sort by 1st column
            if item[0] == before:
                i = i + 1
            else:
                before = item[0]
                i = 1
        else:  # switch is true
            if item[1] == before:
                i = i + 1
            else:
                before = item[1]
                i = 1
        output.write('\t'.join(item) + '\t' + str(i) + '\n')
    res.close()
    output.close()
    
    return outputfile


def uniqFirst2Col(index_sort,indexColumn):
    """
    This function extracts the unique first two columns. The file has at least three columsns, the first two are 
    gene ids, the last one is index. The file should be sorted by 1st then 2nd columns.
    
    * index_sort: filename 
    
    * indexColumn: the index of the indexColumn
    """
    res = open(index_sort,'r') # cho2human.top250.gene.uniq.index.sort.txt
    outputfile = index_sort[:-3] + 'uniqline.txt'
    output = open(outputfile,'w')
    before = res.readline()[:-1].split('\t')
    next(res)
    for line in res:
        item = line[:-1].split('\t')
        if item[0] == before[0] and item[1] == before[1]:
            before[indexColumn] = str(min(int(item[indexColumn]),int(before[indexColumn])))
        else:
            output.write('\t'.join(before) + '\n')
            before = item
    output.write('\t'.join(before) + '\n')
    res.close()
    output.close()

    return outputfile

def nonMappedID(mergeFile,annoGeneID,gene_info):
    """
    This function extracts non mapped IDs in annotation file by extracting 
    the non overlapping genes, which means gene Ids in annoGeneID, but not in
    mergeFile (mapping) file
    
    * mergeFile: file with geneIDs in the 1st column
    
    * annoGeneID: file with geneIDs in the 1st column
    
    * gene_info: file with geneIDs in 1st column, gene name in second column.
    """
    # 1. bulid library
    Info = open(gene_info,'r')
    IDname = {}
    for line in Info:
        item = line[:-1].split('\t')
        IDname[item[0]] = item[1]
    Info.close()
    
    # 2. detect nonOverlap
    file1 = open(mergeFile,'r')
    gene1 = []
    for line in file1:
        item = line[:-1].split('\t')
        gene1.append(item[0])
    file1.close()
    
    file2 = open(annoGeneID,'r')
    nonMapFile = annoGeneID[:-3] + 'nonMap.txt'
    output = open(nonMapFile,'w')
    for line in file2:
        item = line[:-1].split('\t')
        if item[0] not in gene1:
            output.write(item[0]+ '\t'+ IDname[item[0]] + '\n')
    # sort by name
    sortNonMapFile = nonMapFile[:-3] + 'sortbyname.txt'
    file2.close()
    output.close()
    cmd = ('sort -k2 {input} > {output}').format(input=nonMapFile,output=sortNonMapFile)
    subprocess.call(cmd,shell=True)
    return sortNonMapFile
# file1 = uniq1stGene('/data/shangzhong/CHO2Human/2wayBlastPresult/human2cho.top1.gene.uniqline.txt')
# file2 = uniq1stGene('/data/shangzhong/CHO2Human/2wayBlastPresult/cho2human.top1.gene.uniqline.txt')

# mergeFile = '/data/shangzhong/CHO2Human/MergeMapping.txt'
# annoGeneID = '/data/shangzhong/CHO2Human/CHO_geneID.txt'
# gene_info = '/data/shangzhong/CHO2Human/namemapping/141028gene_info.cho.txt'
# res = nonMappedID(mergeFile,annoGeneID,gene_info)

def ConvertIDBetweenOrganisms(org1ID,mapping,output):
    """
    This functions convert IDs in org1 into IDs in org2 based on Homology. mapping file is ID mapping files,
    1st column is IDs in org1, 2nd column is IDs in org2. Mapping means they have the same function.
    
    * org1ID: file containing the id file of org1.
    
    * mapping: file containing ids mapping from org1 to org2.
    
    * output: file containing ids of org2
    """
    # build dictionary
    dic = {}
    res = open(mapping,'r')
    for line in res:
        item = line[:-1].split('\t')
        if ';' in item[1]:
            dic[item[0]] = item[1].split(';')
        else: 
            dic[item[0]] = item[1].split(',')
    res.close()
    # convert
    outputfile = open(output,'w')
    inputfile = open(org1ID,'r')
    humanID = []
    for line in inputfile:
        gene1 = line[:-1]
        try:
            gene2 = dic[gene1]
            humanID.extend(gene2)
        except:
            print gene1, 'has no mapping'
            continue
    humanID = list(set(humanID))
    outputfile.write('\n'.join(humanID) + '\n')
    outputfile.close()
 

def geneID2proteinID(geneIDFile,gene2refseq):
    """
    This function gets all protein ids corresponds to the gene id.
    
    * geneIDFile: string. filename which stores gene ids each line.
    
    * gene2refseq: string. filename which stores the gene ids mapping to protein ids.
    """
    # 1. build dictionary
    dic = {}
    res = open(gene2refseq,'r')
    for line in res:
        item = line[:-1].split('\t')
        if item[0] in dic:
            if item[2] != '-':
                dic[item[0]].append(item[2])
        else:
            dic[item[0]] = [item[2]]
    res.close()
    # 2. output transferred ids
    outfile = geneIDFile[:-3] + 'prID.txt'
    with open(outfile,'w') as output:
        with open(geneIDFile,'r') as inputfile:
            for line in inputfile:
                item = line[:-1]
                output.write('\n'.join(dic[item]) + '\n')
    
    return outfile