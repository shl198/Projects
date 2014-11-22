"""
This file includes functions of all kinds of ID converts, such as gene symbols to gene ids, protein ids to gene ids...
"""
import subprocess
from Bio import Entrez
Entrez.email = "shl198@eng.ucsd.edu"
#===============================================================================
#  This part is for DE analysis
#===============================================================================

def geneSymbol2EntrezID(Dict,output_path,inputpath,sym2ID='yes'):
    import os
    """
    This fuction converts gene symbol to entriz_gene_ID, and ID to symbol
    
    * Dict: Convert file. 1st column is ID, 2nd column is symbol
    
    * sym2ID: if yes: symbol to ID
              if no : ID to symbol
    """
    allFiles = os.listdir(inputpath)
    filelist = [f for f in allFiles if f.endswith('.txt')]
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
        result.close()
        output.close()

#===============================================================================
# This part is more useful for blast results
#===============================================================================

def extract_from_gene2ref(filename,taxID,organism,columnNum=[]):
    """
    This file extracts rows by taxonomy id, and columns. 
    
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
filename = '/data/shangzhong/CHO2Human/namemapping/141028gene_info.gz'
extract_from_gene2ref(filename,10029,'cho',[2,3,9])


def prID2geneIDRemote(prID):
    """
    This function searches latest gene ID using given protein ID
    through biopython. It connets to romote ncbi latest database. 
    Remote search is slow,so only use this if you cannot map gene IDs 
    using local gene2refseq file.
    
    * prID: protein entrez id.
    """
    handle = Entrez.efetch(db='protein',id=str(prID),rettype='gb')
    record = handle.read()
    geneIndex = record.rindex('GeneID:')  # get the last index
    endIndex = geneIndex + 7
    while record[endIndex] != '\"':
        endIndex = endIndex + 1
    geneId = record[geneIndex + 7:endIndex]
    
    return geneId


def prIDMap2geneIDMap(prMapFile,gene2refseq1,gene2refseq2):
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
    uniqFile = outputfile[:-3] + 'uniqline.txt'
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
            item[0] = prID2geneIDRemote(item[0])

        try:
            item[1] = protein_gene_lib[item[1]]
        except:
            # get updated ids.
            item[1] = prID2geneIDRemote(item[1])

        geneMap.write('\t'.join(item) + '\n')
    geneMap.close()
    prMap.close()
    
    # # # get unique id mappings (each line is unique, but genes are not unique)
    cmd = ('sort {input} | uniq > {output}').format(input=outputfile,output=uniqFile)
    subprocess.call(cmd,shell=True)
    
    return uniqFile

def uniq1stGene(prMappingFile):
    """
    This file list unique gene ids in the first column, the second column would list 
    all genes mapping to the gene id in the 1st column.
    
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

def nonOverlapID(mergeFile,annoGeneID,gene_info):
    """
    This function extracts non mapped IDs in annotation file by extracting 
    the non overlapping genes, which means gene Ids in annoGeneID, but not in
    mergeFile file
    
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

# mergeFile = '/data/shangzhong/CHO2Mouse/MergeMapping.txt'
# annoGeneID = '/data/shangzhong/CHO2Mouse/CHO_geneID.txt'
# gene_info = '/data/shangzhong/CHO2Mouse/namemapping/141028gene_info.cho.txt'
# res = nonOverlapID(mergeFile,annoGeneID,gene_info)


