import sys
import pandas as pd
import subprocess,os
sys.path.append('/home/shangzhong/Codes/Pipeline')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
from Modules.f05_IDConvert import extract_from_gene2ref,mRNA_prIDMap2geneIDMap,uniq1stGene,extractIDFromBlastName,mRNA_prID2geneIDRemote,indexUniqline,uniqFirst2Col
from Modules.p04_ParseBlast import extract_blast_ID_map


def blastp2geneMap(blastFiles,organism1,taxID1,organism2,taxID2,gene2refseq,topNum,IDtype='protein'):
    """
    This function transfers blastp results into geneID mapping results
    
    * blastFiles: a list of 2  way blastp tabular result files. eg: [blast1.txt,blast2.txt]
    
    * organism1: string. the 1st organism. eg: 'cho'
    
    * organism2: string. the 2nd organism. eg: 'human'
    
    * gene2refseq: filename. eg: 'gene2refseq'
    
    returns a list of four files. For first 2: each with two columns of gene ID mapping.
    For thre rest 2: gene2refseq files for each organism
    """
    # extract gene, accession, protein mapping
    if IDtype == 'protein':
        columnNum = [2,6,7,16]
    else:
        columnNum = [2,4,5,16]
    org1ref = extract_from_gene2ref(gene2refseq,taxID1,organism1,columnNum)
    org2ref = extract_from_gene2ref(gene2refseq,taxID2,organism2,columnNum)
    result = []
    switches = ['False','True']
    for blast,switch in zip(blastFiles,switches):
        # extract protein id map
        pr_id_map = extract_blast_ID_map(blast,topNum,switch)  # pr_id_map: cho2human.top1.txt
        # protein id mapping to gene id mapping
        gene_id_map = mRNA_prIDMap2geneIDMap(pr_id_map,org1ref,org2ref,IDtype) # cho2human.top1.gene.txt
        # get unique line of mapping
        uniqFile = gene_id_map[:-3] + 'uniqline.txt'   # cho2human.top1.gene.uniqline.txt
        interFile = gene_id_map[:-3] + 'inter.txt'
        # # # get unique id mappings (each line is unique, but genes are not unique)
        if switch == 'True':
            cmd1 = ('awk -F $\'\\t\' \'BEGIN {OFS} {printrow}\' {input} > {output}').format(
                                      OFS='{FS=\"\\t\"; OFS=FS}',printrow='{print $1,$2}',input=gene_id_map,output=interFile)
            cmd2 = ('sort -k2,2 -k1,1 {input} | uniq > {output}').format(input=interFile,output=uniqFile)
        else:
            cmd1 = ('awk -F $\'\\t\' \'BEGIN {OFS} {printrow}\' {input} > {output}').format(
                                      OFS='{FS=\"\\t\"; OFS=FS}',printrow='{print $1,$2}',input=gene_id_map,output=interFile)
            cmd2 = ('sort -k1,1 -k2,2 {input} | uniq > {output}').format(input=interFile,output=uniqFile)
            
        subprocess.call(cmd1,shell=True)
        subprocess.call(cmd2,shell=True)
        
        subprocess.call(('rm {inter}').format(inter=interFile),shell=True)
        # unique gene ID in 1st column
        uniq = uniq1stGene(uniqFile) # cho2human.top1.gene.uniqline.uniq1stgene.txt
        result.append(uniq)
    result.extend([org1ref,org2ref])
    return result
    

def intersectMapping(file1,file2,sep =','):
    """
    This function merges the two gene ID mapping files.
    Only output the overlapped mapping.
    
    * file1: filename. the file has 2 columns: 1st is unique gene ids, 2nd is all genes mapping to 1st one. 
    * file2: same with file1.
    """
    blst1 = {} # stores the final mapping results
    blst2 = {}
    blst = {}
    gene_list = []
    res1 = open(file1,'r'); res2 = open(file2,'r');
    # build blst dictionary
    for line in res1:
        item = line[:-1].split('\t')
        gene_list.append(item[0])
        blst1[item[0]] = item[1].split(',')
    res1.close()
    
    for line in res2:
        item = line[:-1].split('\t')
        gene_list.append(item[0])
        blst2[item[0]] = item[1].split(',')
    res2.close()
    gene_list = list(set(gene_list))
    # get intersect
    for gene in gene_list:
        if gene in blst1 and gene in blst2:
            intersect = set(blst1[gene]).intersection(blst2[gene])
            blst[gene] = intersect
        else:
            blst[gene] = []
        
    
#         if item[0] in blst:
#             intersect = set(blst[item[0]]).intersection(item[1].split(','))
# #             # if there is no intersection, keep the value it already has
# #             if intersect == set():
# #                 continue
# #             else:
# #                 blst[item[0]] = intersect
#             blst[item[0]] = intersect
#         else:
#             blst[item[0]] = [item[1]]
#     res2.close()
    #out put blst
    filename = file1[:-3] + 'final.txt'
    result = open(filename,'w')
    for item in blst:
        result.write(item + '\t' + sep.join(blst[item]) + '\n')
    result.close()
    return filename

# file1 = '/data/shangzhong/CHO2Human/2wayBlastPresult/cho2human.top1.gene.uniqline.uniq1stgene.txt'
# file2 = '/data/shangzhong/CHO2Human/2wayBlastPresult/human2cho.top1.gene.uniqline.uniq1stgene.txt'
# filename = intersectMapping(file1,file2)

#===============================================================================
# The next part is to deal with the gene ids that are non overlapped in the 2 way blastp
#===============================================================================
def DB4unOverlap(unmapGeneIDs,org1ref,org2ref,blastFiles,topPrNum,topGeneNum):
    """
    This function tries to find why the 2wayblastP result don't have overlapped mapping ids.
    It builds a database file that has all gene ids mappings from both sides and then we can check whethe there 
    are some overlapps
    
    * unmapGeneIDs: filename. gene ids that don't have overlapped mapping results.
    
    * or
    """
    # -------------- 1. get the protein ids of unoverlapped gene ids ---------------
#     org1ref = '/data/shangzhong/CHO2Mouse/2wayBlastPresult/141026gene2refseq.cho.txt'
#     org2ref = '/data/shangzhong/CHO2Mouse/2wayBlastPresult/141026gene2refseq.mouse.txt'
#     # -------------- 2. get top 5 protein mappings ---------------------------------
#     blastFiles = ['/data/shangzhong/CHO2Mouse/2wayBlastPresult/all/cho2mouse.txt','/data/shangzhong/CHO2Mouse/2wayBlastPresult/all/mouse2cho.txt']
    indexFile = []
    switches = ['False','True']
    for blast,switch in zip(blastFiles,switches):
        pr_id_map = extract_blast_ID_map(blast,topPrNum,switch)  # pr_id_map: cho2human.top5.txt
        # protein id mapping to gene id mapping
        gene_id_map = mRNA_prIDMap2geneIDMap(pr_id_map,org1ref,org2ref,switch,IDtype='protein') # cho2human.top5.gene.txt
        # delete the consecutive repeated lines, for lines repeated but are seperated by other lines, they will retain
        uniq_id_map = gene_id_map[:-3] + 'uniq.txt'  # cho2human.top5.gene.uniq.txt
        cmd = ('rev {input} | uniq -f 2 | rev > {output}').format(input=gene_id_map,output=uniq_id_map)
        subprocess.call(cmd,shell=True)
        # index the id mapping, this is for one id in org1 has multiple ids in org2 mapped to 
        index_map = indexUniqline(uniq_id_map,switch)  # cho2human.top5.gene.uniq.index.txt
        # sort index file
        sort_map = index_map[:-3] + 'sort.txt'   # cho2human.top5.gene.uniq.index.sort.txt
        if switch == 'False':  # sort based on 1st column, then on 2nd column
            cmd = ('sort -k1,1n -k2,2n {input} > {output}').format(input=index_map,output=sort_map)
        else:
            cmd = ('sort -k2,2n -k1,1n {input} > {output}').format(input=index_map,output=sort_map)
        subprocess.call(cmd,shell=True)
        #sort_map = '/data/shangzhong/CHO2Human/2wayBlastPresult/top5/cho2human.top250.gene.uniq.index.sort.txt'
        # get uniq first two ids.
        uniqline_map = uniqFirst2Col(sort_map)  # cho2human.top5.gene.uniq.index.sort.uniqline.txt
        # sort by 1rs and 3rd columns for cho2human, and 2nd and 3rd columns for human2cho
        sortbyindex_map = uniqline_map[:-3] + 'index.txt' # cho2human.top5.gene.uniq.index.sort.uniqline.index.txt
        if switch == 'False':
            cmd = ('sort -k1,1n -k5,5n {input} > {output}').format(input=uniqline_map,output=sortbyindex_map)
        else:
            cmd = ('sort -k2,2n -k5,5n {input} > {output}').format(input=uniqline_map,output=sortbyindex_map)
        subprocess.call(cmd,shell=True)
        final_index = indexUniqline(sortbyindex_map,switch) # cho2human.top5.gene.uniq.index.sort.uniqline.index.index.txt
        indexFile.append(final_index)
    return indexFile
    # -------------- 3. get unoverlapped proteins into a list ---------------
    # get unonverlapped ids
    unmapGeneIDs = '/data/shangzhong/CHO2Mouse/finalresult/CHO2Mouse_nonOverlap.txt'
    geneIds = {}
    with open(unmapGeneIDs,'r') as inputfile:
        for line in inputfile:
            geneIds[line[:-1]] = [[] for i in range(4)]    
    # -------------- 4. merge two index files into one ---------------------------------
    # indexFile = ['/data/shangzhong/CHO2Mouse/2wayBlastPresult/all/cho2mouse.top250.gene.uniq.index.sort.uniqline.index.index.txt',
    #               '/data/shangzhong/CHO2Mouse/2wayBlastPresult/all/mouse2cho.top250.gene.uniq.index.sort.uniqline.index.index.txt']
    topGeneNum = 5
    res = open(indexFile[0],'r')
    for line in res:
        item = line[:-1].split('\t')
        if item[0] in geneIds:
            if int(item[4]) > topGeneNum:
                continue
            else:
                geneIds[item[0]][1].append(item[1])
                geneIds[item[0]][0].append(item[4])
        else:
            continue
    
    res = open(indexFile[1],'r')
    for line in res:
        item = line[:-1].split('\t')
        if item[0] in geneIds:
            if int(item[4]) > topGeneNum:
                continue
            else:
                geneIds[item[0]][2].append(item[1])
                geneIds[item[0]][3].append(item[4])
        else:
            continue 
            
    outputfile = indexFile[0][:-22] + 'nonmap.txt'  # cho2human.top5.nonmap.txt
    output = open(outputfile,'w')
    for key in geneIds:
        outline = ('{key}\t{cho2human}\t{human2cho}\n-\t{cho2humanIndex}\t{human2choIndex}\n').format(key=key,
                            cho2human=','.join(geneIds[key][1]),human2cho=','.join(geneIds[key][2]),
                            human2choIndex=','.join(geneIds[key][3]),cho2humanIndex=','.join(geneIds[key][0]))
        output.write(outline)
    output.close()
    print 'done'
    
def mapnonOverlap(unmapGeneIDs,indexFile):
# -------------- 3. get unoverlapped proteins into a list ---------------
# get unonverlapped ids

    geneIds = []
    with open(unmapGeneIDs,'r') as inputfile:
        for line in inputfile:
            geneIds.append(line[:-1])
    
    
    df1 = pd.read_table(indexFile[0],names=['org1','org2','per','len','index1','index2'])
    df2 = pd.read_table(indexFile[1],names=['org1','org2','per','len','index1','index2'])
    mappedGene = {}  # stores the founded gene mappings
    for gene in geneIds:  # for each id 
        map1 = df1.loc[(df1['org1'] == int(gene)) & (df1['index1'] <= 4)]
        map2 = df2.loc[(df2['org1'] == int(gene)) & (df2['index1'] <= 4)]
        if map1.values.tolist() == [] or map2.values.tolist() == []:
            continue
        # get gene ids of organism 2 that mapping to organism 1 
        map1_org2IDs = map1.org2.values.tolist()
        for ID in map1_org2IDs:
            try:  # get the index in map2. (human2cho)
                rowNum = map2.org2[map2.org2 == int(ID)].index.tolist()
                if map2.loc[rowNum,'index1'].values <= 3:
                    mappedGene[gene] = [str(ID)]
                    break
            except:
                continue
    
    still_unmap = unmapGeneIDs[:-3] + 'stillunmap.txt'
    unmap_output = open(still_unmap,'w')
    mapped = unmapGeneIDs[:-3] + 'map.txt'
    mapped_output = open(mapped,'w')
    for gene in geneIds:
        if gene in mappedGene:
            mapped_output.write(gene + '\t' + ','.join(mappedGene[gene]) + '\n')
        else:
            unmap_output.write(gene + '\n')
    unmap_output.close()
    mapped_output.close()
    return mapped

def name2geneMap(organism1,taxID1,organism2,taxID2,gene_info):
    """
    This function maps gene ids between two species by full name.
    
    * organism1: string. the 1st organism. eg: 'cho'
    
    * taxID1: taxonomy ID. eg: 10029
    
    * organism2: string. the 2nd organism. eg: 'human'
    
    * taxID2: taxonomy ID. eg: 9606
    
    return: [cho2human.txt, cho_gene_info.txt]. For 'cho2human.txt' file, ID mapping based on gene names. 
    1st column is unique geIDs, 2nd column is all of the genes mapping to 
    """
    org1_lib = {}; org2_lib = {}
    gene1 = []; gene2=[]
    org1_file = extract_from_gene2ref(gene_info,taxID1,organism1,columnNum=[2,9])
    org2_file = extract_from_gene2ref(gene_info,taxID2,organism2,columnNum=[2,9])
    outputFile = org1_file[:-3] + organism1 + '2' + organism2 + '.txt'
    # build libray for organism1
    result = open(org1_file,'r')
    for line in result:
        item = line[:-1].split('\t')
        # if gene name has many ids
        if item[1] in org1_lib:
            org1_lib[item[1]].append(item[0])
        else:
            org1_lib[item[1]] = [item[0]]
        gene1.append(item[1])
    # build library for organism2
    result = open(org2_file,'r')
    for line in result:
        item = line[:-1].split('\t')
        if item[1] in org2_lib:
            org2_lib[item[1]].append(item[0])
        else:
            org2_lib[item[1]] = [item[0]]
        gene2.append(item[1])
    # find intersect gene name
    intersectName = set(gene1).intersection(gene2)
    # from name to gene id, and output
    output = open(outputFile,'w')
    for name in intersectName:
        geneID1 = org1_lib[name]
        geneID2 = org2_lib[name]
        # geneID1 and geneID2 are lists and can have many items in side it.
        for fir in geneID1:
            output.write(fir + '\t' + ';'.join(geneID2) + '\n')
    output.close()
    return [outputFile,org1_file]
    

def inparanoidParse(inparonoid):
    """
    This functions parse results from inparonoid2ID
    * inparonoid: inparonoid file name
    """
    res = open(inparonoid,'r')
    next(res)
    outputFile = inparonoid[:-3] + 'gene.txt'
    output = open(outputFile,'w')
    for line in res:
        line = line[:-1]
        item = line.split('\t')
        output.write(item[1] + '\t' + item[2] + '\n')
    res.close()
    output.close()
    
    return outputFile


def inparanoid2ID(table,org1ref,org2ref):
    """
    This file process the table file got directly from inparonoid
    
    * table: file name from inparonoid
    return a file with 6 columns: organism1_gene_symbol, organism1_geneID, organism2_gene_symbol,organism2_geneID,
                                  confidence_score, cluster number
    
    * org1ref: refseq for organism1, has 3 columns. 1st is gene id, 2nd is protein accession, 3rd is protein id.
    * org2ref: same with org1ref.
    """
    # build dictionary proteinID:geneID
    protein2gene = {}
    protein2symbol = {}
    res1 = open(org1ref,'r')
    for line in res1:
        item = line[:-1].split('\t')
        if item[2] != '\-':
            protein2gene[item[2]] = item[0]
            protein2symbol[item[2]] = item[3]
    res1.close()
    
    res2 = open(org2ref,'r')
    for line in res2:
        item = line[:-1].split('\t')
        if item[2] != '\-':
            protein2gene[item[2]] = item[0]
            protein2symbol[item[2]] = item[3]
    res2.close()
    
    res = open(table,'r')
    outputfile = table[:-3] + 'txt'
    output = open(outputfile,'w')
    next(res)
    for line in res:
        item = line[:-1].split('\t')
        group1 = item[2].split(); group2 = item[3].split()
        name1 = group1[0]; name2 = group2[0];
        prID1 = extractIDFromBlastName(name1)
        prID2 = extractIDFromBlastName(name2)
        try:
            geneID1 = protein2gene[prID1]
            symbol1 = protein2symbol[prID1]
        except:
            geneID1 = mRNA_prID2geneIDRemote(prID1,IDtype='protein')
            symbol1 = 'none'
        try:
            geneID2 = protein2gene[prID2]
            symbol2 = protein2symbol[prID2]
        except:
            geneID2 = mRNA_prID2geneIDRemote(prID2,IDtype='protein')
            symbol2 = 'none'
        
        result = [symbol1,geneID1,geneID2,symbol2]
        output.write('\t'.join(result) + '\n')
    output.close()
    return outputfile
# org1ref = '/data/shangzhong/CHO2Mouse/2wayBlastPresult/141026gene2refseq.cho.txt'
# org2ref = '/data/shangzhong/CHO2Mouse/2wayBlastPresult/141026gene2refseq.mouse.txt'
# table = '/data/shangzhong/CHO2Mouse/inparanoid_4.1/table.refseq68_cho.faa-refseq68_mouse.faa'
# result = inparanoid2ID(table,org1ref,org2ref)


def MergeMapResults(outputFile,*files):
    """
    This function mereges the mapping results.
    
    * outputFile: file stores output
    
    * files: you can list as many files as you want. each file should have two columns. 
             1st column is organism 1, 2nd column is organism 2.
    """
    dic = {}
    for filename in files: # loop for each file
        res = open(filename,'r')
        for line in res:
            item = line[:-1].split('\t')
            # split 2nd column
            if ',' in item[1]:
                split = item[1].split(',')
            else:
                split = item[1].split(';')
            gene = item[0]
            # key overlapped, get intersection
            if gene in dic:
                #continue
                intersect = set(dic[gene]).intersection(split)
                if intersect == set():  # if there is no overlap, then keep the formater ids
                    continue
                else:
                    dic[gene] = list(intersect)
            else:
                try:
                    dic[gene] = [item[1]] # if want to output overlapped genes, replace [item[1]] with split
                    #dic[gene] = split
                except ValueError:
                    print item
    merge = open(outputFile,'w')
    unmapFile = outputFile[:-3] + 'unmap.txt'
    unmap = open(unmapFile,'w')
    for key in dic:
        if dic[key] == ['']:
            unmap.write(key + '\n')
        else:
            merge.write(key + '\t' + ','.join(dic[key]) + '\n')
    merge.close()
    unmap.close()
    
    return [outputFile,unmapFile]


def MergeBasedOnSignificance(outputFile,*files):
    """
    This function merged based on significance. The only add ids from latter files if they don't exist 
    in the former files.
    """
    dic = {}
    for filename in files: # loop for each file
        res = open(filename,'r')
        for line in res:
            item = line[:-1].split('\t')
            # split 2nd column
            if ',' in item[1]:
                split = item[1].split(',')
            else:
                split = item[1].split(';')
            gene = item[0]
            # key overlapped, get intersection
            if gene in dic:
                continue
            else:
                try:
                    dic[gene] = [item[1]] # if want to output overlapped genes, replace [item[1]] with split
                    #dic[gene] = split
                except ValueError:
                    print item
    merge = open(outputFile,'w')
    unmapFile = outputFile[:-3] + 'unmap.txt'
    unmap = open(unmapFile,'w')
    for key in dic:
        if dic[key] == ['']:
            unmap.write(key + '\n')
        else:
            merge.write(key + '\t' + ','.join(dic[key]) + '\n')
    merge.close()
    unmap.close()
    
    return [outputFile,unmapFile]

def extendByName(mapFile,gene_info_cho):
    """
    This function extends mapping by name. Assign the same mappings for ids with the same names.
    
    * mapFile: geneID mapping file
    
    * gene_info_cho: id, name mapping file. 1st column is id, 2nd column is name.
    """
    # 1. build 2 dictionaries for gene_info, name:IDs, ID:name
    nameID = {}
    IDname = {}
    
    gene_info = open(gene_info_cho,'r')
    for line in gene_info:
        item = line[:-1].split('\t')
        IDname[item[0]] = item[1]
        
        if item[1] in nameID:
            # id repeated
            if item[0] in nameID[item[1]]:
                continue
            else:
                nameID[item[1]].append(item[0])
        else:
            nameID[item[1]] = [item[0]]
    gene_info.close()
    # 2. read the file
    dic = {}
    res = open(mapFile,'r')
    for line in res:
        item = line[:-1].split('\t')
        if ',' in item[1]:
            split = item[1].split(',')
        else:
            split = item[1].split(';')
        dic[item[0]] = split
    res.close()
    # 3. extend by name
    outputDic = {}  # stores the final output mappings.
    for key in dic:
        if key not in outputDic: # if the id is not in the outputdic then add to outputDic, otherwise skip the gene
            # id --> name --> ids
            try:
                name = IDname[key]  # is a list
                geneIDs = nameID[name]  # is a list
            except:
                print key,'is not in gene_info file'
                continue
            # get the uniform map for the ids with same names.
            list_item = []  # each element in it maps to the same gene name in another organism.
            for gene in geneIDs:
                if gene in dic:
                    if dic[gene] != ['']:
                        list_item.append(set(dic[gene]))
                    else:
                        continue
                else:
                    continue
            try:
                inter = set.intersection(*list_item)  # get genes that maps to all the ids with same name in another organism.
            except:
                inter = set()
            if inter != set():  # has the uniform mapping
                for gene in geneIDs:
                    outputDic[gene] = list(inter)
            else:
                for gene in geneIDs:
                    try:
                        outputDic[gene] = dic[gene]
                    except:
                        continue
    print 'done'
    final = mapFile[:-3] + 'final.txt'
    output = open(final,'w')
    # write to file
    for key in outputDic:
        output.write(key + '\t' + ';'.join(outputDic[key]) + '\n')
    output.close()
    return final