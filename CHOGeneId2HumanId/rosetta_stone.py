import subprocess,os,sys
sys.path.append('/home/shangzhong/Codes/Pipeline')
from Modules.f05_IDConvert import extract_from_gene2ref,prIDMap2geneIDMap
from Modules.p04_ParseBlast import blastName2ID

def rosetta_pair(group):
    """
    This functions detects rosetta pairs, given a list of lines which have same 2nd column.
    
    * group: a list of blast tabular results. each item in the list is a line of blast tabular result.
            they should have same 2nd reference name.
    """
    rosetta = []
    for i in range(len(group)-1):
        for j in range(i + 1,len(group)):
            if (max(int(group[i][8]),int(group[i][9])) < min(int(group[j][8]),int(group[j][9]))):
                ref = group[i][1]; firstName = group[i][0]; secondName = group[j][0]
                firBegin=group[i][8];firEnd=group[i][9]
                secBegin=group[j][8];secEnd=group[j][9]
                rosetta.append([ref,firstName,secondName,firBegin,firEnd,secBegin,secEnd])
            elif (max(int(group[j][8]),int(group[j][9])) < min(int(group[i][8]),int(group[i][9]))):
                ref = group[i][1]; firstName = group[j][0]; secondName = group[i][0]
                firBegin=group[j][8];firEnd=group[j][9]
                secBegin=group[i][8];secEnd=group[i][9]
                rosetta.append([ref,firstName,secondName,firBegin,firEnd,secBegin,secEnd])
            else:
                continue
    return rosetta

def rosetta_stone(blastFile,organism1,taxID1,organism2,taxID2,gene2refseq):
    """
    This functions detects rosetta stone pairs in the blast result.
    
    * blastFile: fileanme of tabular blast output.
    """
    # # (1) change first 2 columns into ID, extract top hit
    IDblast = blastName2ID(blastFile)
    # # (2) extract cho, human from gene2refseq
    org1ref = extract_from_gene2ref(gene2refseq,taxID1,organism1,columnNum=[2,6,7])
    org2ref = extract_from_gene2ref(gene2refseq,taxID2,organism2,columnNum=[2,6,7])
    # # (3) replace protein ID with gene ID
    IDmap = prIDMap2geneIDMap(IDblast,org1ref,org2ref)
    # # (4) sort file by 2nd column
    sortFile = IDmap[:-3] + 'sort.txt'
    cmd = ('sort -k2,2 -n {input} > {output}').format(input=IDmap,output=sortFile)
    subprocess.call(cmd,shell=True)
    os.remove(IDblast)
    
    # # (5) now come to the parsing rosetta stone pairs
    filename = sortFile
    res = open(filename,'r')
    inter = 'inter.txt'
    interOut = open(inter,'w')
    outputFile = filename[:-16] + 'rosetta.txt'
    output = open(outputFile,'w') 
    # combine lines with same 2 columns.
    id_pair = res.readline()[:-1].split('\t')
    # start from second line
    for line in res:
        line = line[:-1]
        item = line.split('\t')
        # should merge
        if id_pair[0] == item[0] and id_pair[1] == item[1]:
            id_pair[8] = str(min(int(id_pair[8]),int(id_pair[9]),int(item[8]),int(item[9])))
            id_pair[9] = str(max(int(id_pair[8]),int(id_pair[9]),int(item[8]),int(item[9])))
        else:
            interOut.write('\t'.join(id_pair) + '\n')
            id_pair = item
    interOut.write('\t'.join(id_pair) + '\n')
    interOut.close()
    res.close()
    # now interOut has unique first 2 columns of blast tabular results
    res = open(inter,'r')
    group = [res.readline()[:-1].split('\t')] # group stores lines with same reference
    for line in res:
        item = line[:-1].split('\t')
        if item[1] == group[-1][1]:
            group.append(item)
        else:
            # whether to do rosetta pair detection
            if len(group) > 1:
                rosetta = rosetta_pair(group)
                # output to file
                if rosetta != []:
                    for pair in rosetta:
                        output.write('\t'.join(pair) + '\n')
            # there is no rosetta pairs
            group = [item]
    if len(group) > 1:
                rosetta = rosetta_pair(group)
                # output to file
                if rosetta != []:
                    for pair in rosetta:
                        output.write('\t'.join(pair) + '\n')
                        
    output.close()
    res.close()
    return outputFile
blastFile = '/data/shangzhong/CHO2Human/2wayBlastPresult/human2cho.txt'
gene2refseq = '/data/shangzhong/CHO2Human/2wayBlastPresult/141026gene2refseq.gz'
rosetta = rosetta_stone(blastFile,'cho',10029,'human',9606,gene2refseq)

