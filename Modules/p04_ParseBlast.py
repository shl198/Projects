from p03_ParseSam import get_name

def annotateBlast(blastFile,anno_type):
    """
    This function annotates the results got by running all sorts of 
    blast algorithm. It adds the names of reference genomes
    
    * blastFiles: str. blast file name    
    * anno_type: str. 'nucleotide' or 'protein'
    
    output a list of files: [f1.blast.anno.txt,f2.blast.anno.txt,...]
    """
    result = open(blastFile,'r')
    annotation_file = blastFile[:-3] + 'anno.txt'
    output = open(annotation_file,'w')
    header = ['Reference Organism genome','Sequence type','query id',
              'subject id','identity','length','mismatch','gapopen',
              'query start','query end','subject start','subject end',
              'evalue','bitscore']
    output.write('\t'.join(header)+'\n')
    i = 0
    for item in result:
        i = i +1
        print i
        line = item.split('\t')
        [gene_name,gene_type] = get_name(line[1],anno_type)
        line.insert(0,gene_name)
        line.insert(1,gene_type)
        output.write('\t'.join(line))
    result.close()
    output.close()
    

def extract_blast_ID_map(blastFile,topNum = 1,switch='False'):
    """
    This function extracts the top prtein id mapping results from blast tabular 
    files. 
    
    * blastFiles: blast filename
    
    * topNum: an integer indicates how many top hits you want. default is 1
    
    return file with 2 columns of protein IDs.
    """
    id_map = [[]] * 2
    i = 0 # number of tophits
    outputFile = blastFile[:-3] + 'top' + str(topNum) + '.txt'
    output = open(outputFile,'w')
    result = open(blastFile,'r')
    for line in result:
        item = line.split('\t')
        query_index = item[0].index('|ref')
        ref_index = item[1].index('|ref')
        query_id = item[0][3:query_index]
        ref_id = item[1][3:ref_index]
        # output the 1st line
        if id_map[0] == []:
            id_map[0] = query_id
            id_map[1] = ref_id
            i = i + 1
            # switch the column
            if switch == 'True':
                outline = ('{ref}\t{query}\t{iden}\t{len}\n').format(ref=ref_id,query=query_id,
                                                                     iden=item[2],len=item[3])
                output.write(outline)  # changed here
            else:
                outline = ('{query}\t{ref}\t{iden}\t{len}\n').format(ref=ref_id,query=query_id,
                                                                     iden=item[2],len=item[3])
                output.write(outline)
        # if both are equal
        elif (query_id == id_map[0] and ref_id == id_map[1]):
            continue
        # if second item don't equal, count
        elif (query_id == id_map[0] and ref_id != id_map[1]):
            if i == topNum:
                continue
            else:
                i = i + 1
                id_map[0] = query_id
                id_map[1] = ref_id
                if switch == 'True':
                    outline = ('{ref}\t{query}\t{iden}\t{len}\n').format(ref=ref_id,query=query_id,
                                                                     iden=item[2],len=item[3])
                    output.write(outline)
                else:
                    outline = ('{query}\t{ref}\t{iden}\t{len}\n').format(ref=ref_id,query=query_id,
                                                                     iden=item[2],len=item[3])
                    output.write(outline)
        # if first item don't equal
        else:
            i = 1
            id_map[0] = query_id
            id_map[1] = ref_id
            if switch == 'True':
                outline = ('{ref}\t{query}\t{iden}\t{len}\n').format(ref=ref_id,query=query_id,
                                                                     iden=item[2],len=item[3])
                output.write(outline)
            else:
                outline = ('{query}\t{ref}\t{iden}\t{len}\n').format(ref=ref_id,query=query_id,
                                                                     iden=item[2],len=item[3])
                output.write(outline)
    output.close()
    return outputFile

def blastName2ID(blastFile,topNum=1):
    """
    This function changes query name and refernece name into id
    
    * blastFile: tablular blast file output
    
    * topNum: an int indicates the number of top hits
    """
    id_map = [[]] * 2
    i = 0 # number of tophits
    
    outputFile = blastFile[:-3] + 'ID.txt' # store intermediate results
    res = open(blastFile,'r')
    output = open(outputFile,'w')
    # replace the refname and query name with ID only.
    for line in res:
        line = line[:-1]
        item = line.split('\t')
        index = item[0].replace('|','x',1).index('|');item[0] = item[0][3:index]
        index = item[1].replace('|','x',1).index('|');item[1] = item[1][3:index]
        # output first line
        if id_map[0] == []:
            id_map[0] = item[0]
            id_map[1] = item[1]
            i = i + 1
            output.write('\t'.join(item) + '\n')
        # if both are equal
        elif (item[0] == id_map[0] and item[1] == id_map[1]):
            output.write('\t'.join(item) + '\n')
        # if second item don't equal, count
        elif (item[0] == id_map[0] and item[1] != id_map[1]):
            if i == topNum:
                continue
            else:
                i = i + 1
                id_map[1] = item[1]
                output.write('\t'.join(item) + '\n')
        # if first item don't equal
        else:
            i = 1
            id_map[0] = item[0]
            id_map[1] = item[1]
            output.write('\t'.join(item) + '\n')
    output.close()
    return outputFile