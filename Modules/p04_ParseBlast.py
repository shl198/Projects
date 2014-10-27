from p03_ParseSam import get_name

def annotateBlast(blastFiles,anno_type):
    """
    This function annotates the results got by running all sorts of 
    blast algorithm. It adds the names of reference genomes
    
    blastFiles: a list of blast tabular result. [f1.blast.txt,f2.blast.txt,...]
    
    anno_type: 'nucleotide' or 'protein'
    
    output a list of files: [f1.blast.anno.txt,f2.blast.anno.txt,...]
    """
    for blastFile in blastFiles:
        result = open(blastFile,'r')
        annotation_file = blastFile[:-3] + 'anno.txt'
        output = open(annotation_file,'w')
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
        
def protein_ID_map(blastFiles,topNum = 1):
    """
    This function extracts the top mapping results from blast tabular files
    
    blastFiles: a list of blast files. [f1.blast.txt,f2.blast.txt,...]
    
    topNum: an integer indicates how many top hits you want. default is 1
    """
    id_map = [[]] * 2
    i = 0 # number of tophits
    for blast in blastFiles:
        outputFile = blast[:-3] + 'top1.txt'
        output = open(outputFile,'w')
        result = open(blast,'r')
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
                output.write(query_id + '\t' + ref_id + '\n')
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
                    output.write(query_id + '\t' + ref_id + '\n')
            # if first item don't equal
            else:
                i = 1
                id_map[0] = query_id
                id_map[1] = ref_id
                output.write(query_id + '\t' + ref_id + '\n')
        output.close()
blastFiles = ['/data/shangzhong/CHO2Human/Human/homo2cho.txt']
gene_ID_map(blastFiles)

