def ID_Convert(Dict,output_path,inputpath):
    import os
    """
    This fuction convert gene symbol to entriz_gene_ID
    Dict should the be convert file
    """
    allFiles = os.listdir(inputpath)
    filelist = [f for f in allFiles if f.endswith('.txt')]
    #=========== build dictionary ============
    dic = {}
    result = open(Dict,'r')
    for line in result:
        item = line.split(' ')
        dic[item[1][:-1]] = item[0]
    for filename in filelist:
        result = open(inputpath + '/' + filename,'r')
        output = open(output_path + '/' + filename,'w')
        for line in result:
            item = line.split('\t')
            try:
                item[0] = dic[item[0]]
                output.write('\t'.join(item))
            except:
                output.write('\t'.join(item))
        output.close()
