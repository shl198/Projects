import subprocess
def remove(files):
    """
    this function can remove files provided
    """
    if isinstance(files,str):
        subprocess.call('rm {file}'.format(file=files),shell=True)
    if isinstance(files,list):
        cmd = ''
        for f in files:
            cmd = cmd + 'rm {file} & '.format(file=f)
        subprocess.call(cmd[:-3],shell=True)
    

def get_parameters(parFile):
    """
    This function list all parameters for all pipelines.
    And return a dictionary
    """
    res = open(parFile)
    dic = {}
    for line in res:
        if line[0].isalpha():
            item = line.split('\t')
            value = item[1] # combine all values into a single string
            if ',' in value:
                rg = value.split(',')
                rg[-1] = rg[-1][:-1]
                dic[item[0]] = rg
            else:
                dic[item[0]] = item[1][:-1]
        else:
            continue
    if isinstance(dic['readGroup'],str):
                dic['readGroup'] = [dic['readGroup'][:-1]]
    return dic

