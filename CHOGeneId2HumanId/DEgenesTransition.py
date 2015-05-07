import pandas
import os
import ipdb
import numpy
from natsort import natsorted
import subprocess

# mappingFile = '/data/shangzhong/CHO2Human/finalMergeWithmRNA.final.txt'
# ref = pandas.read_csv(mappingFile, delimiter='\t', names=["CHO","HUMAN"], header=None, index_col = 0)
# humanList = []
# for ids in ref["HUMAN"].values:
#     if pandas.isnull(ids):
#         continue
#     ids = str(ids)
#     if ',' in ids:
#         ids = ids.split(",")
#     elif ';' in ids:
#         ids = ids.split(";")
#     if ids.__class__ is list:
#         humanList.extend(ids)
#     else:
#         humanList.extend([ids])
# humanList = list(set(humanList))
# f = open('CHO_human.txt', 'w')
# for item in humanList:
#     f.write("%s\n" % item)
# f.close()


# source = "/data/shangzhong/DE/new"
# deseq_files = [f for f in os.listdir(source) if f.endswith('csv')]
# 
# diff_data = {}            
# for f in deseq_files:
#     diff_data[f] = pandas.read_csv(source + '/' + f, header=0, index_col=0)
# 
# for key in diff_data.keys():
#     mapped = ref.loc[diff_data[key].index]
#     humanList = []
#     for ids in mapped['HUMAN']:
#         if pandas.isnull(ids):
#             print key,'does not have mapping'
#             continue
# 
#         ids = str(ids)
#         if ',' in ids:
#             ids = ids.split(",")
#         elif ';' in ids:
#             ids = ids.split(";")
#         if ids.__class__ is list:
#             humanList.extend(ids)
#         else:
#             humanList.extend([ids])
#     humanList = list(set(humanList))
#     f = open(source + '/' + key + '_human.txt', 'w')
#     for item in humanList:
#         f.write("%s\n" % item)
#     f.close()


source = "/data/shangzhong/DE/GlycoCosmic_gsnap/DE_result"
os.chdir(source)
deseq_files = [f for f in os.listdir(source) if f.endswith('csv')]
mappingFile = '/data/shangzhong/CHO2Mouse/finalMergeWithmRNA.final.txt'
res = open(mappingFile,'r')
dic = {}
for line in res:
    item = line[:-1].split('\t')
    if ';' in item[1]:
        split = item[1].split(';')
    else:
        split = item[1].split(',')
    dic[item[0]] = split
res.close()

for f in deseq_files:
    outFile = f[:-4] + '_mouse.txt'
    output = open(outFile,'w')
    res = open(f,'r')
    firstline = res.readline().split(',')
    output.write('\t'.join(firstline))
    next(res)
    for line in res:
        item = line[:-1].split(',')
        try:
            mapids = dic[str(item[0])]
        except:
            print item[0],'is not in database'
            continue
        for gene in mapids:
            output.write(gene + '\t' + '\t'.join(item[1:]) + '\n')
    output.close()
    # extract the unique ones
    finalresult = outFile[:-3] + 'uniq.txt'
    cmd = ('sort -k1,1n {input} | rev | uniq -f 6 | rev | sort -k 7,7g > {output}').format(
                                            input=outFile,output=finalresult)
    subprocess.call(cmd,shell=True)
    
# sort -k1,1 pgsa_WT_VS_KO_human.txt | rev | uniq -f 6 | rev | sort -k 7,7 > test.txt
    
    
