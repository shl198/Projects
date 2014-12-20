import pandas
import os
import ipdb
import numpy
from natsort import natsorted
filepath = '/data/shangzhong/DE/htseq'
os.chdir(filepath)
mappingFile = '/data/shangzhong/CHO2Mouse/finalMergeWithmRNA.final.txt'
ref = pandas.read_csv(mappingFile, delimiter='\t', names=["CHO","HUMAN"], header=None, index_col = 0)
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


source = "/data/shangzhong/DE/diff_result"
deseq_files = [f for f in os.listdir(source) if f.endswith('csv')]

diff_data = {}            
for f in deseq_files:
    diff_data[f] = pandas.read_csv(source + '/' + f, header=0, index_col=0)

for key in diff_data.keys():
    mapped = ref.loc[diff_data[key].index]
    humanList = []
    for ids in mapped['HUMAN']:
        if pandas.isnull(ids):
            print key,'does not have mapping'
            continue

        ids = str(ids)
        if ',' in ids:
            ids = ids.split(",")
        elif ';' in ids:
            ids = ids.split(";")
        if ids.__class__ is list:
            humanList.extend(ids)
        else:
            humanList.extend([ids])
    humanList = list(set(humanList))
    f = open(source + '/' + key + '_mouse.txt', 'w')
    for item in humanList:
        f.write("%s\n" % item)
    f.close()