import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('ggplot')
#===============================================================================
#                     1. plot glycosytransferase expression level in unit of TPM
#===============================================================================
# file1
ex_file = '/data/shangzhong/CHO_gene_expression_fpkm.txt'
ex_df = pd.read_csv(ex_file,sep='\t',header=0,index_col=0,low_memory=False)
 
gly_g_file = '/data/shangzhong/glycosygenes.txt'
gly_df = pd.read_csv(gly_g_file,sep='\t',header=None,names=['geneid','symbol'])
gly_df['geneid'] = gly_df['geneid'].astype(str)
dic = gly_df.set_index('geneid')['symbol'].to_dict()
genes = dic.keys()
genes = [str(g) for g in genes]
print(genes)

total = ex_df.sum()
ex_df = ex_df.div(total,axis=1)*10**6
ex_df = ex_df.loc[genes].T
ex_df.columns = [dic[x] for x in ex_df.columns.tolist()]
print(ex_df)
ex_df.to_csv('/data/shangzhong/glyco_gene_tpm.csv')
ex_df.index = range(ex_df.shape[0])

# ax = ex_df.plot(title='tpm of glycosytransferase genes in all samples we have')
# ax.set_xlabel('index of samples')
# ax.set_ylabel('tpm')

ax1 = ex_df.plot(title='tpm of glycosytransferase genes in all samples we have')
ax1.set_xlabel('index of samples')
ax1.set_ylabel('tpm')
ax1.set_ylim([0,1])
plt.show()