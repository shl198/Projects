import sys
import subprocess,os
sys.path.append('/home/shangzhong/Codes/Pipeline')
from IdMappingModule import *
#sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # disable buffer
#===============================================================================
#  The following is the final pipeline
#===============================================================================

# # # #========== 1. get 2way blast results ================================
# blastFiles = ['/data/shangzhong/CHO2Mouse/2wayBlastPresult/cho2mouse.txt',
#                '/data/shangzhong/CHO2Mouse/2wayBlastPresult/mouse2cho.txt']
# gene2refseq = '/data/shangzhong/CHO2Mouse/2wayBlastPresult/141026gene2refseq.gz'
# blast_result = blastp2geneMap(blastFiles,'cho','10029','mouse','10090',gene2refseq,topNum=1)
# final_blast = intersectMapping(blast_result[0],blast_result[1],',')
#       
# # # #=========== 2. get name map results ============================
# gene_info = '/data/shangzhong/CHO2Mouse/namemapping/141028gene_info.gz'
# name_result = name2geneMap('cho','10029','mouse','10090',gene_info)
#     
# # # #=========== 3. get inparonoid results ========================== 
# inparonoid = '/data/shangzhong/CHO2Mouse/inparanoid_4.1/cho2mouse.inpara.txt'
# inpara_result = inparanoidParse(inparonoid)
#      
#    
# # #=========== 4. merge mapping files from all sources ==========================
mergemap = '/data/shangzhong/CHO2Mouse/MergeMapping.txt'
final_blast = '/data/shangzhong/CHO2Mouse/2wayBlastPresult/cho2mouse.top1.gene.uniqline.uniq1stgene.final.txt'
name_result = ['/data/shangzhong/CHO2Mouse/namemapping/141028gene_info.cho.cho2mouse.txt','/data/shangzhong/CHO2Mouse/namemapping/141028gene_info.cho.txt']
inpara_result = '/data/shangzhong/CHO2Mouse/inparanoid_4.1/cho2mouse.inpara.gene.txt'
merge = MergeMapResults(mergemap,name_result[0],inpara_result,final_blast)   # merge = [MergeMapping.txt,MergeMapping.unmap.txt]
   
# # #=========== 5. map the unmapped genes ==========================
indexFile = ['/data/shangzhong/CHO2Mouse/2wayBlastPresult/all/cho2mouse.top250.gene.uniq.index.sort.uniqline.index.index.txt',
             '/data/shangzhong/CHO2Mouse/2wayBlastPresult/all/mouse2cho.top250.gene.uniq.index.sort.uniqline.index.index.txt']
mapped = mapnonOverlap(merge[1],indexFile)  # MergeMapping.unmap.map.txt
# # #=========== 6. concatenate MergeMapping and unmap.map ==========================
furtherMerge = mergemap[:-3] + 'furtherMerge.txt'
concat = ('cat {mergemap} {mapped} > {output}').format(mergemap=mergemap,mapped=mapped,output=furtherMerge)  
subprocess.call(concat,shell=True)  # MergeMapping.furtherMerge.txt
# # #=========== 7. extend mapping results by name ==========================
finalbyname = extendByName(furtherMerge,name_result[1])  # MergeMapping.furtherMerge.final.txt
#===============================================================================
#        Process mRNA results
#===============================================================================
# # #=========== 8. get the mRNA mapping result ==========================
blastFiles = ['/data/shangzhong/CHO2Mouse/mRNA/cho2mouse.txt','/data/shangzhong/CHO2Mouse/mRNA/mouse2cho.txt']
organism1 = 'cho'; taxID1 = '10029'; organism2 = 'mouse'; taxID2 = '10090'
gene2refseq = '/data/shangzhong/CHO2Mouse/mRNA/141026gene2refseq.gz'
topNum = 1
mRNAres = blastp2geneMap(blastFiles,organism1,taxID1,organism2,taxID2,gene2refseq,topNum,IDtype='nucleotide')
finalmRNA = intersectMapping(mRNAres[0],mRNAres[1],sep =',')
outputFile = '/data/shangzhong/CHO2Mouse/finalMergeWithmRNA.txt'
finalbyname = '/data/shangzhong/CHO2Mouse/MergeMapping.furtherMerge.final.txt'
finalMerge = MergeBasedOnSignificance(outputFile,finalbyname,finalmRNA) # [finalMergeWithmRNA.txt,finalMergeWithmRNA.unmap.txt]
finalbyname = extendByName(finalMerge[0],name_result[1])  # finalMergeWithmRNA.final.txt