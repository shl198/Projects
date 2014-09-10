source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")
directory = '/data/shangzhong/Diff_express/Esko_project/htseqcount_geneid'
sampleFiles <- grep("htseq_C",list.files(directory),value=TRUE)
sampleCondition <- c(rep('WT',3),rep('test',3))
sampleTable <- data.frame(sampleName=sampleFiles,fileName=sampleFiles,
                          condition=sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=
                directory, design=~condition)
ddsHTSeq$condition <- factor(ddsHTSeq$condition, 
                             levels=c("WT","test"))
dds <- DESeq(ddsHTSeq)
res <- results(dds)
resOrdered <- res[order(res$padj),]
head(resOrdered)
#plot MA
plotMA(res, main="C4_6 VS C2_3", ylim=c(-2,2))
#export differentially expressed genes the threshold is 0.05
result <- resOrdered[complete.cases(resOrdered),]
sig_result <- result[result$padj < 0.05 & (result$log2FoldChange >= 0.58 | result$log2FoldChange <=-0.58),]
write.csv(sig_result,"C4_6VS_C1_3.csv")
write.csv(resOrdered,'test.csv')
#transformation dta
rld <- rlog(dds)
# Heat map
library("RColorBrewer")
library("gplots")
#heat map of distance matrix
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition))
colours = colorRampPalette(rev(brewer.pal(9,"Blues"))) (255)
heatmap.2(mat, trace="none", col = colours,margin=c(10,10))
print(plotPCA(rld, intgroup=c("condition")))
# gene clustering
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)
heatmap.2(assay(rld)[topVarGenes,],scale="row",trace="none",dendrogram="column",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu"))),margin=c(9,9))

#=====================================================================
#============= analyze the mbmec sample ==============================
sampleFiles <- grep("htseq_M",list.files(directory),value=TRUE)
sampleCondition <- c(rep('creN',2),rep('creP',4),rep('Ext',4))
sampleTable <- data.frame(sampleName=sampleFiles,fileName=sampleFiles,
                          condition=sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=
                                         directory, design=~condition)
ddsHTSeq$condition <- factor(ddsHTSeq$condition, 
                             levels=c("Ext","creN","creP"))
background <- assay(ddsHTSeq)
background <- background[rowSums(background[,-1])>0,]
background <- background[apply(background[,-1], 1, function(x) !all(x==0)),]
write.csv(background,'background.csv')

dds <- DESeq(ddsHTSeq)
#================== (1) creN VS creP ========================= 
res <- results(dds, contrast = c("condition","creN","creP"))
resOrdered <- res[order(res$padj),]
result <- resOrdered[complete.cases(resOrdered),]
sig_result <- result[result$padj < 0.01 & (result$log2FoldChange >= 0.7 | result$log2FoldChange <=-0.7),]
write.csv(sig_result,"creN_VS_creP.csv")

#================== (2) crepN VS Ext ========================= 
res <- results(dds, contrast = c("condition","creN","Ext"))
resOrdered <- res[order(res$padj),]
result <- resOrdered[complete.cases(resOrdered),]
sig_result <- result[result$padj < 0.01 & (result$log2FoldChange >= 1 | result$log2FoldChange <=-1),]
write.csv(sig_result,"creN_VS_Ext.csv")
write.csv(result,'test.csv')
#================== (3) crepP VS Ext ========================= 
res <- results(dds, contrast = c("condition","creP","Ext"))
resOrdered <- res[order(res$padj),]
result <- resOrdered[complete.cases(resOrdered),]
sig_result <- result[result$padj < 0.01 & (result$log2FoldChange >= 0.7 | result$log2FoldChange <=-0.7),]
write.csv(sig_result,"creP_VS_Ext.csv")

#transformation
rld <- rlog(dds)
#heat map of distance matrix
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- paste(rld$condition)
colnames(mat) <- rownames(mat)
library("gplots")
library("RColorBrewer")
colours = colorRampPalette(rev(brewer.pal(9,"Blues"))) (255)
heatmap.2(mat, trace="none", col = colours,margin=c(10,10))
print(plotPCA(rld, intgroup=c("condition")))

# gene_clustering
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),50)
heatmap.2(assay(rld)[topVarGenes,],scale="row",trace="none",dendrogram="column",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu"))),margin=c(9,9))

