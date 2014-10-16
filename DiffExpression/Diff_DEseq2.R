#!/usr/bin/env Rscript
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
if (!require("DESeq2")) {
  install.packages("DESeq2", dependencies = TRUE)
  library(DESeq2)
}
#====== 1. set working directory and list files ========
#filePath <- '/data/shangzhong/Diff_express/htseq_GlycoCosmic'
args <- commandArgs(TRUE)
filePath <- args[1]
setwd(filePath)
DE_file <- list.files(filePath,pattern='\\.txt$')
compare <- read.csv("DE_pairs.csv",sep=',') # file decide which samples to compare

#====== 2. 
col <- colnames(compare)
row <- nrow(compare)

for (i in 1:row) {
    # 1st part, list the control sample and test sample
    control <- c()
    test <- c()
    for (j in 1:length(col)) {
        if (grepl('control',col[j],ignore.case=TRUE) & !is.na(compare[i,j])) {
            control <- c(control,DE_file[compare[i,j]])
        } else if (grepl('test',col[j],ignore.case=TRUE) & !is.na(compare[i,j])) {
            test <- c(test,DE_file[compare[i,j]])
        }
    }
    
    # 2nd part, do the DE for one compare
    sampleFiles <- c(control,test)
    sampleCondition <- c(rep('control',length(control)),rep('test',length(test)))
    sampleTable <- data.frame(sampleName=sampleFiles,fileName=sampleFiles,
                              condition=sampleCondition)
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=
                                             filePath, design=~condition)
    ddsHTSeq$condition <- factor(ddsHTSeq$condition, 
                                 levels=c("control","test"))
    dds <- DESeq(ddsHTSeq)
    res <- results(dds)
    resOrdered <- res[order(res$padj),]
    result <- resOrdered[complete.cases(resOrdered),]
    sig_result <- result[result$padj < 0.01 & (result$log2FoldChange >= 0.58 | result$log2FoldChange <=-0.58),]
    outputFile <- paste('Line',toString(i),'.csv',sep="")
    write.csv(sig_result,outputFile)
    # output cluster figure
    pdf(paste(strsplit(outputFile,'\\,')[[1]][1],'.pdf',sep=""))
    rld <- rlog(dds)
    # Heat map
    library("RColorBrewer")
    library("gplots")
    #heat map of distance matrix
    distsRL <- dist(t(assay(rld)))
    mat <- as.matrix(distsRL)
    rownames(mat) <- colnames(mat) <- with(colData(dds),paste(condition,sampleFiles , sep=" : "))
    colours = colorRampPalette(rev(brewer.pal(9,"Blues"))) (255)
    heatmap.2(mat, trace="none", col = colours,margin=c(17,17))
    dev.off()
}