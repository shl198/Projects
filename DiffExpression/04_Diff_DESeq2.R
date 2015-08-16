#!/usr/bin/env Rscript
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(gtools)
if (!require("DESeq2")) {
  install.packages("DESeq2", dependencies = TRUE)
  library(DESeq2)
}
#====== 1. set working directory and list files ========
#filePath <- '/data/shangzhong/Contamination_Result'
args <- commandArgs(TRUE)
filePath <- args[1]
setwd(filePath)
DE_file <- list.files(filePath,pattern='\\.txt$')
DE_file <- mixedsort(DE_file)
compare <- read.csv("04_Diff_pairs.csv",sep=',') # file decide which samples to compare

#====== 2. Do the differential expression ========
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
    sig_result <- result[result$padj < 0.05 & (result$log2FoldChange >= log2(1.5) | result$log2FoldChange <=-log2(1.5)),]
    outputFile <- paste('Line',toString(i),'.csv',sep="")
    write.csv(result,outputFile)
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
