library(gtools)
if (!require("DESeq2")) {
  install.packages("DESeq2", dependencies = TRUE)
  library(DESeq2)
}

filePath <- '/data/shangzhong//DE//Winzeler//pberg_counts'
filePath <- getwd()
setwd(filePath)
DE_file <- list.files(filePath,pattern='\\.txt$')
DE_file <- mixedsort(DE_file)

sampleFiles <- c(DE_file)
sampleCondition <- c(rep('infect2h',3),rep('infecte48h',3),
                     rep('uninfect2h',3),
                     rep('uninfect48h',3))
                     
sampleNames <- sampleFiles

sampleTable <- data.frame(sampleName=sampleNames,filename=sampleFiles,condition=sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=
                                         filePath, design=~condition)

ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
data <- counts(ddsHTSeq,normalized=F)

#=========== Generate heatmap for the samples ==================================
library("RColorBrewer")
library("gplots")rld <- rlogTransformation(ddsHTSeq,blind=TRUE)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colData(rld)$condition
colnames(mat) <- colnames(rld)

hmcol <- colorRampPalette(brewer.pal(9,"Blues"))(255)
heatmap.2(mat,trace='none',col=rev(hmcol),margin=c(11,11))


#=========== Correlation plot for biological replicates  ==========================
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}


pairs(log2(data[,1:3]),lower.panel=panel.smooth, upper.panel=panel.cor,
      main='Human infected at 2h. log2(count)')


