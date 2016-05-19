# stsPlots.2.0.R
#
# Author: Meredith Ashby, Lawrence Lee
# Date: 03/14
#
# Updated to work with current R version and updated libraries
#
# R version 3.0.3 (2014-03-06)
# ggplot2_0.9.3.1    plyr_1.8.1         Hmisc_3.14-3       
# Biostrings_2.30.1  
#
# Description: These are a set of plots generated from the sts.csv file that help in assessing loading of the SMRTcell and pre-mapping quality.

require(ggplot2)
require(methods)
require(reshape2)
require(plyr)

PacBioCols <- c(rgb(156, 19,46,maxColorValue=255),#red 
               rgb(49, 94,161,maxColorValue=255),# dkblue 
               rgb(242,175,50,maxColorValue=255), # yellow
               rgb(142, 186,200,maxColorValue=255), #ltblue
                rgb(0, 0, 0, maxColorValue=255),#black 
               rgb(116, 118, 121,maxColorValue=255),# ltgrey
               rgb(51, 51, 51,maxColorValue=255), #dkgrey 
               rgb(255, 255,255, maxColorValue=255))#white


generateStsPlots <- function(stsCsvPath, pdfOutPath) {

  pdf(pdfOutPath, width=11, height=8.5, onefile=T)
  r <- read.table(stsCsvPath, header=TRUE, sep=",", strip.white=T)
  rn <- subset(r, ZmwType=="SEQUENCING")
  rn$HQReadLength <- rn$HQRegionEnd - rn$HQRegionStart

  
  # PAGE 1:  ReadScore heatmap 
  rn$BinnedReadScore <- round_any(rn$ReadScore, 0.75, floor)

  instructions <- "Confirm an even distribution\nof pass / fail Readscores."
  annotation <- paste("ZMWs passing filter:\n", sum(rn$ReadScore >= 0.75), ' ZMWs \n',format(sum(rn$ReadScore >= 0.75) / length(rn$ReadScore) * 100, digits=4), " %\n", sep="")

  
  d <- qplot(data=rn, x=X, y=Y, geom="tile", fill=factor(BinnedReadScore)) +
      scale_fill_manual(values=c(rgb(156, 19,46,maxColorValue=255), rgb(242,175,50,maxColorValue=255))) +
      labs(title="ReadScore Heatmap : Loading / Alignment Metric", fill="Binned Readscore") +
      coord_cartesian(xlim=c(-175,250)) + theme(legend.position = c(0.1,0.9)) +
      annotate(geom='text', x=240, y=-170, vjust=0, hjust=1,label=annotation) +  
      annotate(geom='text', x=240, y=240, vjust=1, hjust=1,label=instructions)   

  show(d)

  print("plot 1 done")

  # PAGE 2:  Productivity heatmap 
  prod0 <- sum(rn$Productivity==0)
  prod1 <- sum(rn$Productivity==1)
  prod2 <- sum(rn$Productivity==2)

  prod0Pct <- paste("Prod=0 ZMWs: ", prod0, " (",format(prod0 / length(rn$Productivity) * 100, digits=4), "%)", sep="")
  prod1Pct <- paste("Prod=1 ZMWs: ", prod1, " (",format(prod1 / length(rn$Productivity) * 100, digits=4), "%)", sep="")
  prod2Pct <- paste("Prod=2 ZMWs: ", prod2, " (",format(prod2 / length(rn$Productivity) * 100, digits=4), "%)", sep="")
  activePct <- paste("Total Active ZMWs: ",(prod1+prod2), " (", format((prod1 + prod2) / length(rn$Productivity) * 100, digits=4), "%)", sep="")
  instructions <- "Poisson distribution indicates that\nProd=1 ZMWs are maximized at 36.7%\nif 63.2% of total ZMWs are active.\n"
  annotation <- paste(prod0Pct, '\n', prod1Pct, '\n', prod2Pct, '\n', activePct, sep="")
  
  d <- qplot(data=rn, x=X, y=Y, geom="tile", fill=factor(Productivity)) +
      scale_fill_manual(values=c(rgb(156, 19,46,maxColorValue=255), rgb(242,175,50,maxColorValue=255), rgb(49, 94,161,maxColorValue=255))) +
      labs(title="Productivity Heatmap : Loading Metric", fill="Productivity") +
      coord_cartesian(xlim=c(-175,250)) + theme(legend.position = c(0.1,0.9)) +
      annotate(geom='text', x=240, y=min(rn$Y), vjust=0, hjust=1,label=annotation) +
      annotate(geom='text', x=240, y=max(rn$Y), vjust=1, hjust=1,label=instructions)
  
  show(d)
  print("plot 2 done")

  
  # PAGE 3: ReadScore Distribution vs. Productivity

  d <- qplot(data=subset(rn, rn$ReadScore > 0.1), x=ReadScore, geom="bar", fill=factor(Productivity)) +
    labs(title="ReadScore Distribution by Productivity : Sequencing Quality Metric", y="Counts", fill="Productivity") +
    geom_vline(xintercept = 0.75, color="black") +
    facet_grid(Productivity~.) + theme(legend.position = c(0.1,0.9))  +
    scale_fill_manual(values=PacBioCols)
  show(d)
  print("plot 3 done")


  # PAGE 4: HQReadLength boxplot 
  p <- ggplot(rn[rn$ReadScore > 0.1,], aes(x=factor(round_any(HQReadLength, 500)), y=ReadScore, color=as.factor(Productivity))) +
    labs(title="ReadScore Distribution by HQ Readlength : Sequencing Quality Metric", x="Binned HQ Readlength", color="Productivity") +
    theme(axis.text.x=element_text(angle=-90, hjust=0)) +
    geom_boxplot(outlier.size=0) + theme(legend.position = c(0.9,0.15)) +
    scale_colour_manual(values=PacBioCols)
  show(p)
  print("plot 4 done")

  
  # PAGE 5:Comparison of raw vs HQ region only readlength distributions - How much of the sequencing-zmw data is crap?
  hq <- melt(rn, measure=c("NumBases","HQReadLength"), id=c("X","Y","ZmwType","ReadScore"))
  colnames(hq) <- c("X","Y","ZmwType","ReadScore","Raw_vs_HQ_RL","value")
  instructions <- "A large discrepancy between raw and\nHQ read lengths indicates noisy ZMWs.\n\nCompare HQ read length plot to expected\n readlength values for your movie length."

  p <- qplot(data=subset(hq, hq$ReadScore > 0.1), x=value, geom='freqpoly', color=Raw_vs_HQ_RL ) +
    labs(x="Length in Bases", y="Counts", title="Raw vs HQ Region Readlength Distributions : Sequencing Quality Metric", color="Raw vs. HQ Readlength") +
    theme(legend.position = c(0.9,0.9)) +
    scale_colour_manual(values=PacBioCols) +
    annotate(geom='text', x=max(hq$value), y=0, vjust=0, hjust=1,label=instructions)

  show(p)
  print("plot 5 done")

  
  # PAGE 6: What fraction of called bases are in an HQ region?
  rn$HQFraction <- rn$HQReadLength / rn$NumBases
  instructions <- "A large discrepancy between raw and\nHQ readlengths indicates noisy ZMWs."
  p <- qplot(data=subset(rn, ReadScore > 0.1), x=HQFraction, geom='histogram') +
    labs(x="Fraction of Called Bases that are in an HQRegion", y="Counts", title="Per ZMW Sequencing / Noise Ratio") +
    annotate(geom='text', x=0, y=500, vjust=0, hjust=0, label=instructions) +
    scale_colour_manual(values=PacBioCols)

  show(p)
  print("plot 6 done")

  
  # PAGE 7: SNR heatmap
  mb <- melt(rn, measure.var=c("SnrMean_T", "SnrMean_G", "SnrMean_A", "SnrMean_C"), id.var=c("X", "Y", "ZmwType", "ReadScore"))
  q1 <- quantile(subset(mb, ZmwType=="SEQUENCING")$value, 0.01, na.rm=TRUE)
  q99 <- quantile(mb$value, 0.99, na.rm=TRUE)
  mb$Baseline <- pmax(0.9*q1, pmin(mb$value, q99))

  p <- qplot(X,Y, fill=Baseline, geom="tile", data=mb, main = "Mean SNR Heatmap : Alignment Metric") +
    facet_wrap(~ variable) + 
    scale_fill_gradientn(colours = rainbow(7))
  
  show(p)
  print("plot 7 done")


  # PAGE 8: SNR distributions
  instructions <- data.frame(variable=c('SnrMean_T','SnrMean_G','SnrMean_A','SnrMean_C'), Annotation=c('Minimum: 4.0\nOptimum: 4.5-7', '','Minimum: 5.5\nOptimum: 7-10',''))

  p <- qplot(data=subset(mb, ReadScore > 0.1), x=value, geom='freqpoly', color=variable) +
    labs(x="ZMW Mean SNR", y="Counts", title="Per Channel Mean SNR Distribution", color="Channel SNR") +
    facet_wrap(~variable) + coord_cartesian(xlim=c(0,16)) +
    scale_x_continuous(breaks=seq(2,14,2)) +
    theme(legend.position = "none") +
    geom_text(aes(label=Annotation), data=instructions, x=1, y=0, vjust=0, hjust=0, show_guide=FALSE) +
    scale_colour_manual(values=PacBioCols)
  show(p)
  print("plot 8 done")
  
  
  # PAGE 9: IPD histogram : oxygen exclusion test

  prod <- subset(rn,ReadScore >= 0.75 & Productivity == 1 & BaseIpd < 2)
  ipdMean <- mean(prod$BaseIpd)
  interpretation <- ifelse(ipdMean<0.25,'Normal','Slow')
  annotation <- paste('Mean IPD(<0.25 is normal): ', round_any(ipdMean,0.01), 's\nnReads(>22k is normal): ',dim(prod)[1], '\nMode of distribution should\n be to the left of vertical line\n\n',interpretation, '\n',sep="")
  
  p <- qplot(data=prod,x=BaseIpd, y=..count.., binwidth=0.02, theme="bw", geom='freqpoly') +
    labs(title="IPD Histogram : Oxygen Exclusion", x="Inter-pulse Duration", y="Counts") +
    geom_vline(xintercept=0.25,lwd=1,color=ifelse(ipdMean<0.25,'blue','red')) +
#    scale_x_continuous("Base IPD (s)") +
    annotate(geom='text', x=max(prod$BaseIpd), y=0, vjust=0, hjust=1,label=annotation)
  
  show(p)
  print("plot 9 done")

  dev.off()
}

plotStsCsvFolder <- function(stsCsvFolder) {
  outputFolder <- file.path(stsCsvFolder, "Analysis_Results")
  tryCatch(dir.create(outputFolder), simpleWarning = function(e) stop(paste("Unable to create directory: ", outputFolder)))

  csvFiles <- list.files(path = stsCsvFolder, pattern = "*.sts.csv")
  print(csvFiles)

  for (stsCsvFile in csvFiles) {
    print(stsCsvFile)
    pdfOutput <- file.path(outputFolder, paste(strsplit(stsCsvFile, split='\\.csv')[[1]][1], '.pdf', sep=""))
    generateStsPlots(file.path(stsCsvFolder,stsCsvFile), pdfOutput)
  }
}

runStsPlots <- function(args) {

  if (args[1] !="-file" & args[1]!="-folder")
    print("stsPlots requires either the -file or -folder flag  Ex: Rscript --vanilla stsPlots.2.0.R -file yourfile.sts.csv")

  inputType <- args[1]

  if (inputType=="-file") {
    if (length(args)!=2) print("You must supply a file name as a second argument.")
    stsCsvFile <- args[2]
    pdfOutput <- paste(strsplit(stsCsvFile, split='\\.csv')[[1]][1],'.pdf',sep="")
    generateStsPlots(stsCsvFile,pdfOutput)
  }

  if (inputType == "-folder") {
    if (length(args)!=2) print("You must supply a folder path as a second argument")
    stsCsvFolder <- args[2]
    plotStsCsvFolder(stsCsvFolder)
  }
}

args <- commandArgs(trailingOnly = T)
runStsPlots(args)
