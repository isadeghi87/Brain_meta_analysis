#st5a_networkAnalysis.R

##CrossDisorder: step 5a_WGCNA on mega analysis
##This script takes a mega-analysis of gene expression studies, datExpr,
##and runs WGCNA to find modules of co-expressed genes


rm(list=ls())
library(WGCNA)
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")
load("./codes/st4_combinedData.Rdata")



## NETWORK ANALYSIS
## ----------------
multiExpr = vector(mode="list", length=1)
multiExpr[[1]] = list(data=as.data.frame(t(datExpr)))
bsize = 5000
nSets = 1
powers = c(seq(1,9,by=1),seq(10,30,by=2))
enableWGCNAThreads()
allowWGCNAThreads()
n=1

## Compute soft-threshold for mega-analysis


  
  multiExpr[[n]]$softThresh = pickSoftThreshold(data= multiExpr[[n]]$data, 
                                                networkType = "signed", 
                                                corFnc = "bicor",
                                                verbose = 5,
                                                powerVector = powers,
                                                blockSize = bsize)
  
  sft = multiExpr[[n]]$softThresh
  
  pdf("./results/figures/WGCNA/WGCNA.softThresh.pdf")
  par(mfrow=c(1,2))
  n = 1
  plot(sft$fitIndices[,1], 
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
       xlab="Soft Threshold Power",
       ylab="Scale free R^2",
       type="n")
  text(sft$fitIndices[,1],
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels = powers,
       cex = 0.7,
       col="blue",  
       xlab="Soft Threshold Power",
       ylab="Scale free R^2")
  abline(h=0.8, col="black")
  plot(sft$fitIndices[,1], 
       sft$fitIndices[,5],
       xlab = "Soft threshold power",
       ylab = "Mean connectivity",
       type = "n")
  text(sft$fitIndices[,1], 
       sft$fitIndices[,5],
       labels = powers,
       cex = 0.7, 
       col="black")
  dev.off()



wgcna_parameters = list(powers =  10)
wgcna_parameters$minModSize = 40  
wgcna_parameters$minHeight = 0.1
wgcna_parameters$bsize = 25000  ##block size needs to be larger than dim(datExpr)[1]
wgcna_parameters$ds = 2  ##deep split parameter contorls number of modules
wgcna_parameters$networkType = "signed"    ## using signed networks
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = TRUE

n = 1

## ROBUST WGCNA 
## --> we want to generate a "consensus TOM" based on resampling datExpr so as to make 
## --> module definitions robust to outliers
 multiExpr[[n]]$netData = 
  blockwiseModules(datExpr=multiExpr[[n]]$data, 
                   maxBlockSize=wgcna_parameters$bsize,
                   networkType=wgcna_parameters$networkType, 
                   corType = wgcna_parameters$corFnc,
                   power = wgcna_parameters$powers[n], 
                   mergeCutHeight= wgcna_parameters$minHeight, 
                   nThreads=4, 
                   saveTOMFileBase=paste("./results/figures/WGCNA/network_signed", 
                                         wgcna_parameters, "_exprSet", 
                                         as.character(n), sep=""), 
                   saveTOMs=TRUE, 
                   minModuleSize= wgcna_parameters$minModSize,
                   pamStage=wgcna_parameters$pamStage, 
                   reassignThreshold=1e-6, 
                   verbose = 3, 
                   deepSplit=wgcna_parameters$ds)



