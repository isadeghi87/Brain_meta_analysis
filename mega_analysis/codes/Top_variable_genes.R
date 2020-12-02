## Here we compute top variable genes 
rm(list=ls()); options(stringsAsFactors=F)

library(ggplot2);library(sva); 
library(WGCNA);library(dplyr);
library(stringr);library(RColorBrewer)
library(ggpubr);library(edgeR)
library(Rmagic); library(Rtsne);
library(colorspace);library(sva);

### set working directory

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")

# options
enableWGCNAThreads()
allowWGCNAThreads()
condition = TRUE

### load normalized data from each disease

load("./codes/st2_lm_multiExpr.Rdata")
rm(multiExpr)

allMeta$Brain_Lobe = gsub("Basal_ganglia","Basal ganglia",
                          allMeta$Brain_Lobe)
allMeta$Brain_Lobe = as.factor(allMeta$Brain_Lobe)

var = c("Brain_Lobe","Brain_Region","Study","Sex")
for (i in var){allMeta[,i] = factor(allMeta[,i])}
allMeta$Dx = factor(allMeta$Dx, levels = c("CTL","AD","PD","PA","PSP","Scz","ASD","BP","MDD"))

dx = as.numeric(as.factor(allMeta$Dx))
st = as.numeric(as.factor(allMeta$Study))
reg = as.numeric(as.factor(allMeta$Brain_Lobe))

# function for plot 
cols <- c("lightgrey","aquamarine","red","blue","chartreuse",
          "brown","blueviolet","yellow","orange")
mygg = function(data,color,title){
  plt = ggplot(data, aes(x =tSNE1 ,
                         y =tSNE2,
                         color = color))+
    geom_point( alpha =0.7, size = 1)+
    theme_classic2()+
    labs(title = title)
  return(plt)
}


#### filter genes ####

keep = rowSums(cpm(allExp)>0.5)>= 0.3 * ncol(allExp)
datExp = allExp[keep,]

## Compute sample-wise variance
var.per.sample <- apply(datExp, 2, var)
head(var.per.sample)

## Inspect the distribution of sample-wise variance
hist(var.per.sample, breaks=20)

## Compute gene-wise variance
var.per.gene <- apply(datExp, 1, var)

## Inspect the distribution of gene-wise variance
hist(var.per.gene, breaks=100)

## Sort genes per decreasing variance
genes.by.decr.var <- sort(var.per.gene,decreasing=TRUE)

## Print the 5 genes with highest variance
head(genes.by.decr.var)

## Select the 50 top-ranking genes in the list sorted by variance.
## This list of genes will be used below as training variables for
## supervised classification.
top.nb <- 50 ## This number can be changed for testing
genes.selected.by.var <- names(genes.by.decr.var[1:top.nb])

## Check the names of the first selected genes
head(genes.selected.by.var, n=20) 

## Create a data frame to store gene values and ranks 
## for different selection criteria.
gene.ranks <- data.frame(var=var.per.gene)
head(gene.ranks)

## Beware, we rank according to minus variance, because 
## we want to associate the lowest ranks to the highest variances
gene.ranks$var.rank <- rank(-gene.ranks$var, ties.method='random')
head(gene.ranks)

## Check the rank of the 5 genes with highest and lowest variance, resp.
gene.ranks[names(genes.by.decr.var[1:5]),]

## Plot the expression profiles of the two genes with highest variance
g1 <- names(genes.by.decr.var[1])
print(g1)

(g2 <- names(genes.by.decr.var[2]))

x <- as.vector(as.matrix(datExp[g1,]))
y <- as.vector(as.matrix(datExp[g2,]))
plot(x,y,
     col=as.numeric(as.factor(allMeta$Dx)),
     # type='n',
     # panel.first=grid(col='black'), 
     main="2 genes with the highest variance", 
     xlab=paste('gene', g1), 
     ylab=paste('gene', g2))
text(x, y,labels=sample.labels,col=sample.colors,pch=0.5)
legend('bottomleft',
       col=dx,
       legend=levels(allMeta$Dx),pch=1,cex=0.7,bg='white',bty='o')

## Specify the URL of the base for the course
dir.base <- 'http://www.bigre.ulb.ac.be/courses/statistics_bioinformatics'

## Load the general configuration file
source('./codes/config.R')

## Load the general configuration file
source("./codes/util_student_test_multi.R")

## Define a vector indicating whether each sample 
## belongs to the subtype of interest (e.g. "Bo") or not.
sample.labels = as.character(allMeta$Dx)
names(sample.labels) = colnames(datExp)
current.group <- "CTL"
one.group.vs.others<- sample.labels
one.group.vs.others[sample.labels != current.group] <- "other"
print(table(one.group.vs.others))

## Test the mean equality between ctl and all others 
welch.one.group.vs.others <- t.test.multi(datExp, one.group.vs.others)


## Update the gene rank table
test.name <- paste(current.group, '.vs.others.sig', sep='')
gene.ranks[,test.name] <- welch.one.group.vs.others$sig ## Add a column with significances
gene.ranks[,paste(test.name, ".rank", sep="")] <- rank(-welch.one.group.vs.others$sig,
                                                       ties.method='random') ## Add a column with significance ranks

## Do the same Welch test for other groups
grp = c("AD","PSP","PA","PD","ASD","Scz","BP","MDD")
for (current.group in grp) {
  one.group.vs.others<- sample.labels
  one.group.vs.others[sample.labels != current.group] <- "other"
  
  ## Test the mean equality between Bo subtype and all other subtypes 
  welch.one.group.vs.others <- t.test.multi(datExp, one.group.vs.others)
  
  ## Update the gene rank table
  test.name <- paste(current.group, '.vs.others.sig', sep='')
  gene.ranks[,test.name] <- welch.one.group.vs.others$sig
  gene.ranks[,paste(test.name, ".rank", sep="")] <- rank(-welch.one.group.vs.others$sig, ties.method='random')
}


## Check the resulting gene table
head(gene.ranks)
head(gene.ranks[order(gene.ranks$CTL.vs.others.sig.rank),])

## Store the gene rank table in a text file (tab-separated columns)
write.table(gene.ranks, file= "./results/tables/allExpr_gene_ranks.tsv", 
            sep='\t', quote=F, col.names=NA)

## Plot variance against significance of the Welch test for the 4 different groups
plot(gene.ranks[,c("var", "CTL.vs.others.sig","AD.vs.others.sig",
                    "PSP.vs.others.sig",
                   "PA.vs.others.sig","PD.vs.others.sig",
                   "Scz.vs.others.sig","ASD.vs.others.sig",
                   "BP.vs.others.sig","MDD.vs.others.sig")], col="grey")


#### PCA ####

## load the stats library to use the princomp() and prcomp() function
library(stats) 

## Perform the PCA transformation using all genes and seleclted gens
exp = datExp[selected.genes,]

## log trnasform
lcpm = cpm(datExp, log = T)
expr.prcomp <- prcomp(t(lcpm),cor=TRUE)


## Analyze the content of the prcomp result: 
## the result of the method prcomp() is an object 

plot(expr.prcomp, 
     main='Variance  per component', 
     xlab='Component')

## Get the standard deviation and variance per principal component
sd.per.pc <- expr.prcomp$sdev
var.per.pc <- sd.per.pc^2

## Display the percentage of total variance explained by each 
sd.per.pc.percent <- 100 * sd.per.pc/sum(sd.per.pc)
var.per.pc.percent <- 100 * var.per.pc/sum(var.per.pc)
pdf("./results/figures/pca/Variance_perPCA_allExp.pdf")
barplot(var.per.pc.percent[1:10], main='Percent of variance  per component',
        xlab='Component', 
        ylab='Percent variance', 
        col='#BBDDFF')
dev.off()


