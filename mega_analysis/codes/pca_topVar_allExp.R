## pca_topVar_allExp.R
#
rm(list=ls()); options(stringsAsFactors=F)

library(ggplot2);library(sva); 
library(WGCNA);library(dplyr);
library(stringr);library(RColorBrewer)
library(ggpubr);library(edgeR);
library(Rmagic); 
library(Rtsne);
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
geneRank = read.csv("./results/tables/allExpr_gene_ranks.tsv",
                          sep = "\t",row.names = 1)

## Create a list of gene names sorted by decreasing variance
sorted.names <- rownames(geneRank)[order(geneRank$var, decreasing=TRUE)]


## Select the 200 top-ranking genes sorted by decreasing variance
top.variables <- 200
selected.genes <- sorted.names[1:top.variables]
rm(sorted.names)

# pca 
datExp = allExp[selected.genes,]
lcpm = cpm(datExp, log = T)
set.seed(123)
pc = prcomp(t(lcpm))

pc.data = pc$x[,1:10]
dup = duplicated(pc.data)
pc.data = pc.data[!dup,]

# tsne
set.seed(123)
tsn = Rtsne(pc.data, pca = F,theta = 0)
data.tsne = as.data.frame(tsn$Y)
rownames(data.tsne) = rownames(pc.data)
colnames(data.tsne) = c("tSNE1","tSNE2")
data.tsne$Dx = allMeta$Dx[match(rownames(data.tsne),rownames(allMeta))]
data.tsne$study = allMeta$Study[match(rownames(data.tsne),rownames(allMeta))]
data.tsne$Brain_Lobe = allMeta$Brain_Lobe[match(rownames(data.tsne),rownames(allMeta))]


p1  = mygg(data =data.tsne,
           color = data.tsne$Dx, title = "disease")+
  scale_color_manual(values = cols)
p2  = mygg(data =data.tsne,color = data.tsne$study, title = "Study")
p3  = mygg(data =data.tsne,color = data.tsne$Brain_Lobe, title = "Brain lobe")+
  scale_color_manual(values = cols)


plots = ggarrange(p1,p2,p3, ncol = 1, nrow = 3)
ggsave(plot = plots, "./results/figures/pca/pca_topVar_preBatch.pdf",
       width = 7, height = 10) 


#### pca for batch-effect corrected data ####
load("./codes/combined_batch_corrected.Rdata")

exp = datExp.comb[selected.genes,]

#### tsne

set.seed(123)
tsn = Rtsne(t(exp), pca = T,
            theta = 0,
            check_duplicates = F,
            max_iter = 5000)
data.tsne = as.data.frame(tsn$Y)
rownames(data.tsne) = colnames(exp)
colnames(data.tsne) = c("tSNE1","tSNE2")
data.tsne$Dx = allMeta$Dx
data.tsne$study = allMeta$Study
data.tsne$Brain_Lobe = allMeta$Brain_Lobe

p1  = mygg(data =data.tsne,
           color = data.tsne$Dx, title = "disease")+
  scale_color_manual(values = cols)
p2  = mygg(data =data.tsne,color = data.tsne$study, title = "Study")
p3  = mygg(data =data.tsne,color = data.tsne$Brain_Lobe, title = "Brain lobe")+
  scale_color_manual(values = cols)

p = ggarrange(p1,p2,p3,nrow = 3)
ggsave("./results/figures/pca/pca_topVar_corrected.pdf",width = 8,height = 12)
  