##st4_combineDatasets.R
## Here we combine studies together and performs batch effect correction with ComBat to normalize them
#
rm(list=ls()); options(stringsAsFactors=F)

#source("http://bioconductor.org/biocLite.R")
library(ggplot2);library(sva); library(WGCNA);library(dplyr);library(stringr);
library(ggalluvial);library(Biobase); library(ggpubr)
library(Rmagic); library(Rtsne); library(colorspace);library(sva);
library(dendextend);library(SC3);library(pheatmap);library(phateR)
library(foreach);library(doParallel)
### set working directory

setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")

# options
enableWGCNAThreads()
allowWGCNAThreads()
condition = TRUE

### load normalized data from each disease

files = dir("./codes/", pattern = "_metadata")
n=length(files)
multiExpr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("./codes/",files[[i]],sep = ""))
  multiExpr[[i]]$datExpr = datExpr
  
  
  
}

genes = rownames(multiExpr[[1]]$datExpr)
for(i in 2:length(multiExpr)) genes = union(genes, rownames(multiExpr[[i]]$datExpr))

all_datExpr = data.frame(row.names = genes)
all_datMeta = data.frame(matrix(NA, nrow=0, ncol=8));
for(i in 1:length(multiExpr)) {
  # combine datExpr
  all_datExpr = cbind(all_datExpr, 
                      multiExpr[[i]]$datExpr[match(genes, rownames(multiExpr[[i]]$datExpr)),])
  
  # combine datMeta
  datMetaNA = as.data.frame(matrix(NA, ncol=8, nrow=nrow(multiExpr[[i]]$datMeta)))
  idx = pmatch(c("Study", "Subject_ID", "Dx", "Brain_Region", "Age", "Sex", "PMI", "Brain_lobe"),colnames(multiExpr[[i]]$datMeta))
  datMetaNA[,which(!is.na(idx))] = multiExpr[[i]]$datMeta[,na.omit(idx)]
  all_datMeta=rbind(all_datMeta, datMetaNA)
}

## paralle computing 
require(doParallel)
cores <- makeCluster(detectCores(), type='PSOCK')

system <- Sys.info()['sysname']

cl <- NULL
if (system == 'Windows') {
  cl <- makeCluster(getOption('cl.cores', cores))
  registerDoParallel(cl)
  registerDoSEQ()
  on.exit(stopCluster(cl))
} else {
  options('mc.cores' = cores)
  registerDoParallel(cores)
}


##
colnames(all_datMeta) = c("Study", "Subject_ID", "Dx", "Brain_Region", "Age", "Sex", "PMI", "Brain_lobe")

#rename normal to control
all_datMeta$Dx = gsub("Normal", "CTL",all_datMeta$Dx)
all_datMeta$Dx = factor(all_datMeta$Dx, 
                        levels = c("CTL", "AD","PD", "PA","PSP",
                                   "Scz", "ASD", "BP", "MDD"))
all_datMeta$PMI[is.na(all_datMeta$PMI)] = mean(all_datMeta$PMI, na.rm = T)
all_datMeta$Study <- as.factor(all_datMeta$Study)
all_datMeta$Sex = as.character(all_datMeta$Sex)
all_datMeta$Sex[all_datMeta$Sex == "Male"] = "male"
all_datMeta$Sex[all_datMeta$Sex == "Female"] = "female"
all_datMeta$Sex[all_datMeta$Sex == "unknown"] = "male"
all_datMeta$Sex <- as.factor(all_datMeta$Sex)
all_datMeta$Brain_lobe[all_datMeta$Brain_lobe == "DLPFC"] ="Frontal"
all_datMeta$Brain_Region = as.factor(all_datMeta$Brain_Region)
all_datMeta$Brain_lobe <- as.factor(all_datMeta$Brain_lobe)

#--remove duplicate samples
idx <- !duplicated(colnames(all_datExpr)) 
all_datExpr <- all_datExpr[,idx]#---10313 genes
all_datMeta <- all_datMeta[idx,] #4627 Samples, 3276 Subjects


# handle NAs

all_datExpr=  apply(all_datExpr,2,FUN = function(x){
  x[is.na(x)]=mean(x,na.rm=T)
  return(x)
})


#Normalize by Study
rownames(all_datMeta) = colnames(all_datExpr)
all_datExpr = as.matrix(all_datExpr)

datMeta = all_datMeta
# model
mod = model.matrix(~Dx+Sex, data=all_datMeta)
modcombat = model.matrix(~Dx + Sex,data = all_datMeta)

datExpr = ComBat(all_datExpr, batch=all_datMeta$Study,
                 mod=modcombat, prior.plots = F)

#####  Compute distance matrix and hierarchical clustering

if (condition){
  d.usa <- dist(t(datExpr), method = "euc")
  h.usa <- hclust(d.usa, method = "ward.D2")
}

h.usa$labels = colnames(datExpr)

##### Generate coloured dentrogram (indicate datasets)
source("./codes/A2R.R")

###
dend <- as.dendrogram(h.usa)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=9)

# add labels
labels(dend) <- paste(as.character(all_datMeta$Dx)[order.dendrogram(dend2)])

cols <- c("lightgrey","aquamarine","red","blue","chartreuse",
          "brown","blueviolet","yellow","orange")

ggplot(ggDend, labels = T,offset_labels = 15) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")

all_datMeta$color = as.numeric(all_datMeta$Dx)
all_datMeta$color = cols[all_datMeta$color]
all_datMeta$lobe_color = as.numeric(all_datMeta$Brain_lobe)
all_datMeta$lobe_color = cols[all_datMeta$lobe_color]

labels_colors(dend) <- all_datMeta$color[order.dendrogram(dend)]

# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend, hang_height=0.1)

# reduce the size of the labels:

dend <- set(dend, "labels_cex", 0.2) %>% 
  set("leaves_pch", 19) %>%
  set("leaves_col",all_datMeta$lobe_color[order.dendrogram(dend)]) %>% 
  set("leaves_cex",0.2) %>% set("branches_lty", 5) %>% 
  set("branches_lwd", 0.1) 

phyl = as.phylo(h.usa)
tree = ggtree(phyl,layout = "circular")
p = ggtree(h.usa, layout = "circular")
p = p %<+% all_datMeta
p = p + 
  geom_point(aes(color = Dx)) + 
  theme(legend.position = "right")

ggplot(p)+ geom_tree()+theme_tree()
df = all_datMeta[,c("Dx"), drop =F]
p1 <- gheatmap(tree, df)+
  scale_color_manual(values = cols)

# And plot:

pdf("./results/figures/samples_dendrogram2.pdf", width = 4,height = 4)
par(mar= c(0,0,0,0))
circlize_dendrogram(dend,labels = T,)

legend(x = 0.1, y= 0.4, legend = levels(all_datMeta$Dx),
       fill =cols, cex=0.5,border = F,bty = "n", title = "Disease")

legend(x = -0.3, y = 0.4, legend = levels(all_datMeta$Brain_lobe),
       cex=0.5, border = F,pch = 19, title = "Region",
       col = cols[1:7],bty = "n")

dev.off()

##
annot = data.frame(Disease = as.factor(datMeta$Dx),
                   Region = as.factor(datMeta$Brain_lobe))
rownames(annot)= colnames(datExpr)
pheatmap(mat = datExpr,
         show_rownames = F,
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = T,
         treeheight_row = 2,
         treeheight_col = 2,
         annotation_col = annot,
         filename = "./results/figures/heatmap_cluster.pdf"
)
Heatmap(datExpr)


# using tSNE
set.seed(123) # for reproducibility

###  post correction PCA ###
if(condition){
  pc.fit = prcomp(t(datExpr),scale. = T, center = F)
}
### PCA plots ####

dat = as.data.frame(pc.fit$x)
dat$Disease = all_datMeta$Dx
dat$Region = all_datMeta$Brain_lobe

pcDx = ggplot(dat, aes(x =PC1 ,
                       y =PC2,
                       color = Disease))+
  geom_point( alpha =0.8, size = 1)+
  scale_color_manual(values = cols)+
  theme_classic2()

pcRegion = ggplot(dat, aes(x =PC1 ,
                           y =PC2))+
  geom_point(aes(color = Region), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")+
  theme_classic()

pcPlots = ggarrange(pcDx, pcRegion,nrow = 1)


### tSNE ####

tsne <- Rtsne(dat[,1:10], 
              #  perplexity=500,
              #  verbose=TRUE, 
              # max_iter = 1000,
              #  dims = 3,
              pca = F
              # partial_pca = T,
              # eta = 100
)
data = as.data.frame(tsne$Y)
data$Dx = all_datMeta$Dx
data$Study =all_datMeta$Study 
data$region = all_datMeta$Brain_lobe


colnames(data)= c("tSNE1", "tSNE2", "Disease","Study","Region")
postDx = ggplot(data, aes(x =tSNE1 ,
                          y =tSNE2,
                          color = Disease))+
  geom_point( alpha =0.8, size = 1)+
  scale_color_manual(values = cols)+
  theme_classic2()

postregion = ggplot(data, aes(x =tSNE1 ,
                              y =tSNE2))+
  geom_point(aes(color = Region), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")+
  theme_classic()

postplot = ggarrange(postDx, postregion,nrow = 1)
ggsave(plot = postplot, filename = "./results/figures/post_tSNE_allgenes.pdf",
       width = 11, height = 5)

postStudy = ggplot(data, aes(x =tSNE1 ,
                             y =tSNE2))+
  geom_point(aes(color = Study), size = 1)+
  # scale_color_brewer(type = "div",palette = "Set1")+
  labs(x = "tSNE1", y = "tSNE2",
       title = "Post-correction study")+
  theme_bw()

library(ggpubr)
pcPlots = ggarrange(preDx, preregion, preStudy,
                    postDx, postregion, postStudy, 
                    ncol = 3, nrow = 2)

ggsave(filename = "./results/figures/PCA_Plots_combined.pdf",
       plot = pcPlots, width = 20, height = 12)



### PCA using union of DEGs ####
# load data #

load("./codes/DEGs.Rdata")

# filter data for DEGs
id  = unique(sigData$id)
sigExp =  datExpr[id,]
sigPC = prcomp(t(sigExp), scale. = F, center = T)

sig.mds = cmdscale(dist(t(sigExp)), eig = T)
dat = as.data.frame(sig.mds$points)
dat$Disease = datMeta$Dx
dat$Region = datMeta$Brain_lobe

# pca plot

sigDxPC = ggplot(dat, aes(x =V1 ,
                          y =V2,
                          color = Disease))+
  geom_point( alpha =0.8, size = 1)+
  scale_color_manual(values = cols)+
  labs(x = "PC1",
       y = "PC2",
       title = "PCA-union of DEGs")+
  theme_classic2()

sigRegPC = ggplot(dat, aes(x =V1 ,
                           y =V2))+
  geom_point(aes(color = Region), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")+
  labs(x = "PC1",
       y = "PC2",
       title = "PCA-using union of DEGs")+
  theme_classic()



# tsne 
# fist do a PCA to reduce dimensions
sigPC = prcomp(t(sigExp), scale. = F, center = T)
sigTSNE = Rtsne(sigPC$x[,1:10])

dat = as.data.frame(sigTSNE$Y)
colnames(dat) = c("tSNE1","tSNE2")
dat$Disease = datMeta$Dx
dat$Region = datMeta$Brain_lobe

# plot

sigDxTsne = ggplot(dat, aes(x =tSNE1 ,
                            y =tSNE2,
                            color = Disease))+
  geom_point( alpha =0.8, size = 1)+
  scale_color_manual(values = cols)+
  theme_classic2()

sigRegTsne = ggplot(dat, aes(x =tSNE1 ,
                             y =tSNE2))+
  geom_point(aes(color = Region), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")+
  theme_classic()


sigPlots = ggarrange(sigDxPC,sigRegPC,sigDxTsne,sigRegTsne, ncol = 2, nrow = 2)
ggsave(plot= sigPlots, filename = "./results/figures/PCA_allDEGs.pdf",
       width = 14,height = 10,device = "pdf")

## PCA using stringent genes (p < 0.05 & abs(logFC) > 0.5))

stringentGenes = subset(alldx,abs(logFC > 0.5) & p.value < 0.05)
idx = unique(stringentGenes$id)
strExp = datExpr[idx,]

mds = cmdscale(dist(t(strExp)), eig = T)
dat = as.data.frame(mds$points)
dat$Disease = datMeta$Dx
dat$Region = datMeta$Brain_lobe

# pca plot

strDxPC = ggplot(dat, aes(x =V1 ,
                          y =V2,
                          color = Disease))+
  geom_point( alpha =0.8, size = 1)+
  scale_color_manual(values = cols)+
  labs(x = "PC1",
       y = "PC2",
       title = "PCA-top stringent DEGs")+
  theme_classic2()

strRegPC = ggplot(dat, aes(x =V1 ,
                              y =V2))+
  geom_point(aes(color = Region), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")+
  labs(x = "PC1",
       y = "PC2",
       title = "PCA-top stringent DEGs")+
  theme_classic()

# tsne 
# fist do a PCA to reduce dimensions
strPC = prcomp(t(strExp), scale. = F, center = T)
strTSNE = Rtsne(strPC$x[,1:10])

dat = as.data.frame(strTSNE$Y)
colnames(dat) = c("tSNE1","tSNE2")
dat$Disease = datMeta$Dx
dat$Region = datMeta$Brain_lobe

# plot

strDxTsne = ggplot(dat, aes(x =tSNE1 ,
                          y =tSNE2,
                          color = Disease))+
  geom_point( alpha =0.8, size = 1)+
  scale_color_manual(values = cols)+
  labs( title = "tSNE-top stringent DEGs")+
  theme_classic2()

strRegTsne = ggplot(dat, aes(x =tSNE1 ,
                           y =tSNE2))+
  geom_point(aes(color = Region), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")+
  labs(title = "tSNE-top stringent DEGs")+
  theme_classic()

strPlots = ggarrange(strDxPC,strRegPC,strDxTsne, strRegTsne,
                     nrow = 2, ncol = 2)
ggsave(filename = "./results/figures/PCA_stringent_DEGs.pdf",
       width = 14, height = 10)


degData = list(strExp,sigExp)
names(degData) = c("topDEG", "unionDEG")

#####  hierarchical clustering using stringent genes

for (i in names(degData)){
  expData = degData[[i]]
if (condition){
  strDis<- dist(t(expData), method = "euc")
  strHc <- hclust(strDis, method = "ward.D2")
}

str.dend <- as.dendrogram(strHc)

# Color the branches based on the clusters:
str.dend <- color_branches(str.dend, k=9)

# add labels
labels(str.dend) <- paste(as.character(all_datMeta$Dx)[order.dendrogram(str.dend)])
labels_colors(str.dend) <- all_datMeta$color[order.dendrogram(str.dend)]
# We hang the dendrogram a bit:
str.dend <- hang.dendrogram(str.dend, hang_height=20)

# reduce the size of the labels:

str.dend <- set(str.dend, "labels_cex", 0.2) %>% 
  set("leaves_pch", 19) %>%
  set("leaves_col",all_datMeta$lobe_color[order.dendrogram(str.dend)]) %>% 
  set("leaves_cex",0.2) %>% set("branches_lty", 5) %>% 
  set("branches_lwd", 0.01) 


pdf( paste("./results/figures/samples_dendrogram_",i,".pdf",sep = ""),
           width = 4,height = 4)
par(mar= c(0,0,0,0))
circlize_dendrogram(str.dend,labels = T,)

legend(x = 0.1, y= 0.4, legend = levels(all_datMeta$Dx),
       fill =cols, cex=0.5,border = F,bty = "n", title = "Disease")

legend(x = -0.3, y = 0.4, legend = levels(all_datMeta$Brain_lobe),
       cex=0.5, border = F,pch = 19, title = "Region",
       col = cols[1:7],bty = "n")

dev.off()

}

#clustering using phateR package
# phate.test = phate(t(datExpr))
# 
# # runs phate with different parameters
# tree.phate <- phate(t(datExpr), gamma=1 ,
#                     t= 10,knn = 5,
#                     init=phate.test)
# #PLOT
# ggplot(tree.phate, aes(x=PHATE1, y=PHATE2, color=datMeta$Dx)) +
#   geom_point()+
#   scale_color_manual(values = cols)+
#   theme_clean()
# 

### projection score
#computing variance for each gene
datVar= data.frame(matrix(NA, nrow = nrow(datExpr),ncol = 2))

var =apply(datExpr,1,var)
var =as.data.frame(unlist(var)) 

datVar$gene = rownames(var)
datVar$var = var$`unlist(var)` 
datVar = datVar[order(datVar$var, decreasing = T),]

# save 
save(file = "./codes/projectionData.Rdata",datExpr,datVar)
write.table(datVar,file = "./results/tables/GeneVariance.tsv", col.names = T)
write.csv(datExpr,  file = "./results/tables/datExpr.csv")

# filter for top Var
datVar = subset(datVar,var > 1)

#load top variable genes
geneSet = read.table("./results/tables/projectionScore.txt")

#filter data for selected genes
topID = datVar$gene
topVar = datExpr[topID,]
selected = as.table(datExpr[geneSet$V1,])

# PCA 
projPCA = prcomp(t(topVar),scale. = F, center = T)
projPcaData = as.data.frame(projPCA$x)

# tsne
proj.tsne = Rtsne(projPCA$x[,1:10],num_threads = 4,)
dims = as.data.frame(proj.tsne$Y)
colnames(dims) = c("PC1","PC2")
pcaList = list(projPcaData,dims)

for (i in 1:length(pcaList)) {
    pcaList[[i]]$disease  = datMeta$Dx
    pcaList[[i]]$region = datMeta$Brain_lobe
}


plotList = vector("list",2)
names(plotList) = paste(rep(c("disease","region"),each=2),c(1,2),sep = "_")

for (i in 1:2) {
    dat = pcaList[[i]]
    
     p1 = ggplot(dat, aes(x = PC1, y = PC2
                           ))+
     geom_point(size = 1,
                aes(color = disease))+
     scale_color_manual(values = cols)+
     labs(x = ifelse(i == 1, "PC1","tSNE1"),
          y = ifelse(i == 1, "PC2", "tSNE2"))+
       theme_classic()
     
     p2 = ggplot(dat, aes(x = PC1, y = PC2
     ))+
       geom_point(size = 1,
                  aes(color = region))+
       scale_color_manual(values = cols)+
       labs(x = ifelse(i == 1, "PC1","tSNE1"),
            y = ifelse(i == 1, "PC2", "tSNE2"))+
       theme_classic()
     
 p = ggarrange(p1,p2, nrow = 1)
   # name = paste0(var,"_",i)
   # print(name)
   plotList[[i]] = p
    
}

 p = ggarrange(plotlist = plotList,nrow = 2)

ggsave("./results/figures/projectionPCA.pdf", width = 10, height = 8)
 
# run pca usign top genes
top_genes = read.csv("./results/tables/top_genes.csv",col.names = F)
topExpr = datExpr[top_genes$FALSE.,]

set.seed(1) # for reproducibility
topGenes_tsne <- Rtsne(t(topExpr),dims = 2)
top_data = as.data.frame(topGenes_tsne$Y)

ggplot(top_data, aes(x =top_data[,1] ,
                     y =top_data[,2] ))+
  geom_point(aes(color = datMeta$Dx), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")



save.image("./codes/Combined.Rdata")
save(file ="./codes/st4_combinedData.Rdata", datMeta, datExpr, all_datExpr)






