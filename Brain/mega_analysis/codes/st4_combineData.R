##st4_combineDatasets.R
## Here we combine studies together and performs batch effect correction with ComBat to normalize them
#
rm(list=ls()); options(stringsAsFactors=F)

#source("http://bioconductor.org/biocLite.R")
library(ggplot2);library(sva); library(WGCNA);library(dplyr);library(stringr);
library(ggalluvial);library(Biobase); library(ggpubr)
library(Rmagic); library(Rtsne); library(colorspace);library(edgeR);
library(dendextend);library(SC3);library(pheatmap)

### set working directory

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")

# options
enableWGCNAThreads()
allowWGCNAThreads()
condition = TRUE




### load normalized data from each disease

load("./codes/st2_lm_multiExpr.Rdata")

### compute pca using all genes

CPM = cpm(allExp)
pcAll = prcomp(allExp, center = T)

# tsne 
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



##QC Pre-Combat
pdf("./results/figures/PreCorrection_AllCovaritates.pdf", height = 11, width =12.5 )
par(mfrow = c(4,2)
## Covariates
#---Subjects
plot(all_datMeta$Dx, col=c(1:9), main= "Subjects")

#---Age
plot(all_datMeta$Age ~ all_datMeta$Dx, col= c(1:9), 
     main="Age", ylab="year", xlab="")

#---Sex
sex_col = if_else(all_datMeta$Sex=="female", "yellow", "blue")
plot(all_datMeta$Sex ~ all_datMeta$Dx, col=as.numeric(all_datMeta$Sex), 
     main="Sex", ylab="", xlab="")



#---PMI
plot(as.numeric(all_datMeta$PMI) ~ all_datMeta$Dx, col=c(1:9),
     main="PMI", ylab="", xlab="", ylim = c(0,700))

all_datMeta$Study = as.factor(all_datMeta$Study)

age_col = numbers2colors(all_datMeta$Age, blueWhiteRed(100), 
                         signed=F, centered=T, lim=c(min(all_datMeta$Age, na.rm=T),
                                                     max(all_datMeta$Age, na.rm=T)))
pmi_col = numbers2colors(all_datMeta$PMI, blueWhiteRed(100), 
                         signed=F, centered=T, lim=c(min(all_datMeta$PMI, na.rm=T),
                                                     max(all_datMeta$PMI, na.rm=T)))

plot(density(all_datExpr[,1]), xlim=c(-5,15), ylim=c(0, 0.7), col = as.numeric(all_datMeta$Study[1]), 
     xlab="Intensity (log2)", ylab="Density", main="Mega-Analysis: Pre-batch-effect correction")
for(i in 2:dim(all_datExpr)[[2]])
  lines(density(all_datExpr[,i]), xlim=c(0,20), col = as.numeric(all_datMeta$Study[i]))  
legend("topleft", (levels(all_datMeta$Study)), col=c(1:8), pch=16, cex=0.5)


#mds
mds = cmdscale(dist(t(all_datExpr)), eig = T)   
pc1 = mds$eig[1]^2 / sum(mds$eig^2)   
pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(all_datMeta$Study)), 
     pch=20, main="Multidimensional Scaling Plot\nPre-batch-effect correction", 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), 
     ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
legend("bottomleft", (levels(all_datMeta$Study)), col=c(1:8), pch=16, cex=0.5)

dev.off()


#Normalize by Study
rownames(all_datMeta) = colnames(all_datExpr)
all_datExpr = as.matrix(all_datExpr)
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

# pdf("./results/figures/hcluster.pdf", width = 15, height = 8)
# par(mar = c(2,2,2,10))
# A2Rplot(h.usa, 
#         fact.sup = all_datMeta$Study,
#         box = FALSE, 
#         show.labels = F,
#         col.up = "black", 
#         main="", 
#         k=length(levels(factor(all_datMeta$Study))) ) # k = changes the detail/number of subgroups shown.
# 
# ##### Generate coloured dendrogram (indicate groups)
# A2Rplot(h.usa, 
#         fact.sup = datMeta$Dx, 
#         box = F, 
#         show.labels = F, 
#         # col.up = "black",
#         # col.down = rainbow(9),
#         main=" ",knot.pos = "bary",
#         k = length(levels(factor(datMeta$Dx)))) # k = changes the detail/number of subgroups shown.
# 
# dev.off()        
# 
# ggdendrogram(h.usa,labels = T,leaf_labels = F)
######

###
dend <- as.dendrogram(h.usa)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=9)

# add labels
labels(dend) <- paste(as.character(all_datMeta$Dx)[order.dendrogram(dend2)])

cols <- c("lightgrey","aquamarine","red","blue","chartreuse",
          "brown","blueviolet","yellow","orange")

all_datMeta$color = as.numeric(all_datMeta$Dx)
all_datMeta$color = cols[all_datMeta$color]
all_datMeta$lobe_color = as.numeric(all_datMeta$Brain_lobe)
all_datMeta$lobe_color = cols[all_datMeta$lobe_color]

labels_colors(dend) <- all_datMeta$color[order.dendrogram(dend)]

# We shall add the dx to the labels:


# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend, hang_height=0.1)

# reduce the size of the labels:

dend <- set(dend, "labels_cex", 0.5) %>% 
  set("leaves_pch", 19) %>%
  set("leaves_col",all_datMeta$lobe_color[order.dendrogram(dend)]) %>% 
  set("leaves_cex",0.5)



url <- paste0("https://raw.githubusercontent.com/TreeViz/",
                         "metastyle/master/design/viz_targets_exercise/")
              
 info <- read.csv(paste0(url, "tip_data.csv"))
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
pdf("./results/figures/samples_dendrogram2.pdf", width = 10,height = 10)

circlize_dendrogram(dend,labels = T)

legend(x = 0.5, y= 0.5, legend = levels(all_datMeta$Dx),
       fill =cols, cex=1,border = F,bty = "n", title = "Disease")

legend(x = -0.5, y = 0.5, legend = levels(all_datMeta$Brain_lobe),
        cex=1, border = F,pch = 19, title = "Region",
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
#pre batch correction

pca = prcomp(t(all_datExpr), scale. = T,retx = T)
values = as.data.frame(pca$x)

#post batch correction

post.pca = prcomp(t(datExpr), scale. = T, retx = T)
post.values = as.data.frame(post.pca$x)


#3D pca plot 
library(scatterplot3d)

#colors for plot
col = vector()
colors = brewer.pal(9,"Set1")
for (i in 1:length(as.numeric(all_datMeta$Dx))){
  for (j in 1:9){
    if( j == as.numeric(all_datMeta$Dx[i])){
      cl = colors[j]
    }
  }
  col = c(col,cl)
}

#3d plot
pdf("./results/figures/PCA.pdf",
    width = 15, height = 5)
par(mfrow=c(1:2)) 
with(values, {
  s3d <- scatterplot3d(
    x = PC1, 
    y = PC2, 
    z = PC3,
    color = col, 
    pch = 19, 
    type = "h",
    lty.hplot = 3,
    angle = -45,
    scale.y = .75,
    main = "3-D PCA\npre-batch correction",
    xlab = "PC1",
    ylab = "PC2",
    zlab = "PC3")
  
  #add legend
  legend(#location
    "topright",
    inset=.05,
    # suppress legend box, shrink text 50%
    bty="n",
    cex=.7,
    title="Diagnosis",
    levels(all_datMeta$Dx),
    fill=brewer.pal(9,"Set1"))
})

#postcorrection plot
with(post.values, {
  s3d <- scatterplot3d(
    x = PC1, 
    y = PC3, 
    z = PC2,
    color = col, 
    pch = 19, 
    type = "h",
    lty.hplot = 3,
    angle = 45,
    scale.y = .5,
    main = "3-D PCA\npost-batch correction",
    xlab = "PC1",
    ylab = "PC3",
    zlab = "PC2")
  
  #add legend
  legend(#location
    "topright",
    inset=.05,
    # suppress legend box, shrink text 50%
    bty="n",
    cex=.7,
    title="Diagnosis",
    levels(all_datMeta$Dx),
    fill=brewer.pal(9,"Set1"))
})
dev.off()


# using tSne
set.seed(123) # for reproducibility

#pre correction
pretsn = Rtsne(t(values),
               pca = F
               # verbose=TRUE, 
               # max_iter = 10000,
               # dims = 3,
               # partial_pca = T,
               # eta = 2000
)
predata = as.data.frame(pretsn$Y)
predata$Dx = as.character(all_datMeta$Dx)
predata$Study = as.character(all_datMeta$Study)
predata$region = as.character(all_datMeta$Brain_lobe)

#post correction
pc.fit = prcomp(t(datExpr),scale. = T, center = F)
set.seed(123)
tsne <- Rtsne(pc.fit$x[,1:10], 
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

# preDx = ggplot(predata, aes(x =predata[,1] ,
#                             y =predata[,2]))+
#   geom_point(aes(color = Dx), size = 1)+
#   scale_color_brewer(type = "div",palette = "Set1")+
#   labs(x = "tSNE1", y = "tSNE2",
#        title = "Pre-correction disease")+
#   theme_bw()
# 
# preStudy = ggplot(predata, aes(x =predata[,1] ,
#                                y =predata[,2]))+
#   geom_point(aes(color = Study), size = 1)+
#   # scale_color_brewer(type = "div",palette = "Set1")+
#   labs(x = "tSNE1", y = "tSNE2",
#     title = "Pre-correction study")+
#   theme_bw()
# 
# preregion = ggplot(predata, aes(x =predata[,1] ,
#                                 y =predata[,2]))+
#   geom_point(aes(color = region), size = 1)+
#   scale_color_brewer(type = "div",palette = "Set1")+
#   labs(x = "tSNE1", y = "tSNE2",
#     title = "Pre-correction region")+
#   theme_bw()
# 

# groupData = data %>% dplyr::group_by(Dx) %>% select(V1,V2) %>% summarize_all(mean)
colnames(data)= c("tSNE1", "tSNE2", "Disease","Study","Region")
postDx = ggplot(data, aes(x =tSNE1 ,
                    y =tSNE2,
                    color = Disease))+
   geom_point( alpha =0.8, size = 1)+
  # scale_color_brewer(type = "div",palette = "Set1")+
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

ggsave(filename = "./results/figures/PCA_Plots_combined/mega_analysis.pdf",
       plot = pcPlots, width = 20, height = 12)



#clustering using phateR package
phate.test = phate(t(datExpr))

# runs phate with different parameters
tree.phate <- phate(t(datExpr), gamma=1 ,
                    t= 10,knn = 5,
                    init=phate.test)
#PLOT
ggplot(tree.phate, aes(x=PHATE1, y=PHATE2, color=datMeta$Dx)) +
  geom_point()+
  scale_color_manual(values = cols)+
  theme_clean()

#projection score

#computing variance for each gene
datVar= data.frame(matrix(NA, nrow = nrow(datExpr),ncol = 2))
colnames(datVar) = c("gene", "var")

var =apply(datExpr,1,var)
var =as.data.frame(unlist(var)) 

datVar$gene = rownames(var)
datVar$var = var$`unlist(var)` 

write.table(datVar,file = "./results/tables/GeneVariance.tsv", col.names = T)
write.csv(datExpr,  file = "./results/tables/datExpr.csv")

#load top variable genes
geneSet = read.table("./results/tables/projectionScore.txt")

#filter data for selected genes
selected = as.table(datExpr[geneSet$V1,])
proj.tsne = Rtsne(t(selected), dims = 3, num_threads = 4)
dims = as.data.frame(proj.tsne$Y)
rownames(dims) = colnames(selected)
#plot 
ggplot(dims, aes(x = V3, y = V1, color = all_datMeta$Dx))+
  geom_point(size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")

# run pca usign top genes
top_genes = read.csv("./results/tables/top_genes.csv",col.names = F)
topExpr = datExpr[top_genes$FALSE.,]

set.seed(123) # for reproducibility
topGenes_tsne <- Rtsne(t(topExpr),dims = 2)
top_data = as.data.frame(topGenes_tsne$Y)

ggplot(top_data, aes(x =top_data[,1] ,
                 y =top_data[,2] ))+
  geom_point(aes(color = datMeta$Dx), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")


##QC - PostCombat

pdf("./results/figures/PostCorrection_AllCovariats.pdf", height = 12, width = 13)
par(mfrow=c(3,3))
plot(density(datExpr[,1]),xlim=c(0,10), ylim=c(0, 0.5), 
     col = as.numeric(all_datMeta$Study[1]),  xlab="Intensity (log2)", ylab="Density", 
     main="Mega-Analysis: Post-batch-effect correction")
for(i in 2:dim(datExpr)[[2]])
  lines(density(datExpr[,i]), xlim=c(0,10), ylim=c(0,0.5), col = as.numeric(all_datMeta$Study[i]))  
legend("topright", levels(all_datMeta$Study), col=c(1:8), pch=16,cex=0.7)

#MDS Plot
mds = cmdscale(dist(t(datExpr)), eig = T)
pc1 = mds$eig[1]^2 / sum(mds$eig^2)
pc2 = mds$eig[2]^2 / sum(mds$eig^2)
#study
plot(mds$points, col=as.numeric(as.factor(all_datMeta$Study)), 
     pch=20, main="MDS: Study", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), 
     ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
legend("bottomleft", 
       levels(all_datMeta$Study), 
       col=c(1:length(levels(all_datMeta$Study))), 
       pch=16, 
       cex=0.65)

#Dx
plot(mds$points, col=as.numeric(as.factor(all_datMeta$Dx)), 
     pch=16, main="MDS: Diagnosis",  xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), 
     ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
legend("topright", 
       levels(all_datMeta$Dx), 
       col=c(1:length(levels(all_datMeta$Dx))), 
       pch=16, 
       cex=0.6)

#Sex
plot(mds$points, col=sex_col, pch=16, main="MDS - Sex", 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""),
     ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
legend("topright", 
       levels(all_datMeta$Sex), 
       col=sex_col, 
       pch=16, 
       cex=0.8)

#Age
plot(mds$points, col=age_col, pch=16, main="MDS - Age",
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""),
     ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 

#PMI
plot(mds$points, col=pmi_col, pch=16, main="MDS - PMI", 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), 
     ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 



#brain region
plot(mds$points, col=as.numeric(all_datMeta$Brain_Region), main = "MDS-Brain region",
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), 
     ylab = paste("PC2 (", signif(100*pc2,3), "%)", sep = ""))
legend("topright", levels(all_datMeta$Brain_Region), col = c(1:length(as.factor(all_datMeta$Brain_Region))),
       pch = 16, cex = 0.4)

#brain lobe
plot(mds$points, col=as.numeric(all_datMeta$Brain_lobe), main = "MDS-Brain lobes",
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), 
     ylab = paste("PC2 (", signif(100*pc2,3), "%)", sep = ""))
legend("topright", levels(all_datMeta$Brain_lobe), col = c(1:length(as.factor(all_datMeta$Brain_lobe))),
       pch = 16, cex = 0.5)


dev.off()

#---dendogram
pdf("./results/figures/dendogram.pdf", width = 12,height = 5)
# tree = hclust(dist(t(datExpr)), method = "average")
par(mfrow=c(1,1))
plotDendroAndColors(h.usa,
                    cbind(as.numeric(datMeta$Dx),
                          as.numeric(datMeta$Brain_lobe)),
                    groupLabels = c("Disease","Region"),
                    cex.colorLabels=0.8, addGuide = F,
                    cex.dendroLabels=0.15,
                    main="Dendogram of samples")

dev.off()

# a flowchart plot 
datMeta = all_datMeta
library(ggalluvial)
# filtering for disease, lobe and region
df = datMeta %>%
  group_by(Dx, Brain_lobe, Sex)%>%
  count() %>% arrange(desc(n))

df$Dx = gsub("Normal", "Control", df$Dx)
df$Dx = factor(df$Dx, levels = unique(df$Dx))
df$Brain_lobe = factor(df$Brain_lobe, levels = unique(df$Brain_lobe))
df$Sex = factor(df$Sex, levels = unique(df$Sex))
colnames(df) = gsub("Dx", "Disease", colnames(df))
g <- ggplot(df,
       aes(axis1 = Disease,
           axis2 = Brain_lobe,
           axis3 = Sex,
           y = n)) +
  stat_alluvium(aes(fill = Disease), lode.guidance = "forward")+
  stat_stratum()+
  geom_label_repel(stat = "stratum",
            infer.label = TRUE, 
            size =8) +
  scale_x_discrete(limits = c("Disease", "Brain_lobe","Sex"),
                   expand = c(0.00005,0.00005)
                   ) +
  scale_fill_brewer(type = "div",palette = "Set1",aesthetics = "fill")+
  labs(title = "The samples used in this study",
       y = "Number of samples") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(size = 14,colour = "black"),
        text = element_text(size =14, colour = "black"))
g
ggsave(filename ="./results/figures/OverviewData.pdf",
       plot = g,
       device = "pdf",units = "cm",
       width = 30, height = 38, dpi = 300)

#usign all covariates
df2 <- datMeta %>%
  dplyr::group_by(Dx, Study, Sex, Brain_Region, Brain_lobe)%>%
  count()

p <- ggplot(df2,
       aes(axis1 = Dx,
           axis2 = Study,
           axis3 = Sex,
           axis4 = Brain_Region,
           axis5 = Brain_lobe,
           y = n)) +
  geom_alluvium(aes(fill = Dx)) +
  geom_stratum() +
  geom_text(stat = "stratum", 
            label.strata = TRUE, size =3) +
  scale_x_discrete(limits = c("Dx", "Study","Sex", "Brain_Region", "Brain_lobe"),
                   expand = c(0.00005,0.00005)) +
  scale_fill_brewer(type = "div",palette = "Set1",aesthetics = "fill")+
  labs(title = "The overview of the samples types in the study",
       y = "Number") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        text = element_text(size =12))

ggsave(filename ="./results/figures/OverView_allCov.pdf",
       plot = p,
       device = "pdf",units = "cm",
       width = 35, height = 40,dpi = 300)

#correlation between disease and regions
#marge datexpr and datmet
dxCors = as.data.frame(matrix(NA, 
                                  ncol = 8, 
                                  nrow = length(unique(datMeta$Brain_lobe))))
rownames(dxCors) = unique(datMeta$Brain_lobe)
colnames(dxCors) = c("AD","PD","PA","PSP","Scz","ASD","BP","MDD")

#p values
dxPval = dxCors
#calculate linear regression between modules and covariates
for (dx in colnames(dxCors)){
  for (l in rownames(dxCors)){
    mcor =lm(ME.datMeta[,mod] ~ ME.datMeta[,var], data = ME.datMeta )
    moduleCors[mod,var] = summary(mcor)$adj.r.squared
    modulesPval[mod, var] = anova(mcor)$`Pr(>F)`[1]
  }
}
#export 
sampleTable = datMeta %>% group_by(Study, Dx, Brain_lobe, Brain_Region) %>% count()
write.csv("./results/tables/SamplesTable.csv", x = sampleTable)

#Dx numbers
DxSample = datMeta %>% group_by(Dx) %>% count() %>% arrange(desc(n))
StudyDxSubject = datMeta %>% group_by(Study, Dx) %>% count()
study =  datMeta %>% group_by(Study) %>% count()
tablePlot = ggtexttable(DxSample)
regionPlot = 
ggsave(filename = "./results/figures/DxSamples.pdf",plot = tablePlot, width = 3, height = 4)

save.image("./codes/Combined.Rdata")
save(file ="./codes/st4_combinedData.Rdata", datMeta, datExpr, all_datExpr)






