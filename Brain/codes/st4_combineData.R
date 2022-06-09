##st4_combineDatasets.R
## Here we combine studies together and performs batch effect correction with ComBat to normalize them
#
rm(list=ls()); options(stringsAsFactors=F)

#source("http://bioconductor.org/biocLite.R")
library(ggplot2);library(sva); library(WGCNA);library(dplyr);library(stringr);
library(ggalluvial);library(Biobase); library(ggpubr)
library(Rmagic); library(Rtsne)
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")

#load normalized data from each disease
files = dir("./codes/", pattern = "_metadata")
n=length(files)
multiExpr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/",files[[i]],sep = ""))
  multiExpr[[i]]$datExpr = datExpr
  multiExpr[[i]]$datMeta = datMeta
  rm(datExpr,datMeta)
}

genes = rownames(multiExpr[[1]]$datExpr)
for(i in 2:length(multiExpr)) genes = intersect(genes, rownames(multiExpr[[i]]$datExpr))

all_datExpr = data.frame(row.names = genes)
all_datMeta = data.frame(matrix(NA, nrow=0, ncol=8));

for(i in 1:length(multiExpr)) {
  all_datExpr = cbind(all_datExpr, multiExpr[[i]]$datExpr[match(genes, rownames(multiExpr[[i]]$datExpr)),])
  datMetaNA = as.data.frame(matrix(NA, ncol=8, nrow=nrow(multiExpr[[i]]$datMeta)))
  idx = pmatch(c("Study", "Subject_ID", "Dx", "Brain_Region", "Age", "Sex", "PMI", "Brain_lobe"),colnames(multiExpr[[i]]$datMeta))
  datMetaNA[,which(!is.na(idx))] = multiExpr[[i]]$datMeta[,na.omit(idx)]
  all_datMeta=rbind(all_datMeta, datMetaNA)
}

colnames(all_datMeta) = c("Study", "Subject_ID", "Dx", "Brain_Region", "Age", "Sex", "PMI", "Brain_lobe")
all_datMeta$PMI[is.na(all_datMeta$PMI)] = mean(all_datMeta$PMI, na.rm = T)
all_datMeta$Study <- as.factor(all_datMeta$Study)
all_datMeta$Sex = as.character(all_datMeta$Sex)
all_datMeta$Sex[all_datMeta$Sex == "Male"] = "male"
all_datMeta$Sex[all_datMeta$Sex == "Female"] = "female"
all_datMeta$Sex[all_datMeta$Sex == "unknown"] = "male"
all_datMeta$Sex <- as.factor(all_datMeta$Sex)
all_datMeta$Brain_Region = as.factor(all_datMeta$Brain_Region)
all_datMeta$Brain_lobe[all_datMeta$Brain_lobe == "DLPFC"] ="Frontal"
all_datMeta$Brain_lobe <- as.factor(all_datMeta$Brain_lobe)

#--remove duplicate samples
idx <- !duplicated(colnames(all_datExpr)) 
all_datExpr <- all_datExpr[,idx]#---10313 genes
all_datMeta <- all_datMeta[idx,] #4627 Samples, 3276 Subjects

 ##QC Pre-Combat
pdf("./results/figures/PreCorrection_AllCovaritates.pdf", height = 11, width =12.5 )
par(mfrow = c(4,2))
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
mod = model.matrix(~Dx+Sex, data=all_datMeta)
modcombat = model.matrix(~Dx + Sex,data = all_datMeta)
datExpr = ComBat(all_datExpr, batch=all_datMeta$Study,mod=modcombat, prior.plots = F)

#####  Compute distance matrix and hierarchical clustering
d.usa <- dist(t(datExpr), method = "euc")
h.usa <- hclust(d.usa, method = "ward.D2")


##### Generate coloured dentrogram (indicate datasets)
h.usa$labels = colnames(datExpr)
pdf("./results/figures/Study_effects_cluster_datasets_ComBat.pdf", width = 30, height = 15)
par(mar = c(2,2,2,10))
A2Rplot(h.usa, 
        fact.sup = all_datMeta$Study,
        box = FALSE, 
        show.labels = F,
        col.up = "black", 
        main="", 
        k=length(levels(factor(all_datMeta$Study))) ) # k = changes the detail/number of subgroups shown.

##### Generate coloured dentrogram (indicate groups)
A2Rplot(h.usa, 
        fact.sup = all_datMeta$Dx, 
        box = F, 
        show.labels = F, 
        col.up = "black", main=" ",
        k = length(levels(factor(all_datMeta$Dx)))) # k = changes the detail/number of subgroups shown.

dev.off()        

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
tsne <- Rtsne(t(postvalues), 
              # perplexity=50, 
              #  verbose=TRUE, 
              #  max_iter = 5000,
              #  dims = 3,
              pca = F
              # partial_pca = T,
              # eta = 100
              )
data = as.data.frame(tsne$Y)
data$Dx = as.character(all_datMeta$Dx)
data$Study =as.character(all_datMeta$Study) 
data$region = as.character(all_datMeta$Brain_lobe)


preDx = ggplot(predata, aes(x =predata[,1] ,
                 y =predata[,2]))+
  geom_point(aes(color = Dx), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")+
  labs(x = "tSNE1", y = "tSNE2",
    title = "Pre-correction disease")+
  theme_bw()

preStudy = ggplot(predata, aes(x =predata[,1] ,
                               y =predata[,2]))+
  geom_point(aes(color = Study), size = 1)+
  # scale_color_brewer(type = "div",palette = "Set1")+
  labs(x = "tSNE1", y = "tSNE2",
    title = "Pre-correction study")+
  theme_bw()

preregion = ggplot(predata, aes(x =predata[,1] ,
                                y =predata[,2]))+
  geom_point(aes(color = region), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")+
  labs(x = "tSNE1", y = "tSNE2",
    title = "Pre-correction region")+
  theme_bw()

# groupData = data %>% dplyr::group_by(Dx) %>% select(V1,V2) %>% summarize_all(mean)
postDx = ggplot(data, aes(x =data[,1] ,
                    y =data[,2], color = Dx ))+
   geom_point( alpha =0.8, size = 1)+
   scale_color_brewer(type = "div",palette = "Set1")+
  labs(x = "tSNE1", y = "tSNE2",
    title = "Post-correction disease")+
  theme_bw()


postStudy = ggplot(data, aes(x =data[,1] ,
                                y =data[,2]))+
  geom_point(aes(color = Study), size = 1)+
  # scale_color_brewer(type = "div",palette = "Set1")+
  labs(x = "tSNE1", y = "tSNE2",
    title = "Post-correction study")+
  theme_bw()

postregion = ggplot(data, aes(x =data[,1] ,
                                 y =data[,2]))+
  geom_point(aes(color = region), size = 1)+
  scale_color_brewer(type = "div",palette = "Set1")+
  labs(x = "tSNE1", y = "tSNE2",
       title = "Post-correction region")+
  theme_bw()

library(ggpubr)
pcPlots = ggarrange(preDx, preregion, preStudy,
                    postDx, postregion, postStudy, 
          ncol = 3, nrow = 2)

ggsave(filename = "./results/figures/PCA_Plots_combined/mega_analysis.pdf",
       plot = pcPlots, width = 20, height = 12)

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

set.seed(1) # for reproducibility
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

# a flowchart plot 
datMeta = all_datMeta
library(ggalluvial)
# filtering for disease, lobe and region
df = datMeta %>%
  group_by(Dx, Brain_lobe, Brain_Region)%>%
  count() %>% arrange(desc(n))

df$Dx = gsub("Normal", "Control", df$Dx)
df$Dx = factor(df$Dx, levels = unique(df$Dx))
df$Brain_lobe = factor(df$Brain_lobe, levels = unique(df$Brain_lobe))
df$Brain_Region = factor(df$Brain_Region, levels = unique(df$Brain_Region))
colnames(df) = gsub("Dx", "Disease", colnames(df))
g <- ggplot(df,
       aes(axis1 = Disease,
           axis2 = Brain_lobe,
           axis3 = Brain_Region,
           y = n)) +
  stat_alluvium(aes(fill = Disease), lode.guidance = "forward")+
  stat_stratum()+
  geom_label(stat = "stratum",
            label.strata = TRUE, size =4) +
  scale_x_discrete(limits = c("Disease", "Brain_lobe","Brain_Region"),
                   expand = c(0.00005,0.00005)) +
  scale_fill_brewer(type = "div",palette = "Set1",aesthetics = "fill")+
  labs(title = "The source of the samples in the study",
       y = "Number") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        text = element_text(size =14, colour = "black"))
g
ggsave(filename ="./results/figures/OverviewData.pdf",
       plot = g,
       device = "pdf",units = "cm",
       width = 30, height = 60, dpi = 300)

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

save(file ="./codes/st4_combined/mega_analysisData.Rdata", datMeta, datExpr, all_datExpr, all_datExpr)






