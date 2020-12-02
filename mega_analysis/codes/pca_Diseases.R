#### PCA for controls of combined data

library(ggplot2);library(sva); 
library(WGCNA);library(dplyr);
library(stringr);library(RColorBrewer)
library(ggpubr);library(edgeR)
library(Rmagic); library(Rtsne);
library(colorspace);library(sva);

### set working directory

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")


# load data
load("./codes/Combined.Rdata")


#### control samples 
id = which(allMeta$Dx == "CTL")
dxMeta = allMeta[-id,]

## expression data 
dxExp = datExp[,-id]

# corrected exp
dxExp.comb = datExp.comb[,-id]


#### PCA for controls before batch #### 
dx.cpm  = cpm(dxExp, log = T)
pc.dx = prcomp(t(dx.cpm))
pc.dx.res  = pc.dx$x[,1:20]

#### tsne ####
tsn.dx = Rtsne(pc.dx.res,theta = 0,check_duplicates = F, 
                max_iter = 1000)

tsn.dx.res = as.data.frame(tsn.dx$Y)
rownames(tsn.dx.res) = rownames(pc.dx.res)
colnames(tsn.dx.res) = c("tSNE1","tSNE2")
tsn.dx.res$Dx = dxMeta$Dx
tsn.dx.res$study = dxMeta$Study
tsn.dx.res$Brain_Lobe = dxMeta$Brain_Lobe

#### PCA for controls after batch correction#### 
pc.dx.comb = prcomp(t(dxExp.comb))
pc.dx.comb.res  = pc.dx.comb$x[,1:20]

#### tsne ####
tsn.dx.comb = Rtsne(pc.dx.comb.res,theta = 0,check_duplicates = F, 
                     max_iter = 1000)

tsn.dx.comb.res = as.data.frame(tsn.dx.comb$Y)
rownames(tsn.dx.comb.res) = rownames(pc.dx.comb.res)
colnames(tsn.dx.comb.res) = c("tSNE1","tSNE2")
tsn.dx.comb.res$Dx = dxMeta$Dx
tsn.dx.comb.res$study = dxMeta$Study
tsn.dx.comb.res$Brain_Lobe = dxMeta$Brain_Lobe




# colors
cols <- c("lightgrey","aquamarine","red","blue","chartreuse",
          "brown","blueviolet","yellow","orange")

# function of plot
mygg = function(data,color,title){
  plt = ggplot(data, aes(x =tSNE1 ,
                         y =tSNE2,
                         color = color))+
    geom_point( alpha =0.7, size = 1)+
    theme_classic2()+
    labs(title = title)
  return(plt)
}

# plot
p1  = mygg(data =tsn.dx.res,
           color = tsn.dx.res$study,
           title = "Study - pre batch correction")
p2  = mygg(data =tsn.dx.res,color = tsn.dx.res$Brain_Lobe,
           title = "Brain lobe- pre batch correction")+
  scale_color_manual(values = cols)
p3  = mygg(data =tsn.dx.res,color = tsn.dx.res$Dx,
           title = "Disease- pre batch correction")+
  scale_color_manual(values = cols)

p4  = mygg(data =tsn.dx.comb.res,
           color = tsn.dx.comb.res$study,
           title = "Study - batch corrected")
p5 = mygg(data =tsn.dx.comb.res,
          color = tsn.dx.comb.res$Brain_Lobe,
          title = "Brain lobe - batch corrected")+
  scale_color_manual(values = cols)
p6 = mygg(data =tsn.dx.comb.res,
          color = tsn.dx.comb.res$Dx,
          title = "Brain lobe - batch corrected")+
  scale_color_manual(values = cols)

dx.plots = ggarrange(p1,p2,p3,p4,p5,p6, ncol = 3, nrow = 2)

ggsave("./results/figures/pca/pca_Diseases.pdf",
       width = 16, height = 8)


#### mds plot 

plotMDS(t(dxExp.comb), col = as.numeric(as.factor(dxMeta$Dx)))
