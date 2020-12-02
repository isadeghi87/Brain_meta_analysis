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
ctlMeta = allMeta[id,]
dxMeta = allMeta[-id,]

## expression data 
ctlExp = datExp[,id]
dxExp = datExp[,-id]

# corrected exp
ctlExp.comb = datExp.comb[,id]
dxExp.comb = datExp.comb[,-id]


#### PCA for controls before batch #### 
ctl.cpm  = cpm(ctlExp, log = T)
pc.ctl = prcomp(t(ctl.cpm))
pc.ctl.res  = pc.ctl$x[,1:20]

#### tsne ####
tsn.ctl = Rtsne(pc.ctl.res,theta = 0,check_duplicates = F, 
                max_iter = 1000)

tsn.ctl.res = as.data.frame(tsn.ctl$Y)
rownames(tsn.ctl.res) = rownames(pc.ctl.res)
colnames(tsn.ctl.res) = c("tSNE1","tSNE2")
tsn.ctl.res$Dx = ctlMeta$Dx
tsn.ctl.res$study = ctlMeta$Study
tsn.ctl.res$Brain_Lobe = ctlMeta$Brain_Lobe

#### PCA for controls after batch correction#### 
pc.ctl.comb = prcomp(t(ctlExp.comb))
pc.ctl.comb.res  = pc.ctl.comb$x[,1:20]

#### tsne ####
tsn.ctl.comb = Rtsne(pc.ctl.comb.res,theta = 0,check_duplicates = F, 
                max_iter = 1000)

tsn.ctl.comb.res = as.data.frame(tsn.ctl.comb$Y)
rownames(tsn.ctl.comb.res) = rownames(pc.ctl.comb.res)
colnames(tsn.ctl.comb.res) = c("tSNE1","tSNE2")
tsn.ctl.comb.res$Dx = ctlMeta$Dx
tsn.ctl.comb.res$study = ctlMeta$Study
tsn.ctl.comb.res$Brain_Lobe = ctlMeta$Brain_Lobe




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
p1  = mygg(data =tsn.ctl.res,
           color = tsn.ctl.res$study,
           title = "Study - pre batch correction")
p2  = mygg(data =tsn.ctl.res,color = tsn.ctl.res$Brain_Lobe,
           title = "Brain lobe- pre batch correcyion")+
  scale_color_manual(values = cols)

p3  = mygg(data =tsn.ctl.comb.res,
           color = tsn.ctl.comb.res$study,
           title = "Study - batch corrected")
p4 = mygg(data =tsn.ctl.comb.res,
           color = tsn.ctl.comb.res$Brain_Lobe,
           title = "Brain lobe - batch corrected")+
  scale_color_manual(values = cols)

ctl.plots = ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2)

ggsave("./results/figures/pca/pca_CTL.pdf",
       width = 12, height = 8)
