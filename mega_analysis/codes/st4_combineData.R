##st4_combineDatasets.R
## Here we combine studies together and performs batch effect correction with ComBat to normalize them
#
rm(list=ls()); options(stringsAsFactors=F)

library(ggplot2);library(sva); 
library(WGCNA);library(dplyr);
library(stringr);library(RColorBrewer)
library(ggpubr);library(edgeR)
library(Rmagic); library(Rtsne);
library(ComplexHeatmap)
library(colorspace);

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
names(cols) = c("CTL","AD","PD","PA","PSP","Scz","ASD","BP","MDD")
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

lcpm = cpm(datExp, log = T)

# pca 
set.seed(123)
pc = prcomp(t(lcpm))

pc.data = pc$x[,1:10]
dup = duplicated(pc.data)
pc.data = pc.data[!dup,]

# tsne
set.seed(123)
tsn = Rtsne(pc.data, pca = F, dims = 3,theta = 0)
data.tsne = as.data.frame(tsn$Y)
rownames(data.tsne) = rownames(pc.data)
colnames(data.tsne) = c("tSNE1","tSNE2","tSNE3")
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
ggsave(plot = plots, "./results/figures/pca/pca_allGenes.pdf",
       width = 7, height = 10) 



#### correct for batch effect ####
design = model.matrix(~Dx + Brain_Lobe,data = allMeta)
v = voom(datExp,design = design)
exp = as.matrix(v$E)
datExp.comb = ComBat(exp, batch = allMeta$Study,mod = design)


## pca 
set.seed(123)
pc.comb = prcomp(t(datExp.comb))

## Get the standard deviation and variance per principal component

sd.per.pc <- pc.comb$sdev
var.per.pc <- sd.per.pc^2

## Display the percentage of total variance explained by each 

sd.per.pc.percent <- 100 * sd.per.pc/sum(sd.per.pc)
var.per.pc.percent <- 100 * var.per.pc/sum(var.per.pc)

pdf("./results/figures/pca/Variance_perPCA_combat.pdf")

barplot(var.per.pc.percent[1:20], main='Percent of variance  per component',
        xlab='Component', 
        ylab='Percent variance', 
        col='#BBDDFF')
dev.off()


pc.data = pc.comb$x[,1:20]
dup = duplicated(pc.data)
pc.data = pc.data[!dup,]

# tsne

set.seed(123)
# tsne 
thet = 0
perp = c(20,50,100)
int = c(500,1000,2000)

  for ( p in perp){
    for (i in int){
      tsn = Rtsne(pc.data, pca = F,
                      theta = t, perplexity = p, max_iter = i,
                      check_duplicates = F)
      dat= as.data.frame(tsn$Y)
      rownames(dat) = rownames(pc.data)
      colnames(dat) = c("tSNE1","tSNE2")
      dat$Dx = allMeta$Dx[match(rownames(dat),rownames(allMeta))]
      dat$study = allMeta$Study[match(rownames(dat),rownames(allMeta))]
      dat$Brain_Lobe = allMeta$Brain_Lobe[match(rownames(dat),rownames(allMeta))]
      
      
      
      tsDx = mygg(dat,color = dat$Dx, title = paste(t,p,i, sep = "-"))+
        scale_color_manual(values = cols)
      tsStudy = mygg(dat, color = dat$study, paste(t,p,i, sep = "-"))
      tsReg = mygg(dat, dat$Brain_Lobe, paste(t,p,i, sep = "-"))+
        scale_color_manual(values = cols)
      
      plots = ggarrange(tsDx,tsStudy,tsReg, ncol = 1, nrow = 3)
      
      # print(plots)
      ggsave(plot= plots,
             filename = paste("./results/figures/pca/pca_combat",t,p,i,".pdf",sep = ""),
             width = 6, height = 12)
    }
}

tsn = Rtsne(t(datExp.comb), pca = T ,theta = 0,check_duplicates = F)
data.tsne = as.data.frame(tsn$Y)
rownames(data.tsne) = rownames(pc.data)
colnames(data.tsne) = c("tSNE1","tSNE2")
data.tsne$Dx = allMeta$Dx[match(rownames(data.tsne),rownames(allMeta))]
data.tsne$study = allMeta$Study[match(rownames(data.tsne),rownames(allMeta))]
data.tsne$Brain_Lobe = allMeta$Brain_Lobe[match(rownames(data.tsne),rownames(allMeta))]

ggplot(data.tsne, aes(x = tSNE1, y = tSNE2, color  = Dx))+
  geom_point()+
  scale_color_manual(values = cols)
  
p1  = mygg(data =data.tsne,
           color = data.tsne$Dx, title = "disease")+
  scale_color_manual(values = cols)

p2  = mygg(data =data.tsne,color = data.tsne$study, title = "Study")
p3  = mygg(data =data.tsne,color = data.tsne$Brain_Lobe, title = "Brain lobe")+
  scale_color_manual(values = cols)


plots = ggarrange(p1,p2,p3, ncol = 1, nrow = 3)
ggsave(plot = plots, "./results/figures/pca/pca_filtered_combat.pdf",
       width = 7, height = 10) 


## heatmap ####
region_col = rainbow(7)
names(region_col) = levels(allMeta$Brain_Lobe)

col = columnAnnotation(group = factor(allMeta$Dx),
                 col =list(group = cols),
                 show_annotation_name=F)
bottom = columnAnnotation(Region = factor(allMeta$Brain_Lobe),
                          col = list(Region = cols),
                          show_annotation_name=F)

png("./results/figures/combined/Combined_exp_heatmap.png",res = 300)
Heatmap(as.matrix(datExp.comb),
        col = brewer.pal(11,"BrBG"),
        name = "log(CPM)",
        cluster_rows = T,
        use_raster = TRUE,raster_device = "png",
        clustering_distance_rows = "euclidean",
        cluster_columns = T,
        top_annotation = col,
        bottom_annotation = bottom,
        show_row_names = F,
        show_column_names =F,
        border = T,
        # height = unit(15,"cm") ,
        # width = unit(3,"cm"),
        show_row_dend = F,
        show_column_dend = F,
        show_heatmap_legend = T,
        heatmap_legend_param = gpar(fontsize=1)) 


dev.off()



save(datExp.comb,allMeta, file = "./codes/combined_batch_corrected.Rdata")
save(file = "./codes/Combined.Rdata",datExp,datExp.comb,allExp,allMeta)

