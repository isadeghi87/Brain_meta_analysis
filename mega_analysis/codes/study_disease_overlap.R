#### study_disease_overlap.R
#### check correlation between studies and disease results

library(GeneOverlap); library(stringr)
library(tidyr)
library(dplyr)
library(ggpubr)

### set working directory

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/")



files = list.dirs(path = "./studies/",full.names = T, recursive = T)
files = files[grep("tables",files)]

overlaps = vector("list", 5)
names(overlaps) = c("AD","Scz","ASD","BP","MDD")
value_list = overlaps

for (dx in names(overlaps)){
  stu.files = list.files(files,pattern = dx ,full.names = T)
  
  n = length(stu.files)
  studies = vector("list",n)
  studies_names = str_split(basename(stu.files),
                                pattern = paste("_",dx,sep = ""),
                            simplify = T)[,1] 
  # studies
  for (i in 1:n){
    dat = read.csv(stu.files[i], row.names = 1)
    studies[[i]] = dat
  }
  
  # disease specific stats
  sumstat = read.csv(paste("./disease_specific/",dx,
                           "_metaAnalysis/tables/",dx,"_sumstats.csv",
                           sep = ""),row.names = 1)
  
  #filter for significant genes ###
  sumstat_sig = subset(sumstat, p.value <0.05)
  id1 = rownames(sumstat_sig)
  
  
  ## check overlap ##
  overlap = as.data.frame(matrix(NA, 
                                     ncol = n, 
                                     nrow = 4))
  rownames(overlap) = c("p.value","Odd ratio","Jaccard index", "Simpson index")
  colnames(overlap) = studies_names
  for (i in 1:n){
    dat = studies[[i]]
    dat = subset(dat,p.value<0.05)
    id2= rownames(dat)
    go.obj = newGeneOverlap(listA = id1,
                            listB = id2,
                            genome.size = nrow(sumstat))
    go.obj = testGeneOverlap(go.obj)
    
    # extract pvalue
    overlap[1,i] = signif(getPval(go.obj),2)
    
    # extract odd ratio
    overlap[2,i] = signif(getOddsRatio(go.obj),2)
    
    # get jaccard ratio
    overlap[3,i] = signif(getJaccard(go.obj),2)
    
    # simpson index 
    simp = signif(length(intersect(id1,id2))/min(length(id1),length(id2)),2)
    overlap[4,i] = simp
    overlaps[[dx]] = overlap
  
  }
}

for (dx in names(overlaps)) {
  write.csv(overlaps[[dx]],
              file = paste("./mega_analysis/results/tables/",dx,"_studies_overlap.csv",
                           sep = ""))
  
}

## reform data for plotting

data.overlap = data.frame()
table.plots = vector("list",5)
names(table.plots) = names(overlaps)

for (dx in names(overlaps)) {
  dat = as.data.frame(overlaps[[dx]])
  dat["Disease",] = dx
  dat["study",] = colnames(dat)
  dat = t(dat)
  data.overlap = rbind(data.overlap,dat)
  p = ggtexttable(dat,rows = NULL,
                  theme = ttheme(base_style = "mOrangeWhite"))
  table.plots[[dx]] = p
}

rownames(data.overlap) = NULL
for (i in 1:4) {
  data.overlap[,i] = as.numeric(data.overlap[,i])
}

 p1 = ggarrange(plotlist = table.plots, ncol = 1,
                label.y = 1,labels = names(table.plots))
ggsave(p1, filename = "./mega_analysis/results/figures/overlaps/table_overlaps.pdf",
       height = 12,width = 6)
 
## plot

library(ggplot2)
library(ggthemes)
cols <- c("aquamarine","red","blue","chartreuse",
          "brown","blueviolet","yellow","orange")

p_j = ggplot(data.overlap)+
             geom_point(aes(fill = Disease,x = data.overlap$`Jaccard index`,
                            y = reorder(study, -data.overlap$`Jaccard index`), 
                            alpha = p.value <0.05), 
                        # to use a border color pch needs to be > 20
                        colour = "black",
                        size = 5,
                        pch = 21)+
  facet_wrap(~Disease,nrow = 1)+
  scale_fill_manual(values = cols)+
  scale_x_continuous(limits = c(0,1))+
  labs(x = "Jaccard index",
       y = "Study")+
  theme_clean(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45))

p_s = ggplot(data.overlap)+
  geom_point(aes(fill = Disease,x = data.overlap$`Simpson index`,
                 y = reorder(study, -data.overlap$`Jaccard index`), 
                 alpha = p.value <0.05), 
             # to use a border color pch needs to be > 20
             colour = "black",
             size = 5,
             pch = 22)+
  facet_wrap(~Disease,nrow = 1)+
  scale_fill_manual(values = cols)+
  scale_x_continuous(limits = c(0,1))+
  labs(x = "Simpson's index",
       y = "Study")+
  theme_clean(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45))

p = ggarrange(p_j,p_s,ncol = 1,common.legend = T)

ggsave("./mega_analysis/results/figures/overlaps/Disease_study_overlap.pdf",
       width = 13, height = 10)


