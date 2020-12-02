#### study_disease_logFC_correlation.R
#### check correlation between studies and disease results


library(GeneOverlap)
library(tidyr)
library(dplyr)
### set working directory

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/")


### data for studies
files = list.dirs(path = "./studies/",full.names = T, recursive = T)
files = files[grep("tables",files)]

overlaps = vector("list", 5)
names(overlaps) = c("AD","Scz","ASD","BP","MDD")

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
  
  ## check overlap ##
  overlap = as.data.frame(matrix(NA, 
                                 ncol = 2, 
                                 nrow = n))
  colnames(overlap) = c("rho","p.value")
  rownames(overlap) = studies_names
  for (i in 1:n){
    dat = studies[[i]]
    idx= intersect(rownames(dat),rownames(sumstat))
    subdat = sumstat[idx,"logFC",drop = F]
    subdat2 = dat[idx,"logFC",drop = F]
    
    R = cor.test(subdat$logFC, subdat2$logFC, method = "spearman", use = "pairwise.complete.obs")
    rho = cor(subdat$logFC, subdat2$logFC, method="spearman", use="pairwise.complete.obs")
    
    overlap[i,"rho"] = as.numeric(rho)
    overlap[i,"p.value"] = R$p.value
    overlaps[[dx]] = overlap
    
    
  }
}

for (dx in names(overlaps)) {
  write.csv(overlaps[[dx]],
            file = paste("./mega_analysis/results/tables/",dx,"_studies_logFC_correlation.csv",
                         sep = ""))
  
}

## reform data for plotting

data.overlap = data.frame()

for (dx in names(overlaps)) {
  dat = as.data.frame(overlaps[[dx]])
  dat[,"Disease"] = dx
  dat[,"study"] = rownames(dat)
  data.overlap = rbind(data.overlap,dat)
}

rownames(data.overlap) = NULL
for (i in 1:2) {
  data.overlap[,i] = as.numeric(data.overlap[,i])
}
data.overlap$fdr = p.adjust(data.overlap$p.value,"fdr")
data.overlap$p.symbol = ""
data.overlap$p.symbol[data.overlap$fdr <0.05 & data.overlap$rho > 0] = "*"
data.overlap$p.symbol[data.overlap$fdr <0.01 & data.overlap$rho >0] = "**"
data.overlap$p.symbol[data.overlap$fdr <0.001 & data.overlap$rho > 0] = "***"
data.overlap$p.symbol[data.overlap$fdr == 0 & data.overlap$rho > 0 ] = "***"

## plot
library(ggplot2)
library(ggthemes)

# color
cols <- c("aquamarine","red","blue","chartreuse",
          "brown","blueviolet","yellow","orange")

p1 = ggplot(data.overlap, aes(x = rho,
                             y = reorder(study, -rho),
                            fill = Disease,
                            label = p.symbol))+ 
  geom_bar(stat = "identity", color = "black")+
  facet_wrap(~Disease,nrow = 1)+
  scale_fill_manual(values = cols)+
  labs(x = "Spearman's correlation",
       y = "Study")+
  scale_x_continuous(limits = c(min(data.overlap$rho),1), 
                     breaks = seq(0,0.8,by = 0.2))+
  geom_text(color="red",
            size=5,hjust = -0.1,position = "identity")+  
  theme_clean(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45))

ggsave("./mega_analysis/results/figures/overlaps/Disease_study_logFC_correlation.pdf",
       width = 12, height = 6)


