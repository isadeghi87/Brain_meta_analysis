#### study_disease_allComparison.R
#### check correlation between studies and disease results

library(GeneOverlap); library(stringr)
library(tidyr)
library(dplyr)
library(ggpubr)

### set working directory

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/")



files = list.dirs(path = "./studies/",full.names = T, recursive = T)
files = files[grep("tables",files)]

overlaps = vector("list", 8)
names(overlaps) = c("AD","PD","PSP","PA","Scz","ASD","BP","MDD")
value_list = overlaps
genList = list()

for (dx in names(overlaps)){
  # list of study files
  stu.files = list.files(files,pattern = dx ,full.names = T)
  
  n = length(stu.files)
  studies = vector("list",n)
  # name of studies
  studies_names = str_split(basename(stu.files),
                            pattern = paste("_",dx,sep = ""),
                            simplify = T)[,1] 
  # read studies
  for (i in 1:n){
    dat = read.csv(stu.files[i], row.names = 1)
    dat = subset(dat, p.value < 0.05)
    id = rownames(dat)
    studies[[i]] = dat
    n_name = paste(studies_names[[i]],dx, sep = "_")
    genList[[n_name]] = id
    
  }
  
  # disease specific stats
  sumstat = read.csv(paste("./disease_specific/",dx,
                           "_metaAnalysis/tables/",dx,"_sumstats.csv",
                           sep = ""),row.names = 1)
  
  # significant genes
  sumstat_sig = subset(sumstat, p.value<0.05)
  id2 = rownames(sumstat)
  # name of combined disease specific
  id2_name = paste0(dx,"_combined")
  genList[[id2_name]] = id2

}

n = length(genList)
simpson.df = as.data.frame(matrix(ncol = n,nrow = n))
rownames(simpson.df) = colnames(simpson.df) = names(genList)  
jacc.df = simpson.df

for (s in 1:(n)) {
    for(d in 1:(n)){
      idA = genList[[s]]
      idB = genList[[d]]
      simp_val = signif(length(intersect(idA,idB))/min(length(idA),length(idB)),2)
      simpson.df[s,d] = simp_val
      ## check overlap ##
      overlap = as.data.frame(matrix(NA, 
                                     ncol = n, 
                                     nrow = 4))
      rownames(overlap) = c("p.value","Odd ratio","Jaccard index", "Simpson index")
      colnames(overlap) = studies_names
      go.obj = newGeneOverlap(listA = idA,
                              listB = idB,
                              genome.size = length(union(idA,idB)))
      go.obj = testGeneOverlap(go.obj)
        
        # extract pvalue
        # overlap[s,d] = signif(getPval(go.obj),2)
        
        # get jaccard ratio
        jacc.df[s,d] = signif(getJaccard(go.obj),2)
        
        
        
    }
    
}

library(ComplexHeatmap);library(RColorBrewer)

Heatmap(as.matrix(jacc.df),
                        cluster_rows = F,
                        cluster_columns = F,
                        col = brewer.pal(11,"Greens"),
        name = "Jaccard index",column_names_rot = 45)




