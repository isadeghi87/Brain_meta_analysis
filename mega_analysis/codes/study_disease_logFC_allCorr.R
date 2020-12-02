#### study_disease_logFC_allCorr.R
#### check correlation between studies and disease results


library(GeneOverlap)
library(tidyr)
library(dplyr)
### set working directory

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/")


### data for studies
files = list.dirs(path = "./studies/",full.names = T, recursive = T)
files = files[grep("tables",files)]
files = list.files(files,full.names = T,recursive = T)
n = length(files)

# a list of studies
data = list()

for (i in 1:n){
  # read study
    dat = read.csv(files[i], row.names = 1)
    # study name
    stud_names = str_split(basename(files[[i]]),
                           pattern ="_stats",
                           simplify = T)[,1] 
    # save to list
    data[[stud_names]] = dat
}

  
# disease specific stats
dx_files = list.files(path = "disease_specific/",pattern = ".csv",full.names = T,recursive = T)
dx_files = dx_files[!grepl("lobe",dx_files)]
dx_files = dx_files[!grepl("Sex",dx_files)]

n = length(dx_files)
for (i in 1:n){
  # read disease results
  dat = read.csv(dx_files[i], row.names = 1)
  # study name
  dx_name = str_split(basename(dx_files[[i]]),
                         pattern ="_sumstats",
                         simplify = T)[,1] 
  # save to list
  data[[dx_name]] = dat
}


#### intersect of common genes 
n = length(data)
com_genes = intersect(rownames(data[[1]]), rownames(data[[2]]))
for(i in 3:n) {com_genes = intersect(com_genes, rownames(data[[n]]))}

#### keep only com_genes for each dataset

for (i in 1:n) {
  dat = data[[i]]
  dat = dat[com_genes,]
  data[[i]] = dat
  
}

## data frame of all logFC
df = as.data.frame(matrix(NA, ncol = n, nrow = length(com_genes)))
rownames(df) = com_genes
colnames(df) = names(data)

# get logFC 
for (i in 1:n){
  dat = data[[i]]
  df[,i] = dat$logFC

}

## spearman correlation test
cor_res = cor(df,method = "spearman",
               use = "pairwise.complete.obs")

## correlation plot
# corrplot(cor_res)


Heatmap(as.matrix(cor_res),
        cluster_rows = T,
        cluster_columns = T,
        # col = brewer.pal(11,"Greens"),
        name = "correlation",
        column_names_rot = 45)
