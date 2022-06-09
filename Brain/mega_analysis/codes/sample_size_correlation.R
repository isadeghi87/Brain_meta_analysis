### sample_size_correlation.R
### check correlation of transcriptome correlation with sample size
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/")


# libraries

## data frame to collect sample sizes
sample_size = as.data.frame(matrix(nrow = 1, ncol = 8))
dx = c("AD","ASD","Scz","BP", "MDD","PD","PA","PSP")
colnames(sample_size) = dx

# load data from each phenotype to calculate sample size
files = list.files(".",pattern = "metaAnalysis.Rdata",full.names = T,
                   include.dirs = T,recursive = T)
# load data except fot PD, 
for (i in dx[-6]) {
  print(i)
  f = files[grep(i,files)]
  load(f)
  n = ncol(datExp)
  print(n)
  sample_size[,i] = n
  rm(datExp)
}


## sample size for PD is 77 
sample_size[,"PD"] = 77

## load correlation data 
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/Brain/mega_analysis/")

load("./codes/st3_overlap_lm.Rdata")


### add sample size to the table
comparisons = as.data.frame(t(combn(seq(1,length(dx)),2)))
rownames(comparisons) = paste(dx[comparisons[,1]],
                                     dx[comparisons[,2]],
                                  sep = "-")

comparisons$size  = comparisons$V1
for (i in 1:nrow(comparisons)) {
  
comparisons[i,"size"] = sum(
  as.numeric(sample_size[,comparisons[i,1]]),
                        as.numeric(sample_size[comparisons[i,2]]))
}

## add values to matched comparisons in correlation data
barplot$sample_size = barplot$Mean
barplot$sample_size = comparisons$size[match(rownames(comparisons),
                                             rownames(barplot))]

# correlation test ####
cor_res = cor.test(barplot$Mean,barplot$sample_size,
                   method = "spearman",
                   use = "pairwise.complete.obs",
                   exact = T)
### we obtained a rho = -0.096 and p = 0.62

