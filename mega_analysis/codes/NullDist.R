#Step3--calculateNullDistribution
#NullDist.R
seed = commandArgs(trailingOnly=TRUE)
set.seed(seed)

options(stringsAsFactors=F)
#source("http://bioconductor.org/biocLite.R"); biocLite("lme4")
library(lme4); library(corrplot);library(jmuOutlier)
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")

#--##load Compiled expression & metadata for permutation
load("./codes/Permutation_data.RData") 
genes = rownames(multiExpr[[1]]$datExpr)
allmeta = matrix(NA,nrow=length(genes), length(multiExpr))
colnames(allmeta) = c("AD","Scz","BP","MDD","ASD","PD","PSP","PA")
allmeta=  as.data.frame(allmeta)

lmer_apply=function(x, datMeta) {
  if(length(unique(datMeta$Subject_ID)) < nrow(datMeta)) {
    return(summary(lmer(x ~ Dx + Study + (1 | Subject_ID),data=datMeta))$coefficients[2,1])
  } else if(length(unique(datMeta$Study)) > 1) {
    return(summary(lm(x ~ Dx + Study,data=datMeta))$coefficients[2,1])
  }  else {
    return(summary(lm(x ~ Dx, data = datMeta))$coefficients[2,1])
  }
}


for(i in 1:length(multiExpr)) {
  print(i)
  subj = unique(multiExpr[[i]]$datMeta$Subject_ID)
  subj_dx = data.frame(Subject_ID = subj, Dx = multiExpr[[i]]$datMeta$Dx[match(subj, multiExpr[[i]]$datMeta$Subject_ID)])
  subj_dx$Dx = subj_dx$Dx[order(runif(nrow(subj_dx)))] ##Randomly shuffle Dx assignment for each subject
  multiExpr[[i]]$datMeta$Dx = subj_dx$Dx[match(multiExpr[[i]]$datMeta$Subject_ID,subj_dx$Subject_ID)]
  allmeta[,i] = apply(multiExpr[[i]]$datExpr,1,lmer_apply, multiExpr[[i]]$datMeta)
  
}


cor_vec = vector(mode="numeric")
comparisons = t(combn(seq(1,ncol(allmeta)),2))


for(i in 1:nrow(comparisons)) {
  r = cor(allmeta[,comparisons[i,1]], allmeta[,comparisons[i,2]], method = "spearman", use="pairwise.complete.obs")
  print(r)
  cor_vec = c(cor_vec,r)
}
cor_vec
#---permutation p-value
for(i in 1:nrow(comparisons)) {
  r = cor(allmeta[,comparisons[i,1]], allmeta[,comparisons[i,2]], method = "spearman", use="pairwise.complete.obs")
  print(r)
  cor_vec = c(cor_vec,r)
}
cor_vec
perm_vec = vector(mode = "numeric")
for (i in 1:nrow(comparisons)){
  perm.p <- perm.cor.test(allmeta[,comparisons[i,1]], allmeta[,comparisons[i,2]], method = "spearman",
                          num.sim = 100000)$p.value
  perm_vec = c(perm_vec,perm.p)
}

write.table(cor_vec, file="./results/tables/NullDistribution.txt",row.names = F,col.names =F)
write.table(perm_vec, file="./results/tables/perm_pvalue.txt",row.names = F,col.names =F)
save.image("./codes/PermTest.Rdata")
