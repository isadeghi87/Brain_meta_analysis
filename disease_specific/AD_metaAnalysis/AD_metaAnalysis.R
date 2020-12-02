#--- meta analysis for AD diagnosis using batch correction and linear model

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_metaAnalysis/")

library(WGCNA); library(nlme); library(reshape); 
library(ggplot2);library(biomaRt); library(limma); 
library(readr); library(ggplot2); library(purrr); 
library(stringr); library(dplyr);
library(statmod);library(RColorBrewer);
library(sva);library(edgeR)

# options
enableWGCNAThreads()
allowWGCNAThreads()
condition = TRUE

#----AD
load("../../mega_analysis/codes/st2_lm_multiExpr.Rdata")

AD_datMeta = multiExpr[["AD"]]$datMeta 
AD_datExpr <- as.matrix(multiExpr[["AD"]]$datExpr)
rm(multiExpr,allExp,allMeta)

AD_datMeta$Dx = factor(AD_datMeta$Dx)
AD_datMeta$Study = factor(AD_datMeta$Study)
AD_datMeta$Brain_Region = factor(AD_datMeta$Brain_Region)
AD_datMeta$Brain_Lobe = factor(AD_datMeta$Brain_Lobe)

dx = as.numeric(as.factor(AD_datMeta$Dx))
st = as.numeric(as.factor(AD_datMeta$Study))
rg = as.numeric(as.factor(AD_datMeta$Brain_Region))


# check MDS 

mds = cmdscale(dist(t(AD_datExpr)), eig = T,k = 4)   
pc1 = mds$eig[1]^2 / sum(mds$eig^2)   
pc2 = mds$eig[2]^2 / sum(mds$eig^2)


### filter genes
id = rowSums(cpm(AD_datExpr)>0.5)>= 0.3 *ncol(AD_datExpr)
datExp  = AD_datExpr[id,]

design = model.matrix(~Dx, data = AD_datMeta)
v = voom(datExp,design =design)

#### batch correction
datExp.comb = ComBat(dat = as.matrix(v$E),
                     batch = AD_datMeta$Study,
                     mod = design)



pdf("./figures/pca/pca_norm_counts_batch.pdf", height = 12, width = 6)
par(mfrow=c(2,1))

plotMDS.default(datExp.combat,col= st, pch = dx, main = "Batch-corrected",
                dim.plot = c(1,3))
plotMDS(datExp.combat,col= rg, pch = dx, main = "Batch-corrected")
dev.off()


# pca post-batch corrected
pcNorm = prcomp(t(datExp.comb))

# tsne 
tsn.All = Rtsne(pcNorm$x[,1:10],pca = F,theta = 0)

dat= as.data.frame(tsn.All$Y)
colnames(dat) = c("tSNE1","tSNE2")
dat$Disease = AD_datMeta$Dx
dat$Region = AD_datMeta$Brain_Lobe
dat$Study = AD_datMeta$Study

# plot

mygg = function(data,color,title){
  plt = ggplot(data, aes(x =tSNE1 ,
                         y =tSNE2,
                         color = color))+
    geom_point( alpha =1, size = 2)+
    # scale_color_manual(values = cols)+
    theme_classic2()+
    labs(title = title)
  return(plt)
}

tsDx = mygg(dat,color = dat$Disease, title = "batch-corrected")
tsStudy = mygg(dat, color = dat$Study, "batch-corrected")
tsReg = mygg(dat, dat$Region, "batch-corrected")


sigPlots = ggarrange(tsDx,tsStudy,tsReg, ncol = 1, nrow = 3)
ggsave(plot= sigPlots, filename = "./figures/pca/pca_norm_corrected.pdf",
       width = 6, height = 10)

#### differential expression ####

design = model.matrix(~Dx + Study + Sex + Brain_Lobe, data = AD_datMeta)
colnames(design) = gsub("Dx","",colnames(design))


####  Differential gene expression ####
dat = as.matrix(datExp.comb)
meta = matrix(NA, nrow=nrow(dat), ncol=3)
for(i in 1:nrow(dat)) {
  if(i%%100==0) print(i)
  expr = as.numeric(dat[i,])
  tryCatch({
    meta[i,] = summary(lme(expr~ Dx + Study + Sex+Brain_Lobe ,data = AD_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

meta=as.data.frame(meta)
colnames(meta) = c("logFC", "SE", "p.value")
rownames(meta) = rownames(dat)
meta$fdr = p.adjust(meta$p.value, "fdr")
meta= meta[!apply(is.na(meta),1,any),]


#save data
write.csv(file = "./tables/AD_sumstats.csv",x = meta)

#### DGE for brain lobe ####

lobe_stats <- vector("list", length = length(unique(AD_datMeta$Brain_Lobe)))
names(lobe_stats) <- levels(AD_datMeta$Brain_Lobe)
for (lobe in levels(AD_datMeta$Brain_Lobe)){
  idx = which(AD_datMeta$Brain_Lobe == lobe)
  Meta = AD_datMeta[idx,]
  Expr = dat[,idx]
  meta = matrix(NA, nrow=nrow(Expr), ncol=3)
  if (length(unique(Meta$Subject_ID)) < nrow(Meta)){
    for(i in 1:nrow(Expr)) {
      if(i%%100==0) print(i)
      expr = as.numeric(Expr[i,])
      tryCatch({
        meta[i,] = summary(lme(expr~ Dx+Sex ,data = Meta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
      }, error=function(e){})
    }
    meta=as.data.frame(meta)
    colnames(meta) = c("logFC", "SE", "p.value")
    rownames(meta) = rownames(dat)
    meta$fdr = p.adjust(meta$p.value, "fdr")
  } else {
    mod <- model.matrix(~Dx+Sex , data = Meta)
    fit <- lmFit(Expr, mod)
    efit <- eBayes(fit, trend = T, robust = T)
    meta <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
    meta <- meta[,c(1,4,5)]
    colnames(meta) = c("logFC", "p.value", "fdr")
  }
  lobe_stats[[lobe]] = meta
}



lobe_sumstat= do.call("cbind", lobe_stats)
write.csv(lobe_sumstat, "./tables/AD_lobe_sumstats.csv")

save(file = "./AD_metaAnalysis.Rdata",AD_datExpr,AD_datMeta,datExp.comb,datExp)
