#--- meta analysis for PSP using batch correction and linear model
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/PSP_metaAnalysis/")

library(WGCNA); library(nlme); library(reshape); 
library(ggplot2);library(biomaRt); library(limma); 
library(readr); library(ggplot2); library(purrr); 
library(stringr); library(dplyr);
library(statmod);library(RColorBrewer);
library(sva);library(edgeR)

#----PSP
load("../../mega_analysis/codes/st2_lm_multiExpr.Rdata")

PSP_datMeta = multiExpr[["PSP"]]$datMeta 
PSP_datExpr <- as.matrix(multiExpr[["PSP"]]$datExpr)
rm(multiExpr,allExp,allMeta)

PSP_datMeta$Dx = factor(PSP_datMeta$Dx)
PSP_datMeta$Study = factor(PSP_datMeta$Study)
PSP_datMeta$Brain_Region = factor(PSP_datMeta$Brain_Region)
PSP_datMeta$Brain_Lobe = factor(PSP_datMeta$Brain_Lobe)

dx = as.numeric(as.factor(PSP_datMeta$Dx))
st = as.numeric(as.factor(PSP_datMeta$Study))
rg = as.numeric(as.factor(PSP_datMeta$Brain_Region))




pdf("./figures/pca/pca_no_norm.pdf")

plotMDS(log2(PSP_datExpr),col= dx, pch = rg)
dev.off()


### filter genes
id = rowSums(cpm(PSP_datExpr)>0.5)>= 0.3 *ncol(PSP_datExpr)
datExp  = PSP_datExpr[id,]

# # pca for normalized counts
# pcNorm = prcomp(t(datExpr[keep,]),center = T,scale. = T)

# tsne 
tsn.All = Rtsne(pcNorm$x[,1:10],pca = F,)

dat= as.data.frame(tsn.All$Y)
colnames(dat) = c("tSNE1","tSNE2")
dat$Disease = PSP_datMeta$Dx
dat$Region = PSP_datMeta$Brain_Lobe
dat$Study = PSP_datMeta$Study

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

tsDx = mygg(dat,color = dat$Disease, title = "pre-batch")
tsStudy = mygg(dat, color = dat$Study, "pre-batch")
tsReg = mygg(dat, dat$Region, "pre-batch")

plots = ggarrange(tsDx,tsStudy,tsReg, ncol = 3, nrow = 1)
ggsave(plot= plots, filename = "./figures/pca/pca_norm.pdf",
       width = 14, height = 5)


# tsne 
tsn.batch = Rtsne(pcBatch$x[,1:10],pca = F)

dat.batch= as.data.frame(tsn.batch$Y)
colnames(dat.batch) = c("tSNE1","tSNE2")
dat.batch$Disease = PSP_datMeta$Dx
dat.batch$Region = PSP_datMeta$Brain_Lobe
dat.batch$Study = PSP_datMeta$Study

# plot

disPlot = mygg(data = dat.batch,
               color = dat.batch$Disease, title = "corrected")
stPlot = mygg(dat.batch,color = dat.batch$Study, title = "corrected")
regPlot = mygg(dat.batch,color = dat.batch$Region,  title = "corrected")

sigPlots = ggarrange(disPlot,stPlot,regPlot, ncol = 3, nrow = 1)
ggsave(plot= sigPlots, filename = "./figures/pca/pca_norm_corrected.pdf",
       width = 14, height = 5)

#### differential expression ####

design = model.matrix(~Dx + Sex + Brain_Lobe, data = PSP_datMeta)
colnames(design) = gsub("Dx","",colnames(design))
v = voom(datExp, design)


####  Differential gene expression ####
dat = as.matrix(v$E)
meta = matrix(NA, nrow=nrow(dat), ncol=3)
for(i in 1:nrow(dat)) {
  if(i%%100==0) print(i)
  expr = as.numeric(dat[i,])
  tryCatch({
    meta[i,] = summary(lme(expr~ Dx  + Sex +Brain_Lobe ,data = PSP_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

meta=as.data.frame(meta)
colnames(meta) = c("logFC", "SE", "p.value")
rownames(meta) = rownames(dat)
meta$fdr = p.adjust(meta$p.value, "fdr")
meta= meta[!apply(is.na(meta),1,any),]


#save data
write.csv(file = "./tables/PSP_sumstats.csv",x = meta)

#### DGE for brain lobe ####

lobe_stats <- vector("list", length = length(unique(PSP_datMeta$Brain_Lobe)))
names(lobe_stats) <- levels(PSP_datMeta$Brain_Lobe)
for (lobe in levels(PSP_datMeta$Brain_Lobe)){
  idx = which(PSP_datMeta$Brain_Lobe == lobe)
  Meta = PSP_datMeta[idx,]
  Expr = dat[,idx]
  meta = matrix(NA, nrow=nrow(Expr), ncol=3)
  if (length(unique(Meta$Subject_ID)) < nrow(Meta)){
    for(i in 1:nrow(Expr)) {
      if(i%%100==0) print(i)
      expr = as.numeric(Expr[i,])
      tryCatch({
        meta[i,] = summary(lme(expr~ Dx + Sex ,data = Meta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
      }, error=function(e){})
    }
    meta=as.data.frame(meta)
    colnames(meta) = c("logFC", "SE", "p.value")
    rownames(meta) = rownames(dat)
    meta$fdr = p.adjust(meta$p.value, "fdr")
  } else {
    mod <- model.matrix(~Dx+Sex , data = Meta)
    fit <- lmFit(v[,idx], mod)
    efit <- eBayes(fit, trend = T, robust = T)
    meta <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
    meta <- meta[,c(1,4,5)]
    colnames(meta) = c("logFC", "p.value", "fdr")
  }
  lobe_stats[[lobe]] = meta
}



lobe_sumstat= do.call("cbind", lobe_stats)
write.csv(lobe_sumstat, "./tables/PSP_lobe_sumstats.csv")

save(file = "./PSP_metaAnalysis.Rdata",PSP_datExpr,PSP_datMeta,datExp)

