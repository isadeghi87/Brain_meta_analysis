#--- meta analysis for BP using batch correction and linear model
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/BP_metaAnalysis/")

library(WGCNA); library(nlme);
library(ggplot2);library(biomaRt); library(limma); 
library(readr); library(ggplot2); library(purrr); 
library(stringr); library(dplyr);library(Rtsne);library(ggthemes)
library(statmod);library(RColorBrewer);
library(sva);library(edgeR);library(ggpubr)

#----BP
load("../../mega_analysis/codes/st2_lm_multiExpr.Rdata")

BP_datMeta = multiExpr[["BP"]]$datMeta 
BP_datExpr <- as.matrix(multiExpr[["BP"]]$datExpr)
rm(multiExpr,allExp,allMeta)

BP_datMeta$Dx = factor(BP_datMeta$Dx)
BP_datMeta$Study = factor(BP_datMeta$Study)
BP_datMeta$Brain_Region = factor(BP_datMeta$Brain_Region)
BP_datMeta$Brain_Lobe = gsub("Basal_ganglia","Basal ganglia",BP_datMeta$Brain_Lobe)
BP_datMeta$Brain_Lobe = factor(BP_datMeta$Brain_Lobe)

dx = as.numeric(as.factor(BP_datMeta$Dx))
st = as.numeric(as.factor(BP_datMeta$Study))
rg = as.numeric(as.factor(BP_datMeta$Brain_Region))

# log cpm 
lcpm = cpm(BP_datExpr,log = T)

# plot
pdf("./figures/pca/pca_no_norm.pdf", width = 6, height = 10)
par(mfrow=c(2,1))
plotMDS(log2(lcpm),col= st, pch = dx, main = "BP-study")
plotMDS(log2(lcpm),col= rg, pch = dx, main = "BP brain region")
dev.off()


#### filter genes ####
id = rowSums(cpm(BP_datExpr)>0.5)>= 0.3 *ncol(BP_datExpr)
datExp  = BP_datExpr[id,]
design = model.matrix(~Dx, data = BP_datMeta)
v = voom(datExp,design =design)

#### batch correction
datExp.comb = ComBat(dat = as.matrix(v$E),
                     batch = BP_datMeta$Study,
                     mod = design)


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

# # pca for normalized counts
lcpm = cpm(datExp.comb, log = T)
pcNorm = prcomp(t(datExp.comb))

# tsne 
thet = c(0,0.5,1)
perp = c(20,100,200)
int = c(500,1000,2000)

for(t in thet){
  for ( p in perp){
    for (i in int){
tsn.All = Rtsne(pcNorm$x[,1:50], pca = F,
                theta = t, perplexity = p, max_iter = i)
dat= as.data.frame(tsn.All$Y)
colnames(dat) = c("tSNE1","tSNE2")
dat$Disease = BP_datMeta$Dx
dat$Region = BP_datMeta$Brain_Lobe
dat$Study = BP_datMeta$Study


tsDx = mygg(dat,color = dat$Disease, title = paste(t,p,i, sep = "-"))
tsStudy = mygg(dat, color = dat$Study, paste(t,p,i, sep = "-"))
tsReg = mygg(dat, dat$Region, paste(t,p,i, sep = "-"))

plots = ggarrange(tsDx,tsStudy,tsReg, ncol = 3, nrow = 1)

# print(plots)
ggsave(plot= plots,
       filename = paste("./figures/pca/pca_",t,p,i,".pdf",sep = ""),
       width = 14, height = 5)
    }
  }
}


#### differential expression ####
design = model.matrix(~Dx + Sex + Brain_Lobe, data = BP_datMeta)
colnames(design) = gsub("Dx","",colnames(design))


####  Differential gene expression ####
meta = matrix(NA, nrow=nrow(datExp.comb), ncol=3)
for(i in 1:nrow(datExp.comb)) {
  if(i%%100==0) print(i)
  expr = as.numeric(datExp.comb[i,])
  tryCatch({
    meta[i,] = summary(lme(expr~ Dx  + Sex + Brain_Lobe ,data = BP_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

meta=as.data.frame(meta)
colnames(meta) = c("logFC", "SE", "p.value")
rownames(meta) = rownames(datExp.comb)
meta$fdr = p.adjust(meta$p.value, "fdr")
meta= meta[!apply(is.na(meta),1,any),]


#save data
write.csv(file = "./tables/BP_sumstats.csv",x = meta)

#### DGE for brain lobe ####

lobe_stats <- vector("list", length = length(unique(BP_datMeta$Brain_Lobe)))
names(lobe_stats) <- levels(BP_datMeta$Brain_Lobe)
for (lobe in levels(BP_datMeta$Brain_Lobe)){
  idx = which(BP_datMeta$Brain_Lobe == lobe)
  Meta = BP_datMeta[idx,]
  Expr = datExp.comb[,idx]
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
    rownames(meta) = rownames(datExp.comb)
    meta$fdr = p.adjust(meta$p.value, "fdr")
  } else {
    mod <- model.matrix(~Dx+Sex , data = Meta)
    fit <- lmFit(datExp.comb[,idx], mod)
    efit <- eBayes(fit, trend = T, robust = T)
    meta <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
    meta <- meta[,c(1,4,5)]
    colnames(meta) = c("logFC", "p.value", "fdr")
  }
  lobe_stats[[lobe]] = meta
}



lobe_sumstat= do.call("cbind", lobe_stats)
write.csv(lobe_sumstat, "./tables/BP_lobe_sumstats.csv")

save(file = "./BP_metaAnalysis.Rdata",BP_datExpr,BP_datMeta,datExp,datExp.comb)


