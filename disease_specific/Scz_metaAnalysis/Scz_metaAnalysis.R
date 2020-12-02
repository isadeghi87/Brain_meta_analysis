#--- meta analysis for Scz using batch correction and linear model
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/Scz_metaAnalysis/")

library(WGCNA); library(nlme);
library(ggplot2);library(biomaRt); library(limma); 
library(readr); library(ggplot2); library(purrr); 
library(stringr); library(dplyr);library(Rtsne);library(ggthemes)
library(statmod);library(RColorBrewer);
library(sva);library(edgeR);library(ggpubr)

#----Scz
load("../../mega_analysis/codes/st2_lm_multiExpr.Rdata")

Scz_datMeta = multiExpr[["Scz"]]$datMeta 
Scz_datExpr <- as.matrix(multiExpr[["Scz"]]$datExpr)
rm(multiExpr,allExp,allMeta)

Scz_datMeta$Dx = factor(Scz_datMeta$Dx)
Scz_datMeta$Study = factor(Scz_datMeta$Study)
Scz_datMeta$Brain_Region = factor(Scz_datMeta$Brain_Region)
Scz_datMeta$Brain_Lobe = factor(Scz_datMeta$Brain_Lobe)

dx = as.numeric(as.factor(Scz_datMeta$Dx))
st = as.numeric(as.factor(Scz_datMeta$Study))
rg = as.numeric(as.factor(Scz_datMeta$Brain_Region))

# log cpm 
lcpm = cpm(Scz_datExpr,log = T)

# plot
pdf("./figures/pca/pca_no_norm.pdf", width = 6, height = 10)
par(mfrow=c(2,1))
plotMDS(log2(lcpm),col= st, pch = dx, main = "Scz-study")
plotMDS(log2(lcpm),col= rg, pch = dx, main = "Scz brain region")
dev.off()


#### filter genes ####
id = rowSums(cpm(Scz_datExpr)>0.5)>= 0.3 *ncol(Scz_datExpr)
datExp  = Scz_datExpr[id,]
design = model.matrix(~Dx, data = Scz_datMeta)
v = voom(datExp,design =design)

#### batch correction
datExp.comb = ComBat(dat = as.matrix(v$E),
                     batch = Scz_datMeta$Study,
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
pcNorm = prcomp(t(datExp.comb))

# tsne 
theta = 0
perp = 30
int = c(1000,2000, 5000)

    for (i in int){
      tsn.All = Rtsne(pcNorm$x[,1:50], pca = F,
                      theta = theta, perplexity = perp, max_iter = i,
                      check_duplicates = F)
      dat= as.data.frame(tsn.All$Y)
      colnames(dat) = c("tSNE1","tSNE2")
      dat$Disease = Scz_datMeta$Dx
      dat$Region = Scz_datMeta$Brain_Lobe
      dat$Study = Scz_datMeta$Study
      
      
      tsDx = mygg(dat,color = dat$Disease, title = paste(i, sep = "-"))
      tsStudy = mygg(dat, color = dat$Study, paste(i, sep = "-"))
      tsReg = mygg(dat, dat$Region, paste(i, sep = "-"))
      
      plots = ggarrange(tsDx,tsStudy,tsReg, ncol = 3, nrow = 1)
      
      # print(plots)
      ggsave(plot= plots,
             filename = paste("./figures/pca/pca_",i,".pdf",sep = ""),
             width = 14, height = 5)
    }


#### differential expression ####
design = model.matrix(~Dx + Sex + Brain_Lobe, data = Scz_datMeta)
colnames(design) = gsub("Dx","",colnames(design))


####  Differential gene expression ####
meta = matrix(NA, nrow=nrow(datExp.comb), ncol=3)
for(i in 1:nrow(datExp.comb)) {
  if(i%%100==0) print(i)
  expr = as.numeric(datExp.comb[i,])
  tryCatch({
    meta[i,] = summary(lme(expr~ Dx  + Sex + Brain_Lobe ,data = Scz_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

meta=as.data.frame(meta)
colnames(meta) = c("logFC", "SE", "p.value")
rownames(meta) = rownames(datExp.comb)
meta$fdr = p.adjust(meta$p.value, "fdr")
meta= meta[!apply(is.na(meta),1,any),]


#save data
write.csv(file = "./tables/Scz_sumstats.csv",x = meta)

#### DGE for brain lobe ####

lobe_stats <- vector("list", length = length(unique(Scz_datMeta$Brain_Lobe)))
names(lobe_stats) <- levels(Scz_datMeta$Brain_Lobe)
for (lobe in levels(Scz_datMeta$Brain_Lobe)){
  idx = which(Scz_datMeta$Brain_Lobe == lobe)
  Meta = Scz_datMeta[idx,]
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
write.csv(lobe_sumstat, "./tables/Scz_lobe_sumstats.csv")

save(file = "./Scz_metaAnalysis.Rdata",Scz_datExpr,Scz_datMeta,datExp,datExp.comb)

