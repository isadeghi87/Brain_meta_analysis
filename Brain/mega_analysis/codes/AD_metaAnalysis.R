#--- meta analysis for each diagnosis using batch correction and linear model
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/")
library(WGCNA); library(nlme); library(reshape); 
library(ggplot2); library(corrplot); library(biomaRt); library(cqn); library(limma); library(edgeR);
library(Glimma); library(readr); library(ggplot2); library(purrr); library(stringr); library(dplyr);
library(RRHO);library(devtools); library(statmod); library(pheatmap);
library(dendextend); library(cerebroViz); library(RColorBrewer); library(sva);library(ggpubr)

#----AD
#---laod metadata from step 2, compile data 
load("./codes/st2_Permutation_data.RData")
AD_datMeta <- multiExpr[["AD"]]$datMeta  #---1340 samples, 490 unique Subjects
AD_datExpr <- as.matrix(multiExpr[["AD"]]$datExpr)  #---12438 genes 


## Covariates
pdf("./AD_metaAnalysis/figures/AD_Covaritates.pdf", height = 13, width =15 )
par(mfrow = c(5,3))

#---Subjects
plot(AD_datMeta$Dx, col=c("blue", "green"), main= "Subjects", ylim = c(0,1000))

#---Age
A= anova(lm(as.numeric(AD_datMeta$Age) ~ AD_datMeta$Dx))
p= A$"Pr(>F)"[1]
plot(AD_datMeta$Age ~ AD_datMeta$Dx, col= c("blue", "green"), 
     main = paste("Age \np=", signif(p,2), sep=""), ylab="year", xlab="")

#---Sex
A = chisq.test(AD_datMeta$Sex, AD_datMeta$Dx)
p = A$p.value
plot(AD_datMeta$Sex ~ AD_datMeta$Dx, col=ifelse(AD_datMeta$Sex=="female", "green", "blue"), 
     main=paste("Sex \np=", signif(p,2)), ylab="", xlab="")

#---RIN
A = anova(lm(as.numeric(AD_datMeta$RIN) ~ AD_datMeta$Dx))
p = A$`Pr(>F)`[1]
plot(as.numeric(AD_datMeta$RIN) ~ AD_datMeta$Dx, col=c("blue", "green"),
     main=paste("RIN \np=", signif(p,2)), ylab="", xlab="", ylim=c(1,10))

#---PMI
A = anova(lm(as.numeric(AD_datMeta$PMI) ~ AD_datMeta$Dx))
p = A$`Pr(>F)`[1]
plot(as.numeric(AD_datMeta$PMI) ~ AD_datMeta$Dx, col=c("blue", "green"),
     main=paste("PMI \np=", signif(p,2)), ylab="", xlab="", ylim = c(0,700))

#--Brain region
A = chisq.test(AD_datMeta$Brain_Region, AD_datMeta$Dx)
p = A$p.value
plot(AD_datMeta$Dx ~ AD_datMeta$Brain_Region, col = c("blue", "green"),
     main=paste("Brain region \np=", signif(p,2)), ylab="", xlab="")

#--study
AD_datMeta$Study = as.factor(AD_datMeta$Study)
A = chisq.test(AD_datMeta$Study, AD_datMeta$Dx)
p = A$p.value
plot(AD_datMeta$Dx ~ AD_datMeta$Study, col = c("blue", "green"),
     main=paste("Study \np=", signif(p,2)), ylab="", xlab="")

#--density plot
plot(density(AD_datExpr[,1]), xlim=c(-2,12), ylim=c(0, 0.4), col = as.numeric(AD_datMeta$Study[1]), 
     xlab="Intensity (log2)", ylab="Density", main="Expression density plot")
for(i in 2:dim(AD_datExpr)[[2]])
  lines(density(AD_datExpr[,i]), xlim=c(-2,12), col = as.numeric(AD_datMeta$Study[i]))  
legend("topright", (levels(AD_datMeta$Study)), col=c(1:8), pch=16, cex=0.9)

#--PCA plot
mds = cmdscale(dist(t(AD_datExpr)), eig = T)   
pc1 = mds$eig[1]^2 / sum(mds$eig^2)   
pc2 = mds$eig[2]^2 / sum(mds$eig^2)

#mds dx
plot(mds$points, col=c("blue", "green"), 
     pch=20, main="MDS-Diagnosis \nPre-batch-effect correction", 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), 
     ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
legend("topright", levels(AD_datMeta$Dx), col=c("blue", "green"), pch=16, cex=0.9)
#mds study
plot(mds$points, col=as.numeric(as.factor(AD_datMeta$Study)), 
     pch=20, main="MDS-Study \nPre-batch-effect correction", 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), 
     ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
legend("topright", (levels(AD_datMeta$Study)), col=c(1:3), pch=16, cex=0.9)



#Normalize by Study
mod = model.matrix(~Dx+Age+Sex, data=AD_datMeta)
batch = as.factor(AD_datMeta$Study)
combat.Expr = ComBat(AD_datExpr, batch=batch, mod=mod, prior.plots = F)

#--post correction density plot 
plot(density(combat.Expr[,1]), xlim=c(-2,12), ylim=c(0, 0.4), col = as.numeric(AD_datMeta$Study[1]), 
     xlab="Intensity (log2)", ylab="Density", main="Expression density plot\n post-batch-effect correction")
for(i in 2:dim(combat.Expr)[[2]])
  lines(density(combat.Expr[,i]), xlim=c(-2,12), col = as.numeric(AD_datMeta$Study[i]))  
legend("topright", (levels(AD_datMeta$Study)), col=c(1:3), pch=16, cex=0.9)

#--PCA plot
post.mds = cmdscale(dist(t(combat.Expr)), eig = T)   
pc1 = post.mds$eig[1]^2 / sum(post.mds$eig^2)   
pc2 = post.mds$eig[2]^2 / sum(post.mds$eig^2)

#post.mds dx
plot(post.mds$points, col=c("blue", "green"), 
     pch=20, main="MDS-Diagnosis \npost-batch-effect correction", 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), 
     ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
legend("topright", levels(AD_datMeta$Dx), col=c("blue", "green"), pch=16, cex=0.9)
#mds study
plot(post.mds$points, col=as.numeric(as.factor(AD_datMeta$Study)), 
     pch=20, main="MDS-Study \npost-batch-effect correction", 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), 
     ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
legend("topright", (levels(AD_datMeta$Study)), col=c(1:3), pch=16, cex=0.9)

#---dendogram
tree = hclust(dist(t(combat.Expr)), method = "average")
par(mfrow=c(1,1))
sex_col = ifelse(AD_datMeta$Sex == "male", "blue", "green")
age_col = numbers2colors(AD_datMeta$Age, blueWhiteRed(100), 
                         signed=F, centered=T, lim=c(min(AD_datMeta$Age, na.rm=T),
                                                     max(AD_datMeta$Age, na.rm=T)))
pmi_col = numbers2colors(AD_datMeta$PMI, blueWhiteRed(100), 
                         signed=F, centered=T, lim=c(min(AD_datMeta$PMI, na.rm=T),
                                                     max(AD_datMeta$PMI, na.rm=T)))
rin_col = numbers2colors(AD_datMeta$RIN, blueWhiteRed(100), 
                         signed=F, centered=T, lim=c(min(AD_datMeta$RIN, na.rm=T),
                                                     max(AD_datMeta$RIN, na.rm=T)))
plotDendroAndColors(tree, cbind(as.numeric(AD_datMeta$Dx), 
                                as.numeric(AD_datMeta$Study),
                                sex_col, age_col, pmi_col, rin_col), 
                    groupLabels = c("Dx", "Study", "Sex", "Age", "PMI", "RIN"), 
                    cex.colorLabels=0.8, cex.dendroLabels=0.15,
                    main="Dendrogram\npost-batch-effect correction")
dev.off()


#---Differential gene expression

mod <- model.matrix(~Dx + Brain_Region + Sex + PMI + Age + RIN, data = AD_datMeta)
fit <- lmFit(combat.Expr, mod)
efit <- eBayes(fit, trend = T, robust = T)
AD_metaAnalysis <- topTable(efit, coef = 2, number = Inf, sort.by = "none") 
AD_metaAnalysis_sig <- AD_metaAnalysis[AD_metaAnalysis$adj.P.Val<0.05,]
write.csv(AD_metaAnalysis, "./AD_metaAnalysis/tables/AD_meta_sumstats.csv")

#---DGE for brain region
region_stats <- vector("list", length = length(unique(AD_datMeta$Brain_Region)))
names(region_stats) <- levels(AD_datMeta$Brain_Region)
for (region in levels(AD_datMeta$Brain_Region)){
  idx = which(AD_datMeta$Brain_Region == region)
  mod.region = model.matrix(~ Dx + Age + Sex + PMI + RIN, data=AD_datMeta[idx,]) 
  fit.region <- lmFit(combat.Expr[,idx], mod.region)
  efit.region <- eBayes(fit.region, trend = T, robust = T)
  region_stats[[region]] = topTable(efit.region, coef=2, number=Inf, sort.by="none")
}

AD_region_sumstat= do.call("cbind", region_stats)
write.csv(AD_region_sumstat, "./AD_metaAnalysis/tables/AD_region_sumstats.csv")
save.image("./AD_metaAnalysis/AD_metaAnalysis.Rdata")
