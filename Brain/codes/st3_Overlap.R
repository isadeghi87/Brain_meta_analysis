#Step5_Overlap check
#--Overlap.R
rm(list=ls()); options(stringsAsFactors=F)
source("http://bioconductor.org/biocLite.R")
library(ggplot2); library(mada); library(reshape)
library(NMF); library(WGCNA); library(corrplot)
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")


##Load data
ad_meta = read.csv("./results/AD_metaanalysis.csv", row.names=1)
asd_meta = read.csv("./results/ASD_metaanalysis.csv", row.names=1)
scz_meta = read.csv("./results/Scz_metaanalysis.csv", row.names=1)
bp_meta = read.csv("./results/BP_metaanalysis.csv", row.names=1)
mdd_meta = read.csv("./results/MDD_metaanalysis.csv", row.names=1)
pd_meta = read.csv("./results/PD_metaanalysis.csv", row.names=1)
pa_meta = read.csv("./results/PA_metaanalysis.csv", row.names=1)
psp_meta = read.csv("./results/PSP_metaanalysis.csv", row.names=1)


all_genes = intersect(
                      intersect(
                        intersect(
                          intersect(
                            intersect(
                              intersect(
                                intersect(rownames(ad_meta),
                                          rownames(asd_meta)), 
                                rownames(scz_meta)), 
                              rownames(bp_meta)), 
                            rownames(mdd_meta)), 
                          rownames(pd_meta)),
                        rownames(pa_meta)),
                      rownames(psp_meta))

allmeta = matrix(NA,nrow=length(all_genes), 8)
allmeta[,1] = ad_meta$logFC[match(all_genes, rownames(ad_meta))]
allmeta[,2] = asd_meta$logFC[match(all_genes, rownames(asd_meta))]
allmeta[,3] = scz_meta$logFC[match(all_genes, rownames(scz_meta))]
allmeta[,4] = bp_meta$logFC[match(all_genes, rownames(bp_meta))]
allmeta[,5] = mdd_meta$logFC[match(all_genes, rownames(mdd_meta))]
allmeta[,6] = pd_meta$logFC[match(all_genes, rownames(pd_meta))]
allmeta[,7] = pa_meta$logFC[match(all_genes, rownames(pa_meta))]
allmeta[,8] = psp_meta$logFC[match(all_genes, rownames(psp_meta))]

colnames(allmeta) = c("AD", "ASD","Scz", "BP", "MDD", "PD", "PA", "PSP")
rownames(allmeta) = all_genes
allmeta=  as.data.frame(allmeta)
cordat <- cor(allmeta,
    use="pairwise.complete.obs",
    method="spearman")

write.csv(file = "./results/tables/Pairwise_cor.csv", cordat)
pdf(file = "./results/figures/pairwise_cor.pdf", height = 8, width = 11)
corrplot(cordat,
         title = "Spearman's correlation between diseases",
         mar = c(1,1,4,2), 
         addCoef.col = T,
         order = "hclust", 
         hclust.method = "complete", 
         method = "square", 
         tl.col = "black", 
         number.cex = 0.75, tl.cex = 1
          )
dev.off()



#---classify diseases to neurological and psychiatric

#psychiatric
psy_genes <- intersect(intersect(intersect(rownames(scz_meta), rownames(bp_meta)), rownames(asd_meta)), rownames(mdd_meta))

psy_meta = matrix(NA,nrow=length(psy_genes), 4)
psy_meta[,1] = scz_meta$logFC[match(psy_genes, rownames(scz_meta))]
psy_meta[,2] = asd_meta$logFC[match(psy_genes, rownames(asd_meta))]
psy_meta[,3] = bp_meta$logFC[match(psy_genes, rownames(bp_meta))]
psy_meta[,4] = mdd_meta$logFC[match(psy_genes, rownames(mdd_meta))]
colnames(psy_meta) = c("Scz", "ASD", "BP", "MDD")
rownames(psy_meta) = psy_genes
psy_meta=  as.data.frame(psy_meta)

psy_cordat <- cor(psy_meta,
              use="pairwise.complete.obs",
              method="spearman")

write.csv(file = "./results/tables/Psych_cor.csv", psy_cordat)
pdf(file = "./results/figures/psych_pairwise_cor.pdf", height = 8, width = 11)
corrplot(psy_cordat,
         title = "Spearman's correlation between psychiatric disorders",
         mar = c(1,3,4,1), 
         addCoef.col = T,
         order = "hclust", 
         hclust.method = "complete", 
         method = "square", 
         tl.col = "black", 
         number.cex = 0.75, tl.cex = 1)
dev.off()


#neurological
nl_genes <- intersect(intersect(intersect(rownames(ad_meta), rownames(pd_meta)), rownames(pa_meta)), rownames(psp_meta))

nl_meta = matrix(NA,nrow=length(nl_genes), 4)
nl_meta[,1] = ad_meta$logFC[match(nl_genes, rownames(ad_meta))]
nl_meta[,2] = pd_meta$logFC[match(nl_genes, rownames(pd_meta))]
nl_meta[,3] = pa_meta$logFC[match(nl_genes, rownames(pa_meta))]
nl_meta[,4] = psp_meta$logFC[match(nl_genes, rownames(psp_meta))]
colnames(nl_meta) = c("AD", "PD", "PA", "PSP")
rownames(nl_meta) = nl_genes
nl_meta=  as.data.frame(nl_meta)

nl_cordat <- cor(nl_meta,
                  use="pairwise.complete.obs",
                  method="spearman")

write.csv(file = "./results/tables/Neurological_cor.csv", nl_cordat)
pdf(file = "./results/figures/neurol_pairwise_cor.pdf", height = 8, width = 11)
corrplot(nl_cordat,
         title = "Spearman's correlation between neurological diseases",
         mar = c(1,3,4,1), 
         addCoef.col = T,
         order = "hclust", 
         hclust.method = "complete", 
         method = "square", 
         tl.col = "black", 
         number.cex = 0.75, tl.cex = 1)
dev.off()




 #---ASD vs others
dat= data.frame(x=allmeta$ASD, y=allmeta$Scz, Disease="Scz")
dat = rbind(dat, data.frame(x=allmeta$ASD, y=allmeta$MDD, Disease="MDD"))
dat = rbind(dat, data.frame(x=allmeta$ASD, y=allmeta$BP, Disease = "BP"))
dat = rbind(dat, data.frame(x = allmeta$ASD, y = allmeta$AD, Disease = "AD"))
dat = rbind(dat, data.frame(x = allmeta$ASD, y = allmeta$PD, Disease = "PD"))
dat = rbind(dat, data.frame(x = allmeta$ASD, y = allmeta$PA, Disease = "PA"))
dat = rbind(dat, data.frame(x = allmeta$ASD, y = allmeta$PSP, Disease = "PSP"))

asd_corplot <- ggplot(dat,aes(x,y,group=Disease,color=Disease)) +   
  xlim(c(-1.5,1.5)) + ylim(c(-1.5,1.5)) +
  geom_point(alpha=.5,size=2) +   
  geom_smooth(method="lm",fullrange=T) + 
   scale_color_brewer(palette="Dark2")+
  labs(x="ASD log2FC", y="log2FC", title = "ASD vs other diseases") +  
  geom_abline(intercept=0,slope=1,colour=scales::alpha("black",.7),lty=2)+
  theme(plot.title = element_text(size=20, face="bold", hjust = 0.5))
  
asd_corplot

ggsave(filename="./results/figures/ASD_correlation.pdf",asd_corplot,width=11,height=8)

ggplot(allmeta, aes(x = ASD, y = Scz))+
  geom_point()

df <- as.data.frame(allmeta)
for (i in 1:nrow(comparisons)){
  p <-  ggplot( data = df, 
                aes(x = df[,comparisons[i,1]],
                    y = df[,comparisons[i,2]]))+
    geom_point()+
    labs(x = colnames(df[comparisons[i,1]]), y = colnames(df[comparisons[i,2]]),
         title = paste(colnames(df[comparisons[i,1]]), colnames(df[comparisons[i,2]]), sep = "-"))+
    geom_smooth(method = "lm", fullrange=T)
  print(p)
}

#Load null
null = data.frame(read.delim("./results/tables/NullDistribution.txt",head=F))
null = sort(null$V1)
null = cbind(null, data.frame(prob=1-abs(2*seq(1:length(null))/length(null)-1)))
hist(null$null, 50, main="Null Distribution\n40,000 permutations", xlab="Spearman's Rho")

#Make Barplot
comparisons = t(combn(seq(1,ncol(allmeta)),2))
barplot = data.frame(Mean = NA, SEM=NA, p.fdr=NA)
for (i in 1:dim(comparisons)[1]) {
  x = comparisons[i,1]
  y = comparisons[i,2]
  R = cor.test(allmeta[,x], allmeta[,y], method = "spearman", use = "pairwise.complete.obs")
  rho =cor(allmeta[,x], allmeta[,y], method="spearman", use="pairwise.complete.obs")
  sem = (tanh(atanh(rho + 1.96/sqrt(nrow(allmeta)-3))) - tanh(atanh(rho - 1.96/sqrt(nrow(allmeta)-3))))/3.92
  barplot[i,] = c(rho, sem, R$p.value)
  rownames(barplot)[i] = paste(colnames(allmeta)[x],colnames(allmeta)[y],sep="-")
}

barplot$p.fdr = p.adjust(barplot$p.fdr, method="fdr")
barplot$p.bootstrap = null$prob[findInterval(barplot$Mean, null$null)]
barplot$p.symbol = ""
barplot$p.symbol[barplot$p.fdr<0.05] = "*"
barplot$p.symbol[barplot$p.fdr<0.01] = "**"
barplot$p.symbol[barplot$p.fdr<0.001] = "***"
barplot$Comparison = rownames(barplot)

barplot_comp <- ggplot(barplot,
                       aes(x = reorder(Comparison, -Mean), 
                           y=Mean, label=p.symbol,
                            fill = Comparison)) +  
  geom_bar(stat="identity",width=0.75) +
  geom_errorbar(aes(ymin=(Mean - SEM), ymax=(Mean + SEM)), position=position_dodge(width=0.8), width=0.25,size=0.25) +   
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) +   
  labs(x="", y=expression(paste("Transcriptome correlation (", rho, ")", sep="")), 
       title = "Transcriptome correlation across diseases") +     	
  theme(
    axis.text.x=element_text(angle=50, size=12, hjust=1),
    axis.text.y=element_text(size=10, vjust=0.5), 
    legend.position = "none", 
    plot.title = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.x = element_text(size=10, vjust=-0.35, face="bold"),
    axis.title.y = element_text(size=12, vjust=0.5, face = "bold"),
    plot.margin=unit(c(2,2,1,2),"mm")
  ) + 
  geom_text(color="red",size=4,aes(y=Mean+ sign(Mean)*SEM + sign(Mean)*.02))+
  scale_y_continuous(breaks = seq(-1,1, by= 0.2))

barplot_comp

ggsave("./results/figures/comparisons_barplot.pdf", barplot_comp,height=8, width=11)


#---remove non-significant genes


##Calculate the slope of transcriptome overlap using principle components regression
pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  return(list(slope,intercept))
}

dat= data.frame(AD=allmeta$AD, 
                Scz = allmeta$Scz,
                ASD=allmeta$ASD,
                BP=allmeta$BP,
                MDD=allmeta$MDD,
                PD = allmeta$PD,
                PA = allmeta$PA,
                PSP = allmeta$PSP)
dat2= melt(dat,id=1)
dat2$value = as.numeric(dat2$value)

fit.dat <- vector(mode = "list", length = 8)
fit.ad = pcreg(allmeta$AD)
fit.asd = pcreg(allmeta$SCZ, allmeta$ASD)
fit.bd = pcreg(allmeta$SCZ, allmeta$BD)
fit.mdd=pcreg(allmeta$SCZ, allmeta$MDD)

dat2$variable = as.character(dat2$variable)
dat2$variable = gsub("ASD", paste("ASD, slope=", signif(fit.asd[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("BD", paste("BD, slope=", signif(fit.bd[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("MDD", paste("MDD, slope=", signif(fit.mdd[[1]],2), sep=""), dat2$variable)


TxSlope.array=ggplot(dat2,aes(x=SCZ,y=value,color=variable)) + geom_point(alpha=.5) + 
  geom_abline(slope=1, lty=2) + xlim(-1,1) + ylim(-1,1) + 
  geom_abline(slope=fit.asd[[1]], intercept = fit.asd[[2]], color="#F8766D") + 
  geom_abline(slope=fit.bd[[1]], intercept = fit.bd[[2]], color="#00BA38") + 
  geom_abline(slope=fit.mdd[[1]], intercept = fit.mdd[[2]], color="#619CFF") +
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) 

ggsave(TxSlope.array, filename="./results/figures/Manuscript/Fig2B.pdf",width=5,height=5)

TxSlope.array


dat= data.frame(SCZ=rnaseq.scz$log2FC, ASD=rnaseq.asd$log2FC, BD=rnaseq.bd$log2FC)
dat2= melt(dat,id=1)
dat2$value = as.numeric(dat2$value)

fit.asd = pcreg(rnaseq.scz$log2FC, rnaseq.asd$log2FC)
fit.bd = pcreg(rnaseq.scz$log2FC, rnaseq.bd$log2FC)


dat2$variable = as.character(dat2$variable)
dat2$variable = gsub("ASD", paste("ASD, slope=", signif(fit.asd[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("BD", paste("BD, slope=", signif(fit.bd[[1]],2), sep=""), dat2$variable)

TxSlope.rnaseq=ggplot(dat2,aes(x=SCZ,y=value,color=variable)) + geom_point(alpha=.5) + 
  geom_abline(slope=1, lty=2) + xlim(-1,1) + ylim(-1,1) + 
  geom_abline(slope=fit.asd[[1]], intercept = fit.asd[[2]], color="#F8766D") + 
  geom_abline(slope=fit.bd[[1]], intercept = fit.bd[[2]], color="#00BA38") +
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) 






##Plot dendrogram fo the top Genes
all_genes = intersect(intersect(intersect(intersect(rownames(asd_meta), rownames(scz_meta)), rownames(bd_meta)), rownames(mdd_meta)), rownames(aad_meta))
gene.symbols= asd_meta$symbol[match(all_genes, rownames(asd_meta))]

all_beta = matrix(NA,nrow=length(all_genes), 5)
all_beta[,1] = asd_meta$beta[match(all_genes, rownames(asd_meta))]
all_beta[,2] = scz_meta$beta[match(all_genes, rownames(scz_meta))]
all_beta[,3] = bd_meta$beta[match(all_genes, rownames(bd_meta))]
all_beta[,4] = mdd_meta$beta[match(all_genes, rownames(mdd_meta))]
all_beta[,5] = aad_meta$beta[match(all_genes, rownames(aad_meta))]
all_beta=as.data.frame(all_beta)

all_pvals = matrix(NA,nrow=length(all_genes), 5)
all_pvals[,1] = asd_meta$fdr[match(all_genes, rownames(asd_meta))]
all_pvals[,2] = scz_meta$fdr[match(all_genes, rownames(scz_meta))]
all_pvals[,3] = bd_meta$fdr[match(all_genes, rownames(bd_meta))]
all_pvals[,4] = mdd_meta$fdr[match(all_genes, rownames(mdd_meta))]
all_pvals[,5] = aad_meta$fdr[match(all_genes, rownames(aad_meta))]
all_pvals=  as.data.frame(all_pvals)

colnames(all_beta) = colnames(all_pvals) = c("ASD", "SCZ", "BD", "MDD", "AAD")


rowsums=rowSums(all_beta)
idx = order(rowsums,decreasing = T)[1:25]
idx = c(idx, order(rowsums)[1:25])
gene.symbols[idx]
mat.plot = as.matrix(-sign(all_beta[idx,]) * log10(all_pvals[idx,]))

rownames(mat.plot) = gene.symbols[idx]
colnames(mat.plot) = c("ASD", "SCZ", "BD", "MDD", "AAD")

textMat = signif(all_beta,2)
#textMat[all_pvals.fdr>0.05] = ''
textMat = textMat[idx,]







pdf("./results/figures/Manuscript/FigS2.pdf", width=5,height=8)
aheatmap(mat.plot,color=blueWhiteRed(100)[90:0],cexRow=.7, cexCol=1,fontsize=8,border_color="grey",scale="none",treeheight = 10, txt=textMat,Colv=F)
dev.off()


## Compile Table S1 - Transcriptome Signatures for Each Disease
union_genes = union(union(union(union(union(rownames(asd_meta), rownames(scz_meta)), rownames(bd_meta)), rownames(ibd_meta)), rownames(mdd_meta)), rownames(aad_meta))
all_disorders = vector(mode="list",length=6)
all_disorders[[1]] = asd_meta
all_disorders[[2]] = scz_meta
all_disorders[[3]] = bd_meta
all_disorders[[4]] = mdd_meta
all_disorders[[5]] = aad_meta
all_disorders[[6]] = ibd_meta

tableS1 = as.data.frame(matrix(NA,nrow=length(union_genes), 18))
for(i in 1:6) {
  all_disorders[[i]] = all_disorders[[i]][match(union_genes, rownames(all_disorders[[i]])),]
  tableS1[,3*i-2] = all_disorders[[i]]$beta
  tableS1[,3*i-1] = all_disorders[[i]]$p
  tableS1[,3*i] = all_disorders[[i]]$fdr
}

i = 1
for(dx in c("ASD", "SCZ", "BD", "MDD", "AAD", "IBD")) {
  for(var in c("beta_log2FC", "P.value", "FDR")) {
    colnames(tableS1)[i] = paste(dx,var,sep=".")
    i=i+1
  }
}


library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",  host="feb2014.archive.ensembl.org") 
#f=listFilters(ensembl); a =listAttributes(ensembl)
bm = getBM(attributes = c("ensembl_gene_id", "external_gene_id", "hgnc_symbol", "entrezgene", "chromosome_name", "start_position", "end_position"),
           filters="ensembl_gene_id", values=union_genes, mart=ensembl)
bm = bm[match(union_genes, bm$ensembl_gene_id),]
tableS1 = cbind(bm, tableS1)
write.csv(file="./results/tables/Manuscript/TableS1 - Microarray Meta-Analyses.csv", tableS1)
save.image("./codes/st3_overlap.Rdata")
