source("http://bioconductor.org/biocLite.R")
library(ggplot2); library(mada); library(reshape); library(pheatmap); library(RColorBrewer)
library(NMF); library(WGCNA); library(corrplot); library(purrr); 
library(dplyr); library(biomaRt);library(RRHO); library(ggplot2); library(venn)
library(circlize);library(ComplexHeatmap)
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")


##Load data

ad_meta = read.csv("./AD_metaAnalysis/AD_mainFocus_stats.csv", row.names=1)
asd_meta = read.csv("./ASD_metaAnalysis/tables/ASD_meta_sumstats.csv", row.names=1)
scz_meta = read.csv("./Scz_metaAnalysis/tables/Scz_meta_sumstats.csv", row.names=1)
bp_meta = read.csv("./BP_metaAnalysis/tables/BP_meta_sumstats.csv", row.names=1)
mdd_meta = read.csv("./MDD_metaAnalysis/tables/MDD_meta_sumstats.csv", row.names=1)
pa_meta = read.csv("./PA_metaAnalysis/PA_mainFocus_stats.csv", row.names=1)
psp_meta = read.csv("./PSP_metaAnalysis/PSP_mainFocus_stats.csv", row.names=1)


all_genes = intersect(
  intersect(
    intersect(
      intersect(
        intersect(
          intersect(rownames(ad_meta),
                  rownames(asd_meta)),
        rownames(scz_meta)),
      rownames(bp_meta)),
    rownames(mdd_meta)),
  rownames(pa_meta)),
  rownames(psp_meta)) #---10226 genes

allmeta = matrix(NA,nrow=length(all_genes), 7)
allmeta[,1] = ad_meta$logFC[match(all_genes, rownames(ad_meta))]
allmeta[,2] = asd_meta$logFC[match(all_genes, rownames(asd_meta))]
allmeta[,3] = scz_meta$logFC[match(all_genes, rownames(scz_meta))]
allmeta[,4] = bp_meta$logFC[match(all_genes, rownames(bp_meta))]
allmeta[,5] = mdd_meta$logFC[match(all_genes, rownames(mdd_meta))]
allmeta[,6] = pa_meta$logFC[match(all_genes, rownames(pa_meta))]
allmeta[,7] = psp_meta$logFC[match(all_genes, rownames(psp_meta))]
colnames(allmeta) = c("AD", "ASD","Scz", "BP", "MDD",  "PA", "PSP")

rownames(allmeta) = all_genes
allmeta=  as.data.frame(allmeta)

#---correlation 
cordat <- cor(allmeta,
              use="pairwise.complete.obs",
              method="spearman")
corrplot(cordat)
pdf(file = "./results/figures/pairwise_MainFocus_correlation.pdf", height = 8, width = 11)
corrplot(cordat,
         title = "Transcriptome correlation between diseases (main cause)",
         mar = c(3,3,4,2), 
         addCoef.col = T,
         order = "hclust",type = "upper", 
         hclust.method = "complete", 
         method = "circle", 
         tl.col = "black", 
         number.cex = 1.5, 
         tl.cex = 1.3
)
dev.off()


save.image(file = "./codes/st3_overlap_MainFocus.Rdata")
