#--- meta analysis for PD using batch correction and linear model
setwd("/nfs/users2/rg/isadeghi/data/combined/disease_specific/")
library(WGCNA); library(nlme); library(reshape); 
library(ggplot2); library(corrplot); library(biomaRt); library(limma); library(edgeR);
library(readr); library(ggplot2); library(purrr); library(stringr); library(dplyr);
library(RRHO);library(devtools); library(statmod); library(pheatmap);
library(dendextend); library(cerebroViz); library(RColorBrewer); library(sva);
library(ggpubr); library(venn); library(ggridges); library(gridExtra) 


#----PD
load("./codes/st2_lm_multiExpr.Rdata")
datMeta <-
  lm_multiExpr[["PD"]]$datMeta #---73 samples, 73 Subjects
datExpr <-
  as.matrix(lm_multiExpr[["PD"]]$datExpr)#---17252 genes

datMeta$Dx = gsub("Normal","Control", datMeta$Dx)
datMeta$Sex <- as.factor(datMeta$Sex)
datMeta$Dx <- as.factor(datMeta$Dx)
datMeta$Study <- as.factor(datMeta$Study)
datMeta$Subject_ID <- as.factor(datMeta$Subject_ID)
datMeta$Brain_Region <- as.factor(datMeta$Brain_Region)
datMeta$Brain_lobe <- as.factor(datMeta$Brain_Region)


## Covariates
pdf(
  "./PD_metaAnalysis/figures/PD_Covaritates.pdf",
  height = 6,
  width = 11.5
)
par(mfrow = c(2, 3))

#---Subjects
plot(datMeta$Dx,
     col = c("black", "red"),
     main = "Subjects")

#---Age
A = anova(lm(as.numeric(datMeta$Age) ~ datMeta$Dx))
p = A$"Pr(>F)"[1]
plot(
  datMeta$Age ~ datMeta$Dx,
  col = c("black", "red"),
  main = paste("Age \np=", signif(p, 2), sep = ""),
  ylab = "year",
  xlab = ""
)

#---PMI
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   
plot(datMeta$PMI ~ datMeta$Dx, col=c("black", "red"), ylim =c(0,20),
     main=paste("PMI \np=", signif(p,2)), ylab="", xlab="")

#---RIN
A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   
plot(as.numeric(datMeta$RIN) ~ datMeta$Dx, col=c("black", "red"),
     main=paste("RIN \np=", signif(p,2)), ylab="", xlab="")

#---Sex
# A = chisq.test(datMeta$Sex, datMeta$Dx)
# p = A$p.value
# plot(
#   datMeta$Sex ~ datMeta$Dx,
#   col = c("red", "black"),
#   main = paste("Sex \np=", signif(p, 2)),
#   ylab = "",
#   xlab = ""
# )




#--density plot
plot(
  density(datExpr[, 1]),
  xlim = c(-2, 12),
  ylim = c(0, 0.4),
  col = as.numeric(datMeta$Dx[1]),
  xlab = "Intensity (log2)",
  ylab = "Density",
  main = "Expression density plot"
)
for (i in 2:dim(datExpr)[[2]])
  lines(
    density(datExpr[, i]),
    xlim = c(-2, 12),
    col = as.numeric(datMeta$Dx[i])
  )
legend(
  "topright",
  (levels(datMeta$Dx)),
  col = c(1:length(levels(datMeta$Dx))),
  pch = 16,
  cex = 0.7
)

#--PCA plot
mds = cmdscale(dist(t(datExpr)), eig = T)
pc1 = mds$eig[1] ^ 2 / sum(mds$eig ^ 2)
pc2 = mds$eig[2] ^ 2 / sum(mds$eig ^ 2)

#mds dx
plot(
  mds$points,
  col = c("black", "red"),
  pch = 20,
  main = "MDS-Diagnosis",
  xlab = paste("PC1 (", signif(100 * pc1, 3), "%)", sep = ""),
  ylab = paste("PC2 (", signif(100 * pc2, 3), "%)", sep = "")
)

legend(
  "topright",
  levels(datMeta$Dx),
  col = c("black", "red"),
  pch = 16,
  cex = 0.9
)


## Remove outliers based on network connectivity z-scores
normadj <- (0.5 + 0.5 * bicor(datExpr)) ^ 2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku - mean(ku)) / sqrt(var(ku))

# #z-score plot
# plot(
#   1:length(z.ku),
#   z.ku,
#   col = c("black" , "red"),
#   pch = 19,
#   xlab = "",
#   ylab = "Z-score",
#   main = "Outlier detection"
# )
# legend(
#   "bottomright",
#   legend = levels(datMeta$Dx),
#   col = c("black", "red"),
#   pch = 19,
#   cex = .7
# )
# abline(h = -2, lty = 2)
# outliers = (z.ku < -2)
# table(outliers) #-- outliers
# datExpr = datExpr[, !outliers]
# datMeta = datMeta[!outliers, ]



#---dendogram
# tree = hclust(dist(t(datExpr)), method = "average")
# par(mfrow = c(1, 1))
# dx_col = ifelse(datMeta$Dx == "Normal", "black", "red")
# sex_col = ifelse(datMeta$Sex == "male", "black", "red")
# age_col = numbers2colors(
#   datMeta$Age,
#   blueWhitered(100),
#   signed = F,
#   centered  = T,
#   lim = c(min(datMeta$Age, na.rm = T),
#           max(datMeta$Age, na.rm =
#                 T))
# )
# plotDendroAndColors(
#   tree,
#   cbind(
#     dx_col,
#     age_col
#   ),
#   groupLabels = c("Dx", "Age"),
#   cex.colorLabels = 0.8,
#   cex.dendroLabels = 0.15,
#   main = "Dendrogram\npost-batch-effect correction"
# )
dev.off()




#---------------------------------------------Differential gene expression

mod <-
  model.matrix( ~ Dx + PMI + RIN + Age, data = datMeta)
fit <-
  lmFit(
    datExpr,
    mod
  )

efit <- eBayes(fit, trend = T, robust = T)
PD_metaAnalysis <-
  topTable(efit,
           coef = 2,
           number = Inf,
           sort.by = "none")
PD_metaAnalysis_sig <-  PD_metaAnalysis[PD_metaAnalysis$P.Value < 0.05, ] #---4470 genes
PD_metaAnalysis_sig$symbol = datProbes$external_gene_name[match(rownames(PD_metaAnalysis_sig),
                                                                datProbes$ensembl_gene_id_version)]
up <- PD_metaAnalysis_sig[PD_metaAnalysis_sig$logFC>0,]# 1994 genes
down <- PD_metaAnalysis_sig[PD_metaAnalysis_sig$logFC<0,]#2476  genes
PD_metaAnalysis$symbol = datProbes$external_gene_name[match(rownames(PD_metaAnalysis),
                                                            datProbes$ensembl_gene_id_version)]
write.csv(PD_metaAnalysis,
          "./PD_metaAnalysis/tables/PD_meta_sumstats.csv")
write.csv(PD_metaAnalysis_sig,
          "./PD_metaAnalysis/tables/PD_meta_sig_sumstats.csv")

#--number of significant genes
sig.genes <- as.data.frame(matrix(NA, nrow = 2, ncol = 1))
rownames(sig.genes) <- c("Up", "Down")
sig.genes[1,1] <- PD_metaAnalysis_sig[PD_metaAnalysis_sig$logFC>0,] %>%
  nrow()
sig.genes[2,1] <- PD_metaAnalysis_sig[PD_metaAnalysis_sig$logFC<0,] %>%
  nrow()
colnames(sig.genes) = "Frontal"
write.table(sig.genes, "./PD_metaAnalysis/tables/PD_sig_geneNumber.csv", row.names = T)
sigGeneTable = ggtexttable(t(sig.genes), theme = ttheme("mBlueWhite"))
sigGeneTable
ggsave(filename = "./PD_metaAnalysis/figures/sigGenesNumber.pdf", width = 4, height = 5)



##------Annotating gene IDs---###

getinfo <-
  c(
    "ensembl_gene_id_version",
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "band",
    "gene_biotype",
    "percentage_gene_gc_content"
  )
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
datProbes <- getBM(
  attributes = getinfo,
  filters = c("ensembl_gene_id_version"),
  values = rownames(datExpr),
  mart = mart
)
datProbes <-
  datProbes[match(rownames(datExpr), datProbes$ensembl_gene_id_version), ]
datProbes$length = datProbes$end_position - datProbes$start_position
to_keep = !is.na(datProbes$length)
datProbes = datProbes[to_keep, ]


#----check number of biotypes
PD_metaAnalysis_sig$biotype <-
  datProbes$gene_biotype[match(rownames(PD_metaAnalysis_sig),
                               datProbes$ensembl_gene_id_version)]
n_bio <- PD_metaAnalysis_sig %>% group_by(biotype) %>%
  count() %>% arrange(desc(n)) %>% filter(n >= 10) #---filter n >10

#plot
p <- ggplot(data = n_bio, aes(x = reorder(biotype, n),
                              y = n)) +
  ylim(0,4000)+
  geom_bar(stat = "identity", fill = "cornflowerblue") +
  coord_flip() +
  geom_text(aes(label = n, hjust = 1)) +
  labs(x = "", y = "", title = "Number of gene biotypes in PD")

p
ggsave(filename = "./PD_metaAnalysis/figures/PD_biotypes.pdf", height = 4, width = 9, device = "pdf")



save(file = "./codes/PD_metadata.Rdata",datMeta,datExpr)
save.image("./PD_metaAnalysis/PD_metaAnalysis.Rdata")
