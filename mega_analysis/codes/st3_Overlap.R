#Step3_Overlap.R
# here we compare transcriptome correlation for protein coding 
# and non-coding RNAs across diseases

options(stringsAsFactors=F)
library(ggplot2); library(reshape); library(RColorBrewer)
library(NMF); library(WGCNA); library(corrplot); library(purrr); 
library(dplyr); library(biomaRt);library(RRHO);  library(venn)
library(ComplexHeatmap); library(ggrepel); library(ggthemes)
library(ggpubr);library(igraph);library(circlize);library(esc)


# set wd
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")


#### Load data ####

ad_meta = read.csv("../disease_specific/AD_metaAnalysis/tables/AD_sumstats.csv", row.names=1)
asd_meta = read.csv("../disease_specific/ASD_metaAnalysis/tables/ASD_sumstats.csv", row.names=1)
scz_meta = read.csv("../disease_specific/Scz_metaAnalysis/tables/Scz_sumstats.csv", row.names=1)
bp_meta = read.csv("../disease_specific/BP_metaAnalysis/tables/BP_sumstats.csv", row.names=1)
mdd_meta = read.csv("../disease_specific/MDD_metaAnalysis/tables/MDD_sumstats.csv", row.names=1)
pd_meta = read.csv("../disease_specific/PD_metaAnalysis/tables/pj283498_PD_sumstats.csv", row.names=1)
pa_meta = read.csv("../disease_specific/PA_metaAnalysis/tables/PA_sumstats.csv", row.names=1)
psp_meta = read.csv("../disease_specific/PSP_metaAnalysis/tables/PSP_sumstats.csv", row.names=1)


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
  rownames(psp_meta)) #---15819 genes

### save intersect genes 
write.table(all_genes,
            file = "./results/tables/intersect_genes.tsv",
              append = F, col.names = F,row.names = F)

# annotate gene ids

getinfo <- c( "ensembl_gene_id",
              "external_gene_name",
              "gene_biotype")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
datProbes = getBM(
  attributes = getinfo,
  filters = "ensembl_gene_id",
  values = all_genes,
  mart = mart)
rownames(datProbes) = datProbes$ensembl_gene_id

## check consistency ####
all(all_genes == datProbes$ensembl_gene_id)

## data frame to save logFC of diseases
allmeta = matrix(NA,nrow=length(all_genes), 8)
allmeta[,1] = ad_meta$logFC[match(all_genes, rownames(ad_meta))]
allmeta[,2] = pd_meta$logFC[match(all_genes, rownames(pd_meta))]
allmeta[,3] = pa_meta$logFC[match(all_genes, rownames(pa_meta))]
allmeta[,4] = psp_meta$logFC[match(all_genes, rownames(psp_meta))]
allmeta[,5] = scz_meta$logFC[match(all_genes, rownames(scz_meta))]
allmeta[,6] = asd_meta$logFC[match(all_genes, rownames(asd_meta))]
allmeta[,7] = bp_meta$logFC[match(all_genes, rownames(bp_meta))]
allmeta[,8] = mdd_meta$logFC[match(all_genes, rownames(mdd_meta))]

colnames(allmeta) = c("AD", "PD", "PA", "PSP","Scz","ASD","BP","MDD")

rownames(allmeta) = all_genes
allmeta=  as.data.frame(allmeta)

## save logFC data
saveRDS("./results/tables/shared_logFCs.rds",object = allmeta)

## corrplot for allmeta ####
pdf(file = "./results/figures/overlaps/All_pairwise_correlation.pdf", 
    height = 6, width = 5)
corDat = cor(allmeta,method = "spearman",use = "pairwise.complete.obs")
corrplot(corDat,
         mar = c(0,0,3,0),
         addCoef.col = T,
         diag = F,
         order = "hclust",
         hclust.method = "complete",
         method = "circle", 
         tl.col = "black", 
         type = "upper",
         number.cex = 0.9,
         tl.cex = 0.9,
         main = paste("Transcriptome correlation (",
                      nrow(allmeta)," genes)")
)
dev.off()

## divide data based on proteing coding and non-coding genes
geneType.cor = vector("list",3) 
names(geneType.cor) = c("All","protein_coding","non-coding RNA")

for (n in names(geneType.cor)){
  if(n == "All"){
    dat = allmeta
  } else if (n == "protein_coding"){
  id = which(datProbes$gene_biotype == n) 
  dat = allmeta[id,]
  } else {
    id = which(datProbes$gene_biotype != "protein_coding")
    dat = allmeta[id,]
  }
##pvalue matrix
corPval = matrix(NA,ncol = 8,nrow = 8)
colnames(corPval)=colnames(dat)
rownames(corPval) = colnames(dat)

#### Make a barplot ####
comparisons = t(combn(seq(1,ncol(dat)),2))
barplot = data.frame(Mean = NA, SEM=NA, p.fdr=NA, Cohen_d = NA)
for (i in 1:dim(comparisons)[1]) {
  x = comparisons[i,1]
  y = comparisons[i,2]
  R = cor.test(dat[,x], dat[,y], method = "spearman", use = "pairwise.complete.obs")
  rho =cor(dat[,x], dat[,y], method="spearman", use="pairwise.complete.obs")
  sem = (tanh(atanh(rho + 1.96/sqrt(nrow(dat)-3))) - tanh(atanh(rho - 1.96/sqrt(nrow(dat)-3))))/3.92
  cohen = esc_rpb(r = rho, p = R$p.value, grp1n = nrow(dat), grp2n = nrow(dat),es.type = "d") 
  barplot[i,] = c(rho, sem, R$p.value, cohen)
  corPval[x,y] = R$p.value; corPval[y,x]=R$p.value
  corPval[x,x] = 0; corPval[y,y]=0
  rownames(barplot)[i] = paste(colnames(dat)[x],colnames(dat)[y],sep="-")
}

barplot$p.fdr = p.adjust(barplot$p.fdr, method="fdr")
barplot$p.symbol = ""
barplot$p.symbol[barplot$p.fdr <0.05] = "*"
barplot$p.symbol[barplot$p.fdr <0.01] = "**"
barplot$p.symbol[barplot$p.fdr <0.001] = "***"
barplot$Comparison = rownames(barplot)

# save to list
geneType.cor[[n]] = barplot


#color
my_col = if_else(barplot$Mean<0,"indianred1", "seagreen")

symbol = ifelse(barplot$Mean >0.1, barplot$p.symbol, "")

#---barplot usign ggplot
barplot_comp <- ggplot(barplot,
                       aes(x = reorder(Comparison, -Mean), 
                           y= Mean, 
                           label = symbol)) +  
  geom_bar(aes(fill = my_col),
           stat="identity",
           width=0.75) +
  geom_errorbar(aes( ymin = (Mean - SEM), 
                     ymax = (Mean + SEM)),
                position = position_dodge( width = 0.8),
                width = 0.25, size = 0.25) +   
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +   
  labs(x="", y=expression(paste("Transcriptome correlation (", rho, ")", sep="")), 
       title = paste0("Using ",n, " genes (n = ", nrow(dat),")")) +     	
  theme_clean()+
  theme(
    axis.text.x=element_text(angle=50, size=10, hjust=1, colour = "black"),
    axis.text.y=element_text(size=10, vjust=0.5), 
    legend.position = "none", 
    plot.title = element_text(size=18, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=12, vjust=0.5, face = "bold", color = "black"),
    plot.margin=unit(c(2,2,1,2),"mm")
  ) + 
  geom_text(color="red",size=4,aes(y=Mean+ sign(Mean)*SEM + sign(Mean)*.02))+
  scale_y_continuous(breaks = seq(-1,1, by= 0.2))

barplot_comp

ggsave(filename = paste("./results/figures/overlaps/",n,"_comparisons_barplot.pdf",sep = ""),
       barplot_comp,
       height=5, 
       width=8.5)

}

#save correlations
write.csv(x = barplot, file = "./results/tables/Correlations_table.csv")



# custom function for scatterplot
my_scatter <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) + 
    geom_point(alpha = 0.5,
               color = "cornflowerblue") + 
    geom_smooth(method= lm,
                color = "red", 
                se=TRUE, ...)
}


#### create scatterplot matrix ####
library(GGally)

my_custom_cor_color <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {
  # assign('data',data,envir = .GlobalEnv)
  # get the x and y data to use the other code
  x <- eval_data_col(data,mapping$x)
  y <- eval_data_col(data,mapping$y)
  
  ct <- cor.test(x,y,method = "spearman", use = "pairwise.complete.obs")
  
  r <- unname(ct$estimate)
  # r = signif(r,2)
  rt <- format(r, digits=2)[1]
  tt <- as.character(rt)
  
  # plot the cor value
  p <- ggally_text(
    label = tt, 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = 5,
    color=color,
    ...
  ) +
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    )
  
  corColors <- RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")[2:10]
  
  if (r <= -0.4) {
    corCol <- corColors[1]
  } else if (r <= -0.3) {
    corCol <- corColors[2]
  }  else if (r <= -0.2) {
    corCol <- corColors[3]
  } else if (r <= 0) {
    corCol <- corColors[5]
  }  else if (r <= 0.2) {
    corCol <- corColors[6]
  } else if (r <= 0.4) {
    corCol <- corColors[7]
  } else if(r<= 0.6){
    corCol <- corColors[8]
  } else {
    corCol <- corColors[9]
  }
  p <- p + theme(
    panel.background = element_rect(fill = corCol)
  )
  
  p
}

###costum pvalue
my_pval <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {
  x <- eval_data_col(data,mapping$x)
  y <- eval_data_col(data,mapping$y)
  
  ct <- cor.test(x,y,method = "spearman", use = "pairwise.complete.obs")
  
  pval <- unname(ct$p.value)
  pval = signif(pval,1)
  pval <- format(pval, digits=1)[1]
  pval <- as.character(pval)
  
  # plot the p value
  p <- ggally_text(
    label = pval, 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = 5,
    color="black",
    ...
  ) +
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    )
  p <- p + theme(panel.background = element_rect(fill = "lightgreen"))
  p
}

scplot =ggpairs(allmeta,
                lower=list(continuous = my_pval),
                upper = list(continuous = my_custom_cor_color),
                diag = NULL
) +
  labs(title = "") +
  theme(axis.text  = element_blank(), axis.ticks = element_blank())

ggsave(filename = "./results/figures/overlaps/CrossDisease_scatterplotMat.pdf",
       height = 5,
       width = 8,
       plot = scplot)


#### MDS ####
adj = cor(allmeta, use = "pairwise.complete.obs", method = "spearman")
mds = cmdscale(dist(t(allmeta)), eig = T);   
mds$points[,1] = -1* mds$points[,1]
g1 <- graph.adjacency(as.matrix(adj),
                      mode="undirected",
                      weighted=T,
                      diag=FALSE)
layoutFR <- mds$points
col_fun = colorRamp2(c(min(adj),0,max(adj)), c("orange","white","cornflowerblue"))

#color of the edges bases on correlation data
g.keep =  delete.edges(g1, which(abs(E(g1)$weight) <0.05))
g.keep = delete_edges(g.keep, which(-1/E(g.keep)$weight >=5))
edgecolors = numbers2colors(E(g.keep)$weight,
                            colors = blueWhiteRed(1000, gamma=1)[850:150],
                            signed=T,
                            centered=T,
                            lim=c(-1,1))
edgeCol = ifelse(E(g.keep)$weight >0, "cornflowerblue", "indianred")

#label each vertex as disease
label <- rownames(mds$points)
my_layout = layout.mds(g.keep)
g.keep$layout = my_layout

#keep only strong correlations
source("./codes/ig2ggplot2.R")
pdf("./results/figures/overlaps/diseases_mds.pdf",height = 7, width = 10, useDingbats = FALSE)
plot.igraph(g.keep,
            vertex.label.dist=0, 
            vertex.size=18,
            vertex.label.color="white",
            vertex.label.family = "sans",
            vertex.label.cex=0.9,
            vertex.color = "black",
            layout = layoutFR,
            # edge.color = edgecolors,
            edge.color = edgeCol,
            edge.label.cex= 0.8,edge.label.font= 2,
            edge.label.color = "black",
            edge.label = signif(E(g.keep)$weight,2),
            edge.width =20*abs(E(g.keep)$weight),
            axes = T,
            xlim = c(-1,1),
            xlab= "Dim1",
            ylab= "Dim2",
            asp=0.8)

# title(main="MDS of diseases based on transcriptome similarity", cex.main = 1)
dev.off()


#### RRHO plots ####

rrho.list <- list()
for ( i in 1:ncol(allmeta)){
  df <- data.frame(gene = rownames(allmeta),
                   logfc = allmeta[,i])
  rrho.list[[i]] = df
}

names(rrho.list) <- colnames(allmeta)
comparison <- t(combn(seq(1,ncol(allmeta)),2))

jet.colors  <- colorRampPalette(
  c("#00007F", "blue", "#007FFF", "cyan", 
    "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"));
layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE));

mat <- as.matrix(cbind(c(rep(0,6),28),
                       c(rep(0,5),26,27),
                       c(rep(0,4),23:25),
                       c(rep(0,3),19:22),
                       c(rep(0,2),14:18),
                       c(0,8:13),
                       c(1:7)))

mat

# png("./results/figures/RRHO_allComparisons.png", width = 6, height = 4,units = "cm",res = 300)
pdf("./results/figures/overlaps/RRHO_allComparisons.pdf", width = 6, height = 4)
layout(mat)
#layout.show(28)
par(mar = c(0.2, 0.2, 0, 0))
for (i in 1:length(comparison)){
  list1 = comparison[i,1]
  list2 = comparison[i,2]
  R.obj <- RRHO(list1 = rrho.list[[list1]], 
                list2 = rrho.list[[list2]],
                BY = T,
                alternative = "enrichment",
                plots = T, 
                log10.ind = T,
                outputdir = "./results/figures/overlaps/RRHO/", 
                labels = c(names(rrho.list[list1]), names(rrho.list[list2]))
  )
  image(R.obj$hypermat,
        xlab = names(rrho.list[list1]),
        ylab=names(rrho.list[list2]),
        col=jet.colors(100), 
        axes=FALSE,
        main="")
}
dev.off()
