#Step3_Overlap check using lm data 
 options(stringsAsFactors=F)
source("http://bioconductor.org/biocLite.R")
library(ggplot2); library(mada); library(reshape); library(pheatmap); library(RColorBrewer)
library(NMF); library(WGCNA); library(corrplot); library(purrr); 
library(dplyr); library(biomaRt);library(RRHO);  library(venn)
 library(circlize);library(ComplexHeatmap); library(ggrepel); library(ggpubr)
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")


##Load data

ad_meta = read.csv("../disease_specific/AD_metaAnalysis/tables/AD_meta_sumstats.csv", row.names=1)
asd_meta = read.csv("../disease_specific/ASD_metaAnalysis/tables/ASD_meta_sumstats.csv", row.names=1)
scz_meta = read.csv("../disease_specific/Scz_metaAnalysis/tables/Scz_meta_sumstats.csv", row.names=1)
bp_meta = read.csv("../disease_specific/BP_metaAnalysis/tables/BP_meta_sumstats.csv", row.names=1)
mdd_meta = read.csv("../disease_specific/MDD_metaAnalysis/tables/MDD_meta_sumstats.csv", row.names=1)
pd_meta = read.csv("../disease_specific/PD_metaAnalysis/tables/PD_meta_sumstats.csv", row.names=1)
pd_meta= rename_(pd_meta, "p.value" = "P.Value", "fdr" = "adj.P.Val") 
pa_meta = read.csv("../disease_specific/PA_metaAnalysis/tables/PA_meta_sumstats.csv", row.names=1)
psp_meta = read.csv("../disease_specific/PSP_metaAnalysis/tables/PSP_meta_sumstats.csv", row.names=1)


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
  rownames(psp_meta)) #---10315 genes

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

#---correlation 
cordat <- cor(allmeta,
              use="pairwise.complete.obs",
              method="spearman")

write.csv(file = "./results/tables/Pairwise_cor_lm.csv", cordat)
pdf(file = "./results/figures/pairwise_correlation.pdf", height = 8, width = 11)
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

#---heatmap for comparing cross disease expression

pheatmap(t(allmeta), 
         color = redWhiteGreen(1000,1),
         scale = "none",
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = F,
         main = "Cross-disease differential gene expression",
         cellheight  = 15, 
         treeheight_col = 20, 
         treeheight_row = 20, 
         legend_breaks = c(-1.5, -1, 0, 1, 1.5, max(allmeta)),
         filename = "./results/figures/CrossDiseaseHeatmap.pdf",
         height = 4, 
         width = 8,
         fontsize_row = 8,
         legend_labels = c("-1.5", "-1", "0", "1","1.5", "logFC\n")
         )

dev.off()
#-------------------------------------RRHO plots
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

png("./results/figures/RRHO_allComparisons.png", width = 6, height = 4,units = "cm",res = 300)
pdf("./results/figures/RRHO_allComparisons.pdf", width = 6, height = 4)
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
                plots = F, 
                log10.ind = T,
                outputdir = "./results/figures/RRHO/all/", 
                labels = c(names(rrho.list[list1]), names(rrho.list[list2]))
  )
image(R.obj$hypermat,
      xlab = names(rrho.list[list1]),
                     ylab=names(rrho.list[list2]),
                     col=jet.colors(100), 
                     axes=FALSE, main="")
}
dev.off()

#---#--------------------------------------------------circos plot
getinfo <-
  c("ensembl_gene_id_version",
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "band",
    "gene_biotype")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
bm <- getBM(
  attributes = getinfo,
  filters = c("ensembl_gene_id_version"),
  values = rownames(allmeta),
  mart = mart)

bm<- bm[match(rownames(allmeta), bm$ensembl_gene_id_version),]
circDat <- cbind(allmeta, bm)%>% dplyr::select(ensembl_gene_id, chromosome_name, start_position,
                                             everything())
circDat <- circDat[,1:11]
rownames(circDat)<- NULL
circDat <- circDat %>% arrange(chromosome_name, start_position)
for (i in 1:nrow(circDat)){
  circDat[i,"chromosome_name"]= gsub(circDat[i,"chromosome_name"],
                                     paste("chr", circDat[i,"chromosome_name"], sep = ""),
                                     circDat[i,"chromosome_name"])
}


#Make Barplot
comparisons = t(combn(seq(1,ncol(allmeta)),2))
barplot = data.frame(Mean = NA, SEM=NA, p.fdr=NA, Cohen_d = NA)
for (i in 1:dim(comparisons)[1]) {
  x = comparisons[i,1]
  y = comparisons[i,2]
  R = cor.test(allmeta[,x], allmeta[,y], method = "spearman", use = "pairwise.complete.obs")
  rho =cor(allmeta[,x], allmeta[,y], method="spearman", use="pairwise.complete.obs")
  sem = (tanh(atanh(rho + 1.96/sqrt(nrow(allmeta)-3))) - tanh(atanh(rho - 1.96/sqrt(nrow(allmeta)-3))))/3.92
  cohen = esc_rpb(r = rho, p = R$p.value, grp1n = nrow(allmeta), grp2n = nrow(allmeta),es.type = "d") 
  barplot[i,] = c(rho, sem, R$p.value, cohen)
  rownames(barplot)[i] = paste(colnames(allmeta)[x],colnames(allmeta)[y],sep="-")
}

barplot$p.fdr = p.adjust(barplot$p.fdr, method="fdr")
barplot$p.symbol = ""
barplot$p.symbol[barplot$p.fdr <0.05] = "*"
barplot$p.symbol[barplot$p.fdr <0.01] = "**"
barplot$p.symbol[barplot$p.fdr <0.001] = "***"
barplot$Comparison = rownames(barplot)

#save correlations
write.csv(x = barplot, file = "./results/tables/Correlations_table.csv")

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
       title = "Transcriptome correlation across diseases") +     	
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

ggsave("./results/figures/comparisons_barplot.pdf",
       barplot_comp,
       height=5, 
       width=8.5)


#---create a scatterplot matrix
# custom function for density plot
my_density <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) + 
    geom_density(alpha = 0.5,
                 fill = "green", ...)
}

# custom function for scatterplot
my_scatter <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) + 
    geom_point(alpha = 0.5,
               color = "cornflowerblue") + 
    geom_smooth(method= lm,
                color = "red", 
                se=TRUE, ...)
}


# create scatterplot matrix
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
    size = 7,
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

ggpairs(allmeta,
        lower=list(continuous = my_scatter),
        upper = list(continuous = my_custom_cor_color),
        diag = list(continuous = my_density)) +
  labs(title = "Cross-Disease Transcriptome overlap") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        axis.text.x.bottom = element_text(size = 6),
        axis.text.y.left = element_text(size = 6),
        strip.text.x = element_text(size = 17),
        strip.text.y = element_text(size = 15))

grid.echo()
f2a = grid.grab()
saveRDS(object = f2a,file = "./codes/Fig2a_scplot.rds")
ggsave(filename = "./results/figures/CrossDisease_scatterplotMat.pdf",
       height = 5,
       width = 8,
       device = "pdf",
       plot = scPlot)

#---------------------MDS
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
#add legend to the graph
# lgd_links = Legend(at = c( -1, 0, 1), col_fun = col_fun, 
                   # title_position = "topleft", title = "rho score")
# lgd_etl = Legend(at = c("Neuropsychiatric disorder",
                    #        "Neurological disorder"),border = c("black", "black"),
                    # legend_gp = gpar(fill= c("tomato2", "green")),
                    # title_position = "topleft",
                    # title = "Type")
# lgd_list_vertical = packLegend(lgd_links, lgd_etl)

#label each vertex as disease
label <- rownames(mds$points)
# label_col = if_else(label %in% c("AD", "PD", "PSP", "PA"), "green","tomato2")

my_layout = layout.mds(g.keep)
g.keep$layout = my_l
#keep only stron correlations
source("./codes/ig2ggplot2.R")
pdf("./results/figures/diseases_mds.pdf",height = 7, width = 10, useDingbats = FALSE)
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
            asp=0.8)
# title(main="MDS of diseases based on transcriptome similarity", cex.main = 1)
grid.echo()
fig2c = grid.grab()
fig2c = as.ggplot(fig2c)

saveRDS(object = fig2c, file ="./codes/fig2c.rds")
# draw(lgd_list_vertical, x = unit(23, "cm"), y = unit(12, "cm"), just = c("top","right"))
dev.off()

#--------------------------------------Plot dendrogram fo the top Genes

##------Annotating gene IDs

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
  values = rownames(mdd_meta),
  mart = mart
)
datProbes <-  datProbes[match(all_genes, datProbes$ensembl_gene_id_version), ]
gene.symbols= datProbes$external_gene_name

#--logFCs
all_logfc = matrix(NA,nrow=length(all_genes), 8)
all_logfc[,1] = ad_meta$logFC[match(all_genes, rownames(ad_meta))]
all_logfc[,2] = asd_meta$logFC[match(all_genes, rownames(asd_meta))]
all_logfc[,3] = scz_meta$logFC[match(all_genes, rownames(scz_meta))]
all_logfc[,4] = bp_meta$logFC[match(all_genes, rownames(bp_meta))]
all_logfc[,5] = mdd_meta$logFC[match(all_genes, rownames(mdd_meta))]
all_logfc[,6] = pd_meta$logFC[match(all_genes, rownames(pd_meta))]
all_logfc[,7] = pa_meta$logFC[match(all_genes, rownames(pa_meta))]
all_logfc[,8] = psp_meta$logFC[match(all_genes, rownames(psp_meta))]

#--P-values
all_pvals = matrix(NA,nrow=length(all_genes), 8)
all_pvals[,1] = ad_meta$p.value[match(all_genes, rownames(ad_meta))]
all_pvals[,2] = asd_meta$p.value[match(all_genes, rownames(asd_meta))]
all_pvals[,3] = scz_meta$p.value[match(all_genes, rownames(scz_meta))]
all_pvals[,4] = bp_meta$p.value[match(all_genes, rownames(bp_meta))]
all_pvals[,5] = mdd_meta$p.value[match(all_genes, rownames(mdd_meta))]
all_pvals[,6] = pd_meta$p.value[match(all_genes, rownames(pd_meta))]
all_pvals[,7] = pa_meta$p.value[match(all_genes, rownames(pa_meta))]
all_pvals[,8] = psp_meta$p.value[match(all_genes, rownames(psp_meta))]

all_pvals=  as.data.frame(all_pvals)

colnames(all_pvals) =  colnames(all_logfc) = c("AD", "ASD","Scz", "BP", "MDD", "PD", "PA", "PSP")

#top differentially expressed genes with using mean
rowsums=rowSums(all_logfc)
idx = order(rowsums,decreasing = T)[1:25]
idx = c(idx, order(rowsums)[1:25])
gene.symbols[idx]
mat.plot = as.matrix(-log10(all_pvals[idx,]))

rownames(mat.plot) = gene.symbols[idx]


textMat = signif(all_logfc,2)
#textMat[all_pvals.fdr>0.05] = ''
textMat = textMat[idx,]
rownames(textMat) = rownames(mat.plot)

#plot the heatmap
pheatmap(textMat,
  color = redWhiteGreen(1000),
  cluster_cols = T,
  cluster_rows = T,
  fontsize_col = 12,
  fontsize_number = 9,
  fontsize_row = 8,
  treeheight_row = 10,
  cellwidth = 25,
  angle_col = 45,
  show_rownames = T,
  show_colnames = T,
  main = "Top DEGs across all diseases",
  legend_breaks = c(-1,-0.5,0,0.5,1, max(textMat)),
  legend_labels = c("-1","-0.5","0","0.5","1","logFC\n"),
  treeheight_col = 10, 
  display_numbers = signif(mat.plot,1),
  filename = "./results/figures/topGenes_heatmap.pdf",
  height = 6,
  width = 6
)

#checking with median
rowmed=rowMedians(all_logfc)
id = order(rowmed,decreasing = T)[1:25]
id = c(id, order(rowmed)[1:25])

#save gene ids for pca
top_genes = rownames(allmeta)[id]
write.csv(top_genes,file = "./results/tables/top_genes.csv", col.names = F)

#
gene.symbols[id]
med.plot = as.matrix( -log10(all_pvals[id,]))
rownames(med.plot) = gene.symbols[id]

labelMat = signif(all_logfc,2)
labelMat = labelMat[id,]
rownames(labelMat) = rownames(med.plot)

#plot the heatmap
pheatmap(labelMat,
         color = redWhiteGreen(1000),
         cluster_cols = T,
         cluster_rows = T,
         fontsize_col = 12,
         fontsize_number = 9,
         fontsize_row = 8,
         treeheight_row = 10,
         cellwidth = 25,
         angle_col = 45,
         show_rownames = T,
         show_colnames = T,
        main = "Top DEGs across all diseases",
          legend_breaks = c(-1,-0.5,0,0.5,1, max(labelMat)),
         legend_labels = c("-1","-0.5","0","0.5","1","logFC\n"),
         treeheight_col = 10, 
         display_numbers = signif(med.plot,1),
        filename = "./results/figures/topGenesMedian_heatmap.pdf",
         height = 6,
         width = 6
)

#check significant genes

#remove pvalue for further visualization
regGenes=regGenes[,1:3]

#divide by up and downregulated
reglist = vector("list", 2)
names(reglist) = c("up", "down")
for (i in names(reglist)){
  if (i == "down"){
    dat = regGenes[regGenes$logFC<0,]
  } else {
    dat = regGenes[regGenes$logFC>0,]
  }
  dat = spread(dat,key = Dx, value = logFC)
  
  #classify for neurological and psychiatric
  # neur = dat[dat$Dx %in% c("AD", "PD", "PSP", "PA"),]
  # psych = dat[dat$Dx %in% c("ASD", "Scz", "MDD", "BP"),]
  
  neur = na.omit(dat[,c("gene","AD", "PD", "PSP", "PA")])
  psych = na.omit(dat[,c(c("gene","ASD", "Scz", "MDD", "BP"))])
  unionGene = union(psych$gene,neur$gene)
  
  dat = allmeta[match(unionGene, rownames(allmeta)),]
  # 
  rownames(dat) =attr$external_gene_name[match(unionGene, attr$ensembl_gene_id_version)]
  reglist[[i]] = dat
}

#plot the heatmap
library(ComplexHeatmap)
pdf("./results/figures/GenesEtiology.pdf",
    width = 18,height = 6)
Heatmap(t(sigdata),
        col = redWhiteGreen(1000,gamma = 0.5),
        name = "logFC",
        cluster_rows = T,
        clustering_distance_rows = "euclidean",
        cluster_columns = F,
        right_annotation = rowAnnotation(Etiology = c("neurological",
                                                      "neuropsychiatric", 
                                                      "neuropsychiatric",
                                                      "neuropsychiatric",
                                                      "neuropsychiatric",
                                                      "neurological",
                                                      "neurological",
                                                      "neurological"),
                                         col =list(Etiology = c("neurological" = "blue",
                                                                "neuropsychiatric" = "black"))),
        top_annotation = columnAnnotation(
          Class = c(rep("up",50), rep("down",108)),
          col = list(Class = c("up" = "green","down" = "red" ))),
        show_row_names = T,
        show_column_names =T,
        row_split = 2,
        width = unit(35,"cm") ,
        height = unit(10,"cm"),
        # column_names_rot = 90,
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 15),
        show_heatmap_legend = T,
        column_title = "Dysregulated genes in each etiology")


dev.off()
## Compile Table of Transcriptome Signatures for Each Disease
union_genes = 
  union(
  union(
    union(
      union(
        union(
          union(
            union(rownames(ad_meta), 
                  rownames(asd_meta)), 
            rownames(scz_meta)), 
          rownames(bp_meta)),
        rownames(mdd_meta)),
      rownames(pd_meta)),
    rownames(pa_meta)),
  rownames(psp_meta))

all_disorders = vector(mode="list",length=8)
names(all_disorders) <- c("AD","ASD","Scz", "BP", "MDD", "PD", "PA", "PSP")
all_disorders[[1]] = ad_meta
all_disorders[[2]] = asd_meta
all_disorders[[3]] = scz_meta
all_disorders[[4]] = bp_meta
all_disorders[[5]] = mdd_meta
all_disorders[[6]] = pd_meta
all_disorders[[7]] = pa_meta
all_disorders[[8]] = psp_meta

gene_table = as.data.frame(matrix(NA,nrow=length(union_genes), 24))
for(i in 1:8) {
  all_disorders[[i]] = all_disorders[[i]][match(union_genes, rownames(all_disorders[[i]])),]
  gene_table[,3*i-2] = all_disorders[[i]]$logFC
  gene_table[,3*i-1] = all_disorders[[i]]$p.value
  gene_table[,3*i] = all_disorders[[i]]$fdr
}

i = 1
for(dx in c("AD", "ASD", "Scz", "BP", "MDD", "PD", "PA", "PSP")) {
  for(var in c("logFC", "P.Value", "adj.P.Val")) {
    colnames(gene_table)[i] = paste(dx,var,sep=".")
    i=i+1
  }
}



getinfo <-
  c(
    "ensembl_gene_id_version",
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "band",
    "gene_biotype"
  )
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
attr <- getBM(
  attributes = getinfo,
  filters = c("ensembl_gene_id_version"),
  values = union_genes,
  mart = mart
)
attr = attr[match(union_genes, attr$ensembl_gene_id_version),]
gene_table = cbind(attr, gene_table)

write.csv(file="./results/tables/disease_signatures.csv", gene_table)

#########################
# Disease-specific DEGs #
#########################

# list of plots
listPlots = vector("list",8)
names(listPlots) = unique(colnames(all_logfc))

#table of number of disease-specific DEGs
DEG.table = as.data.frame(matrix(NA, nrow = 8, ncol = 1))
rownames(DEG.table) = unique(colnames(all_logfc))
colnames(DEG.table) = "DEGs (P < 0.05)"

#frequency of significant genes
freq.gene = vector()
for (i in names(all_disorders)){
  dat = as.data.frame(all_disorders[[i]])
  dat$gene = attr$external_gene_name[match(rownames(dat), attr$ensembl_gene_id_version)]
#filter for significant genes (p<0.05)
  sig= subset(dat, p.value <0.05)
  DEG.table[i,1] = nrow(sig)
  #order top DEGs 
  subdata = order(abs(sig$logFC), decreasing = T)
  subdata = subdata[1:10]
  sub = sig[subdata,]
  
  freq.gene = c(na.omit(freq.gene),sub$gene)
  #color 
  col = ifelse(-log10(dat$p.value)>1.3 , "black", "grey")
  top.col = ifelse(sub$logFC >0, "green", "red")
  
  #plot
  p = ggplot(data = dat , aes(x = logFC,
                            y = -log10(p.value)))+
    geom_point(size = 0.7, color =col)+
    geom_point(data =sub, aes(x =logFC, 
                              y = -log10(p.value)), 
               color = top.col,
               size = 1.5)+
    geom_label_repel(data = sub,
                     aes(x = logFC,
                         y = -log10(p.value),
                         label = gene), 
                     color = "blue",
                     vjust = -1,
                     size = 1.8)+
    geom_hline(yintercept = -log10(0.05), color = "red", lty = 2)+
    labs(y = "-log10(p.value)",
      title = i)+
    theme_bw()+
    theme(text = element_text(size = 10 ,colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 12))
  #save the plot in the list
  listPlots[[i]] = p
}

freq.gene = as.data.frame(table(freq.gene))
freq.gene %<>% arrange(desc(Freq))

p = ggarrange(plotlist = listPlots, ncol = 4, nrow = 2)
saveRDS(p,file = "./codes/fig3a.rds")
ggsave(filename = "./results/figures/Dis_Specific_DEG.pdf", width = 12, height = 6, plot = p)
write.table(DEG.table, "./results/tables/Number_DEGs.csv",col.names = T)

#plot the DEGs tagle
texTable =ggtexttable(t(DEG.table), theme = ttheme("mGreen"))
ggsave(plot = texTable,
       filename = "./results/figures/Dis_DEGs_table.pdf",
       width = 6, height = 2)

# install packages for text mining
library(ggwordcloud)
#data
set.seed(42)
ggplot(data = freq.gene,
       aes( label = freq.gene,
            size = Freq,
            color =Freq)) +
  geom_text_wordcloud_area(perc_step = 0.2,
                           shape = "square",eccentricity = 1,
                           grid_size = 8, show.legend = T) +
  scale_size_area(max_size = 10)
 # scale_size_continuous(breaks = c(1,2,3),
 #              labels = c("1","2","3"))



#-------------------Disease Specific Gene Enrichment
library(gProfileR)
dis_list = list(ad_meta,asd_meta,scz_meta, bp_meta, mdd_meta, pd_meta, psp_meta, pa_meta)
names(dis_list) = c("AD", "ASD", "Scz", "BP", "MDD", "PD", "PSP", "PA")
GO_results= data.frame()
plot_dat = data.frame()

for (i in 1:length(dis_list)){
  dis_list[[i]]$gene_id = attr$ensembl_gene_id[match(rownames(dis_list[[i]]), attr$ensembl_gene_id_version)]
 }

for (i in names(dis_list)){
  query = dis_list[[i]]$gene_id
  go = gprofiler(query, 
                 organism="hsapiens", 
                 ordered_query = F, 
                 significant = T,
                 exclude_iea = F, 
                 region_query = F,
                 max_p_value = 1,
                 correction_method = "fdr", 
                 custom_bg = attr$ensembl_gene_id,
                 max_set_size = 2000,
                 hier_filtering = "moderate", 
                 domain_size = "annotated", 
                 numeric_ns = "", 
                 include_graph = F,
                 src_filter = c("GO", "KEGG"))
  
  go = go[order(go$p.value),]
  GO_results = rbind(GO_results, cbind(i, go))
  toPlot <- go %>% dplyr::select(term.name, p.value)
  toPlot <- toPlot[1:2,]
  plot_dat <- rbind(plot_dat, cbind(i,toPlot))
}
write.csv(GO_results, file="./results/tables/DiseaseSpecificGeneEnrichment.csv")


#---plot for each disease
plotGO = data.frame()
for(dx in names(dis_list)) {
  idx = which(GO_results$i==dx) 
  plotGO = rbind(GO_results[idx[10:1],],plotGO)  
}
plotGO = plotGO[!apply(is.na(plotGO),1,any),]
p <-ggplot(data =plotGO,
           aes(x = reorder(term.name, -p.value),
               y =-log10(p.value),
               fill = i))+
  facet_wrap(~i)+
  geom_bar(stat = "identity", width = 0.8, position = "dodge")+
  coord_flip()+
  geom_abline(slope = 0, intercept = -log10(0.01),lty = 2)+
  scale_color_brewer(palette = "Set1")+
  labs(x = "", 
       y = "-log10(FDR)",
       title = "Disease Specific Pathway Enrichment")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        text = element_text(size = 14),
        legend.position = "none")

p



ggsave(filename = "./results/figures/DiseaseEnrichment.pdf", 
       height = 10,
       width = 11,
       device = "pdf", plot = p)




goDx =ggplot(plotGO, 
             aes(y = reorder(term.name, -log10(p.value)),
                   x= i,
                   color = i))+
  geom_point(aes(size =-log10(p.value)))+
  scale_color_brewer(type = "div", palette = "Set1")+
  labs(x = "", 
       y = "",
       title = "Disease Specific Pathway Enrichment")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        text = element_text(size = 14),
        axis.text.x.bottom = element_text(angle = 45, vjust = 0.9, hjust = 0.8)
        )+
  guides(color = guide_legend(title = "Disease"))

goDx
ggsave(filename = "./results/figures/DiseasePathwayEnrichment.pdf", 
       height = 6,
       width = 8, 
       device = "pdf", 
       plot = goDx)

#--Top 2 pathways per disease
colnames(plot_dat) = c("Disease", "term.name", "p.value")
plot_dat$Disease= as.factor(plot_dat$Disease)
p <- ggplot(data = plot_dat,
            aes(x = reorder(Disease, -p.value),
                y =-log10(p.value),
                fill = term.name,
                group = term.name))+
  geom_bar(stat = "identity", width = 0.5, position = "dodge")+
  coord_flip()+
  geom_abline(slope = 0, intercept = -log10(0.01),lty = 2)+
  scale_fill_manual(values=sort(unique(plot_dat$Disease)))+
  labs(x = "", 
       y = "-log10(FDR)",
       title = "Top pathway enrichment for each disease")+
  guides(fill  = guide_legend(title = "Term name"))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        text = element_text(size = 14))
p

ggsave(filename = "./results/figures/DiseaseEnrichment_topGO.pdf", 
       height = 3,
       width = 5, 
       device = "pdf", 
       plot = p)


#save the data 
save.image("./codes/st3_overlap_lm.Rdata")
