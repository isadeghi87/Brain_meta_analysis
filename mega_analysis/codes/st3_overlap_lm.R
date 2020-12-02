#Step3_Overlap check using lm data 
options(stringsAsFactors=F)
source("http://bioconductor.org/biocLite.R")
library(ggplot2); library(mada); library(reshape); library(pheatmap); library(RColorBrewer)
library(NMF); library(WGCNA); library(corrplot); library(purrr); 
library(dplyr); library(biomaRt);library(RRHO);  library(venn)
library(circlize);library(ComplexHeatmap); library(ggrepel); library(ggpubr)
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")


#### Load data ####

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

##pvalue matrix
corpval = matrix(NA,ncol = 8,nrow = 8)
colnames(corpval)=colnames(allmeta)
rownames(corpval) = colnames(allmeta)

#### Make a barplot ####
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
  corpval[x,y] = R$p.value;corpval[y,x]=R$p.value
  corpval[x,x] = 0; corpval[y,y]=0
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

#### correlation #### 
cordat <- cor(allmeta,
              use="pairwise.complete.obs",
              method="spearman")

write.csv(file = "./results/tables/Pairwise_cor_lm.csv", cordat)
pdf(file = "./results/figures/pairwise_correlation.pdf", 
    height = 4, width = 4)
corrplot(cordat,
         # title = "Spearman's correlation between diseases",
         mar = c(0,0,1,0),
         addCoef.col = T,
         diag = F,
         # order = "hclust", 
         # hclust.method = "complete", 
         method = "circle", 
         tl.col = "black", 
         type = "upper",
         number.cex = 0.9,
         tl.cex = 0.8
)
dev.off()

rplot = ggcorr(data = allmeta,cor_matrix = cordat,low = "red",
               mid = "white",high = "blue",
               midpoint = 0,geom = "circle",
               label = T,label_color = "black",min_size = 3,max_size = 12,
               name =expression(paste("correlation (",rho,")",sep = "")),
               label_round = 2,legend.size = 7,label_alpha = T)

ggsave(filename = "./results/figures/Pairwise_correlationGGplot.pdf",
       width = 5,height = 4)

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

#### circos plot ####
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

grid.echo()
f2a = grid.grab()
saveRDS(object = f2a,file = "./codes/Fig2a_scplot.rds")
ggsave(filename = "./results/figures/CrossDisease_scatterplotMat.pdf",
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
g.keep$layout = my_layout

#keep only stron correlations
source("./codes/ig2ggplot2.R")
pdf("./results/figures/diseases_mds2.pdf",height = 7, width = 10, useDingbats = FALSE)
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
grid.echo()
fig2c = grid.grab()
fig2c = as.ggplot(fig2c)

saveRDS(object = fig2c, file ="./codes/fig2c.rds")
# draw(lgd_list_vertical, x = unit(23, "cm"), y = unit(12, "cm"), just = c("top","right"))
dev.off()

#### Plot dendrogram fo the top Genes ####

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

#### Compile Table of Transcriptome Signatures for Each Disease####
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

#### Disease-specific DEGs ####

# list of plots
listPlots = vector("list",8)
names(listPlots) = unique(colnames(all_logfc))

#table of number of disease-specific DEGs
DEG.table = as.data.frame(matrix(NA, nrow = 8, ncol = 1))
rownames(DEG.table) =  c("AD","PD", "PA","PSP","Scz","ASD","BP","MDD")
colnames(DEG.table) = "DEGs"

#### frequency of significant genes ####
alldx = as.data.frame(matrix(NA,ncol = 6))
colnames(alldx) = c("logFC", "p.value","fdr", "disease","id" ,"gene")
freq.gene = vector()
subtext = alldx

for (i in names(all_disorders)){
  dat = as.data.frame(all_disorders[[i]])
  dat = na.omit(dat[,c("logFC", "p.value","fdr")])
  dat$disease = i #add disease name
  dat$id = rownames(dat)
  dat$gene = attr$external_gene_name[match(rownames(dat), attr$ensembl_gene_id_version)]
  rownames(dat) = NULL
  
  #check if duplicated
  id = !duplicated(dat)
  dat = dat[id,]
  alldx = rbind(alldx,dat)
  
  ###filter for significant genes (p<0.05)
  sig= subset(dat, p.value <0.05 & abs(logFC)>0.5 )
  DEG.table[i,1] = nrow(sig)
  ###order top DEGs 
  subdata = order(abs(sig$logFC), decreasing = T)
  subdata = subdata[1:10]
  sub = sig[subdata,]
  subtext = rbind(subtext,sub)
  freq.gene = c(na.omit(freq.gene),sub$gene)
}
alldx = na.omit(alldx)

#reorder based NDD and NPDs
alldx$disease = factor(alldx$disease, levels = c("AD", "PD", "PA", "PSP",
                                                 "Scz", "ASD", "BP", "MDD"))
subtext = subtext[-1,]
subtext$disease = factor(subtext$disease, levels = c("AD", "PD", "PA", "PSP",
                                                     "Scz", "ASD", "BP", "MDD"))
#color 
col = ifelse(-log10(alldx$p.value)>1.3 & abs(alldx$logFC)>0.5 , "green", "grey")
top.col = ifelse(subtext$logFC >0, "green", "red")

#### volcano plot ####
volc.p = ggplot(data = alldx , aes(x = logFC,
                                   y = -log10(p.value)))+
  geom_point(size = 0.5, color =col)+
  geom_text_repel(data = subtext,
                  aes(x = logFC,
                      y = -log10(p.value),
                      label = gene), 
                  color = "black",
                  # vjust = ,
                  size = 3)+
  geom_hline(yintercept = -log10(0.05), color = "red", lty = 2)+
  labs(y = "-log10(p-value)", x = expression(~log[2]~FC))+
  facet_wrap(~disease, ncol = 4, scales = "free_y")+
    theme_clean()+
  theme(text = element_text(size = 8 ,colour = "black"),
        strip.text.x = element_text(size = 14))

# save plot
ggsave(plot = volc.p, filename = "./results/figures/Dis_Specific_DEG.pdf.pdf",
       width = 12, height = 6)

saveRDS(p,file = "./codes/fig3a.rds")

write.table(DEG.table, "./results/tables/Number_DEGs.csv",col.names = T)

#plot the DEGs table
texTable =ggtexttable(t(DEG.table), theme = ttheme("mGreen"))
ggsave(plot = texTable,
       filename = "./results/figures/Dis_DEGs_table.pdf",
       width = 6, height = 2)


#top differentially expressed genes with using sum
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

# intersect of significant genes
sigData = subset(alldx,p.value < 0.05)

# check which genes are significant in all diseases
n = as.data.frame(table(sigData$gene))
n = n[order(n$Freq,decreasing = T),]## 

#now look at each etiology
neur = c("AD","PD","PA","PSP")
ps = c("Scz","ASD","BP","MDD")

sigData = sigData[,c("logFC","gene", "id","disease")]
sigData = sigData[order(abs(sigData$logFC),decreasing = T),]
for (d in unique(sigData$disease)){
  for (g in unique(sigData$gene)){
    sub = subset(sigData, disease == d & gene == g)
    if (nrow(sub) > 1){
      sigData[sigData$disease==d & sigData$gene ==g,][2,] = 0
      sigData = na.omit(sigData)
    }
  }
  
}

sigData = spread(sigData, key = disease, value = logFC)

# we need logFC for all genes 
logFCs= gene_table %>%  select(ensembl_gene_id_version, contains("logFC"))
colnames(logFCs)= gsub(".logFC","",colnames(logFCs))

neurDat = na.omit(sigData[,c("id","gene","AD", "PD", "PSP", "PA")])
psychDat = na.omit(sigData[,c("id","gene","ASD", "Scz", "MDD", "BP")])

unionGene = union(psychDat$gene,neurDat$gene)
unionId = union( psychDat$id, neurDat$id)

save(file = "./codes/DEGs.Rdata",unionGene,unionId, sigData,alldx)

dat = logFCs[match(unionId,logFCs$ensembl_gene_id_version),]
dat[is.na(dat)==T] = 0
dat = dat %>% select(ensembl_gene_id_version,AD,PD,PA,PSP,Scz,ASD,BP,MDD)

# display only top genes
top = as.character(n$Var1[n$Freq > 3])

#plot the heatmap

library(ComplexHeatmap)

pdf("./results/figures/GenesEtiology.pdf",width = 4.5,height =7)

right = rowAnnotation(gene = anno_mark(at =match(top,unionGene),
                                       labels = top, 
                                       labels_gp = gpar(fontsize=7)))

Heatmap(dat[,2:9],
        col = brewer.pal(11,"BrBG"),
        name = "logFC",
        cluster_rows = T,
        row_split = 2,
        clustering_distance_rows = "euclidean",
        cluster_columns = T,
        top_annotation = columnAnnotation(group = rep(c("NDD","NPD"),each=4),
                                          col =list(group = c("NDD" = "blue",
                                                              "NPD" = "black")),
                                          show_annotation_name=F),
        left_annotation = rowAnnotation(
          group = c(rep("NPD",29), rep("NDD",338)),
          col = list(group = c("NPD" = "black","NDD" = "blue" )),
          show_legend=F, show_annotation_name=F
        ),
        right_annotation = right,
        show_row_names = F,
        show_column_names =T,
        column_names_rot = 45,
        row_labels = unionGene,
        column_split = 2,
        border = T,
        height = unit(15,"cm") ,
        width = unit(3,"cm"),
        #heatmap_width = unit(15,"cm"),
        #heatmap_height = unit(5,"cm"),
        # column_names_rot = 90,
        row_names_gp = gpar(fontsize = 2),
        column_names_gp = gpar(fontsize = 8),
        show_row_dend = F,
        show_column_dend = F,
        show_heatmap_legend = T,
        heatmap_legend_param = gpar(fontsize=1)) 


dev.off()


#### Disease Specific Gene Enrichment ####
library(gProfileR)
dis_list = list(ad_meta,asd_meta,scz_meta, bp_meta, mdd_meta, pd_meta, psp_meta, pa_meta)
names(dis_list) = c("AD", "ASD", "Scz", "BP", "MDD", "PD", "PSP", "PA")

for (i in 1:length(dis_list)){
  dis_list[[i]]$gene_id = attr$ensembl_gene_id[match(rownames(dis_list[[i]]), attr$ensembl_gene_id_version)]
}
GO_results= data.frame()
plot_dat = data.frame()

for (i in names(dis_list)){
  sig = dis_list[[i]]
  sig = subset(sig,p.value <0.05)
  query = sig$gene_id
  go = gprofiler(query, 
                 organism="hsapiens", 
                 ordered_query = F, 
                 significant = T,
                 exclude_iea = F, 
                 region_query = F,
                 max_p_value = 1,
                 correction_method = "fdr",
                 # custom_bg = attr$ensembl_gene_id,
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


goDx =ggplot(plotGO, 
             aes(y = reorder(term.name, -log10(p.value)),
                 x= i,
                 color = i))+
  geom_point(aes(size =-log10(p.value)))+
  scale_color_brewer(type = "div", palette = "Set1")+
  labs(x = "", 
       y = "",
       title = "")+
  theme_minimal()+
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(angle = 45,
                                          vjust = 0.9,
                                          hjust = 0.8,
                                          size = 8),
        axis.text.y.left = element_text(size = 5)
  )+
  guides(color = "none",
         size = guide_legend(title = "-log10(FDR)"))


goDx
ggsave(filename = "./results/figures/DiseasePathwayEnrichment.pdf", 
       height = 6,
       width = 8, 
       device = "pdf", 
       plot = goDx)

#### Heatmap plot for GO ####
#refine some terms
plotGO$term.name = str_replace_all(plotGO$term.name,
                                   c("SRP-dependent cotranslational protein targeting to membrane" =
                                       "SRP-dependent targeting",
                                     "nuclear-transcribed mRNA catabolic process, nonsense-mediated decay"=
                                       "mRNA process",
                                     "protein modification by small protein conjugation"=
                                       "protein modification",
                                     "cell morphogenesis involved in neuron differentiation" =
                                       "neuron differentiation" 
                                   ))

matGo = plotGO %>% select(i,term.name,p.value) %>% 
  pivot_wider(names_from = i,values_from = p.value)
matGo = as.data.frame(matGo)
rownames(matGo)= matGo$term.name
matGo$term.name = NULL
matGo[is.na(matGo) == TRUE]=1
pheatmap(-log10(matGo),
         color = redblue(1000)[500:1000],
         show_rownames = T,
         show_colnames = T,
         angle_col = 45,
         fontsize_col = 5,
         cellwidth = 8,
         cellheight = 5,
         legend_breaks = c(0,10,20,30,40,max(-log10(matGo))),
         legend_labels = c("0","10","20","30","40","-log10(FDR)\n"),
         fontsize = 5,
         fontsize_row = 5,
         cluster_cols = T,
         treeheight_col = 1,
         cluster_rows = F,
         filename = "./results/figures/GO_heatmap_disease.pdf",
         width = 4,height = 5)

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

saveRDS(attr,file = "./codes/attributes.rds")
#save the data 
save.image("./codes/st3_overlap_lm.Rdata")
