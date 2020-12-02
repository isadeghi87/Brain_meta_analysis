#st5d_networkVisualization.R

rm(list=ls())
options(stringsAsFactors = F)
#source("http://bioconductor.org/biocLite.R"); biocLite("igraph")
library(WGCNA);library(ggplot2); library(reshape); 
library(igraph); library(RColorBrewer); 
library(WGCNA); library(corrplot);library(ggthemes)
library(dplyr)

#set working dir
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")

#load final network
load("./codes/finalNetwork.RData")

#### Make module eigengene-MDS plot ####
eigmat = MEs$eigengenes
colnames(eigmat) = gsub("E","",colnames(eigmat))
eigmat = eigmat[,-1]
adj = bicor(eigmat)
cols = cols[-1]
#### correlation between modules ####

pdf(file = "./results/figures/network_analysis/ModulesCorrelationPlot.pdf",
    height = 10, 
    width = 10)
corrplot(signif(adj,1),
         col = blueWhiteRed(n = 1000, gamma = 1)[1000:1],
         # main = "Gene modules pairwise correlation",
         mar = c(0,0,0,0),
         method = "ellipse",
         type = "upper",
         addCoef.col = T,
         order = "hclust",
         hclust.method = "complete",
         tl.cex = 1.5,
         tl.srt = 45,
         is.corr = T,
         diag = F,
         tl.col = "black",
         number.cex = 1.2)
dev.off()


#### MDS ####

mds = cmdscale(dist(t(eigmat)), eig = T);   
mds$points[,1] = -1* mds$points[,1]
g1 <- graph.adjacency(as.matrix(adj),
                      mode="undirected",
                      weighted=T,
                      diag=FALSE)
layoutFR <- mds$points
edgecolors = numbers2colors(E(g1)$weight,
                            colors = blueWhiteRed(100, gamma=1)[100:1],
                            signed=T, 
                            centered=T,
                            lim=c(-1,1))

## remove low weigth edges
pdf("./results/figures/network_analysis/ModuleNetworkGraph.pdf",
    height = 7,
    width = 9,
    useDingbats = FALSE)
par(mar=c(0,0,0,0))
plot.igraph(g1,
            vertex.label.dist=1.3,
            vertex.size=10,
            vertex.label.color="black",
            vertex.label.family = "sans",
            vertex.label.cex=0.9,
            vertex.color = cols,
            layout = layoutFR,
            edge.color = edgecolors,
            edge.width =2,
            asp=1, 
            # main="Module Eigengene MDS"
)

dev.off()

## Part 2) hub genes #### 
mod_gene = data.frame(module = mods,
                       row.names = rownames(datExpr))

cons_kme = signedKME(t(datExpr), MEs$eigengenes)
cons_kme = cons_kme[,-which(colnames(cons_kme)=="kME0")]
colnames(cons_kme) <- gsub("kME", "M", colnames(cons_kme))

maxsize = 10  #plot top 20 hub genes for each module
gene_idx = order(cons_kme[which(colors=="turquoise"),1], decreasing = T)[1:maxsize]

for(i in c(2:9)){
  gene_idx = c(gene_idx, order(cons_kme[,i], decreasing = T)[1:maxsize])
}

hubGenes = character()
hubGenes.kme = numeric()

#---hub for significant modules
for(col in colnames(cons_kme))  {
  modgenes = rownames(datExpr)[which((mods) == col)]
  kmes = cons_kme[modgenes, col]
  top_hubs = modgenes[order(kmes, decreasing=T)[1:maxsize]]
  top_hubs.kme = kmes[order(kmes, decreasing=T)[1:maxsize]]
  hubGenes = c(hubGenes,top_hubs)
  hubGenes.kme = c(hubGenes.kme, top_hubs.kme)
}
gene_idx = match(hubGenes,rownames(datExpr))

#### save hubgenes
saveRDS(object = hubGenes,"./results/tables/hubGenes.rds")

## calculate adjacency matrix 
adjMat = bicor(t(datExpr))
keepgenes = rownames(cons_kme)[gene_idx]

# keep hubgenes
adjMat = adjMat[gene_idx,gene_idx]

## keep edges with corr > 0.65
topcors=0.85^16
adjMat[adjMat< topcors] = 0

geneSymbols = datProbes$external_gene_name[match(keepgenes, datProbes$ensembl_gene_id)]
g1 <- graph.adjacency(as.matrix(adjMat),
                      mode="undirected",
                      weighted=T,
                      diag=FALSE)

mds = cmdscale(dist(t(adjMat)), eig = T)
layoutFR = mds$points
mdsedgecolors = numbers2colors(E(g1)$weight, 
                               colors = blueWhiteRed(100, gamma=1),
                               signed=T,
                               centered=T,
                               lim=c(-1,1))

# remove weigth edges < 0.7
g.copy <- delete.edges(g1, which(E(g1)$weight <0.7))

pdf("./results/figures/network_analysis/TopHubGenes.pdf", width=7, height=7)
par(mar=c(0,0,0,0))
set.seed(12)
plot.igraph(g.copy, 
            vertex.label = geneSymbols,
            vertex.label.dist=0.6, 
            edge.width=0.1,
            vertex.size=4, 
            vertex.frame.color="black",
            vertex.label.color="black",
            vertex.color = colors[gene_idx],
            vertex.label.cex=0.4, 
            layout=layout.fruchterman.reingold(g1),
            edge.color="grey")
title(main = "", cex.main = 1)

dev.off()

saveRDS(object = datProbes,"./results/tables/attributes.rds")
save.image("./codes/St5d_networkVis.Rdata")
