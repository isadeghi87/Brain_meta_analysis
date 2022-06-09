#st5d_networkVisualization.R

rm(list=ls())
options(stringsAsFactors = F)
#source("http://bioconductor.org/biocLite.R"); biocLite("igraph")
library(WGCNA);library(ggplot2); library(reshape); 
library(igraph); library(RColorBrewer); library(WGCNA); library(corrplot)

#set working dir
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")

#load final network
load("./codes/finalNetwork.RData")

#### Make module eigengene-MDS plot ####
eigmat = MEs$eigengenes
colnames(eigmat) <- gsub("ME","", colnames(eigmat))
rownames(eigmat) <- gsub("ME", "", rownames(eigmat))
idx = which(colnames(eigmat) == "grey")
eigmat = eigmat[,-idx]
adj = bicor(eigmat)

#### correlation between modules ####
  pdf(file = "./results/figures/ModulesCorrelationPlot.pdf",
      height = 8, 
      width = 9)
  corrplot(signif(adj,1),
           col = blueWhiteRed(n = 1000, gamma = 1)[1000:1],
           # main = "Gene modules pairwise correlation",
           mar = c(0,0,3,0),
           method = "ellipse",
           type = "upper",
           addCoef.col = T,
           order = "hclust",
           hclust.method = "complete",
           tl.cex = 1.5,
           is.corr = T,
           diag = F,
           tl.col = "black",
           number.cex = 1.4)
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
  
  label <- rownames(mds$points)
  
  pdf("./results/figures/ModuleNetworkGraph.pdf",
      height = 8,
      width = 11,
      useDingbats = FALSE)
  plot.igraph(g1,
              vertex.label.dist=1.5, 
              vertex.size=10,
              vertex.label.color="black",
              vertex.label.family = "sans",
              vertex.label.cex=0.9,
              vertex.color = label,
              layout = layoutFR,
              edge.color = edgecolors,
              edge.width =2,
              asp=1, 
              # main="Module Eigengene MDS"
              )
  
  dev.off()
  
  #### correlation and MDS graph ####
  #-------------------------
  pdf(file = "./results/figures/ModulesCorrelationNetwork.pdf",
      height = 8, 
      width = 11)
  par(mfrow= c(1,2))
  
  corrplot(adj,
           col = blueWhiteRed(n = 100, gamma = 1)[100:1],
           main = "Gene modules pairwise correlation",
           mar = c(1,1,3,1),
           method = "color",
           type = "upper",
           addCoef.col = T,
           hclust.method = "complete",
           tl.cex = 0.9,
           tl.col = "black",
           number.cex = 0.5)
  plot.igraph(g1,
              vertex.label.dist=1.5, 
              vertex.size=10,
              vertex.label.color="black",
              vertex.label.family = "sans",
              vertex.label.cex=0.9,
              vertex.color = label,
              layout=layoutFR,
              edge.color=edgecolors,
              edge.width=2,asp=1, main="Module Eigengene MDS")
  
  
  dev.off()
  
  
  
  ## Part 2) 
  ## Make Module MDS plot
  ## --------------------
  moduleColors = colors
  modColors = data.frame(color= moduleColors,
                         row.names = rownames(datExpr))
  
  cons_kme = signedKME(t(datExpr), MEs$eigengenes)
  cons_kme = cons_kme[,-which(colnames(cons_kme)=="kMEgrey")]
  colnames(cons_kme) <- gsub("kME", "", colnames(cons_kme))
  
  maxsize = 20  #plot top 20 hub genes for each module
  gene_idx = order(cons_kme[which(colors=="black"),1], decreasing = T)[1:maxsize]
  
  for(i in c(2:10)){
    gene_idx = c(gene_idx, order(cons_kme[,i], decreasing = T)[1:maxsize])
  }
  
  hubGenes = character()
  hubGenes.kme = numeric()
  
  #---hub for significant modules
  for(col in colnames(cons_kme))  {
    modgenes = rownames(datExpr)[which((colors) == col)]
    kmes = cons_kme[modgenes, col]
    top_hubs = modgenes[order(kmes, decreasing=T)[1:maxsize]]
    top_hubs.kme = kmes[order(kmes, decreasing=T)[1:maxsize]]
    hubGenes = c(hubGenes,top_hubs)
    hubGenes.kme = c(hubGenes.kme, top_hubs.kme)
  }
  gene_idx = match(hubGenes,rownames(datExpr))
  
  adjMat = bicor(t(datExpr))
  keepgenes = rownames(cons_kme)[gene_idx]
  adjMat = adjMat[gene_idx,gene_idx]
  topcors=0.65^9
  adjMat[adjMat< topcors]=0
  
  geneSymbols = datProbes$external_gene_name[match(keepgenes, datProbes$ensembl_gene_id_version)]
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
  
  
  g.copy <- delete.edges(g1, which(E(g1)$weight <0.7))
  
  pdf("./results/figures/TopHubGenes.pdf", width=5, height=5)
  par(mar=c(0,0,1,0))
  set.seed(12)
  plot.igraph(g.copy, 
              vertex.label = geneSymbols,
              vertex.label.dist=0.6, 
              edge.width=0.1,
              vertex.size=5, 
              vertex.frame.color="black",
              vertex.label.color="black",
              vertex.color = colors[gene_idx],
              vertex.label.cex=0.25, 
              layout=layout.fruchterman.reingold(g1),
              edge.color="grey")
  title(main = "", cex.main = 1)
  
  dev.off()
  
  
saveRDS(hubGenes,file = "./codes/hubgenes.rds")

#### cell types modules ####

  group = colors[gene_idx]
  cellGroups = vector("list", length = 8)
  names(cellGroups) = c("Oligo", "Astro", "Bergmann", "Schwann", "Neuron", "Pyramidal", "GABAergic",
                        "NPCs")
  for (i in c("Oligo", "Schwann")){
    cellGroups[[i]] = which(group =="brown")
  }
  
  for (i in c("Astro", "Bergmann", "NPCs")){
    cellGroups[[i]] = which(group =="blue")
  }
  
  cellGroups[["Neuron"]] = which(group %in% c("magenta", "turquoise") )
  cellGroups[["Pyramidal"]] = which(group =="magenta")
  cellGroups[["GABAergic"]] = which(group =="turquoise")
  
  layoutfr = layout_with_fr(g1)
  

  save.image("./codes/St5d_networkVis.Rdata")
  