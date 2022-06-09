#st5e_GeneEnrichment.R

rm(list=ls()); options(stringsAsFactors = F)
library(pSI); library(gProfileR); library(gplots); library(biomaRt); library(WGCNA)
library(ggpubr); library(pheatmap); library(dplyr); library(magrittr)
library(GeneOverlap); library(factoextra);library(FactoMineR)

setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")
#---load final network
load("./codes/finalNetwork.RData")

eigmat = MEs$eigengenes; colnames(eigmat) = gsub("ME","",colnames(eigmat))
kME = signedKME(t(datExpr), MEs$eigengenes); colnames(kME) = gsub("kME", "", colnames(kME))

#rename normal to control
datMeta$Dx = gsub("Normal","Control",datMeta$Dx)
datMeta$Dx = factor(datMeta$Dx, levels = c("Control", "AD", "PD",
                                           "PA","PSP",
                                           "Scz","ASD","BP","MDD"))
#Write module genes for GO elite
for(c in unique(colors)) {
  mod = datProbes$ensembl_gene_id_version[colors==c]
  write.table(file = paste("./results/tables/GeneSet/module_", c, ".txt", sep=""),
              data.frame(Gene=mod, Id="En"), row.names = F, quote=F)
}
write.table(file="./results/tables/GeneSet/background.txt", 
            data.frame(Gene=datProbes$ensembl_gene_id_version, Id="En"), row.names=F, quote=F)



## Calculate GO enrichment of Disease-associated modules
rownames(datExpr) <- datProbes$ensembl_gene_id[match(rownames(datExpr), datProbes$ensembl_gene_id_version)]
resultsGO = data.frame()
ModulePlots <- list()
plotGO = data.frame()
for (m in unique(colors)) {
  me_name = paste("ME", m, sep="")
  i = which(m == unique(colors))  
  
  #-----------------------------------------GO Enrichment
  query = rownames(datExpr)[colors==m]
  go = gprofiler(query, 
                 organism="hsapiens", 
                 ordered_query = F, 
                 significant = T,
                 exclude_iea = F, 
                 region_query = F,
                 max_p_value = 1,
                 correction_method = "fdr", 
                # custom_bg = rownames(datExpr),
                 max_set_size = 2000,
                 hier_filtering = "moderate", 
                 domain_size = "annotated", 
                 numeric_ns = "", 
                 include_graph = F,
                 src_filter = c("GO", "KEGG"))
  go = go[order(go$p.value),]
  plotdat <- data.frame(term.name=go$term.name[1:5], p.value=go$p.value[1:5])
  plotdat <- na.omit(plotdat)
  plotdat$Module = m
  plotGO= rbind(plotGO, plotdat)
 #--plot
  p <- ggplot(data = plotdat,
         aes(x = reorder(term.name, -p.value),
                            y =-log10(p.value)))+
    geom_bar(stat = "identity", fill= m)+
    coord_flip()+
    geom_abline(slope = 0, intercept = -log10(0.01),lty = 2)+
    labs(x = "", 
         y = "-log10(FDR)",
         title = m)+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          text = element_text(size = 14))
 ModulePlots[[m]]= p
  resultsGO = rbind(resultsGO, cbind(m, go))
}

ModulePlots[["grey"]]=NULL
plotGO = plotGO[!grepl("grey",plotGO$Module),]
GOplots <- ggarrange(plotlist = ModulePlots, 
                     ncol = 2,
                     nrow = 7)
GOplots
ggsave(plot = GOplots, 
       filename = "./results/figures/ModulesGO.pdf", 
       height = 14, width =24, device = "pdf")

allGO <- ggplot(data = plotGO,
       aes(x = reorder(term.name, -p.value),
           y =-log10(p.value),
           group = Module,
           fill = Module))+
  geom_bar(stat = "identity", width = 0.5)+
  coord_flip()+
  facet_wrap(~Module,nrow = 1 )+
  geom_abline(slope = 0, intercept = -log10(0.01),lty = 2)+
  scale_fill_manual(values=sort(unique(plotGO$Module)))+
  labs(x = "", 
       y = "-log10(FDR)",
       title = "Pathway Enrichment for the modules")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        text = element_text(size = 14),
        legend.position = "none")

#save plot and data
ggsave(plot = allGO,
       filename = "./results/figures/ModulesGO2.pdf",
       height = 14, width =24, device = "pdf")
write.csv(file="./results/tables/GeneSet/GO_gprofilerResults.csv", resultsGO)

go_2d =ggplot(plotGO, 
       aes(y = term.name,
           x= Module,
           group =Module))+
  geom_point(aes(size =-log10(p.value)),color = plotGO$Module)+
  scale_color_manual(values = plotGO$Module)+
  labs(x = "", 
       y = "",
       title = "")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        text = element_text(size = 14),
        axis.text.x.bottom = element_text(angle = 45, vjust = 0.9, hjust = 0.8)
  )+
  guides(size = guide_legend(title = "-log10(FDR)"))
ggsave(filename = "./results/figures/GO_modules_2D.pdf",
       width = 8.5,height = 8,plot = go_2d)

## plot GO as heatmap
plotGO$term.name = str_replace_all(plotGO$term.name,
                                   c("SRP-dependent cotranslational protein targeting to membrane" =
                                       "SRP-dependent targeting",
                                     "cotranslational protein targeting to membrane"=
                                       "protein targenting",
                                     "regulation of intracellular signal transduction"=
                                       "regulation of signaling"))
library(tidyr)
matGo = plotGO  %>% 
  pivot_wider(names_from = Module,values_from = p.value)
matGo = as.data.frame(matGo)
rownames(matGo)= matGo$term.name
matGo$term.name = NULL
matGo[is.na(matGo) == TRUE]=1
annot_mod = data.frame(module = colnames(matGo))
rownames(annot_mod)= colnames(matGo)
annot_col = list(module = c(yellow = "yellow",
                            pink = "pink",
                            brown = "brown",
                            blue = "blue",
                            purple = "purple",
                            magenta = "magenta",
                            red = "red",
                            black = "black",
                            turquoise = "turquoise",
                            green="green"))
pheatmap(-log10(matGo),
         color = redblue(1000)[500:0],
         show_rownames = T,
         show_colnames = T,
         angle_col = 45,
         fontsize_col = 5,
         cellwidth = 8,
         cellheight = 5,
         legend_breaks = c(0,20,40,60,max(-log10(matGo))),
         legend_labels = c("0","20","40","60","-log10(FDR)\n"),
         fontsize = 5,
         fontsize_row = 5,
         cluster_cols = F,
         treeheight_col = 1,
         cluster_rows = T,
         treeheight_row = 1,
         annotation_col = annot_mod,
         annotation_colors = annot_col,
         annotation_legend = F,
         annotation_names_col = F,
         filename = "./results/figures/GO_heatmap_modules.pdf",
         width = 3.5,height = 4)
#-----------------------Pick top 2 pathways per module
resultsGOsub2 = data.frame()
for(c in names(ModulePlots)) {
  idx = which(resultsGO$m==c) 
  resultsGOsub2 = rbind(resultsGO[idx[2:1],],resultsGOsub2)  
}
resultsGOsub2 = resultsGOsub2[!grepl("grey",resultsGOsub2$m),]
resultsGOsub2$m = factor(resultsGOsub2$m, levels = unique(resultsGOsub2$m))


pdf("./results/figures/PathwayEnrichment.pdf",width = 8, height=6)
par(oma=c(0,5,0,3), mar=c(5,7,1,0))
bp = barplot(-log10(resultsGOsub2$p.value), 
             main="Top Pathway Enrichments of Modules",
             horiz=T, yaxt='n', 
             col=as.character(resultsGOsub2$m),
             xlab='-log10(FDR)',cex.main=1, cex.axis = .7, border=NA)
axis(2,at=bp,labels=resultsGOsub2$term.name,tick=FALSE,las=2,cex.axis=0.6);
abline(v=-log10(0.01), col="black", lwd=2,lty=2)
dev.off()

#plot as 2D

resultsGOsub2$term.name = factor(resultsGOsub2$term.name, levels = unique(resultsGOsub2$term.name))
goMod =ggplot(resultsGOsub2, 
              aes(y = term.name,
                  x= m,
                  group =m))+
  geom_point(aes(size =-log10(p.value)), color = resultsGOsub2$m)+
  scale_color_manual(values = resultsGOsub2$m)+
  labs(x = "", 
       y = "",
       title = "Modules Top Pathway Enrichments")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        text = element_text(size = 14),
        axis.text.x.bottom = element_text(angle = 45, vjust = 0.9, hjust = 0.8)
  )+
  guides(color = guide_legend(title = "Module"))

goMod
ggsave(filename = "./results/figures/ModulesPathwayEnrichment.pdf", 
       height = 6,
       width = 10, 
       device = "pdf", 
       plot = goMod)
# correspondence analysis for modules and pathways
library("FactoMineR");library("factoextra")

#CA 
library(tidyr)
resultsGOsub = data.frame()

for(c in names(ModulePlots)) {
  idx = which(resultsGO$m==c) 
  resultsGOsub = rbind(resultsGO[idx[3:1],],resultsGOsub)  
}
resultsGOsub = resultsGOsub[!grepl("grey",resultsGOsub$m),]
resultsGOsub$m = factor(resultsGOsub$m, levels = unique(resultsGOsub$m))
corres.dat = resultsGOsub[,c(1,4,13)]
corresDat <- spread(data = corres.dat,key = m, value = p.value)
rownames(corresDat) = corresDat$term.name
corresDat = corresDat[,-1]
corresDat[is.na(corresDat)] =1
corresDat = -log10(corresDat)
ca <- CA(corresDat,  ncp = 5, 
                graph = TRUE)
factoextra::fviz_ca_biplot(ca,
title =  "correspondence analysis of modules pathway enrichment",
                            # select.row = list(contrib =12),
                           # alpha.row = "contrib",
                           # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                           repel = T)
#
datrow = as.data.frame(ca$row$coord)
datrow$label = rownames(datrow)
datcol = as.data.frame(ca$col$coord)
datcol$label = rownames(datcol)

#save a plot as ggplot
library(ggrepel)
ca_plot = ggplot()+
  geom_point(data = datcol,
             aes(x = datcol[,1], y = datcol[,2]),color = "red")+
  geom_label_repel(data = datcol, aes(label = label, 
                                      x = datcol[,1], 
                                      datcol[,2]), 
                   color = unique(datcol$label),
                   size =3,
                   nudge_x= 0.5,
                   nudge_y = 1.3,
                   # hjust= 1 ,
                   segment.colour = "black")+
  geom_point(data = datrow,
             aes(x = datrow[,1], y = datrow[,2]),color = "blue")+
  geom_label_repel(data = datrow, aes(label = label, 
                                      x = datrow[,1],
                                      y = datrow[,2]), 
                   color = "blue",  size = 1.75)+
  labs(x = "Dim1", y= "Dim2", title = "Correspondense analysis of GO for modules")+
  scale_x_continuous(limits = c(-3,3))+
  scale_y_continuous(limits = c(-3,3))+
  geom_hline(yintercept = 0,lty=2, size = 0.5)+
  geom_vline(xintercept = 0, lty=2, size = 0.5)+
  theme_pubr()+
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
ca_plot
ggsave(filename = "./results/figures/Corresp_Analysis.pdf", 
       plot = ca_plot, width = 9, height = 5)
  


##Download human transcriptome data from:
##Zhang, Y. et al. Purification and Characterization of Progenitor and Mature Human Astrocytes Reveals Transcriptional and Functional Differences with Mouse. Neuron 89, 37–53 (2016).
##Supplemental Table 3, "Human data only" sheet

  zhang.datExpr = read.csv("./results/tables/Zhang_mmc3datExpr.csv",skip=3,nrow=23223,head=F) 
  zhang.datMeta = data.frame(row.names=1:41,t(read.csv("./results/tables/Zhang_mmc3datExpr.csv",nrow=3,head=F)[,-1]))
  zhang.datMeta$CellType = NA
  zhang.datProbes = data.frame(symbol=zhang.datExpr$V1)
  zhang.datExpr = zhang.datExpr[,-1]
  colnames(zhang.datMeta) = c("X1", "Age", "Gender","CellType")
  zhang.datMeta$CellType[15:26] = "Astrocyte"
  zhang.datMeta$CellType[27] = "Neuron"
  zhang.datMeta$CellType[28:32] = "Oligo"
  zhang.datMeta$CellType[33:35] = "Microglia"
  zhang.datMeta$CellType[36:37] = "Endothelial"
  zhang.datMeta$CellType[38:41] = "WholeCortex"
  
  zhang.datExpr2 = data.frame(matrix(NA, 
                                     nrow=nrow(zhang.datExpr), 
                                     ncol=5))
  colnames(zhang.datExpr2)=  c("Neuron", "Astrocyte", "Oligo", "Microglia","Endothelial")
  zhang.datExpr2$Neuron = zhang.datExpr[,which(zhang.datMeta$CellType=="Neuron")]
  for(cell in colnames(zhang.datExpr2)[2:5]) {
    zhang.datExpr2[,cell] = apply(zhang.datExpr[,which(zhang.datMeta$CellType==cell)],1,mean)  
  }
  
  #----Annotating gene IDs
  getinfo <-
    c("ensembl_gene_id_version",
      "ensembl_gene_id",
      "external_gene_name")
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset = "hsapiens_gene_ensembl",
                  host = "apr2018.archive.ensembl.org") ## Gencode v28
   bm <- getBM(
    attributes = getinfo,
    filters = "hgnc_symbol",
    values = zhang.datProbes$symbol,
    mart = mart)
 
   ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", 
                     host="feb2014.archive.ensembl.org") 
   bm = getBM(attributes = c("ensembl_gene_id", "external_gene_id"),
              filters = "hgnc_symbol", values = zhang.datProbes$symbol, mart=ensembl)
   zhang.datProbes = data.frame(symbol=zhang.datProbes$symbol, 
                                ensg=bm$ensembl_gene_id[match(zhang.datProbes$symbol, bm$external_gene_id)])
   cr =collapseRows(zhang.datExpr2, 
                     rowGroup = zhang.datProbes$ensg,
                   rowID=1:nrow(zhang.datExpr2))
  zhang.datExpr = cr$datETcollapsed
  rm(zhang.datExpr2, zhang.datProbes)


# Calculate Cell-Type specificity of modules
#Zhang using pSI
set.seed(100)
pSI.output = specificity.index(pSI.in=zhang.datExpr,
                               bts=100,
                               p_max=.1, 
                               e_min=0.3); 
ps.count <- pSI.count(pSI.output)
write.table(ps.count, file = "./results/tables/pSI_specificGeneCount.csv")
cell.p.zhang = matrix(NA, 11,5)
rownames(cell.p.zhang) = unique(colors)
colnames(cell.p.zhang) = colnames(pSI.output)


for(mod in rownames(cell.p.zhang)) {
  f = fisher.iteration(pSI.output, rownames(datExpr)[colors==mod],p.adjust = F)
  cell.p.zhang[mod,] = f$`0.05 - nominal`
}

#remove grey module
cell.p.zhang = cell.p.zhang[-3,]
cell.p.zhang.fdr = p.adjust(cell.p.zhang,"fdr")
dim(cell.p.zhang.fdr) = dim(cell.p.zhang); dimnames(cell.p.zhang.fdr) = dimnames(cell.p.zhang);

dendro.col = as.dendrogram(hclust(as.dist(1-bicor(zhang.datExpr)), method="average"))
dendro.row= as.dendrogram(hclust(
  as.dist(1-bicor(
    eigmat[,-7])),method="average"))

# #---plot heatmap 

# pdf("./results/figures/CellTypesEnrichmentModules.pdf", 
#     width = 10, height = 10)
# heatmap.2(-log10(cell.p.zhang.fdr),
#           col = redWhiteGreen(1000,1)[500:1000],
#           scale = "none",
#           trace = "none",
#           cexRow = 1.3, 
#           cexCol = 1.3, 
#           density.info = "none",
#           # colsep = 0:7,
#           # rowsep = 0:10,
#           sepcolor ="blue",
#           sepwidth=c(0.02,0.02),
#           srtCol = 45,
#           offsetRow = 0,
#           offsetCol = -0.5,
#           Rowv = dendro.row,
#           Colv = dendro.col,
#           key = T, 
#           key.xlab = "-log10(FDR)",
#           cellnote = signif(cell.p.zhang.fdr,1),
#           notecex=1, 
#           notecol = "black",
#           RowSideColors = rownames(cell.p.zhang.fdr),
#           main="Cell-type specific enrichment")
# dev.off()
pheatmap(-log10(cell.p.zhang.fdr),
         # color = brewer.pal(9,"Greens")[0:8],
         color = redWhiteGreen(1000,1)[500:1000],
         cluster_cols = T,
         cluster_rows = T,
        fontsize_col = 10,
         fontsize_number = 10,
         fontsize_row = 10,
         treeheight_row = 5,
          cellwidth = 40,
        cellheight = 30,
         angle_col = 45,
         show_rownames = T,
         show_colnames = T,
         # main = "Cell type-specific gene modules",
         legend_breaks = c(0,50,100,150, max(-log10(cell.p.zhang.fdr))),
         legend_labels = c(0,50,100,150,"-log10 (FDR)\n\n"),
         treeheight_col = 10, 
         display_numbers = signif(cell.p.zhang.fdr,1),
         filename = "./results/figures/CellTypesEnrichmentModules.pdf",
         height = 6,
         width = 6
)



out = cell.p.zhang.fdr
colnames(out) = paste(colnames(out), ".FDR",sep="")
write.csv(out,file="./results/tables/CellTypeEnrichmentFDR.csv")

#check cell types with PanglaoDB database
all_genes = data.frame(row.names = rownames(datExpr), color = colors )
PangDB <- read.csv("~/github/PanglaoDB-master/PanglaoDB_markers_11_Sep_2019.csv", sep ='\t')
PangDB <- na.omit(PangDB)
PangDB <- PangDB[PangDB$organ == "Brain" & PangDB$species %in% c("Hs", "Mm Hs"),]
hPang <- PangDB[PangDB$species == "Hs",]
PangDB <- PangDB[PangDB$cell.type != "Immature neurons",]
PangDB$cell.type <- factor(PangDB$cell.type)

cellGene  = list()
for (cell in levels(PangDB$cell.type)){
cellGene[[cell]] = PangDB[PangDB$cell.type == cell,] %>% dplyr::select(official.gene.symbol)
}
names(cellGene) = levels(PangDB$cell.type)

ModuleGenes <- as.data.frame(all_genes)
ModuleGenes$symbol = datProbes$external_gene_name[match(rownames(ModuleGenes),
                                                                 datProbes$ensembl_gene_id)]
colnames(ModuleGenes) = c("module", "symbol")
ModuleGenes$module = as.factor(ModuleGenes$module)
ModuleGenelist = list()

for (mod in levels(ModuleGenes$module)){
  ModuleGenelist[[mod]] = ModuleGenes[ModuleGenes$module == mod,]
  }

cell_mod = as.data.frame(matrix(NA, ncol = length(ModuleGenelist), nrow = length(cellGene)))
rownames(cell_mod) = names(cellGene);colnames(cell_mod) = names(ModuleGenelist)
cell_mod_ODD = cell_mod
for (mod in 1:length(ModuleGenelist)){
  for ( cell in 1: length(cellGene)){
    go.obj = newGeneOverlap(listA = ModuleGenelist[[mod]]$symbol,
                        listB = cellGene[[cell]]$official.gene.symbol,
                        genome.size = nrow(all_genes))
    go.obj = testGeneOverlap(go.obj)
  cell_mod[cell,mod]= getPval(go.obj)
  cell_mod_ODD[cell,mod] = getOddsRatio(go.obj)
  }
}

# remove grey
cell_mod = cell_mod[,-5]
cell_mod_ODD = cell_mod_ODD[,-5]

mat.plot = as.matrix(-sign(cell_mod_ODD) * log10(cell_mod))

#annotation colors
annot_mod = data.frame(module = colnames(cell_mod))
rownames(annot_mod)= colnames(cell_mod)
annot_col = list(module = c(yellow = "yellow",
                            pink = "pink",
                            brown = "brown",
                            blue = "blue",
                            purple = "purple",
                            magenta = "magenta",
                            red = "red",
                            black = "black",
                            turquoise = "turquoise",
                            green="green"))
pheatmap(-log10(cell_mod),
         color = redWhiteGreen(2000,1)[1000:2000],
         cluster_cols = T,
         cluster_rows = T,
         annotation_col = annot_mod,
         annotation_colors = annot_col,
         annotation_legend = F,
         annotation_names_col = F,
         fontsize_col = 4,
         fontsize = 4,
         fontsize_number = 3,
         fontsize_row = 4,
         treeheight_col = 1, 
         treeheight_row = 1,
         cellwidth = 10,
         cellheight = 5,
         angle_col = 45,
         show_rownames = T,
         show_colnames = T,
         legend_breaks = c(0,5,10,15, max(-log10(cell_mod))),
         legend_labels = c(0,5,10,15,"-log10 (FDR)\n"),
         display_numbers = signif(cell_mod,1),
         filename = "./results/figures/CellModule_PangDB.pdf",
         height = 3,
         width = 3.1,
)


## Calculate TFBS enrichment of Disease-associated modules
resultsTFBS = data.frame()
for (m in unique(colors)) {
 
  ## GO Enrichment
  query = rownames(datExpr)[colors==m]
  go = gprofiler(query, organism="hsapiens", 
                 ordered_query = F,
                 significant = T,
                 exclude_iea = F,
                 region_query = F,
                 max_p_value = 1,
                 correction_method = "fdr",
                 custom_bg = rownames(datExpr),
                 hier_filtering = "moderate",
                 domain_size = "annotated",
                 numeric_ns = "",
                 include_graph = F,
                 src_filter = c("TF"))
  go = go[order(go$p.value),]
  
  resultsTFBS = rbind(resultsTFBS, cbind(m, go[1:min(nrow(go),20),]))
}

write.csv(resultsTFBS,file="./results/tables/TranscriptionFactors.csv")

#plot the top TFBs
topTF = data.frame()
for(c in unique(resultsTFBS$m)) {
  idx = which(resultsTFBS$m == c) 
  print(idx)
  topTF = rbind(resultsTFBS[idx[2:1],],topTF)  
}
topTF = na.omit(topTF)
library(stringr)
labels = str_extract_all(topTF$term.name, pattern = "Factor.*[A-Z]")

pdf("./results/figures/TFBsPlot.pdf",width = 5, height=4)
par(oma=c(1,5,0,1), mar=c(4,4,0,0))
bp = barplot(-log10(topTF$p.value), 
             # main="Top TF binding sites for modules",
             horiz=T, yaxt='n', 
             col=as.character(topTF$m),cex.lab=0.5,
             xlab='-log10(FDR)', cex.axis = .4, border=NA)
axis(2,at=bp,labels=labels,tick=FALSE,las=2,cex.axis=0.4);
abline(v=-log10(0.05), col="black", lwd=2,lty=2)
dev.off()

## Identify hub genes transcription factors
getinfo <-
  c("ensembl_gene_id_version",
    "ensembl_gene_id",
    "external_gene_name")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
Probes <- getBM(
  attributes = getinfo,
  filters = "ensembl_gene_id_version",
  values = rownames(kME),
  mart = mart)
rownames(Probes) <- Probes$ensembl_gene_id

rownames(kME) = rownames(datExpr)
hubGenes = data.frame(Module=NA, Gene=NA, Symbol=NA,Rank=NA)
for (m in unique(colors)) {
  mod = rownames(datExpr)[colors==m]
  mod = mod[order(kME[mod,m],decreasing = T)[1:20]]
  sym = Probes[mod,"external_gene_name"]
  hubGenes = rbind(hubGenes, data.frame(Module=m, Gene=mod, Symbol = sym, Rank=1:20))
}
hubGenes = hubGenes[-1,]



a=listAttributes(mart);f=listFilters(mart)
bm = getBM(attributes = c("ensembl_gene_id_version",
                          "ensembl_gene_id",
                          "external_gene_name",
                          "go_id",
                          "go_linkage_type",
                          "goslim_goa_accession",
                          "goslim_goa_description"), 
           filters = "ensembl_gene_id", 
           values = hubGenes$Gene, 
           mart=mart)
bm = bm[grep("transcription factor",bm$goslim_goa_description),]
hubGenes$TF = bm$goslim_goa_description[match(hubGenes$Gene, bm$ensembl_gene_id)]

write.csv(hubGenes,file="./results/tables/HubGeneTFs.csv", row.names = F)


## eRNA overlap with co-expression modules
## Download eRNA modules from http://www.nature.com/neuro/journal/v18/n8/extref/nn.4063-S12.xlsx
eRNAnetwork=  read.csv("./results/tables/eRNA_Yao.etal.csv")
source("./codes/Fisher_overlap.R")
table.p = matrix(NA, nrow=length(unique(colors)), ncol=19)
rownames(table.p) = unique(colors); 
colnames(table.p) = paste("M", 1:19,sep="")
table.or = table.p
hgnc = Probes$external_gene_name
for (m in unique(colors)) {
  for(e in colnames(table.p)) {
    f = ORA(hgnc[colors==m],eRNAnetwork$Symbol[eRNAnetwork$Module.Label==e], hgnc, eRNAnetwork$Symbol)
    table.or[m,e] = as.numeric(f[[1]])
    table.p[m,e] = as.numeric(f[[2]])
  }
}

table.p.fdr = p.adjust(table.p,"fdr")
dim(table.p.fdr) = dim(table.p); dimnames(table.p.fdr) = dimnames(table.p)
table.p.fdr[table.or<1] = 1

#remove grey module
to_plot = table.p.fdr[-3,]
colnames(to_plot) <- paste("eModule",1:19,sep = "")


pdf("./results/figures/BrainEnhancerEnrichment.pdf",width = 14, height = 12 )
par(margin(2,2,2,2))

pheatmap(-log10(to_plot),
         color = blueWhiteRed(n = 2000)[1000:0],
         cluster_cols = F,
         cluster_rows = T,
         annotation_row = annot_mod,
         annotation_colors = annot_col,
         annotation_legend = F,
         annotation_names_row = F,
         fontsize = 4,
         fontsize_number = 5,
         fontsize_row = 5,
         fontsize_col = 4,
         treeheight_row = 1,
         cellwidth = 16,
         cellheight = 12,
         angle_col = 45,
         show_rownames = T,
         show_colnames = T,
         # main = "Enrichment of Brain Enhancers Modules",
         legend_breaks = c(0,50,100,150,200,250, max(-log10(to_plot))),
         legend_labels = c(0,50,100,150,200,250,"-log10 (FDR)\n"),
         treeheight_col = 1, 
         display_numbers = signif(to_plot,1),
         filename = "./results/figures/BrainEnhancerEnrichment.pdf",
         height = 2.2,
         width = 5.5
)



#Calculate overlap with Winden modules
#Download winden modules from http://msb.embopress.org/content/5/1/291.long#sec-23
winden= read.csv("./results/tables/winden.csv")
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="feb2014.archive.ensembl.org")
a=listAttributes(ensembl);f=listFilters(ensembl)
bm1= getBM(attributes = c("affy_moe430a","ensembl_gene_id"),
           filters = "affy_moe430a", values = winden$Affymetrix.ID, mart=ensembl)
bm2 = getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"), mart=ensembl)
bm = merge(bm1,bm2)
idx = match(winden$Affymetrix.ID, bm$affy_moe430a)
winden$ensg= bm$hsapiens_homolog_ensembl_gene[idx]

table.p = matrix(NA, nrow=length(unique(colors)), ncol=length(table(winden$Module.Assigned)))
rownames(table.p) = unique(colors); 
nm <-na.omit(unique(winden$Module.Assigned))[1:17]
colnames(table.p) = unique(nm)
table.or = table.p

for (m in unique(colors)) {
  for(e in colnames(table.p)) {
    f = ORA(datProbes$ensembl_gene_id[colors==m],winden$ensg[winden$Module.Assigned==e],
            datProbes$ensembl_gene_id, winden$ensg)
    table.or[m,e] = as.numeric(f[[1]])
    table.p[m,e] = as.numeric(f[[2]])
  }
}
table.p.fdr = p.adjust(table.p,"fdr")
dim(table.p.fdr) = dim(table.p); dimnames(table.p.fdr) = dimnames(table.p)
table.p.fdr[table.or<1] = 1

to_plot = data.frame( Module = rownames(table.p.fdr),
                      Mitochondria="Synaptic",
                      Enrichment = (table.or[,"turquoise"]),
                      FDR = table.p.fdr[,"turquoise"])
to_plot=rbind(to_plot,
              data.frame(Module=rownames(table.p.fdr), 
                         Mitochondria="Non-Synaptic",
                         Enrichment=(table.or[,"blue"]),
                         FDR = table.p.fdr[,"blue"]))
to_plot = to_plot[to_plot$Module != "grey",]
to_plot$Module = factor(to_plot$Module)
to_plot$Mitochondria = factor(to_plot$Mitochondria)
to_plot$symbol=""
to_plot$symbol[to_plot$FDR<0.05]="*"
to_plot$symbol[to_plot$FDR<0.01]="**"
to_plot$symbol[to_plot$FDR<0.001]="***"


mitoPlot = ggplot(to_plot, 
              aes(x=Module, 
                  y=Enrichment,
                  fill=Mitochondria,
                  alpha=FDR<0.05)) +
  coord_flip() +
  geom_bar(stat="identity",
           position=position_dodge(),
           width = 0.5)+ 
  geom_hline(yintercept = 1,lty=2) +
  ylab("Enrichment") + 
  xlab("Module") +
  scale_fill_manual(values = c("springgreen","red"))+
  # ggtitle("Mitochondrial Module Enrichment") + 
  theme(plot.title = element_text(hjust=0.5, face = "bold"))+
  theme_bw()
mitoPlot


ggsave(mitoPlot, file="./results/figures/Module_mitochondria.pdf",
       width = 6, height = 3) 
write.csv(to_plot,file="./results/tables/TableS2-Winden.csv",row.names = F)


#---motifs for TFs
# library(RcisTarget)
# 
# genelist <- datProbes$external_gene_name
# geneLists <- list(geneSetName=genelist)
# featherURL <- "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather" 
# download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
# # Search space: 10k bp around TSS - HUMAN
# motifRankings <- importRankings("hg19-tss-centered-10kb-7species.mc9nr.feather")
# # Load the annotation to human transcription factors
# data(motifAnnotations_hgnc)
# 
# motifEnrichmentTable_wGenes <- cisTarget(geneLists, 
#                                          motifRankings,
#                                          motifAnnot=motifAnnotations_hgnc)
# motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
# library(DT)
# datatable(motifEnrichmentTable_wGenes_wLogo[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
#           escape = FALSE, # To show the logo
#           filter="top", options=list(pageLength=5))
# signifMotifNames <- motifEnrichmentTable_wGenes$motif[1:3]
# incidenceMatrix <- getSignificantGenes(geneLists$geneSetName, 
#                                        motifRankings,
#                                        signifRankingNames=signifMotifNames,
#                                        plotCurve=TRUE, maxRank=5000-20, 
#                                        genesFormat="incidMatrix",
#                                        method="aprox")$incidMatrix
# 
# 
# library(reshape2)
# # # edges <- melt(incidenceMatrix)
# edges <- edges[which(edges[,3]==1),1:2]
# colnames(edges) <- c("from","to")
# library(visNetwork)
# motifs <- unique(as.character(edges[,1]))
# genes <- unique(as.character(edges[,2]))
# nodes <- data.frame(id=c(motifs, genes),   
#                     label=c(motifs, genes),    
#                     title=c(motifs, genes), # tooltip 
#                     shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
#                     color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
# visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
#                                         nodesIdSelection = TRUE)

#- BRETIGEA package for cell type in brain
library(BRETIGEA)
exp = datExpr; rownames(exp) = NULL; rownames(exp) = datProbes$external_gene_name
brCells = brainCells(inputMat = exp,nMarker = 100,method = "PCA")
cells = findCells(inputMat = exp, markers = markers_df_brain, nMarker = 100, method = "PCA")

annot = data.frame(disease = as.factor(datMeta$Dx),
                   Region = as.factor(datMeta$Brain_lobe))
rownames(annot)= rownames(brCells)
pdf("./results/figures/celltype_samples.pdf")
pheatmap(as.matrix(brCells),
         scale = "column",
         # color = redblue(n = 9),
         show_rownames = F, show_colnames = T,
         annotation_row = annot,
         cluster_cols = F,cluster_rows = F)
dev.off()

br <- melt(brCells, id= 1:6)
colnames(br) = c("sample", "cell_type", "Specificity")

celldata = data.frame()
for (lobe in levels(datMeta$Brain_lobe)){
idx = which(datMeta$Brain_lobe == lobe)
br_lobe = brCells[idx,]
br_lobe = melt(br_lobe, id =1:6)
br_lobe$lobe = lobe
celldata= rbind(celldata, br_lobe)
}

colnames(celldata) = c("sample", "cell_type", "Specificity", "lobe")
celldata$disease = celldata$sample
celldata$disease = datMeta$Dx[match(celldata$sample,
                                    rownames(datMeta$Dx))]
p= ggplot(celldata, 
         aes(x= disease,
            y = Specificity, 
           color = cell_type))+
labs(x="Cell type",
    y = "Enrichment",
    title = "Cell types enrichment in brain lobes")+
 geom_jitter()+
geom_boxplot(width = 0.2,
            outlier.size = 2)+
  facet_wrap(~lobe, ncol =1)+
  theme(text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5))+
  geom_abline(intercept = 0, slope = 0, lty = 2)
p
ggsave(filename = "./results/figures/Lobes_celltypes.pdf" ,
       plot = p,
       height = 10, 
       width = 4)



save.image("./codes/GO.Rdata")

