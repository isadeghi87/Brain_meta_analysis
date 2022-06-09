

setwd("../Desktop/Iman/Brain_meta/Brain/mega_analysis/")

# libraries
library(nlme);library(dplyr);library(ComplexHeatmap);library(corrplot)
library(RColorBrewer); library(pheatmap)

#### load data ####

load("./codes/finalNetwork.RData")
load("./codes/Combined.Rdata")
logFCs = readRDS(file = "./codes/shared_logFCs.rds")  

# clean data
datMeta$Dx = gsub("Normal","CTL", datMeta$Dx)
datMeta$Dx = factor(datMeta$Dx, levels = c("CTL", "AD","PD","PA",
                                           "PSP","Scz","ASD","BP","MDD"))

# filter only cell-type specific modules

modColors = data.frame(color= colors,
                       row.names = rownames(datExpr))
cellMod = c("yellow","brown","purple","magenta","turquoise") 
modColors = subset(modColors, colors %in% cellMod)
modColors$color = factor(modColors$color, levels = cellMod)

datExpr = datExpr[rownames(modColors),]

# compute DE 
meta = matrix(NA, nrow=nrow(datExpr), ncol=8)
colnames(meta) = c("AD","PD","PA","PSP","Scz","ASD","BP","MDD")

for (dx in colnames(meta)) {
  idx = c(dx,"CTL")
  metaData = datMeta[datMeta$Dx %in% idx,]
  id = match(rownames(metaData),colnames(datExpr))
  expData =  datExpr[,id]

for(i in 1:nrow(expData)) {
  if(i%%100==0) print(i)
  expr = as.numeric(expData[i,])
  tryCatch({
    beta = summary(lme(expr~ Dx +  Study , data = metaData, random=~1|Subject_ID))$tTable[2,1]
    meta[i,dx] = beta
    }, error=function(e){})
}
}

meta=as.data.frame(meta)
rownames(meta) = rownames(datExpr)

#### check correlation of cell types across diseases ####

geneCelltype = modColors
geneCelltype$cell_type = geneCelltype$color
geneCelltype$cell_type = str_replace_all(geneCelltype$cell_type, c("yellow" = "Astrocyte",
                                                           "brown" = "Oligodendrocyte",
                                                           "purple"= "Microglia",
                                                           "magenta"= "Neuron",
                                                           "turquoise" = "Neuron"))

geneCelltype= cbind(geneCelltype,meta)

#### compute correlation #####

corrList = vector("list",4)
names(corrList) = unique(geneCelltype$cell_type)

for (i in names(corrList)){
  dat = geneCelltype[geneCelltype$cell_type == i,3:ncol(geneCelltype)]
  corr = cor(dat, use = "pairwise.complete.obs",method = "spearman")
  corrList[[i]] = corr
  
} 

pdf("./results/figures/cellTypeCorrletions.pdf", width = 7, height = 7)
par(mfrow=c(2,2))

for (i in names(corrList)) {
  corrplot(corrList[[i]],
           col = brewer.pal(8,"PuOr"),
           title = i,
           mar = c(0,0,2,1),
           tl.srt = 45,
           addCoef.col = T,
           diag = T,
           order = "FPC", 
           hclust.method = "complete", 
           method = "shade", 
           tl.col = "black",
           is.corr = T,
           number.cex = 0.8,
           tl.cex = 0.8)
}
dev.off()

# load hub genes
hubGenes = readRDS("./codes/hubgenes.rds")
id = intersect(hubGenes, rownames(geneCelltype))
hubExp = geneCelltype[id,]

# load attributes
attr = readRDS("./codes/attributes.rds")

hubExp$symbol = attr$external_gene_name[match(rownames(hubExp),attr$ensembl_gene_id_version)]
colnames(hubExp) = gsub("color","module",
                        colnames(hubExp))

hubExp$module = factor(hubExp$module, levels = c("yellow",
                                                 "brown",
                                                 "purple",
                                                 "magenta",
                                                 "turquoise"))
hubExp$cell_type = factor(hubExp$cell_type, levels = c("Astrocyte",
                                                       "Oligodendrocyte",
                                                       "Microglia",
                                                       "Neuron"))
hubExp = hubExp[order(hubExp$module),]
annot_col = list(module = c(yellow ="yellow",
                            brown = "brown",
                            purple = "purple",
                            magenta = "magenta",
                            turquoise = "turquoise"),
                 `cell_type` = c(Astrocyte = "yellow",
                                 Oligodendrocyte = "brown",
                                 Microglia = "purple",
                                 Neuron = "black"))

####heatmap of hub genes expression####

dat= hubExp[,3:10]

pdf("./results/figures/hubgeneHeatmap.pdf",
    width = 4.2, height = 8)

par(mar= c(0,2,1,0))
pheatmap(as.matrix(dat),
         color =redWhiteGreen(1000,gamma = 0.5)[0:750],
         show_rownames = T,
         # scale = "none",
         show_colnames = T,
          legend_breaks =c(-2,-1,0,1,max(dat)),
          legend_labels = c("-2","-1","0","1","beta\n"),
         fontsize_row = 5,
         fontsize_col = 7,
         fontsize = 7,
         angle_col = 45,
         annotation_row = hubExp[,c(1,2)],
         cluster_rows = F,
         cluster_cols = T,
         gaps_row = c(20,40,60),
         cellwidth = 14,
         treeheight_col = 1,
         cellheight = 5,
         annotation_colors = annot_col,
         labels_row = hubExp$symbol)

dev.off()

save(file = "./codes/cellTypeCorrelations.Rdata", meta,metaData,datExpr,datMeta,
     colors,modColors,hubExp)
