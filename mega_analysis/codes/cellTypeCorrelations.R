## cellTypeCorrelations.R
### correlation of cell type modules across diseases


setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")

# libraries
library(nlme);library(dplyr);library(ComplexHeatmap);library(corrplot)
library(RColorBrewer); library(pheatmap)

#### load data ####

load("./codes/finalNetwork.RData")
# load("./codes/Combined.Rdata")
# logFCs = readRDS(file = "./results/tables/shared_logFCs.rds")  


# filter only cell-type specific modules

modColors = data.frame(module = mods,
                       color = colors,
                       row.names = rownames(datExpr))
cellMod = c("M2","M5","M8","M10","M13","M14") 
modColors = subset(modColors, module %in% cellMod)
modColors$color = factor(modColors$color, levels = unique())

exp = datExpr[rownames(modColors),]

# compute DE 
meta = matrix(NA, nrow=nrow(exp), ncol=8)
colnames(meta) = c("AD","PD","PA","PSP","Scz","ASD","BP","MDD")

for (dx in colnames(meta)) {
  idx = c(dx,"CTL")
  metaData = allMeta[allMeta$Dx %in% idx,]
  id = match(rownames(metaData),colnames(exp))
  expData =  exp[,id]

for(i in 1:nrow(expData)) {
  if(i%%100==0) print(i)
  expr = as.numeric(expData[i,])
  tryCatch({
    beta = summary(lme(expr~ Dx +  Study , data = metaData, random=~1|Subject_ID))$tTable[2,1]
    meta[i,dx] = beta
    }, error=function(e){})
}
}

meta = as.data.frame(meta)
rownames(meta) = rownames(exp)

#### check correlation of cell types across diseases ####

geneCelltype = modColors
geneCelltype$cell_type = geneCelltype$module
geneCelltype$cell_type = str_replace_all(geneCelltype$cell_type,
                                         c("M2" = "Oligodendrocyte",
                                           "M5" = "Astrocyte",
                                           "M8"= "Neuron",
                                           "M10"= "Microglia",
                                           "M13" = "Neuron",
                                           "M14" = "Endothelial"))

geneCelltype = cbind(geneCelltype, meta)

#### compute correlation #####

corrList = vector("list",5)
names(corrList) = unique(geneCelltype$cell_type)

for (i in names(corrList)){
  dat = geneCelltype[geneCelltype$cell_type == i,3:ncol(geneCelltype)]
  corr = cor(dat[,-1], use = "pairwise.complete.obs",method = "spearman")
  corrList[[i]] = corr
  
} 

pdf("./results/figures/network_analysis/cellTypeCorrletions.pdf", width = 10, height = 7)
par(mfrow=c(2,3))

for (i in names(corrList)) {
  corrplot(corrList[[i]],
           col = brewer.pal(8,"PuOr"),
           title = i,
           mar = c(1,0,2,1),
           tl.srt = 45,
           addCoef.col = T,
           diag = T,
           order = "FPC", 
           hclust.method = "complete", 
           method = "shade", 
           tl.col = "black",
           is.corr = T,
           number.cex = 1,
           tl.cex = 0.8)
}
dev.off()

# load hub genes
hubGenes = readRDS("./results/tables//hubgenes.rds")
id = intersect(hubGenes, rownames(geneCelltype))
hubExp = geneCelltype[id,]

# load attributes
attr = readRDS("./results/tables/attributes.rds")

hubExp$symbol = attr$external_gene_name[match(rownames(hubExp),attr$ensembl_gene_id)]

hubExp$module = factor(hubExp$module, levels = c("M2",
                                                 "M5",
                                                 "M8",
                                                 "M13",
                                                 "M10",
                                                 "M14"))
hubExp$cell_type = factor(hubExp$cell_type, levels = c( "Oligodendrocyte",
                                                        "Astrocyte",
                                                        "Neuron", 
                                                        "Microglia",
                                                       "Endothelial"))
hubExp = hubExp[order(hubExp$module),]

## color for modules
mod_col = unique(hubExp$color)
names(mod_col) = levels(hubExp$module)

## color for cell types

annot_col = list(module = mod_col,
                 cell_type = c(Oligodendrocyte = "blue",
                                 Astrocyte = "green",
                                 Neuron = "salmon",
                                 Microglia = "purple",
                                 Endothelial = "cyan"))

####heatmap of hub genes expression####

dat = hubExp[,4:11]

pdf("./results/figures/network_analysis/hubGeneHeatmap.pdf",
    width = 4.5, height = 7)

par(mar= c(0,0,0,0))
pheatmap(as.matrix(dat),
         color =divergex_hcl(1000,palette = "PRGn")[100:1000],
         show_rownames = T,
         show_colnames = T, 
          legend_breaks =c(-0.5,0,0.5,1,max(dat)),
          legend_labels = c("-0.5","0","0.5","1","beta\n"),
         fontsize_row = 5,
         fontsize_col = 7,
         border_color = "grey60",
         fontsize = 7,
         angle_col = 45,
         annotation_row = hubExp[,c(1,3)],
         cluster_rows = F,
         cluster_cols = T,
         gaps_row = c(10,20,40,50),
         cellwidth = 14,
         cellheight = 6.5,
         treeheight_col = 1,
         annotation_colors = annot_col,
         labels_row = hubExp$symbol)

dev.off()

save(file = "./codes/cellTypeCorrelations.Rdata", meta,metaData,datExpr,allMeta,
     colors,modColors,hubExp)
