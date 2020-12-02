#5b_networkParam.R
#---Here we compute the dendogram and network parameters for the modules

rm(list=ls())
library(WGCNA); library(biomaRt);library(ggpubr); 
library(ComplexHeatmap); library(pheatmap);library(RColorBrewer)
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")
load("./codes/Combined_for_network.Rdata")

enableWGCNAThreads()
allowWGCNAThreads()
condition =TRUE

#load TOM comptued by rWGCNA
load("./results/figures/WGCNA/network_signed20_exprSet1-block.1.RData")

if(condition){
  geneTree = hclust(1-as.dist(TOM), method="average")
}

# Iterate WGCNA parameters for robustness -- this takes a while
colors = vector(mode="list")
labels = vector(mode="list")
if (condition){
  for (pam in c(FALSE,TRUE)) {
    for (minModSize in c(50,100, 200)) {
      for (dthresh in c(0.1, 0.2)) {
        for(ds in c(0:4)) { 
          print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
          
          tree = cutreeHybrid(
            dendro = geneTree,
            minClusterSize = minModSize,
            pamStage = pam, 
            cutHeight = 0.999,
            deepSplit = ds,
            distM = as.matrix(1-as.dist(TOM)))
          merged = mergeCloseModules(
            exprData = t(datExpr), 
            colors = tree$labels, 
            cutHeight=dthresh)
          colors = cbind(colors, labels2colors(merged$colors))
          labels = c(labels, 
                     paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
        }
      }
    }
  }
}


save(file="./results/figures/WGCNA//WGCNA_diffParams.Rdata", geneTree, colors, labels)

for (i in 1:60){
  l = length(unique(colors[,i]))
  print(c(i,l))
}

if(condition){
pdf("./results/figures/WGCNA/WGCNA_diffParams.pdf", width=10, height=14)

par(mar= c(2,5,2,2))
plotDendroAndColors(geneTree,
                    colors,
                    addGuide=T,
                    dendroLabels=F)

plotDendroAndColors(geneTree, 
                    colors,
                    groupLabels = labels,
                    addGuide= TRUE,
                    dendroLabels=FALSE,
                    main="Dendrogram",
                    cex.colorLabels=0.5)
dev.off()

colors2 = cbind(colors[,24], colors)
labels = c("Final Modules", labels)

pdf("./results/figures/WGCNA/WGCNA_finalModule.pdf", width=10, height=14)
plotDendroAndColors(geneTree,
                    colors2,
                    groupLabels = labels,
                    addGuide=T,
                    dendroLabels=F,
                    cex.colorLabels=0.5)
dev.off()

}



# Finalized Parameters MMS=200,DCOR=0.1,PAM=FALSE, ds= 0
# --------------------
# Parameters to Use: "DS=4,MMS=50,DCOR=0.1,PAM=FALSE" ->14 modules
wgcna_parameters = list(powers =  20)
wgcna_parameters$minModSize = 200
wgcna_parameters$minHeight = 0.1
wgcna_parameters$bsize = 20000  ##block size needs to be larger than dim(datExpr)[1]
wgcna_parameters$ds = 3  ##deep split parameter contorls number of modules
wgcna_parameters$networkType = "signed"    ## using signed networks
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = FALSE

if(condition){
  tree = cutreeHybrid(
    dendro = geneTree,
    minClusterSize = wgcna_parameters$minModSize,
    pamStage = wgcna_parameters$pamStage,
    cutHeight = 0.999, 
    deepSplit = wgcna_parameters$ds,
    distM=as.matrix(1-as.dist(TOM)))

  merged = mergeCloseModules(exprData= t(datExpr),
                             colors = tree$labels,
                             cutHeight = wgcna_parameters$minHeight)
}

# collect colors 
colors = labels2colors(merged$colors)

# rename modules
mods = paste0("M",merged$colors)
names(mods) = colors
# obtain number of genes per module
colTable = as.data.frame(table(mods))

# plot the table
colplot = ggtexttable(colTable,cols = c("Module", "Number of genes"), theme = ttheme("mBlueWhite"))
ggsave(filename = "./results/figures/network_analysis/ModulesGenesNumber.pdf", 
       plot = colplot,
       height = 5, width = 4)
write.table(table(colors), 
            file = "./results/figures/WGCNA/ModuleColorsNumber.csv", row.names = F, col.names = T )

# colors
cols = table(names(mods)) %>% sort(decreasing = T) %>% names
colTable = colTable[order(colTable$Freq,decreasing = T),]
colTable$color = names(mods)[match(colTable$mods,mods)]
colTable$mods = factor(colTable$mods,levels = colTable$mods)
colTable$color = factor(colTable$color,levels = colTable$color)

mod_p = ggplot(colTable, aes( x = mods,
                              y = Freq,
                              fill = mods))+
  geom_point(
             colour= "black",
             stroke=1,
             size = 4,shape = 21)+
  scale_fill_manual(values = levels(colTable$color))+
  scale_x_discrete()+
  scale_y_log10()+
  labs(y = "module size",
       x = "module")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 8, color = "black"),
        axis.text.x.bottom  =element_text(angle = 45, hjust = 1))
mod_p

ggsave("./results/figures/network_analysis/ModuleSize.pdf", 
       width = 5, height = 3,
       mod_p)
pdf(file = "./results/figures/WGCNA/WGCNA_finalModules.pdf", 
    height = 3, width = 8)

plotDendroAndColors(geneTree,
                    colors,
                    rowTextAlignment = "center", 
                    groupLabels = "Module",
                    main = "",
                    cex.colorLabels = 0.5,
                    addGuide = T,
                    dendroLabels = F)

dev.off()
MEs = moduleEigengenes(expr = t(datExpr),
                       merged$colors, 
                       softPower = wgcna_parameters$powers,
                       scale = F)
kMEtable = signedKME(t(datExpr),
                     MEs$eigengenes)
tableSup = data.frame(kMEtable)

######check correlation with metadata
MEs0 = as.data.frame(MEs$eigengenes)
# MEs0 = orderMEs(MEs0)
MEs0 = MEs0[,colnames(MEs0) != "ME0"]
colnames(MEs0) = gsub("ME","M",colnames(MEs0))

# merge datmeta and modules eigengene
ME.datMeta = cbind(allMeta, MEs0)

# a df of correlations
moduleCors = as.data.frame(matrix(NA, 
                                  ncol = ncol(allMeta), 
                                  nrow = ncol(MEs0)))
rownames(moduleCors) = colnames(MEs0);
colnames(moduleCors) = colnames(allMeta)

# exlcude sample and subjects from correlation 
moduleCors = moduleCors[ ,-c(1:2)]
#p values
modulesPval = moduleCors

#calculate linear regression between modules and covariates
for (mod in rownames(moduleCors)){
  for (var in colnames(moduleCors)){
    mcor =lm(ME.datMeta[,mod] ~ME.datMeta[,var], data = ME.datMeta )
    moduleCors[mod,var] = summary(mcor)$adj.r.squared
    modulesPval[mod, var] = anova(mcor)$`Pr(>F)`[1]
  }
}

# rownames(moduleCors) = gsub("ME","", rownames(moduleCors))
###heatmap 
moduleCors = as.matrix(moduleCors)
moduleCors[moduleCors <0]=0
modulesPval = as.matrix(modulesPval)
pval = signif(-log10(modulesPval),1)
symbol = pval
symbol[moduleCors >0.2] ="**"
symbol[moduleCors <0.2 & moduleCors >0.1 ]="*"
symbol[moduleCors <0.1] = ""

#color
library(RColorBrewer)
col_rnorm = brewer.pal(9,name = "Blues")
column = c("Study","Disease", "Age", "Sex","Brain region","Brain lobe")

annot_mod = data.frame(module = rownames(moduleCors))
rownames(annot_mod)= rownames(moduleCors)
mod_col = unique(names(mods))
names(mod_col) = unique(mods)
annot_col = list(module = mod_col)

pdf("./results/figures/network_analysis/Modules_covariates_heatmap.pdf", 
    width = 4,height = 4)

pheatmap(signif(moduleCors,2),
         color = blueWhiteRed(1000)[500:0],
         # annotation_row = annot_mod,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "correlation",
         # annotation_colors = annot_col,
         clustering_method = "single",
         annotation_legend = F,
         treeheight_row = 1,
         treeheight_col = 1,
         labels_col = column,
         angle_col = 45,
         fontsize = 7,
         fontsize_col = 7, 
         fontsize_row = 6,
         cellwidth = 20,
         cellheight = 15,
         legend_breaks = c(0,0.1,0.2,0.3,max(moduleCors)),
         legend_labels =c("0","0.1","0,2","0.3", "adj. R square\n")
         
)
dev.off()

## correlation of each gene with covariates
exp_meta = cbind(allMeta,t(datExpr))

# data frame to collect results
cor_res = as.data.frame(matrix(NA, ncol = 5,
                               nrow = nrow(datExpr)))
colnames(cor_res) = c("Study","Dx","Age","Sex","Brain_Lobe" )
rownames(cor_res)= rownames(datExpr)
p_res = cor_res

# correlation for each gene and covariate
for(cov in colnames(cor_res)){
  for(id in rownames(datExpr)){
    
    mcor =lm(exp_meta[,id] ~ exp_meta[,cov], data = exp_meta )
    cor_res[id,cov] = summary(mcor)$adj.r.squared
    p_res[id, cov] = anova(mcor)$`Pr(>F)`[1]  }
}

cor_res= t(cor_res)
cor_res= as.data.frame(cor_res)
# add colors
cor_res$color = colors

#----Annotating gene IDs
getinfo <- c( "ensembl_gene_id",
              "external_gene_name",
              "gene_biotype")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
datProbes <- getBM(
  attributes = getinfo,
  filters = c("ensembl_gene_id"),
  values = rownames(tableSup),
  mart = mart)

tableSup = cbind(datProbes, 
                data.frame(Module.Color=colors,
                           Module.name = paste0("M",merged$colors)),
                tableSup)

write.csv(file="./results/tables/Modules_kMEtable.csv", tableSup)
save(file="./codes/finalNetwork.RData",
     datExpr,
     allMeta,
     datProbes,
     geneTree,
     wgcna_parameters,
     colors,
     cols,mods,
     MEs,
     kMEtable)
