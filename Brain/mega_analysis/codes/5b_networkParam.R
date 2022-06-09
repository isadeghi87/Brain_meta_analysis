#5b_networkParam.R
#---Here we compute the dendogram and network parameters for the modules

rm(list=ls())
library(WGCNA); library(biomaRt);library(ggpubr); library(ComplexHeatmap)
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")
load("./codes/st4_combinedData.Rdata")

enableWGCNAThreads()
allowWGCNAThreads()
condition =TRUE
#load TOM comptued by rWGCNA
load("./results/figures/WGCNA/network_signed10_exprSet1-block.1.RData")
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
  
  
 
  colors2 = cbind(colors[,15], colors)
  labels = c("Final Modules", labels)
  
 
  plotDendroAndColors(geneTree,
                      colors,
                      groupLabels = labels,
                      addGuide=T,
                      dendroLabels=F,
                      cex.colorLabels=0.5)
  dev.off()




# Finalized Parameters
# --------------------
# Parameters to Use: "DS=4,MMS=50,DCOR=0.1,PAM=FALSE" ->14 modules
wgcna_parameters = list(powers =  10)
wgcna_parameters$minModSize = 50
wgcna_parameters$minHeight = 0.1
wgcna_parameters$bsize = 25000  ##block size needs to be larger than dim(datExpr)[1]
wgcna_parameters$ds = 1  ##deep split parameter contorls number of modules
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
}

if(condition){
merged = mergeCloseModules(exprData= t(datExpr),
                           colors = tree$labels,
                           cutHeight = wgcna_parameters$minHeight)
}
colors = labels2colors(merged$colors)
colTable = as.data.frame(table(colors))
colplot =ggtexttable(colTable,cols = c("Module", "Number of genes"), theme = ttheme("mBlueWhite"))
ggsave(filename = "./results/figures/ModulesGenesNumber.pdf", plot = colplot,height = 5, width = 4)
write.table(table(colors), 
            file = "./results/figures/WGCNA/ModuleColorsNumber.csv", row.names = F, col.names = T )

mods =ggplot(colTable, aes( x = reorder(colors, -Freq), y = Freq))+
  geom_point(aes(fill = colors),
             colour= "black",
             stroke=1,
             size = 4,shape = 21)+
  scale_fill_manual(values = as.character(colTable$colors))+
  scale_y_log10()+
  labs(y = "module size",
       x = "module")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 8, color = "black"),
        axis.text.x.bottom  =element_text(angle = 45, hjust = 1))

ggsave("./results/figures/ModuleSize.pdf", 
       width = 5, height = 3,
       mods)
pdf(file = "./results/figures/WGCNA/WGCNA_finalModule.pdf", 
    height = 3, width = 8)
plotDendroAndColors(geneTree,
                    colors,rowTextAlignment = "center", 
                    groupLabels = "Module",
                    main = "",
                    cex.colorLabels = 0.5,
                    addGuide = T,
                    dendroLabels = F)

dev.off()
MEs = moduleEigengenes(expr = t(datExpr),
                       colors, 
                       softPower = wgcna_parameters$powers)
kMEtable = signedKME(t(datExpr),
                     MEs$eigengenes)
tableS1 = data.frame(kMEtable)
colnames(tableS1) = paste0("CD", 1:11, ".", colnames(tableS1))

######check correlation with metadata
MEs0 = as.data.frame(MEs$eigengenes)
MEs0 = orderMEs(MEs0)
MEs0 = MEs0[,colnames(MEs0) != "MEgrey"]

#merge datmeta and modules eigengene
ME.datMeta = cbind(datMeta, MEs0)

# a df of correlations
moduleCors = as.data.frame(matrix(NA, 
                                  ncol = ncol(datMeta), 
                                  nrow = ncol(MEs0)))
rownames(moduleCors) = colnames(MEs0);
colnames(moduleCors) = colnames(datMeta)
moduleCors = moduleCors[,colnames(moduleCors) != "Subject_ID"]

#p values
modulesPval = moduleCors
#calculate linear regression between modules and covariates
for (mod in rownames(moduleCors)){
  for (var in colnames(moduleCors)){
    mcor =lm(ME.datMeta[,mod] ~ ME.datMeta[,var], data = ME.datMeta )
    moduleCors[mod,var] = summary(mcor)$adj.r.squared
    modulesPval[mod, var] = anova(mcor)$`Pr(>F)`[1]
    }
}
rownames(moduleCors) = gsub("ME","", rownames(moduleCors))

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
column = c("Study","Disease","Brain region", "Age", "Sex","PMI","Brain lobe")

annot_mod = data.frame(module = rownames(moduleCors))
rownames(annot_mod)= rownames(moduleCors)
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

pdf("./results/figures/Modules_covariates_heatmap.pdf", 
    width = 4,height = 3
    )

pheatmap(signif(moduleCors,2),
         color = bluered(1000)[500:0],
         annotation_row = annot_mod,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "correlation",
         annotation_colors = annot_col,
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
#----Annotating gene IDs
getinfo <-
  c("ensembl_gene_id_version",
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "band",
    "gene_biotype")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
datProbes <- getBM(
  attributes = getinfo,
  filters = c("ensembl_gene_id_version"),
  values = rownames(tableS1),
  mart = mart
)
datProbes <-  datProbes[match(rownames(tableS1), datProbes$ensembl_gene_id_version), ]

tableS1 = cbind(datProbes, 
                data.frame(Module.Color=colors,
                           Module.name = paste0("CD",merged$colors)),
                tableS1)

write.csv(file="./results/tables/Modules_kMEtable.csv", tableS1)
save(file="./codes/finalNetwork.RData",
     datExpr,
     datMeta,
     datProbes,
     geneTree,
     colors,
     wgcna_parameters,
     colors,
     MEs,
     kMEtable)
