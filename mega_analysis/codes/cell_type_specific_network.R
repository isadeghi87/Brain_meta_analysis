## cell_type_specific_network.R
rm(list=ls()); options(stringsAsFactors = F)
library(pSI); library(gProfileR); library(gplots); library(biomaRt); library(WGCNA)
library(ggpubr); library(pheatmap); library(dplyr); library(magrittr)
library(GeneOverlap); library(factoextra);library(FactoMineR)

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")

#---load final network
load("./codes/finalNetwork.RData")

eigmat = MEs$eigengenes
colnames(eigmat) = gsub("ME","M",colnames(eigmat))
kME = signedKME(t(datExpr), MEs$eigengenes)
colnames(kME) = gsub("kME", "M", colnames(kME))

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
getinfo <- c("ensembl_gene_id",
             "external_gene_name")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
bm <- getBM(
  attributes = getinfo,
  filters = "hgnc_symbol",
  values = zhang.datProbes$symbol,
  mart = mart)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host = "apr2018.archive.ensembl.org") 
bm = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
           filters = "hgnc_symbol", values = zhang.datProbes$symbol, mart=ensembl)
zhang.datProbes = data.frame(symbol=zhang.datProbes$symbol, 
                             ensg=bm$ensembl_gene_id[match(zhang.datProbes$symbol, bm$external_gene_name)])
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
cell.p.zhang = matrix(NA, nrow = length(unique(mods)),
                      ncol = 5)
rownames(cell.p.zhang) = unique(mods) %>% sort
colnames(cell.p.zhang) = colnames(pSI.output)


for(mod in rownames(cell.p.zhang)) {
  f = fisher.iteration(pSI.output, rownames(datExpr)[mods==mod],p.adjust = F)
  cell.p.zhang[mod,] = f$`0.05 - nominal`
}

#remove grey module
cell.p.zhang = cell.p.zhang[-1,]

# fdr correction
cell.p.zhang.fdr = p.adjust(cell.p.zhang,"fdr")
dim(cell.p.zhang.fdr) = dim(cell.p.zhang); dimnames(cell.p.zhang.fdr) = dimnames(cell.p.zhang);

dendro.col = as.dendrogram(hclust(as.dist(1-bicor(zhang.datExpr)), method="average"))
dendro.row= as.dendrogram(hclust(
  as.dist(1-bicor(
    eigmat[,-1])),method="average"))


#annotation of modules
annot_mod = data.frame(module = unique(mods)[-2])
rownames(annot_mod) = unique(mods)[-2]

# color for modules
mod_col = unique(names(mods))[-2]
names(mod_col) = unique(mods)[-2]
annot_col = list(module = mod_col)

# plot heatmap
pdf(file= "./results/figures/network_analysis/CellTypesEnrichmentModules.pdf",
    height = 6,
    width = 6)
par(mar=c(0,0,0,0))
pheatmap(-log10(cell.p.zhang.fdr),
         # color = brewer.pal(9,"Greens")[0:8],
         color = redWhiteGreen(1000,1)[500:1000],
         annotation_row = annot_mod,
         annotation_colors = annot_col,
         annotation_legend = F,
         cluster_cols = T,
         cluster_rows = T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         fontsize_col = 10,
         fontsize_number = 10,
         fontsize_row = 10,
         treeheight_row = 5,
         cellwidth = 45,
         cellheight = 20,
         angle_col = 45,
         show_rownames = T,
         show_colnames = T,
         legend_breaks = c(0,100,200, max(-log10(cell.p.zhang.fdr))),
         legend_labels = c(0,100,200,"-log10 (FDR)\n"),
         treeheight_col = 1, 
         display_numbers = signif(cell.p.zhang.fdr,2))

dev.off()

out = cell.p.zhang.fdr
colnames(out) = paste(colnames(out), ".FDR",sep="")
write.csv(out,file="./results/tables/CellTypeEnrichmentFDR.csv")

#check cell types with PanglaoDB database
all_genes = data.frame(row.names = rownames(datExpr),
                       color = colors,
                       module = mods)
PangDB <- read.csv("./results/tables/PanglaoDB_markers_11_Sep_2019.csv", sep ='\t')
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
# ModuleGenes$module = as.factor(ModuleGenes$module)
ModuleGenelist = list()
for (mod in unique(ModuleGenes$module)){
  ModuleGenelist[[mod]] = ModuleGenes[ModuleGenes$module == mod,]
}

cell_mod = as.data.frame(matrix(NA, ncol = length(ModuleGenelist), nrow = length(cellGene)))
rownames(cell_mod) = names(cellGene);
colnames(cell_mod) = names(ModuleGenelist)
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
id = which(colnames(cell_mod) == "M0")
cell_mod = cell_mod[,-id]
cell_mod_ODD = cell_mod_ODD[,-id]

mat.plot = as.matrix(-sign(cell_mod_ODD) * log10(cell_mod))


# plot
pheatmap(-log10(cell_mod),
         color = redWhiteGreen(2000,1)[1000:2000],
         cluster_cols = T,
         cluster_rows = T,border_color = "black",
         # clustering_distance_rows = "euclidean",
         # clustering_distance_cols = "euclidean",
         annotation_col = annot_mod,
         annotation_colors = annot_col,
         annotation_legend = F,
         annotation_names_col = F,
         fontsize_col = 7,
         fontsize = 5,
         fontsize_number = 6,
         fontsize_row = 6,
         treeheight_col = 1, 
         treeheight_row = 1,
         cellwidth = 20,
         cellheight = 10,
         angle_col = 45,
         show_rownames = T,
         show_colnames = T,
         legend_breaks = c(0,10,20,30, max(-log10(cell_mod))),
         legend_labels = c(0,10,20,30,"-log10 (FDR)\n"),
         display_numbers = signif(cell_mod,1),
         filename = "./results/figures/network_analysis/CellModule_PangDB.pdf",
         height = 5.5,
         width = 7,
)


## Calculate TFBS enrichment of Disease-associated modules ####
resultsTFBS = data.frame()
for (m in unique(mods)) {
  
  ## GO Enrichment
  query = rownames(datExpr)[mods==m]
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

# extract tf and motifs from term names
topTF$term.name = str_extract_all(topTF$term.name, pattern = "Factor.*[A-Z]")
topTF$color = names(mods)[match(topTF$m,mods)]
id = which(topTF$m == "M0")
topTF = topTF[-id,]
labels = topTF$term.name

# plot for tf and motifs
pdf("./results/figures/network_analysis/TFBsPlot.pdf",width = 5, height=5)
par(oma=c(3,5,0,1), mar=c(4,5,0,0))
bp = barplot(-log10(topTF$p.value), 
             # main="Top TF binding sites for modules",
             horiz=T, yaxt='n', 
             col = as.character(topTF$color),
             cex.lab=0.5,
             xlab='-log10(FDR)', cex.axis = .4, border=NA)
axis(2,at=bp,labels = labels,tick=FALSE,las=2,cex.axis=0.4);
abline(v=-log10(0.05), col="black", lwd=2,lty=2)
dev.off()


## Identify hub genes transcription factors ####
rownames(datProbes) = datProbes$ensembl_gene_id
hubGenes = data.frame(Module=NA, Gene=NA, Symbol=NA,Rank=NA)
for (m in unique(mods)) {
  mod = rownames(datExpr)[mods == m]
  mod = mod[order(kME[mod,m],decreasing = T)[1:20]]
  sym = datProbes[mod,"external_gene_name"]
  hubGenes = rbind(hubGenes, data.frame(Module=m, Gene=mod, Symbol = sym, Rank=1:20))
}

hubGenes = hubGenes[-1,]

## retrieve attributes for tf
a=listAttributes(mart);f=listFilters(mart)
bm = getBM(attributes = c("ensembl_gene_id",
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

## eRNA overlap with co-expression modules ####
## Download eRNA modules from http://www.nature.com/neuro/journal/v18/n8/extref/nn.4063-S12.xlsx
eRNAnetwork=  read.csv("./results/tables/eRNA_Yao.etal.csv")
colnames(eRNAnetwork)[1] = "Symbol"
source("./codes/Fisher_overlap.R")
moduleNames = unique(mods)[-2]
table.p = matrix(NA, nrow=length(moduleNames), ncol=19)
rownames(table.p) = moduleNames; 
colnames(table.p) = paste("M", 1:19,sep="")
table.or = table.p
hgnc = datProbes$external_gene_name
for (m in moduleNames) {
  for(e in colnames(table.p)) {
    f = ORA(hgnc[mods==m],eRNAnetwork$Symbol[eRNAnetwork$Module.Label==e], hgnc, eRNAnetwork$Symbol)
    table.or[m,e] = as.numeric(f[[1]])
    table.p[m,e] = as.numeric(f[[2]])
  }
}

table.p.fdr = p.adjust(table.p,"fdr")
dim(table.p.fdr) = dim(table.p); dimnames(table.p.fdr) = dimnames(table.p)
table.p.fdr[table.p.fdr>0.05]=1

#remove grey module
colnames(table.p.fdr) <- paste("eModule",1:19,sep = "")

text = signif(-log10(table.p.fdr),1)

# plot
pdf("./results/figures/network_analysis/BrainEnhancerEnrichment.pdf",
    width = 5.5,
    height = 3 )
par(margin(0,0,0,0))
pheatmap(-log10(table.p.fdr),
         color = blueWhiteRed(n = 2000)[1000:0],
         cluster_cols = T,
         cluster_rows = T,
         annotation_row = annot_mod,
         annotation_colors = annot_col,
         annotation_legend = F,
         annotation_names_row = F,
         fontsize = 7,
         fontsize_number = 5,
         fontsize_row = 7,
         # fontsize_col = 4,
         treeheight_row = 1,
         cellwidth = 15,
         cellheight = 10,
         angle_col = 45,
         show_rownames = T,
         show_colnames = T,
         # main = "Enrichment of Brain Enhancers Modules",
         legend_breaks = c(0,50,100,150,200,250, max(-log10(table.p.fdr))),
         legend_labels = c(0,50,100,150,200,250,"-log10 (FDR)\n"),
         treeheight_col = 1, 
         display_numbers = text,
         # filename = "./results/figures/BrainEnhancerEnrichment.pdf",
         # height = 4,
         # width = 6
)

dev.off()  

#Calculate mitochondrial overlap with Winden modules
#Download winden modules from http://msb.embopress.org/content/5/1/291.long#sec-23
winden= read.csv("./results/tables/winden.csv")
ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="mmusculus_gene_ensembl",
                  host="feb2014.archive.ensembl.org")
a=listAttributes(ensembl);f=listFilters(ensembl)
bm1= getBM(attributes = c("affy_moe430a","ensembl_gene_id"),
           filters = "affy_moe430a", values = winden$Affymetrix.ID, mart=ensembl)
bm2 = getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"), mart=ensembl)
bm = merge(bm1,bm2)
idx = match(winden$Affymetrix.ID, bm$affy_moe430a)
winden$ensg= bm$hsapiens_homolog_ensembl_gene[idx]

table.p = matrix(NA, nrow=length(moduleNames),
                 ncol=length(table(winden$Module.Assigned)))
rownames(table.p) = moduleNames 
nm = na.omit(unique(winden$Module.Assigned))[1:17]
colnames(table.p) = unique(nm)
table.or = table.p

for (m in moduleNames) {
  for(e in colnames(table.p)) {
    f = ORA(datProbes$ensembl_gene_id[mods==m],winden$ensg[winden$Module.Assigned==e],
            datProbes$ensembl_gene_id, winden$ensg)
    table.or[m,e] = as.numeric(f[[1]])
    table.p[m,e] = as.numeric(f[[2]])
  }
}
table.p.fdr = p.adjust(table.p,"fdr")
dim(table.p.fdr) = dim(table.p); dimnames(table.p.fdr) = dimnames(table.p)
table.p.fdr[table.or<1] = 1

#plot 
to_plot = data.frame( Module = rownames(table.p.fdr),
                      Mitochondria="Synaptic",
                      Enrichment = (table.or[,"turquoise"]),
                      FDR = table.p.fdr[,"turquoise"])
to_plot =    rbind(to_plot,
                   data.frame(Module=rownames(table.p.fdr), 
                              Mitochondria="Non-Synaptic",
                              Enrichment=(table.or[,"blue"]),
                              FDR = table.p.fdr[,"blue"]))
to_plot$color = names(mods)[match(to_plot$Module,mods)]
to_plot$Module = factor(to_plot$Module)
to_plot$Mitochondria = factor(to_plot$Mitochondria)
to_plot$symbol=""
to_plot$symbol[to_plot$FDR<0.05]="*"
to_plot$symbol[to_plot$FDR<0.01]="**"
to_plot$symbol[to_plot$FDR<0.001]="***"


mitoPlot = ggplot(to_plot, 
                  aes(x = reorder(Module,Enrichment), 
                      y = Enrichment,
                      fill = Mitochondria,
                      alpha = FDR<0.05)) +
  coord_flip() +
  geom_bar(stat="identity",
           color = "black",
           position=position_dodge(),
           width = 0.7)+ 
  geom_hline(yintercept = 1,lty=2) +
  ylab("Enrichment") + 
  xlab("Module") +
  scale_fill_manual(values = c("black","yellow"))+
  # ggtitle("Mitochondrial Module Enrichment") + 
  theme(plot.title = element_text(hjust=0.5, face = "bold"))+
  theme_bw(base_size = 10)
mitoPlot


ggsave(mitoPlot, file="./results/figures/network_analysis/Module_mitochondria.pdf",
       width = 6, height = 4) 
write.csv(to_plot,file="./results/tables/TableS2-Winden.csv",row.names = F)

save.image("./codes/GO.Rdata")

