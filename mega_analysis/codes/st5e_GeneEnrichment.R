#st5e_GeneEnrichment.R

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

#Write module genes for GO elite
for(c in unique(mods)) {
  mod = datProbes$ensembl_gene_id[mods==c]
  write.table(file = paste("./results/tables/GeneSet/module_", c, ".txt", sep=""),
              data.frame(Gene=mod, Id="En"), row.names = F, quote=F)
}
write.table(file="./results/tables/GeneSet/background.txt", 
            data.frame(Gene=datProbes$ensembl_gene_id, Id="En"), row.names=F, quote=F)



## Calculate GO enrichment of Disease-associated modules
resultsGO = data.frame()
ModulePlots <- list()
plotGO = data.frame()
for (m in unique(mods)) {
  i = which(m == unique(mods))  
  
  #-----------------------------------------GO Enrichment
  query = rownames(datExpr)[mods==m]
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
  plotdat <- data.frame(term.name=go$term.name[1:5], 
                        p.value=go$p.value[1:5])
  plotdat <- na.omit(plotdat)
  plotdat$Module = m
  plotGO= rbind(plotGO, plotdat)
 #--plot
  p <- ggplot(data = plotdat,
         aes(x = reorder(term.name, -p.value),
                            y =-log10(p.value)))+
    geom_bar(stat = "identity", fill= unique(names(mods)[mods ==m]))+
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

ModulePlots[["M0"]]=NULL
#  data for plot
plotGO = plotGO[!grepl("M0",plotGO$Module),]
# arrange plots together
GOplots <- ggarrange(plotlist = ModulePlots, 
                     ncol = 2)
GOplots
ggsave(plot = GOplots, 
       filename = "./results/figures/ModulesGO.pdf", 
       height = 14, width =24, device = "pdf")

write.csv(file="./results/tables/GO_gprofilerResults.csv", resultsGO)


## plot GO as heatmap
library(tidyr)

# create a matrix for heatmap
matGo = plotGO  %>% 
  pivot_wider(names_from = Module,values_from = p.value)

# shorten the length of term names to 25 char
matGo$term.name = substring(matGo$term.name, 1, 30)
matGo = as.data.frame(matGo)
rownames(matGo)= matGo$term.name
matGo$term.name = NULL
matGo[is.na(matGo) == TRUE] = 1

# annotate colors
annot_mod = data.frame(module = colnames(matGo))
rownames(annot_mod)= colnames(matGo)
mod_col = unique(names(mods))[-2]
names(mod_col) = unique(mods)[-2] 
annot_col = list(module = mod_col)

matGO = orderMEs(matGo)
# heatmap of GO 
pheatmap(-log10(matGo),
         color = redblue(1000)[500:0],
         show_rownames = T,
         show_colnames = T,
         angle_col = 45,
         cellwidth = 10,
         cellheight = 7,
         legend_breaks = c(0,20,40,60,max(-log10(matGo))),
         legend_labels = c("0","20","40","60","-log10(FDR)\n"),
         fontsize = 6,
         fontsize_col = 6,
         fontsize_row = 6,
         treeheight_col = 1,
         treeheight_row = 1,
         cluster_rows = T,
         cluster_cols = F,
         # clustering_method = "ward.D",
         annotation_col = annot_mod,
         annotation_colors = annot_col,
         annotation_legend = F,
         annotation_names_col = F,
         filename = "./results/figures/network_analysis/GO_heatmap_modules.pdf",
         width = 5,height = 7)

