#Step3_Overlap check using lm data 
options(stringsAsFactors=F)
library(ggplot2); library(RColorBrewer); library(purrr); 
library(cowplot);library(clusterProfiler);library(org.Hs.eg.db)
library(gProfileR); library(enrichplot);library(DOSE)
library(pheatmap); library(WGCNA);library(colorspace)

# set wd
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")

#### Load data ####

ad_meta = read.csv("../disease_specific/AD_metaAnalysis/tables/AD_sumstats.csv", row.names=1)
asd_meta = read.csv("../disease_specific/ASD_metaAnalysis/tables/ASD_sumstats.csv", row.names=1)
scz_meta = read.csv("../disease_specific/Scz_metaAnalysis/tables/Scz_sumstats.csv", row.names=1)
bp_meta = read.csv("../disease_specific/BP_metaAnalysis/tables/BP_sumstats.csv", row.names=1)
mdd_meta = read.csv("../disease_specific/MDD_metaAnalysis/tables/MDD_sumstats.csv", row.names=1)
pd_meta = read.csv("../disease_specific/PD_metaAnalysis/tables/pj283498_PD_sumstats.csv", row.names=1)
pa_meta = read.csv("../disease_specific/PA_metaAnalysis/tables/PA_sumstats.csv", row.names=1)
psp_meta = read.csv("../disease_specific/PSP_metaAnalysis/tables/PSP_sumstats.csv", row.names=1)


#### Disease Specific Gene Enrichment ####
dis_list = list(ad_meta, pd_meta, pa_meta, psp_meta,scz_meta,asd_meta, bp_meta, mdd_meta)
names(dis_list) = c("AD", "PD", "PA", "PSP", "Scz","ASD", "BP", "MDD")

#---gene enrichment for each disease usign clusterprofiler

go_plots = vector("list",8)
names(go_plots) = c("AD","PD","PA","PSP","Scz","ASD","BP","MDD")

for (i in names(dis_list)){
  sig = dis_list[[i]]
  sig = subset(sig,p.value <0.05 & abs(logFC)>0.5)
  query = rownames(sig)
  
  # convert to symbol
  ids <- bitr(query, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  ids = ids$ENTREZID
  edo <- enrichDO(ids,)
  ## convert gene ID to Symbol
  edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  p <- cnetplot(edox, categorySize="pvalue",
                colorEdge = T,
                node_label= "category")+
    labs(title = i, size = "-log10 (p-value)")+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          text = element_text(size = 16))


  go_plots[[i]] = p
  
}

plots = ggarrange(plotlist = go_plots, 
          ncol=2, 
          nrow = 4
          )

ggsave("./results/figures/GO/Disease_GO.pdf",
       width = 14, height =12, plot = plots)


### usign gprofiler ####

GO_results= data.frame()
plot_dat = data.frame()

for (i in names(dis_list)){
  sig = dis_list[[i]]
  sig = subset(sig,p.value <0.05 & abs(logFC)>0.5)
  query = rownames(sig)
  go = gost(query, 
                 organism="hsapiens", 
                 ordered_query = F, 
                 significant = T,
                 exclude_iea = F, 
                 # region_query = F,
                 # max_p_value = 1,
            user_threshold = 0.05,
                 correction_method = "fdr",
                 # custom_bg = attr$ensembl_gene_id,
                 # max_set_size = 20000,
                 # hier_filtering = "moderate", 
                 domain_scope  = "annotated", 
                 numeric_ns = "", 
                 # include_graph = F,
                 sources  = c("GO", "KEGG"))
  
  go = as.data.frame(go$result)
  go = go[order(go$p_value),]
  GO_results = rbind(GO_results, cbind(i, go))
  toPlot <- go %>% dplyr::select(term_name, p_value)
  toPlot <- toPlot[1:10,]
  plot_dat <- rbind(plot_dat, cbind(i,toPlot))
}

## plot GO as heatmap
library(tidyr)

# create a matrix for heatmap
matGo = plot_dat  %>% 
  pivot_wider(names_from = i,values_from = p_value)

matGo = GO_results  %>% dplyr::select(term.name,i,p.value)
matGo = subset(matGo, p.value <0.0001)
matGo =  matGo %>% pivot_wider(names_from = i,values_from = p.value)

# shorten the length of term names to 25 char
matGo$term_name = substring(matGo$term_name, 1, 30)
matGo = as.data.frame(matGo)
rownames(matGo)= matGo$term_name
matGo$term_name = NULL
matGo[is.na(matGo) == TRUE] = 1

# heatmap of GO 
pheatmap(-log10(matGo),
         color = rev(sequential_hcl(3000, palette = "Purples")),
         show_rownames = T,
         show_colnames = T,
         angle_col = 45,
         cellwidth = 10,
         cellheight = 7,
         # legend_breaks = c(0,10,20,max(-log10(matGo))),
         # legend_labels = c("0","10","20","-log10(FDR)\n"),
         fontsize = 6,
         fontsize_col = 6,
         fontsize_row = 6,
         treeheight_col = 1,
         treeheight_row = 1,
         cluster_rows = T,
         cluster_cols = T,
         clustering_method = "mcquitty",
         annotation_legend = F,
         annotation_names_col = F,
         filename = "./results/figures/GO/disease_GO_heatmap.pdf",
         width = 4, height = 6)

write.csv(x =  GO_results, file="./results/tables/Disease_GO_results.csv" )

