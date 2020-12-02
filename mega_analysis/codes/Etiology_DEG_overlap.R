#Etiology_DEGs_overlap.R check using lm data 
options(stringsAsFactors=F)
library(ggplot2); library(mada); library(reshape); library(pheatmap); library(RColorBrewer)
library(NMF); library(WGCNA); library(corrplot); library(purrr); 
library(dplyr); library(biomaRt);library(RRHO);  library(venn)
library(circlize);library(ComplexHeatmap); library(ggrepel); library(ggpubr)


# set wd
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")

#### Load data ####

ad_meta = read.csv("../disease_specific/AD_metaAnalysis/tables/AD_sumstats.csv", row.names=1)
asd_meta = read.csv("../disease_specific/ASD_metaAnalysis/tables/ASD_sumstats.csv", row.names=1)
scz_meta = read.csv("../disease_specific/Scz_metaAnalysis/tables/Scz_sumstats.csv", row.names=1)
bp_meta = read.csv("../disease_specific/BP_metaAnalysis/tables/BP_sumstats.csv", row.names=1)
mdd_meta = read.csv("../disease_specific/MDD_metaAnalysis/tables/MDD_sumstats.csv", row.names=1)
pd_meta = read.csv("../disease_specific/PD_metaAnalysis/tables/pj283498_PD_sumstats.csv", row.names=1)
pd_meta= rename_(pd_meta, "p.value" = "P.Value", "fdr" = "adj.P.Val") 
pa_meta = read.csv("../disease_specific/PA_metaAnalysis/tables/PA_sumstats.csv", row.names=1)
psp_meta = read.csv("../disease_specific/PSP_metaAnalysis/tables/PSP_sumstats.csv", row.names=1)

## union of genes for annotating 
union_genes = 
  union(
    union(
      union(
        union(
          union(
            union(
              union(rownames(ad_meta), 
                    rownames(pd_meta)), 
              rownames(pa_meta)), 
            rownames(psp_meta)),
          rownames(scz_meta)),
        rownames(asd_meta)),
      rownames(bp_meta)),
    rownames(mdd_meta))


all_disorders = vector(mode="list",length=8)
names(all_disorders) <- c("AD","PD","PA", "PSP", "Scz", "ASD", "MDD", "BP")
all_disorders[[1]] = ad_meta
all_disorders[[2]] = pd_meta
all_disorders[[3]] = pa_meta
all_disorders[[4]] = psp_meta
all_disorders[[5]] = scz_meta
all_disorders[[6]] = asd_meta
all_disorders[[7]] = mdd_meta
all_disorders[[8]] = bp_meta

gene_table = as.data.frame(matrix(NA,nrow=length(union_genes), 24))
for(i in 1:8) {
  all_disorders[[i]] = all_disorders[[i]][match(union_genes, rownames(all_disorders[[i]])),]
  gene_table[,3*i-2] = all_disorders[[i]]$logFC
  gene_table[,3*i-1] = all_disorders[[i]]$p.value
  gene_table[,3*i] = all_disorders[[i]]$fdr
}

i = 1
for(dx in c("AD","PD","PA", "PSP", "Scz", "ASD", "MDD", "BP")) {
  for(var in c("logFC", "p.value", "fdr")) {
    colnames(gene_table)[i] = paste(dx,var,sep=".")
    i=i+1
  }
}

rownames(gene_table) = union_genes
write.csv(file="./results/tables/disease_signatures.csv", gene_table)

getinfo <- c("ensembl_gene_id",
             "external_gene_name",
             "gene_biotype")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")
attr <- getBM(  attributes = getinfo,
                filters = c("ensembl_gene_id"),
                values = union_genes,
                mart = mart)
attr = na.omit(attr[match(union_genes, attr$ensembl_gene_id),])
gene_table = gene_table[match(attr$ensembl_gene_id,rownames(gene_table)),]


#### frequency of significant genes ####
alldx = as.data.frame(matrix(NA,ncol = 6))
colnames(alldx) = c("logFC", "p.value","fdr", "disease","id" ,"gene")

for (i in names(all_disorders)){
  dat = as.data.frame(all_disorders[[i]])
  dat = na.omit(dat[,c("logFC", "p.value","fdr")])
  dat$disease = i #add disease name
  dat$id = rownames(dat)
  dat$gene = attr$external_gene_name[match(rownames(dat), attr$ensembl_gene_id)]
  rownames(dat) = NULL
  
  #check if duplicated
  id = !duplicated(dat)
  dat = dat[id,]
  alldx = rbind(alldx,dat)
  
}
alldx = na.omit(alldx)

#reorder based NDD and NPDs
alldx$disease = factor(alldx$disease, levels = c("AD", "PD", "PA", "PSP",
                                                 "Scz", "ASD", "BP", "MDD"))

# intersect of significant genes
sigData = subset(alldx, p.value < 0.05 & abs(logFC) > 0.5)

# check which genes are significant in all diseases
freq_genes = as.data.frame(table(sigData$id)) %>% arrange(desc(Freq))
top_freq_genes = subset(freq_genes,Freq>1)
top_freq_id = top_freq_genes$Var1

# we need logFC for all genes 
logFCs= gene_table %>%  dplyr::select(contains("logFC"))
colnames(logFCs)= gsub(".logFC","",colnames(logFCs))


dat = logFCs[top_freq_id,]
dat[is.na(dat)==T] = 0
dat = dat %>% 
  dplyr::select(AD, PD, PA, PSP,
         Scz, ASD, BP, MDD)
# display only top genes
top = as.character(top_freq_genes$Var1[top_freq_genes$Freq > 4])
top_symb = attr$external_gene_name[match(top,attr$ensembl_gene_id)]
View(alldx[alldx$id %in% top_freq_id,])

#plot the heatmap

library(ComplexHeatmap)


right = rowAnnotation(gene = anno_mark(at =match(top,top_freq_id),
                                       labels = top_symb, 
                                       labels_gp = gpar(fontsize=10)))

pdf("./results/figures/differential_expression/GenesEtiology.pdf",width = 4.5,height =7)
Heatmap(as.matrix(dat),
        col = brewer.pal(11,"BrBG"),
        name = "logFC",
        # cluster_rows = T,
        # row_split = 2,
        cluster_columns = F,
        top_annotation = columnAnnotation(group = rep(c("NDD","NPD"),each=4),
                                          col =list(group = c("NDD" = "blue",
                                                              "NPD" = "black")),
                                          show_annotation_name=F),
        right_annotation = right,
        show_row_names = F,
        show_column_names =T,
        column_names_rot = 45,
        row_labels = top_freq_id,
        border = T,
        height = unit(15,"cm") ,
        width = unit(5,"cm"),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 10),
        show_row_dend = F,
        show_column_dend = F,
        show_heatmap_legend = T,
        heatmap_legend_param = gpar(fontsize=1)) 


dev.off()

