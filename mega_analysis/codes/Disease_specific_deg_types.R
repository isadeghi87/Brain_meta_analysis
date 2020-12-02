# Disease_specific_deg_types.R
# Here we want to see percentage of gene types in 
# differentially expressed genes

options(stringsAsFactors=F)
library(ggplot2); library(RColorBrewer); library(purrr); 
library(dplyr);library(tidyverse);

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

## annotate 
getinfo <- c( "ensembl_gene_id",
              "external_gene_name",
              "gene_biotype")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl") ## Gencode v28
datProbes <- getBM(
  attributes = getinfo,
  filters = "ensembl_gene_id",
  values = union_genes,
  mart = mart)
datProbes <-  datProbes[match(union_genes, datProbes$ensembl_gene_id), ]
datProbes = na.omit(datProbes)
rownames(datProbes) = datProbes$ensembl_gene_id


#### Disease list ####
dis_list = list(ad_meta, pd_meta, pa_meta, psp_meta, scz_meta,asd_meta,bp_meta,mdd_meta)
names(dis_list) = c("AD", "PD", "PA", "PSP", "Scz","ASD", "BP", "MDD")

#--- table of DEGs
degTable = as.data.frame(matrix(ncol = 8,nrow = 4))
rownames(degTable) = c("Total","protein_coding","lncRNA","other")
colnames(degTable) = c("AD","PD","PA","PSP","Scz","ASD","BP","MDD")

for (i in names(dis_list)){
  sig = dis_list[[i]]
  sig = subset(sig, p.value <0.05 & abs(logFC)>0.5)
  id = match(rownames(sig), datProbes$ensembl_gene_id)
  probe = datProbes[id,]
  degTable["Total",i]  = nrow(sig)
  
  ## based on gene type
  for (n in c("protein_coding","lncRNA")){
    n_sig = nrow(probe[probe$gene_biotype == n,])
  degTable[n,i]  = n_sig
    
  }
}

degTable[4,] = degTable[1,] - (degTable[2,] + degTable[3,])
degTable$gene_type = rownames(degTable)

## save table 
write.csv(degTable,"./results/tables/Dis_specific_DEGs.csv")

table.plot = ggtexttable(degTable[,1:8],theme = ttheme(base_style = "mRed"))
ggsave("./results/figures/differential_expression/Disease_DEGs_table.pdf",
       width =5,height = 3)


### a barplot of DEGs ####
# spread data 

datDeg = pivot_longer(degTable[-1,],
             cols = 1:8,
             names_to = "Disease",values_to = "n")
datDeg$gene_type = gsub("protein_coding","protein coding",datDeg$gene_type)
datDeg$gene_type = factor(datDeg$gene_type,
                          levels = c("protein coding","lncRNA","other"))
## bar plot 
bp =  ggplot(as.data.frame(datDeg), aes( y = n, x = reorder(Disease,-n) ,
                          fill = gene_type))+
  geom_bar(width = 0.9, stat = "identity",position = "dodge")+
    scale_fill_brewer(palette = "Set1")+
  geom_text(aes(label = n), position = position_dodge(width = 0.8),vjust = -0.5)+
  labs(x = "", y = "Number of genes", fill = "gene type")+
  theme_clean(base_size = 14)+
  theme(axis.text.x.bottom = element_text(angle = 45,hjust = 0.8))
  
  
ggsave("./results/figures/differential_expression/Disease_DEG_types.pdf",
       width = 9, height =5, plot = bp)


## volcano plot for DEGS ####

alldx = as.data.frame(matrix(NA,ncol = 5))
colnames(alldx) = c("logFC", "p.value","fdr", "disease","id" )
subtext = alldx

for (i in names(dis_list)){
  dat = as.data.frame(dis_list[[i]])
  dat = na.omit(dat[,c("logFC", "p.value","fdr")])
  dat$disease = i #add disease name
  dat$id = rownames(dat)
  rownames(dat) = NULL
  
  #check if duplicated
  alldx = rbind(alldx,dat)
  
  ###filter for significant genes (p<0.05)
sig = subset(dat,p.value <0.05& abs(logFC)>0.5)
    ###order top DEGs 
  subdata = order(abs(sig$logFC), decreasing = T)
  subdata = subdata[1:10]
  sub = sig[subdata,]
  subtext = rbind(subtext,sub)
}

alldx = alldx[-1,]
alldx$gene = datProbes$external_gene_name[match(alldx$id, datProbes$ensembl_gene_id)]
subtext = subtext[-1,]
subtext$gene = datProbes$external_gene_name[match(subtext$id, datProbes$ensembl_gene_id)]

#reorder based NDD and NPDs
alldx$disease = factor(alldx$disease, levels = c("AD", "PD", "PA", "PSP",
                                                 "Scz", "ASD", "BP", "MDD"))
subtext$disease = factor(subtext$disease, levels = c("AD", "PD", "PA", "PSP",
                                                     "Scz", "ASD", "BP", "MDD"))
#color 
sig_up = -log10(alldx$p.value)>1.3 & alldx$logFC >0.5
sig_down = -log10(alldx$p.value)>1.3 & alldx$logFC < - 0.5

col = ifelse( sig_up , "green", ifelse(sig_down,"red","grey"))
top.col = ifelse(subtext$logFC >0, "green", "red")

#### volcano plot ####
volc.p = ggplot(data = alldx , aes(x = logFC,
                                   y = -log10(p.value)))+
  geom_point(size = 0.5, color =col)+
  geom_text_repel(data = subtext,
                  aes(x = logFC,
                      y = -log10(p.value),
                      label = gene), 
                  color = "black",
                  # vjust = ,
                  size = 3)+
  geom_hline(yintercept = -log10(0.05), color = "red", lty = 2)+
  labs(y = "-log10(p-value)", x = expression(~log[2]~FC))+
  facet_wrap(~disease, ncol = 4, scales = "free_y")+
  theme_clean()+
  theme(text = element_text(size = 8 ,colour = "black"),
        strip.text.x = element_text(size = 14))

# save plot
ggsave(plot = volc.p, filename = "./results/figures/differential_expression/Dis_Specific_DEG.pdf",
       width = 12, height = 6)

## disease - deg
sig_genes = subset(alldx, p.value <0.05 & abs(logFC)>0.5)
uni_genes = data.frame(table(sig_genes$id)) %>% subset(Freq == 1)
id = match(uni_genes$Var1,sig_genes$id)
specifc_deg = sig_genes[id,]

# table of specific genes for each disease
spec_table = data.frame(table(specifc_deg$disease))
colnames(spec_table) = c("Disease","Specific DEGs")
g = ggtexttable(t(spec_table),theme = ttheme(base_style = "mBlue",base_size = 10) )
g
# save plot
ggsave("./results/figures/differential_expression/Disease_unique_degs.pdf",
       plot = g, width = 5, height = 1)
