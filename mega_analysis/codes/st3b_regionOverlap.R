#---check brain lobes transcriptome overlap between diseases
#st3b_regionOverlap.R

rm(list=ls()); options(stringsAsFactors=F)
source("http://bioconductor.org/biocLite.R")
library(ggplot2); library(mada); library(reshape)
library(NMF); library(WGCNA); library(corrplot); library(GGally);library(dplyr)
library(biomaRt); library(gProfileR);library(circlize); library(ggridges)
library(ComplexHeatmap); library(ggplotify);library(gridGraphics);library(grid)
library(colorspace); library(ggExtra);library(ggpubr)

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")
plot = TRUE

#---load Data for lobes
ad_meta = read.csv("../disease_specific/AD_metaAnalysis/tables/AD_lobe_sumstats.csv", row.names=1)
asd_meta = read.csv("../disease_specific/ASD_metaAnalysis/tables/ASD_lobe_sumstats.csv", row.names=1)
scz_meta = read.csv("../disease_specific/Scz_metaAnalysis/tables/Scz_lobe_sumstats.csv", row.names=1)
bp_meta = read.csv("../disease_specific/BP_metaAnalysis/tables/BP_lobe_sumstats.csv", row.names=1)
mdd_meta = read.csv("../disease_specific//MDD_metaAnalysis/tables/MDD_lobe_sumstats.csv", row.names=1)
pd_meta = read.csv("../disease_specific//PD_metaAnalysis/tables/pj283498_PD_sumstats.csv", row.names=1)
colnames(pd_meta) = paste("PD","Frontal",colnames(pd_meta),sep =".")#rename colnames
pa_meta = read.csv("../disease_specific/PA_metaAnalysis/tables/PA_lobe_sumstats.csv", row.names=1)
psp_meta = read.csv("../disease_specific/PSP_metaAnalysis/tables/PSP_lobe_sumstats.csv", row.names=1)

#--keep overlap genes
all_genes = intersect(
  intersect(
    intersect(
      intersect(
        intersect(
          intersect(
            intersect(rownames(ad_meta),
                      rownames(asd_meta)), 
            rownames(scz_meta)), 
          rownames(bp_meta)), 
        rownames(mdd_meta)), 
      rownames(pd_meta)),
    rownames(pa_meta)),
  rownames(psp_meta)) #---12438 genes

#---keep all genes from all disease
union_genes = 
  union(
    union(
      union(
        union(
          union(
            union(
              union(rownames(ad_meta), 
                    rownames(asd_meta)), 
              rownames(scz_meta)), 
            rownames(bp_meta)),
          rownames(mdd_meta)),
        rownames(pd_meta)),
      rownames(pa_meta)),
    rownames(psp_meta))

#---annotate
getinfo <-
  c(
    "ensembl_gene_id",
    "external_gene_name",
    "gene_biotype"
  )
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl") ## Gencode v28
attr <- getBM(
  attributes = getinfo,
  filters = c("ensembl_gene_id"),
  values = union_genes,
  mart = mart
)
attr = attr[match(union_genes, attr$ensembl_gene_id),]


alldata = list(ad_meta,pd_meta,pa_meta,psp_meta,scz_meta,asd_meta, bp_meta,mdd_meta)
names(alldata) = c("AD","PD","PA","PSP","Scz","ASD","BP","MDD")

#gather all data in one df
allExpr = data.frame()
plot = T
if (plot){
  for (dx in names(alldata)){
    for(region in c("Cerebellum", "Temporal",
                    "Frontal", "Limbic","Occipital","Basal.ganglia", "Insular"
                    )){
      dat = alldata[[dx]]
      dat = dat[,grepl(region,colnames(dat))]
      
      if (ncol(dat)> 0){
      dat = dat[,grepl("logFC|p.value|P.Value", colnames(dat)),] %>%
      as.data.frame()
      dat$gene = rownames(dat)
      rownames(dat) = NULL
      dat$Dx = dx
      dat$Region = region
      colnames(dat) = c("logFC", "p.value","gene", "Disease", "Region")
      allExpr = rbind(na.omit(allExpr),dat)
      }
    }
      }
}

allExpr$symbol = attr$external_gene_name[match(allExpr$gene, attr$ensembl_gene_id_version)]
allExpr$Region = factor(allExpr$Region, levels = c("Cerebellum", "Temporal",
                                                   "Frontal", "Limbic","Occipital","Basal.ganglia", 
                                                   "Insular"))
##save df for more analysis
saveRDS(object = allExpr,file = "./codes/regions_allExpr.rds")


#----------------------------------------Basal ganglia 
bas.gang <- list()
for (i in names(alldata)){
  dat = t(alldata[[i]])
  dat = dat[grepl("Basal.*logFC",rownames(dat)),]
  dat = as.data.frame(dat)
  bas.gang[[i]]=dat
  if (nrow(bas.gang[[i]])==0){bas.gang[[i]]=NULL}
}

#intersect of genes
bg_genes <- intersect(rownames(bas.gang[[1]]),
                      intersect(rownames(bas.gang[[2]]),rownames(bas.gang[[3]])))

bglogfc <- as.data.frame(matrix(NA, nrow = length(bg_genes), ncol = length(bas.gang)))
colnames(bglogfc) <- names(bas.gang)

for (i in 1:length(bas.gang)){
  dat <- bas.gang[[i]]
  dat <- dat[match(bg_genes, rownames(dat)),]
  bglogfc[,i] <- dat
}
rm(bas.gang)
rownames(bglogfc)<- bg_genes

#----Cerebellum 
cer <- list()
for (i in names(alldata)){
  dat = t(alldata[[i]])
  dat = dat[grepl("Cerebellum.*logFC",rownames(dat)),]
  dat = as.data.frame(dat)
  cer[[i]]=dat
  if (nrow(cer[[i]])==0){cer[[i]]=NULL}
}

n = length(cer)
cer_genes <- intersect(rownames(cer[[1]]),rownames(cer[[2]]))
for (i in 1:n){cer_genes=intersect(cer_genes,rownames(cer[[i]]))}

cerlogfc <- as.data.frame(matrix(NA, nrow = length(cer_genes), ncol = n))
colnames(cerlogfc) <- names(cer)

for (i in 1:length(cer)){
  dat <- cer[[i]]
  dat <- dat[match(cer_genes, rownames(dat)),]
  cerlogfc[,i] <- dat
}

rm(cer)
rownames(cerlogfc)<- cer_genes

#----------------------------------------------Frontal 
fr = list()
for (i in names(alldata)){
  dat = t(alldata[[i]])
  dat = dat[grepl("Frontal.*logFC",rownames(dat)),]
  dat = as.data.frame(dat)
  fr[[i]]=dat
  if (nrow(fr[[i]])==0){fr[[i]]=NULL}
}
# fr[["PD"]] = alldata$PD
n = length(fr)
fr_genes <- intersect(rownames(fr[[1]]),rownames(fr[[2]]))
for (i in 1:n){fr_genes = intersect(fr_genes, rownames(fr[[i]]))}

frlogfc <- as.data.frame(matrix(NA, nrow = length(fr_genes), ncol = length(fr)))
colnames(frlogfc) <- names(fr)

for (i in 1:length(fr)){
  dat <- fr[[i]]
  dat <- dat[match(fr_genes, rownames(dat)),]
  frlogfc[,i] <- dat
}
rm(fr)
rownames(frlogfc) = fr_genes

#---Insular (MDD)
ins =data.frame(MDD= mdd_meta$Insular.logFC )
rownames(ins) = rownames(mdd_meta)

#-----------------------------------limbic 
limb <- list()
for (i in names(alldata)){
  dat = t(alldata[[i]])
  dat = dat[grepl("Limb.*logFC",rownames(dat)),]
  dat = as.data.frame(dat)
  limb[[i]]=dat
  if (nrow(limb[[i]])==0){limb[[i]]=NULL}
}

n = length(limb)
lmb_genes <- intersect(rownames(limb[[1]]),rownames(limb[[2]]))
for (i in 1:n){lmb_genes = intersect(lmb_genes, rownames(limb[[i]]))}
limblogfc <- as.data.frame(matrix(NA, nrow = length(lmb_genes), ncol = length(limb)))
colnames(limblogfc) <- names(limb)

for (i in 1:length(limb)){
  dat <- limb[[i]]
  dat <- dat[match(lmb_genes, rownames(dat)),]
  limblogfc[,i] <- dat
}
rownames(limblogfc) = lmb_genes



#-------------------------------------Temporal 

tmp <- list()
for (i in names(alldata)){
  dat = t(alldata[[i]])
  dat = dat[grepl("Temp.*logFC",rownames(dat)),]
  dat = as.data.frame(dat)
  tmp[[i]]=dat
  if (nrow(tmp[[i]])==0){tmp[[i]]=NULL}
}
n = length(tmp)
tmp_genes <- intersect(rownames(tmp[[1]]),rownames(tmp[[2]]))
for (i in 1:n){tmp_genes = intersect(tmp_genes, rownames(tmp[[i]]))}
tmplogfc <- as.data.frame(matrix(NA, nrow = length(tmp_genes), ncol = length(tmp)))
colnames(tmplogfc) <- names(tmp)

for (i in 1:length(tmp)){
  dat <- tmp[[i]]
  dat <- dat[match(tmp_genes, rownames(dat)),]
  tmplogfc[,i] <- dat
}
rm(tmp)
rownames(tmplogfc) = tmp_genes

#-------------------------------Occipital (only ASD)

occlog =data.frame(ASD= asd_meta$Occipital.logFC )
rownames(occlog) = rownames(asd_meta)

#list of logfc
logList = list(bglogfc, cerlogfc, frlogfc, limblogfc, tmplogfc,ins,occlog)
names(logList) = c("Basal_ganglia", "Cerebellum", "Frontal", "Limbic", "Temporal", "Insular", "Occipital")

regionGenes = union(rownames(logList[[1]]), rownames(logList[[2]]))
for (i in 3:length(logList)){
  regionGenes = union(regionGenes, rownames(logList[[i]]))
}

Nm = vector()
for (i in names(logList)){
  for (j in names(logList[[i]])){
    newname = paste(j,i,sep = "_")
    Nm= c(Nm, newname)
  }
}
allregions = as.data.frame(matrix(NA,nrow = length(regionGenes), ncol = length(Nm)))
rownames(allregions) = regionGenes;colnames(allregions) = Nm

logList2 = logList
#rename colnames
for (i in names(logList2)){
    colnames(logList2[[i]]) = paste(colnames(logList2[[i]]),i,sep = "_")
}

#merge data
for (nm in names(logList2)){
  dat = logList2[[nm]]
  for (r in rownames(dat)){
    for (c in colnames(dat)){
      allregions[r,c]= dat[r,c]
    }
  }
}

#order based on Dx
allregions = allregions[order(colnames(allregions))]

reg = names(logList)
dx_col = c("AD" = "black",
           "ASD" = "blue",
           "Scz" = "green",
           "BP" = "brown",
           "MDD" = "red",
           "PD" = "salmon",
           "PA" = "turquoise",
           "PSP" = "purple")

col = c(rep("black",3),
        rep("black",3),
        rep("blue",5),
        rep("brown",3),
        rep("red",4),
        rep("turquoise",2),
        "salmon",
        rep("purple",2),
        rep("green",4))
gcor =ggcorr(allregions,method = "complete.obs",
       digits = 1,
       color=col,
       name = "rho score",
       low = "red",
       mid= "white",
       high = "blue",
       label = T,
       label_size = 4,
       label_color = "black",
       label_alpha = T,
        layout.exp = 3,legend.size = 10,
       hjust = 0.9,
       size = 4)+
  labs(x ="", y= "", title = "Correlation of cortical regions across diseases")+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))
gcor
ggsave(filename = "./results/figures/overlaps/AllregionsCorr.pdf", width = 10, height = 10)


###
corrlist = list()
for ( i in c("Basal_ganglia","Cerebellum","Frontal","Limbic",       
             "Temporal")){
gc =ggcorr(logList[[i]],
             method = "pairwise",
             digits = 1,
             # color=col,
             name = expression(paste(rho, " values")),
             low = "red",
             mid= " white",
             high = "blue",
             label = T,
             label_size = 4,
             label_color = "black",
             label_alpha = F,
           legend.size = 7,
           size = 3)+
  labs(x ="", y= "", 
       title = i)+
  theme(plot.title = element_text(size = 8,hjust = 0.5))
corrlist[[i]]= gc 
}

pc = ggarrange(plotlist = corrlist,
               common.legend = T, legend = "bottom")+
  theme_cleveland()
ggsave(plot = pc, filename = "./results/figures/overlaps/RegionsCorrelation.pdf",
       width = 8, height = 5.5)

#------------------------------------Gene enrichment 
GO_results= data.frame()
plot_dat = data.frame()

for (i in names(logList)){
query = rownames(logList[[i]])
go = gprofiler(query, 
               organism="hsapiens", 
               ordered_query = F, 
               significant = T,
               exclude_iea = F, 
               region_query = F,
               max_p_value = 1,
               correction_method = "fdr", 
               custom_bg = attr$ensembl_gene_id,
               max_set_size = 2000,
               hier_filtering = "moderate", 
               domain_size = "annotated", 
               numeric_ns = "", 
               include_graph = F,
               src_filter = c("GO", "KEGG"))

go = go[order(go$p.value),]
GO_results = rbind(GO_results, cbind(i, go))
toPlot <- go %>% dplyr::select(term.name, p.value)
toPlot <- toPlot[1:2,]
plot_dat <- rbind(plot_dat, cbind(i,toPlot))
}
write.csv(GO_results, file="./results/tables/BrainLobesGeneEnrichment.csv")


#---plot for each region
plotGO = data.frame()
for(lobe in names(logList)) {
  idx = which(GO_results$i==lobe) 
  plotGO = rbind(GO_results[idx[10:1],],plotGO)  
}

p <-ggplot(data =plotGO,
       aes(x = reorder(term.name, -p.value),
           y =-log10(p.value),
           fill = i))+
  facet_wrap(~i)+
  geom_bar(stat = "identity", width = 0.8, position = "dodge")+
  coord_flip()+
  geom_abline(slope = 0, intercept = -log10(0.01),lty = 2)+
  scale_color_brewer(palette = "Set1")+
  labs(x = "", 
       y = "-log10(FDR)",
       title = "Brain lobes pathway enrichment")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        text = element_text(size = 14),
        legend.position = "none")

p
ggsave(filename = "./results/figures/BrainLobes_PathwayEnrichment.pdf", 
       height = 12,
       width = 14, device = "pdf", plot = p)

#plot it as a 2D 
golobe =ggplot(plotGO, aes(y = reorder(term.name,-log10(p.value)),
                         x= i,
                         color = i))+
  geom_point(aes(size =-log10(p.value)))+
  scale_color_brewer(type = "div", palette = "Set1")+
  labs(x = "", 
       y = "",
       title = "Brain Regions Pathway Enrichment")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        text = element_text(size = 14),
          axis.text.x.bottom = element_text(angle = 45, vjust = 0.9, hjust = 1)
  )+
  guides(color = guide_legend(title = "Lobe"))
golobe

ggsave(filename = "./results/figures/BrainlobesPathwayEnrichment.pdf", 
       height = 6,
       width = 9, 
       device = "pdf", 
       plot = golobe)



#--Top 2 pathways per region
colnames(plot_dat) = c("Brain.lobe", "term.name", "p.value")
p <- ggplot(data = plot_dat,
            aes(x = reorder(Brain.lobe, -p.value),
                y =-log10(p.value),
                fill = term.name,
                group = term.name))+
  geom_bar(stat = "identity", width = 0.5, position = "dodge")+
  coord_flip()+
  geom_abline(slope = 0, intercept = -log10(0.01),lty = 2)+
  scale_fill_manual(values=sort(unique(plot_dat$Brain.lobe)))+
  labs(x = "", 
       y = "-log10(FDR)",
       title = "Top pathway enrichment for brain lobes ")+
  guides(fill  = guide_legend(title = "Term name"))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        text = element_text(size = 14))
p

ggsave(filename = "./results/figures/BrainLobes_topGO.pdf", 
       height = 3,
       width = 6, device = "pdf", plot = p)

#---------------------------------------------RRHO
library(RRHO)
jet.colors  <- colorRampPalette(
  c("#00007F", "blue", "#007FFF", "cyan", 
    "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

#Basal ganglia
rrho.list <- list()
dat = logList[["Basal_ganglia"]]
for ( i in 1:ncol(dat)){
  df <- data.frame(gene = rownames(dat),
                   logfc = dat[,i])
  rrho.list[[i]] = df
}

names(rrho.list) <- colnames(dat)
comparison <- t(combn(seq(1,ncol(dat)),2))#Scz-BP-MDD
mat <- as.matrix(cbind(c(0,3),
                       c(1:2)))
# pdf("./results/figures/RRHO_BasalGanglia.pdf",
#     height = 3, width = 3)
png("./results/figures/RRHO_BasalGanglia.png",
    height = 3, width = 3,units = "cm", res = 300)
layout(mat)
#layout.show(28)
par(mar = c(0.2, 0.2, 0, 0))
for (i in 1:nrow(comparison)){  
R.obj <- RRHO(list1 = rrho.list[[1]], 
                list2 = rrho.list[[2]],
                BY = T,
                alternative = "enrichment",
                plots = T, 
                log10.ind = T,
                outputdir = "./results/figures/RRHO/lobe/basal_ganglia/", 
                labels = c(names(rrho.list[1]), names(rrho.list[2]))
  )
image(R.obj$hypermat,
        xlab =names(rrho.list[1]),
        ylab=names(rrho.list[2]),
        col=jet.colors(100), 
        axes=FALSE, main="")
}
dev.off()


#--Cerebellum
rrho.list <- list()
dat = logList[["Cerebellum"]]
for ( i in 1:ncol(dat)){
  df <- data.frame(gene = rownames(dat),
                   logfc = dat[,i])
  rrho.list[[i]] = df
}

names(rrho.list) <- colnames(dat)
comparison <- t(combn(seq(1,ncol(dat)),2))#AD-ASD-PA-PSP
mat <- as.matrix(cbind(c(0,0,6),
                       c(0,4,5),
                       c(1:3)))


# pdf("./results/figures/RRHO_Cerebllum.pdf",
#     height = 3, width = 3)
png("./results/figures/RRHO_Cerebellum.png",
    height = 3, width = 3,units = "cm", res = 300)

layout(mat)
#layout.show(28)
par(mar = c(0.2, 0.2, 0, 0))
for (i in 1:nrow(comparison)){
  list1 = comparison[i,1]
  list2 = comparison[i,2]
  R.obj <- RRHO(list1 = rrho.list[[list1]], 
                list2 = rrho.list[[list2]],
                BY = T,
                alternative = "enrichment",
                plots = T, 
                log10.ind = T,
                outputdir = "./results/figures/RRHO/lobe/cerebellum/", 
                labels = c(names(rrho.list[list1]), names(rrho.list[list2]))
  )
  image(R.obj$hypermat,
        xlab = ifelse(i %in% c(1,8,14,19,23,26,28),
                      names(rrho.list[list1]),""),
        ylab=ifelse(i %in% c(1,8,14,19,23,26,28), names(rrho.list[list2]), ""),
        col=jet.colors(100), 
        axes=FALSE, main="")
}
dev.off()


#---Frontal lobe
rrho.list <- list()
dat = logList[["Frontal"]]
for ( i in 1:ncol(dat)){
  df <- data.frame(gene = rownames(dat),
                   logfc = dat[,i])
  rrho.list[[i]] = df
}

names(rrho.list) <- colnames(dat)
comparison <- t(combn(seq(1,ncol(dat)),2))#AD-ASD-Scz-BP-MDD-PD
mat <- as.matrix(cbind(c(rep(0,4),15),
                       c(rep(0,3),13:14),
             c(rep(0,2),10:12),
             c(0,6:9),
             c(1:5)))



# pdf("./results/figures/RRHO_Frontal.pdf",
#     height = 3, width = 3)

png("./results/figures/RRHO_Frontal.png",
    height = 3, width = 3,units = "cm", res = 300)
layout(mat)
#layout.show(15)
par(mar = c(0.2, 0.2, 0, 0))
for (i in 1:nrow(comparison)){
  list1 = comparison[i,1]
  list2 = comparison[i,2]
  R.obj <- RRHO(list1 = rrho.list[[list1]], 
                list2 = rrho.list[[list2]],
                BY = T,
                alternative = "enrichment",
                plots = T, 
                log10.ind = T,
                outputdir = "./results/figures/RRHO/lobe/frontal/", 
                labels = c(names(rrho.list[list1]), names(rrho.list[list2]))
  )
  image(R.obj$hypermat,
        xlab = names(rrho.list[list1]),
        ylab=names(rrho.list[list2]),
        col=jet.colors(100), 
        axes=FALSE, main="")
}
dev.off()



#Limbic
rrho.list <- list()
dat = logList[["Limbic"]]
for ( i in 1:ncol(dat)){
  df <- data.frame(gene = rownames(dat),
                   logfc = dat[,i])
  rrho.list[[i]] = df
}

names(rrho.list) <- colnames(dat)
comparison <- t(combn(seq(1,ncol(dat)),2))#ASD-Scz-BP-MDD
mat <- as.matrix(cbind(c(0,0,6),
             c(0,4,5),
             c(1:3)))


# pdf("./results/figures/RRHO_Limbic.pdf",
#     height = 3, width = 3)
png("./results/figures/RRHO_Limbic.png",
    height = 3, width = 3,units = "cm", res = 300)

layout(mat)
#layout.show(6)
par(mar = c(0.2, 0.2, 0, 0))
for (i in 1:nrow(comparison)){
  list1 = comparison[i,1]
  list2 = comparison[i,2]
  R.obj <- RRHO(list1 = rrho.list[[list1]], 
                list2 = rrho.list[[list2]],
                BY = T,
                alternative = "enrichment",
                plots = T, 
                log10.ind = T,
                outputdir = "./results/figures/RRHO/lobe/limbic/", 
                labels = c(names(rrho.list[list1]), names(rrho.list[list2]))
  )
  image(R.obj$hypermat,
        xlab = ifelse(i %in% c(1,4,6),
                      names(rrho.list[list1]),""),
        ylab=ifelse(i %in% c(1,4,6), names(rrho.list[list2]), ""),
        col=jet.colors(100), 
        axes=FALSE, main="")
}
dev.off()


#Temporal
rrho.list <- list()
dat = logList[["Temporal"]]
for ( i in 1:ncol(dat)){
  df <- data.frame(gene = rownames(dat),
                   logfc = dat[,i])
  rrho.list[[i]] = df
}

names(rrho.list) <- colnames(dat)
comparison <- t(combn(seq(1,ncol(dat)),2))#AD-ASD-Scz-PA-PSP
mat <- as.matrix(cbind(c(rep(0,3),10),
                       c(rep(0,2),8:9),
                       c(0,5:7),
                       c(1:4)))

# pdf("./results/figures/RRHO_Temporal.pdf",
#     height = 3, width = 3)
png("./results/figures/RRHO_Temporal.png",
    height = 3, width = 3,units = "cm", res = 300)

layout(mat)
#layout.show()
par(mar = c(0.2, 0.2, 0, 0))
for (i in 1:nrow(comparison)){
  list1 = comparison[i,1]
  list2 = comparison[i,2]
  R.obj <- RRHO(list1 = rrho.list[[list1]], 
                list2 = rrho.list[[list2]],
                BY = T,
                alternative = "enrichment",
                plots = T, 
                log10.ind = T,
                outputdir = "./results/figures/RRHO/lobe/temporal/", 
                labels = c(names(rrho.list[list1]), names(rrho.list[list2]))
  )
  image(R.obj$hypermat,
        xlab = ifelse(i %in% c(1,6,10,13,15),
                      names(rrho.list[list1]),""),
        ylab=ifelse(i %in% c(1,6,10,13,15), names(rrho.list[list2]), ""),
        col=jet.colors(100), 
        axes=FALSE, main="")
}
dev.off()


##circplot


#meging cor data from all regions as data.frame

df = data.frame()
for (i in names(logList)){
  cordat = cor(logList[[i]], use = "pairwise.complete.obs")
  diag(cordat)=0
  data = data.frame(from = rep(rownames(cordat), times = ncol(cordat)),
                    to = rep(colnames(cordat), each = nrow(cordat)),
                    value= as.vector(cordat),
                    stringsAsFactors = FALSE)
  data$region = i
  df = rbind(df, data)
}

#remove self correlation
idx = which(df$from == df$to)
idx = c(idx, which(duplicated(df)))
df= df[-idx,]
df=df[order(df$region),]

#remove trivial values
df = df[!duplicated(df$value),]
df = df[df$value>0.1 | df$value< -0.1,]
df = df %>% group_by(from,region) %>% arrange()

#set colorss
col_fun = colorRamp2(c(min(df$value),0,max(df$value)), c("red","white","blue"))
region_col = c("Cerebellum" = "black",
               "Basal_ganglia" = "blue",
               "Temporal"= "red",
               "Frontal" = "green",
               "Limbic" = "brown")



pdf("./results/figures/overlaps/CorrelationCircos.pdf", width = 6, height = 7)
par(mar=c(3,2,0,0))
circos.par(start.degree = -1, gap.after =10,
           cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)
cdm = chordDiagram(df[1:3],
                   grid.col = dx_col, 
             col = col_fun,
             link.lwd = 1, 
             link.lty = 1, 
             link.border = "black", 
             # link.sort = T,
             link.decreasing = T,
             annotationTrack = c("grid", "name"),
             preAllocateTracks = list(track.height=0.01))
# title("Correlation of cortical subregions across diseases")

for(i in seq_len(nrow(cdm))) {
    circos.rect(cdm[i, "x1"], -uy(1, "mm"), 
                cdm[i, "x1"] - abs(cdm[i, "value1"]), -uy(2, "mm"), 
                col = region_col[df$region[i]], border =region_col[df$region[i]] ,
                sector.index = cdm$rn[i], track.index = 3)
    circos.rect(cdm[i, "x2"], -uy(1, "mm"),
                cdm[i, "x2"] - abs(cdm[i, "value2"]), -uy(2, "mm"),
                col = region_col[df$region[i]], border =region_col[df$region[i]] ,
                sector.index = cdm$cn[i], track.index = 3)
}


#add legend
lgd_links = Legend(at = c( -1, 0, 1), col_fun = col_fun, 
                   title_position = "topleft",
                   title = "rho score",
                   direction = "horizontal")
lgd_region = Legend(at = c("Cerebellum",
                           "Basal ganglia",
                           "Temporal",
                           "Frontal" ,
                           "Limbic"),
                    border = c("black", "blue","red","green","brown"),
                    legend_gp = gpar(fill= c("black", "blue","red","green","brown")),
                    title_position = "topleft",
                    title = "Regions", direction = "horizontal")

lgd_list_vertical = packLegend(lgd_links, lgd_region,direction = "horizontal")
draw(lgd_list_vertical, x = unit(2.5, "cm"), y = unit(2.5, "cm"), just = "topleft")

dev.off()



save.image("./codes/st3b_regionOverlap.Rdata")
