#---check brain lobes transcriptome overlap between diseases
#st3b_regionOverlap.R

rm(list=ls()); options(stringsAsFactors=F)
source("http://bioconductor.org/biocLite.R")
library(ggplot2); library(mada); library(reshape)
library(NMF); library(WGCNA); library(corrplot); library(GGally);library(dplyr)
library(biomaRt); library(gProfileR);library(circlize); library(ggridges)
library(ComplexHeatmap); library(ggplotify);library(gridGraphics);library(grid)
library(colorspace); library(ggExtra);library(ggpubr)

setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")
plot = TRUE

#---load Data for lobes
ad_meta = read.csv("../disease_specific/AD_metaAnalysis/tables/AD_lobe_sumstats.csv", row.names=1)
asd_meta = read.csv("../disease_specific/ASD_metaAnalysis/tables/ASD_lobe_sumstats.csv", row.names=1)
scz_meta = read.csv("../disease_specific/Scz_metaAnalysis/tables/Scz_lobe_sumstats.csv", row.names=1)
bp_meta = read.csv("../disease_specific/BP_metaAnalysis/tables/BP_lobe_sumstats.csv", row.names=1)
mdd_meta = read.csv("../disease_specific//MDD_metaAnalysis/tables/MDD_lobe_sumstats.csv", row.names=1)
pd_meta = read.csv("../disease_specific//PD_metaAnalysis/tables/PD_meta_sumstats.csv", row.names=1)
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
    "ensembl_gene_id_version",
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "band",
    "gene_biotype"
  )
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
attr <- getBM(
  attributes = getinfo,
  filters = c("ensembl_gene_id_version"),
  values = union_genes,
  mart = mart
)
attr = attr[match(union_genes, attr$ensembl_gene_id_version),]


alldata = list(ad_meta,asd_meta,scz_meta,bp_meta,mdd_meta,pd_meta,pa_meta,psp_meta)
names(alldata) = c("AD","ASD","Scz","BP","MDD","PD","PA","PSP")

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
#density of logFC for regions
library(forcats)
gp = ggplot(data = allExpr, 
       aes(x = logFC, y= fct_rev(Disease), fill = Region))+
  geom_density_ridges(color = "white", rel_min_height=0.001)+
  geom_vline(xintercept = 0, lty = 2 , color = "red", size = 0.1)+
  facet_wrap(~Region, nrow = 1)+
  theme_ridges(center_axis_labels = T, grid = F)+
  theme_economist()+
  theme(text = element_text(size = 8, colour = "black"),
        legend.position = "none",
        axis.text = element_text(size = 8))+
     labs(y = "", x = "logFC")+
  scale_x_continuous(breaks = c(-1,0,1), limits = c(-1,1))

saveRDS(gp, file = "./codes/fig3c_Densityplot.rds")

ggsave(filename = "./results/figures/RegionsDensity.pdf", plot = gp,
       width = 8, height = 4)

#################################
#####Filter significant genes####
#################################
#A table of DEGs
DEGs = as.data.frame(matrix(NA, nrow = 8, ncol = 7))
rownames(DEGs) = levels(as.factor(unique(allExpr$Disease))); colnames(DEGs) = unique(allExpr$Region)

#### top DEGs####

subdat = data.frame()
top2 = data.frame()
for (dx in unique(allExpr$Disease)){
  for(region in unique(allExpr$Region)){
    dat = subset(allExpr, Disease == dx & Region ==region)
    dat = subset(dat , p.value < 0.05 & abs(logFC)> 0.5)#significant
    DEGs[dx,region] = nrow(dat) #number of DEGs
    idx = order(abs(dat$logFC), decreasing = T)[1:10]#top DEGs
    id = order(abs(dat$logFC), decreasing = T)[1:2]
    sub = dat[idx,]
    top = dat[id,]
    subdat = na.omit(rbind(subdat,sub))
    top2 = na.omit(rbind(top2,top))
      }
}

#save DEG table

DEGs[DEGs == 0] = "-"
write.table(x = DEGs, file = "./results/tables/Regions_DEGs.txt",row.names = T)

#ggtable

DEG.table = ggtexttable(DEGs,theme = ttheme("mBlue",base_size = 7))
saveRDS(DEG.table,"./codes/fig3b_DEGs.rds")
ggsave(plot = DEG.table, 
       file = "./results/figures/Regions_DEGs.pdf",
       width = 5,
       height = 3)

####volcano plot####

library(ggrepel)
pointCol = ifelse(-log10(allExpr$p.value)>1.3 & abs(allExpr$logFC)> 0.5 , "black", "grey")
g = ggplot(data = as.data.frame(allExpr),
       aes(x = logFC, y = -log10(p.value)))+
  geom_point(color = pointCol , size = 0.2)+
  geom_text_repel(data = top2,
                   aes(label = symbol),
                    color = ifelse(top2$logFC >0, "green", "indianred3"),
                   size = 3, 
                  show.legend = F)+
  facet_grid(Region~Disease,
             scales = "free_y",
             )+
  theme_clean()+
  theme(text = element_text(colour = "black", 
                            size = 12),
        strip.text.y  = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.key = element_rect(size = 6))

ggsave(plot= g , file = "./results/figures/RegionsDEGs.pdf",
       width = 10, height = 8)


#### word cloud to find the most frequent genes####
library(ggwordcloud)
freq = as.data.frame(table(subdat$symbol)) %>%  arrange(desc(Freq))
colnames(freq) = c("gene", "Frequency")

#remove frequency <3
freq = subset(freq, Frequency >2)
# 
# set.seed(122)
# g = ggplot(data = freq, aes(label = gene, 
#                         size = Frequency,
#                         color = Frequency))+
#   geom_text_wordcloud_area(shape = "cardioid" ,
#                            perc_step = 0.1,
#                            grid_margin = 2,
#                            show.legend = T )+
#   scale_size_area(max_size = 15)+
#   theme_minimal()+
#   scale_color_gradient(name = "Frequency", 
#                        labels = c("2", "3", "4"),
#                        breaks = c(2,3,4),
#                        low = "orange",
#                        space = "green",
#                        high = "blue")+
#   guides(size = FALSE)+
#   ggtitle("")+
#   theme(legend.key.size = unit(2,"mm"))


# ggsave(filename = "./results/figures/Regions_wordcloud.png",
#        width = 4, height = 4, plot = g, dpi = 300)

#filter only frequent genes

rownames(subdat)= NULL
id = unique(as.character(freq$gene))
top.freq = subdat[subdat$symbol %in% id,]
top.freq$frequency = top.freq$p.value
top.freq$frequency = freq$Frequency[match(top.freq$symbol,
                                            freq$gene)]
top.freq = top.freq[order(top.freq$frequency,decreasing = T),]

# plot

freq.plot = ggplot(top.freq, aes(y = reorder(symbol, frequency), x = Region))+
  geom_point(aes(color = logFC),size = 4)+
  # geom_text_repel(parse = T,aes(label =, symbol,
  #                     #size = frequency,
  #                     color = logFC))+
  # scale_size_continuous(name = "Frequency",
  #                      labels = c( "3", "4", "5"),
  #                      breaks = c(3,4,5))+
  scale_color_continuous_diverging(palette = "Red-Green",
                                   rev = F)+
  facet_wrap(~Disease, nrow = 1)+
  labs(x = "",
       y = "")+
 theme_pubclean(base_size = 12)+
  theme(text = element_text(size = 12),
        axis.text.x.bottom = element_text(size = 7, angle = 45, hjust = 1))


  
  ggsave(plot = freq.plot, 
         filename = "./results/figures/RegionsWordcloud2.pdf",
         width = 10,
         height = 4)
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
rownames(bglogfc)<- attr$ensembl_gene_id[match(bg_genes,attr$ensembl_gene_id_version)]


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
rownames(cerlogfc)<- attr$ensembl_gene_id[match(cer_genes,attr$ensembl_gene_id_version)]


#----------------------------------------------Frontal 
fr = list()
for (i in names(alldata)){
  dat = t(alldata[[i]])
  dat = dat[grepl("Frontal.*logFC",rownames(dat)),]
  dat = as.data.frame(dat)
  fr[[i]]=dat
  if (nrow(fr[[i]])==0){fr[[i]]=NULL}
}
fr[["PD"]] = alldata$PD
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
rownames(frlogfc)<- attr$ensembl_gene_id[match(fr_genes,attr$ensembl_gene_id_version)]


#---Insular (MDD)
ins =data.frame(MDD= mdd_meta$Insular.logFC )
rownames(ins) = attr$ensembl_gene_id[match(rownames(mdd_meta), attr$ensembl_gene_id_version)]

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
rownames(limblogfc) <- attr$ensembl_gene_id[match(lmb_genes, attr$ensembl_gene_id_version)]




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
rownames(tmplogfc) = attr$ensembl_gene_id[match(tmp_genes, attr$ensembl_gene_id_version)]


#-------------------------------Occipital (only ASD)

occlog =data.frame(ASD= asd_meta$Occipital.logFC )
rownames(occlog) = attr$ensembl_gene_id[match(rownames(asd_meta), attr$ensembl_gene_id_version)]

#-------------------------------------------plot scatterplot matrix 
#custom scatter corr plot
my_custom_cor_color <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {
  # assign('data',data,envir = .GlobalEnv)
  # get the x and y data to use the other code
  x <- eval_data_col(data,mapping$x)
  y <- eval_data_col(data,mapping$y)
  
  ct <- cor.test(x,y)
  
  r <- unname(ct$estimate)
  rt <- format(r, digits=1)[1]
  tt <- as.character(rt)
  
  # plot the cor value
  p <- ggally_text(
    label = tt, 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = 7,
    color=color,
    ...
  ) +
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    )
  
  corColors <- RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")[2:10]
  
  if (r <= -0.4) {
    corCol <- corColors[1]
  } else if (r <= -0.3) {
    corCol <- corColors[2]
  }  else if (r <= -0.2) {
    corCol <- corColors[3]
  } else if (r <= 0) {
    corCol <- corColors[5]
  }  else if (r <= 0.2) {
    corCol <- corColors[6]
  } else if (r <= 0.4) {
    corCol <- corColors[7]
  } else if(r<= 0.6){
    corCol <- corColors[8]
  } else {
    corCol <- corColors[9]
  }
  p <- p + theme(
    panel.background = element_rect(fill = corCol)
  )
  
  p
}

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

#custom 
my_density <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) + 
    geom_density(alpha = 0.5,
                 fill =  "green"  , ...)
}

# custom function for scatterplot
my_scatter <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) + 
    geom_point(alpha = 0.5,
               color = "black") + 
    geom_smooth(method= lm,
                color = "red", 
                se=F, ...)
}

allplot =ggpairs(allregions,
        # lower=list(continuous = my_scatter),
        # diag = list(continuous = my_density),
        upper = list(continuous = my_custom_cor_color))+
  labs(title ="All regions comparison") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        # text = element_text(size = 9, face = "bold"),
        axis.text.x.bottom = element_text(size = 6),
        axis.text.y.left = element_text(size = 6),
        strip.text.x = element_text(size = 5),
        strip.text.y = element_text(size = 5))
# print(p)
ggsave(filename = paste("./results/figures/",names(logList[i]),"_TxOverlap.pdf",sep = ""),
       height =4,
       width =6,
       device = "pdf",
       plot = p)
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
ggsave(filename = "./results/figures/AllregionsCorr.pdf", width = 10, height = 10)

corrplot(corr = cor(allregions,
                    use = "pairwise.complete.obs"),
         method = "square",
         type = "upper",tl.col = col,
         # col = colorRamp2(breaks = c(-1,0,1),
         #                  colors = c("red","white","blue")),
         diag = T )
#plot for each region individually
for ( i in 1:length(logList)){
  my_density <- function(data, mapping, ...){
    ggplot(data = data, mapping = mapping) + 
      geom_density(alpha = 0.5,
                   fill =  i  , ...)
  }
  
  # custom function for scatterplot
  my_scatter <- function(data, mapping, ...){
    ggplot(data = data, mapping = mapping) + 
      geom_point(alpha = 0.5,
                 color = "black") + 
      geom_smooth(method= lm,
                  color = "red", 
                  se=F, ...)
  }
  
  # create scatterplot matrix
  p<- ggpairs(logList[[i]],
              lower=list(continuous = my_scatter),
              diag = list(continuous = my_density),
              upper = list(continuous = my_custom_cor_color))+
    labs(title =names(logList[i])) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          # text = element_text(size = 9, face = "bold"),
          axis.text.x.bottom = element_text(size = 6),
          axis.text.y.left = element_text(size = 6),
          strip.text.x = element_text(size = 17),
          strip.text.y = element_text(size = 15))
   # print(p)
  ggsave(filename = paste("./results/figures/",names(logList[i]),"_TxOverlap.pdf",sep = ""),
        height =4,
        width =6,
        device = "pdf",
        plot = p)
}

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
ggsave(plot = pc, filename = "./results/figures/RegionsCorrelation.pdf",
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
df = df[df$value>0.2 | df$value< -0.2,]

#set colorss
col_fun = colorRamp2(c(min(df$value),0,max(df$value)), c("red","white","blue"))
region_col = c("Cerebellum" = "black",
               "Basal_ganglia" = "blue",
               "Temporal"= "red",
               "Frontal" = "green",
               "Limbic" = "brown")


dev.off()

pdf("./results/figures/CorrelationCircos.pdf", width = 6, height = 7)
par(mar=c(3,2,4,1))
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
title("Correlation of cortical subregions across diseases")

# circos.track(track.index = 2, panel.fun = function(x, y) {
#     circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, col = "white",
#                 font = 1, cex = 1, adj = c(0.5, 0.5), niceFacing = TRUE)
#     xlim = CELL_META$xlim
#      # circos.axis(major.at = cdm$x1, 
#      #             labels = df$region,
#      #             labels.cex = 0.5,
#      #             labels.facing = "clockwise")
# }, bg.border = NA)

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
                   title_position = "topleft", title = "rho score")
lgd_region = Legend(at = c("Cerebellum",
                           "Basal_ganglia",
                           "Temporal",
                           "Frontal" ,
                           "Limbic"),border = c("black", "blue","red","green","brown"),
                                 legend_gp = gpar(fill= c("black", "blue","red","green","brown")),
                    title_position = "topleft",
                                 title = "Regions")

lgd_list_vertical = packLegend(lgd_links, lgd_region)

grid.echo()
ribon = as.ggplot(grid.grab())
dev.off()

draw(lgd_list_vertical, x = unit(2.5, "cm"), y = unit(2.5, "cm"), just = "topleft")
grid.echo()
lgd = as.ggplot(grid.grab())

fig2d = ribon + lgd


dev.off()



save.image("./codes/st3b_regionOverlap.Rdata")
