## Regions specific DEGs

library(ggplot2); library(RColorBrewer); library(purrr); 
library(ggridges);library(ggthemes)


# set wd
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")

# options 
options(stringsAsFactors=F)
# colors
cols <- c("aquamarine","red","blue","chartreuse",
          "brown","blueviolet","yellow","orange")
names(cols) = c("AD","PD","PA","PSP","Scz","ASD","BP","MDD")

#---load Data for lobes
ad_meta = read.csv("../disease_specific/AD_metaAnalysis/tables/AD_lobe_sumstats.csv", row.names=1)
asd_meta = read.csv("../disease_specific/ASD_metaAnalysis/tables/ASD_lobe_sumstats.csv", row.names=1)
scz_meta = read.csv("../disease_specific/Scz_metaAnalysis/tables/Scz_lobe_sumstats.csv", row.names=1)
bp_meta = read.csv("../disease_specific/BP_metaAnalysis/tables/BP_lobe_sumstats.csv", row.names=1)
mdd_meta = read.csv("../disease_specific//MDD_metaAnalysis/tables/MDD_lobe_sumstats.csv", row.names=1)
pd_meta = read.csv("../disease_specific/PD_metaAnalysis/tables/pj283498_PD_sumstats.csv", row.names=1)
colnames(pd_meta) = paste("PD","Frontal",colnames(pd_meta),sep =".") #rename colnames
pa_meta = read.csv("../disease_specific/PA_metaAnalysis/tables/PA_lobe_sumstats.csv", row.names=1)
psp_meta = read.csv("../disease_specific/PSP_metaAnalysis/tables/PSP_lobe_sumstats.csv", row.names=1)

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

getinfo <- c("ensembl_gene_id",
             "external_gene_name",
             "gene_biotype")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")
attr <- getBM(
  attributes = getinfo,
  filters = c("ensembl_gene_id"),
  values = union_genes,
  mart = mart)



alldata = list(ad_meta, pd_meta, pa_meta, psp_meta,
               scz_meta,asd_meta,bp_meta,mdd_meta)
names(alldata) = c("AD","PD","PA","PSP","Scz","ASD","BP","MDD")

## gather all data in one df

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

allExpr$symbol = attr$external_gene_name[match(allExpr$gene, attr$ensembl_gene_id)]
allExpr$Region = factor(allExpr$Region, 
                        levels = c("Cerebellum", "Temporal",
                                   "Frontal", "Limbic","Occipital","Basal.ganglia", 
                                   "Insular"))
allExpr$Disease = factor(allExpr$Disease, levels = names(alldata))

##save df for more analysis
saveRDS(object = allExpr,file = "./codes/regions_allExpr.rds")

#density of logFC for regions

library(forcats)
gp = ggplot(data = allExpr, 
            aes(x = logFC, y= fct_rev(Disease), fill = Disease))+
  geom_density_ridges(color = "white", rel_min_height=0.001)+
  geom_vline(xintercept = 0, lty = 2 , color = "red", size = 0.1)+
  facet_wrap(~Region, nrow = 1)+
  scale_fill_manual(values = cols)+
  theme_ridges(center_axis_labels = T, grid = F)+
  theme_economist()+
  theme(text = element_text(size = 8, colour = "black"),
        legend.position = "none",
        axis.text = element_text(size = 8))+
  labs(y = "", x = "logFC")+
  scale_x_continuous(breaks = c(-1,0,1), limits = c(-1,1))


ggsave(filename = "./results/figures/differential_expression/RegionsDensity.pdf", plot = gp,
       width = 8, height = 4)

##### Filter significant genes####

#A table of DEGs
DEGs = as.data.frame(matrix(NA, nrow = 8, ncol = 7))
rownames(DEGs) = names(alldata); 
colnames(DEGs) = unique(allExpr$Region)

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
ggsave(plot = DEG.table, 
       file = "./results/figures/differential_expression/Regions_DEGs.pdf",
       width = 5,
       height = 3)

####volcano plot####

library(ggrepel)

#color 
sig_up = -log10(allExpr$p.value)>1.3 & allExpr$logFC >0.5
sig_down = -log10(allExpr$p.value)>1.3 & allExpr$logFC < - 0.5

col = ifelse( sig_up , "green", ifelse(sig_down,"red","grey"))

g = ggplot(data = as.data.frame(allExpr),
           aes(x = logFC, y = -log10(p.value)))+
  geom_point(color = col , size = 0.2)+
  geom_text_repel(data = top2,
                  aes(label = symbol),
                  color = "black",
                  size = 3, 
                  show.legend = F)+
  facet_grid(Region~Disease,
             scales = "free_y",
  )+
  theme_bw()+
  theme(text = element_text(colour = "black", 
                            size = 12),
        strip.text.y  = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.key = element_rect(size = 6))

ggsave(plot= g , file = "./results/figures/differential_expression/Regions_deg_volcano.pdf",
       width = 12, height = 9)


#### word cloud to find the most frequent genes####
freq = as.data.frame(table(subdat$symbol)) %>%  arrange(desc(Freq))
colnames(freq) = c("gene", "Frequency")

#remove frequency <3
freq = subset(freq, Frequency >2)

rownames(subdat)= NULL
id = unique(as.character(freq$gene))
top.freq = subdat[subdat$symbol %in% id,]
top.freq$frequency = top.freq$p.value
top.freq$frequency = freq$Frequency[match(top.freq$symbol,
                                          freq$gene)]
top.freq = top.freq[order(top.freq$frequency,decreasing = T),]

# plot

freq.plot = ggplot(top.freq, 
                   aes(y = reorder(symbol, frequency), x = Region))+
  geom_point(aes(color = logFC),size = 4)+
  scale_color_continuous_diverging(palette = "Red-green")+
  facet_wrap(~Disease, nrow = 1)+
  labs(x = "",
       y = "")+
  theme_pubclean(base_size = 12)+
  theme(text = element_text(size = 12),
        axis.text.x.bottom = element_text(size = 7, angle = 45, hjust = 1))



ggsave(plot = freq.plot, 
       filename = "./results/figures/differential_expression/Regions_overlapping_genes.pdf",
       width = 10,
       height = 4)
