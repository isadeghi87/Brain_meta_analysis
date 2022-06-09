

# libraries
library(nlme);library(ggplot2);library(doParallel);library(foreach)
library(WGCNA); library(ggthemes); library(ggrepel)

# load data
load("./codes/st4_combinedData.Rdata")

# options
condition = TRUE
enableWGCNAThreads()
allowWGCNAThreads()

cores <- makeCluster(detectCores(), type='PSOCK')

system <- Sys.info()['sysname']

cl <- NULL
if (system == 'Windows') {
  cl <- makeCluster(getOption('cl.cores', cores))
  registerDoParallel(cl)
  registerDoSEQ()
  on.exit(stopCluster(cl))
} else {
  options('mc.cores' = cores)
  registerDoParallel(cores)
}

# Dx factors
datMeta$Dis = as.character(datMeta$Dx)
datMeta$Dis[datMeta$Dis != "CTL"] = "case"
datMeta$Dis = factor(datMeta$Dis, levels = c("CTL","case"))

# DGE between all diseases and controls

dge = matrix(NA, nrow=nrow(datExpr), ncol=3)
for(i in 1:nrow(datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(datExpr[i,])
  tryCatch({
    dge[i,] = summary(lme(expr~ Dis + Study + Sex + Brain_lobe ,data = datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

dge=as.data.frame(dge)
colnames(dge) = c("logFC", "SE", "p.value")
rownames(dge) = rownames(datExpr)
dge$fdr = p.adjust(dge$p.value, "fdr")
#dge= dge[!apply(is.na(dge),1,any),]

# significant
sigData = subset(dge,fdr < 0.05 & abs(logFC) > 0.5)

# load attributes
attr = readRDS("./codes/attributes.rds")
sigData$symbol = attr$external_gene_name[match(
  rownames(sigData),attr$ensembl_gene_id_version
)]

# order top logFC 
sigData = sigData[order(abs(sigData$logFC), decreasing = T),]
topGenes = sigData[1:10,]

# plot
pointCol = ifelse(dge$p.value< 0.05 & abs(dge$logFC) > 0.5, "cornflowerblue","grey")
textCol = ifelse(topGenes$logFC < 0, "red","green")

p = ggplot(dge, aes( x = logFC, y = -log10(p.value)))+
  #geom_point(color = "grey")+
  geom_point(aes(x = logFC, y = -log10(p.value)),
             color = pointCol)+
  geom_text_repel(data = topGenes,
                  aes(label = symbol), 
                  size = 4,
                  color = textCol,
                  show.legend = F
                  )+
  theme_classic(base_size = 14)

ggsave("./results/figures/DxCTL_topGenes.pdf", p, width = 6, height = 4)
# save data
save.image("./codes/DxCTL_DGE.RData")

 