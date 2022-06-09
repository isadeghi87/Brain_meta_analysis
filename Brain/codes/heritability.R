##########
## Heritability analysis

#packages
library(DESeq2); library(dplyr); library(magrittr)
library(HeritSeq); library(ggplot2); library(biomaRt)


## Provide normalized data 
# this contains a matrix of 
#expression data( rows are gene and columns are samples) and Meta data (datMeta)

load("./codes/st4_combinedData.Rdata")
#gene ids
attr = readRDS(file = "./codes/attributes.rds")

## heritability:
result.vst <- fitComputeVPC.lmer(CountMatrix = datExpr, Strains =datMeta$Dx,
                                 test = TRUE)
## Extract parameters
vpc.vst <- as.data.frame(result.vst[[1]])
vpc.vst$gene = rownames(vpc.vst)
vpc.vst$id = attr$external_gene_name[match(vpc.vst$gene,attr$ensembl_gene_id_version)]

## Extract p-values
vpc.vst$p.value <- result.vst[[2]]
vpc.vst$FDR = p.adjust(vpc.vst$p.value,method = "fdr")
pvals <- p.adjust(pval.vst$`P-value`,method = "fdr")

#median
vpc.vst$expression = rowMedians(datExpr)

#sort by highest heritability
vpc.vst = vpc.vst %>% arrange(desc(LMM))

##plot 
g =ggplot(vpc.vst[1:20,], aes(x = reorder(id,-LMM), 
                           y = LMM))+
  geom_point(color = 3, size = 4)+
  labs(x = "gene", y = "heritability")+
  # geom_text(aes(label = signif(LMM,2)),vjust = -0.5)+
  theme_bw()+
   theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45,hjust = 1))
  

ggsave(filename = "./results/figures/Heritability.pdf",plot = g,
       width = 5, height = 3)

#heritability for diseases
vpcs = as.data.frame(matrix(NA,nrow = nrow(datExpr),ncol = 8))
rownames(vpcs) = rownames(datExpr)
colnames(vpcs) = c("AD", "PD", "PA", "PSP", "ASD", "Scz", "MDD", "BP")
pv = vpcs

for ( i in c("AD", "PD", "PA", "PSP", "ASD", "Scz", "MDD", "BP")) {
  m = which(datMeta$Dx %in% c(i,"Normal"))
  meta = datMeta[m,]
  expr = datExpr[,m]
  h = fitComputeVPC.lmer(CountMatrix = expr, 
                         Strains =meta$Dx,
                         test = TRUE)
  vpcs[,i] = h[[1]]
  pv[,i] = h[[2]]
}

vpcs$gene = attr$external_gene_name[match(rownames(vpcs),attr$ensembl_gene_id_version)]
her = melt(vpcs)

m = as.data.frame(rowMedians(t(vpcs[,1:8])))
m$disease = colnames(vpcs)[1:8]
sd= vector()
se= vector()
for (i in 1:8) {
  s = sd(vpcs[,i])
  sd= c(sd,s)
  e=stdErr(vpcs[,i])[1]
  se = c(se,e)
}
m$sd = sd
m$se = m$SD/sqrt(10313)
colnames(m)= c("Heritability", "Disease", "SD","SE")
p = ggplot(m,aes(x= reorder(Disease,-Heritability),y = Heritability, 
               color= Disease))+
  geom_point(show.legend = F, size= 3)+
  geom_errorbar(aes( ymin = (Heritability - SD), 
                     ymax = (Heritability + SD)),
                position = position_dodge( width = 0.8),
                width = 0.25, size = 0.25)+
  theme_bw()+
  labs(y= bquote("Heritability"~(h^2)),
       x = "")+
  theme(legend.position = "none",
        axis.text.x.bottom = element_text(angle = 45, hjust = 0.8))

#save plot
ggsave(p, filename = "./results/figures/heritability_disease.pdf",
       width = 4, height = 3)
#save data
save.image("./codes/Heritability.Rdata")
