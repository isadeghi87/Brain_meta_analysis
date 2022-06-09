#step 2------combine datasets per Diagnosis
#compile_data.R
#------------------------------------------#
library(WGCNA); library(nlme); library(reshape); library(corrplot);
library(ggplot2); library(corrplot); library(biomaRt); library(cqn); library(limma)
library(readr); library(ggplot2); library(purrr); library(stringr); library(dplyr);
library(devtools); library(statmod); library(pheatmap);
library(RColorBrewer)
#-----------------------------------------------#
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/")



#--load normalized AD data
AD.files = dir("./codes/", pattern = "_AD_")
n=length(AD.files)
AD_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/",AD.files[[i]],sep = ""))
  AD_expr[[i]]$datExpr = datExpr
  AD_expr[[i]]$datMeta = datMeta
  AD_expr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}

#---rename colnames of datmeta
AD_expr[[3]]$datMeta$Study <- "MSBB"

colnames(AD_expr[[3]]$datMeta) <- gsub("individualIdentifier", "Subject_ID", colnames(AD_expr[[3]]$datMeta))

genes= intersect(rownames(AD_expr[[1]]$datExpr), rownames(AD_expr[[2]]$datExpr))
genes = intersect(genes, rownames(AD_expr[[3]]$datExpr))
for(i in 1:n) AD_expr[[i]]$datExpr = AD_expr[[i]]$datExpr[match(genes,rownames(AD_expr[[i]]$datExpr)),]

AD_datExpr = data.frame(row.names = genes)
AD_datMeta = data.frame(Subject_ID=NA, Study=NA, Dx=NA, PMI=NA, RIN=NA, Age=NA, Sex=NA, Brain_Region=NA)

for (i in 1:n) {
  AD_datExpr = cbind(AD_datExpr, AD_expr[[i]]$datExpr)
  AD_datMeta = rbind(AD_datMeta, AD_expr[[i]]$datMeta[,c("Subject_ID", "Study", "Dx", "PMI", 
                                                         "RIN", "Age", "Sex", "Brain_Region")])
}

AD_datMeta = AD_datMeta[-1,]

keep = (AD_datMeta$Dx %in% c("Normal", "AD"))
AD_datExpr = AD_datExpr[,keep]
AD_datMeta = AD_datMeta[keep,]
AD_datMeta$Dx = factor(AD_datMeta$Dx, levels = c("Normal", "AD"))
AD_datMeta$Subject_ID = as.factor(AD_datMeta$Subject_ID)
AD_datMeta$Brain_Region = as.factor(AD_datMeta$Brain_Region)
AD_datMeta$PMI[is.na(AD_datMeta$PMI)] = mean(AD_datMeta$PMI, na.rm = T)
AD_datMeta$RIN[is.na(AD_datMeta$RIN)] = mean(AD_datMeta$RIN, na.rm = T)
AD_datMeta$Sex = as.factor((AD_datMeta$Sex))

ad_meta = matrix(NA, nrow=nrow(AD_datExpr), ncol=3)
for(i in 1:nrow(AD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(AD_datExpr[i,])
  tryCatch({
    ad_meta[i,] = summary(lme(expr~ Dx + Study,data = AD_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}
ad_meta=as.data.frame(ad_meta)
colnames(ad_meta) = c("logFC", "SE", "p.value")
rownames(ad_meta) = genes
ad_meta$fdr = p.adjust(ad_meta$p.value, "fdr")
ad_meta$symbol=AD_expr[[1]]$datProbes$external_gene_name[match(genes, AD_expr[[1]]$datProbes$ensembl_gene_id_version)]
ad_meta$biotype <- AD_expr[[1]]$datProbes$gene_biotype[match(genes, AD_expr[[1]]$datProbes$ensembl_gene_id_version)]
sig.ad_meta <-ad_meta[ad_meta$p.value<0.05,] 
write.csv(file="./results/AD_metaanalysis.csv",ad_meta)
write.csv(file="./results/AD_sig_metaanalysis.csv",sig.ad_meta)
#---------------------------------------------------------------------------------------------------------------------------------------#

#--normalized Scz data
Scz.files = dir("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/", pattern = "_Scz_")
n=length(Scz.files)
Scz_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/",Scz.files[[i]],sep = ""))
  Scz_expr[[i]]$datExpr = datExpr
  Scz_expr[[i]]$datMeta = datMeta
  Scz_expr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}


genes= intersect(rownames(Scz_expr[[1]]$datExpr), rownames(Scz_expr[[2]]$datExpr))
for(i in 3:n) genes = intersect(genes, rownames(Scz_expr[[i]]$datExpr))
for(i in 1:n) Scz_expr[[i]]$datExpr = Scz_expr[[i]]$datExpr[match(genes,rownames(Scz_expr[[i]]$datExpr)),]



#---merge datexpr and datmeta
Scz_datExpr = data.frame(row.names = genes)
Scz_datMeta = data.frame(Subject_ID=NA, Study=NA, Dx=NA, Age=NA, Sex=NA, Brain_Region=NA)

for (i in 1:n) {
  Scz_datExpr = cbind(Scz_datExpr, Scz_expr[[i]]$datExpr)
  Scz_datMeta = rbind(Scz_datMeta, Scz_expr[[i]]$datMeta[,c("Subject_ID", "Study", "Dx", "Age", "Sex", "Brain_Region")])
}

Scz_datMeta = Scz_datMeta[-1,]

keep = (Scz_datMeta$Dx %in% c("Normal", "Scz"))
Scz_datExpr = Scz_datExpr[,keep]
Scz_datMeta = Scz_datMeta[keep,]
Scz_datMeta$Dx = factor(Scz_datMeta$Dx, levels = c("Normal", "Scz"))
Scz_datMeta$Subject_ID = as.factor(Scz_datMeta$Subject_ID)
Scz_datMeta$Brain_Region = as.factor(Scz_datMeta$Brain_Region)

Scz_meta = matrix(NA, nrow=nrow(Scz_datExpr), ncol=3)
for(i in 1:nrow(Scz_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(Scz_datExpr[i,])
  tryCatch({
    Scz_meta[i,] = summary(lme(expr~ Dx + Study ,data = Scz_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

Scz_meta=as.data.frame(Scz_meta)
colnames(Scz_meta) <- c("logFC", "SE", "p.value")
rownames(Scz_meta) = genes
Scz_meta$fdr = p.adjust(Scz_meta$p.value, "fdr")
Scz_meta$symbol=Scz_expr[[1]]$datProbes$external_gene_name[match(genes, Scz_expr[[1]]$datProbes$ensembl_gene_id_version)]
Scz_meta$biotype <- Scz_expr[[1]]$datProbes$gene_biotype[match(genes, Scz_expr[[1]]$datProbes$ensembl_gene_id_version)]
sig.Scz_meta <-Scz_meta[Scz_meta$p.value<0.05,] 
write.csv(file="./results/Scz_metaanalysis.csv",Scz_meta)
write.csv(file="./results/Scz_sig_metaanalysis.csv",sig.Scz_meta)
#-------------------------------------------------------------------------------------------------------------------------------#
#--normalized BP data
BP.files = dir("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/", pattern = "_BP_")
n=length(BP.files)
BP_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/",BP.files[[i]],sep = ""))
  BP_expr[[i]]$datExpr = datExpr
  BP_expr[[i]]$datMeta = datMeta
  BP_expr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}



genes= intersect(rownames(BP_expr[[1]]$datExpr), rownames(BP_expr[[2]]$datExpr))
for(i in 3:n) genes = intersect(genes, rownames(BP_expr[[i]]$datExpr))
for(i in 1:n) BP_expr[[i]]$datExpr = BP_expr[[i]]$datExpr[match(genes,rownames(BP_expr[[i]]$datExpr)),]

BP_datExpr = data.frame(row.names = genes)
BP_datMeta = data.frame(Subject_ID=NA, Study=NA, Dx=NA, Age=NA, Sex=NA, Brain_Region=NA)

for (i in 1:n) {
  BP_datExpr = cbind(BP_datExpr, BP_expr[[i]]$datExpr)
  BP_datMeta = rbind(BP_datMeta, BP_expr[[i]]$datMeta[,c("Subject_ID", "Study", "Dx", "Age", "Sex", "Brain_Region")])
}

BP_datMeta = BP_datMeta[-1,]

keep = (BP_datMeta$Dx %in% c("Normal", "BP"))
BP_datExpr = BP_datExpr[,keep]
BP_datMeta = BP_datMeta[keep,]
BP_datMeta$Dx = factor(BP_datMeta$Dx, levels = c("Normal", "BP"))
BP_datMeta$Subject_ID = as.factor(BP_datMeta$Subject_ID)

BP_meta = matrix(NA, nrow=nrow(BP_datExpr), ncol=3)
for(i in 1:nrow(BP_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(BP_datExpr[i,])
  tryCatch({
    BP_meta[i,] = summary(lme(expr~ Dx + Study,data = BP_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

BP_meta=as.data.frame(BP_meta)
colnames(BP_meta) = c("logFC", "SE", "p.value")
rownames(BP_meta) = genes
BP_meta$fdr = p.adjust(BP_meta$p.value, "fdr")
BP_meta$symbol=BP_expr[[1]]$datProbes$external_gene_name[match(genes, BP_expr[[1]]$datProbes$ensembl_gene_id_version)]
BP_meta$biotype <- BP_expr[[1]]$datProbes$gene_biotype[match(genes, BP_expr[[1]]$datProbes$ensembl_gene_id_version)]
sig.BP_meta <-BP_meta[BP_meta$p.value<0.05,] 
write.csv(file="./results/BP_metaanalysis.csv",BP_meta)
write.csv(file="./results/BP_sig_metaanalysis.csv",sig.BP_meta)

#---------------------------------------------------------------------------------------------------------------------------#

#--normalized MDD data
MDD.files = dir("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/", pattern = "_MDD_")
n=length(MDD.files)
MDD_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/",MDD.files[[i]],sep = ""))
  MDD_expr[[i]]$datExpr = datExpr
  MDD_expr[[i]]$datMeta = datMeta
  MDD_expr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}


#---keep common genes only
genes= intersect(rownames(MDD_expr[[1]]$datExpr), rownames(MDD_expr[[2]]$datExpr))
for(i in 1:n) MDD_expr[[i]]$datExpr = MDD_expr[[i]]$datExpr[match(genes,rownames(MDD_expr[[i]]$datExpr)),]

MDD_datExpr = data.frame(row.names = genes)
MDD_datMeta = data.frame(Subject_ID=NA, Study=NA, Dx=NA, PMI=NA, Age=NA, Sex=NA, Brain_Region=NA)

for (i in 1:n) {
  MDD_datExpr = cbind(MDD_datExpr, MDD_expr[[i]]$datExpr)
  MDD_datMeta = rbind(MDD_datMeta, MDD_expr[[i]]$datMeta[,c("Subject_ID", "Study", "Dx", "PMI", 
                                                            "Age", "Sex", "Brain_Region")])
}

MDD_datMeta = MDD_datMeta[-1,]
MDD_datMeta$Brain_Region[MDD_datMeta$Brain_Region=="BA9"]="DLPFC"
keep = (MDD_datMeta$Dx %in% c("Normal", "MDD"))
MDD_datExpr = MDD_datExpr[,keep]
MDD_datMeta = MDD_datMeta[keep,]
MDD_datMeta$Dx = factor(MDD_datMeta$Dx,levels = c("Normal", "MDD"))
MDD_datMeta$Subject_ID = as.factor(MDD_datMeta$Subject_ID)
MDD_datMeta$Brain_Region = as.factor(MDD_datMeta$Brain_Region)
MDD_datMeta$Sex = as.factor(MDD_datMeta$Sex)

MDD_meta = matrix(NA, nrow=nrow(MDD_datExpr), ncol=3)
for(i in 1:nrow(MDD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(MDD_datExpr[i,])
  tryCatch({
    MDD_meta[i,] = summary(lme(expr~ Dx + Study+ Age+
                                 PMI+ Sex ,
                               data = MDD_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

MDD_meta=as.data.frame(MDD_meta)
colnames(MDD_meta) = c("logFC", "SE", "p.value")
rownames(MDD_meta) = genes
MDD_meta$fdr = p.adjust(MDD_meta$p.value, "fdr")
MDD_meta$symbol=MDD_expr[[1]]$datProbes$external_gene_name[match(genes, MDD_expr[[1]]$datProbes$ensembl_gene_id_version)]
MDD_meta$biotype <- MDD_expr[[1]]$datProbes$gene_biotype[match(genes, MDD_expr[[1]]$datProbes$ensembl_gene_id_version)]
sig.MDD_meta <-MDD_meta[MDD_meta$p.value<0.05,] 
write.csv(file="./results/MDD_metaanalysis.csv",MDD_meta)
write.csv(file="./results/MDD_sig_metaanalysis.csv",sig.MDD_meta)
#MDD_meta= bd_meta[!apply(is.na(bd_meta),1,any),]
#--------------------------------------------------------------------------------------------------------------------#

#--normalized ASD data
ASD.files = dir("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/", pattern = "_ASD_")
n=length(ASD.files)
ASD_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/",ASD.files[[i]],sep = ""))
  ASD_expr[[i]]$datExpr = datExpr
  ASD_expr[[i]]$datMeta = datMeta
  ASD_expr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}



#---keep common genes only
genes= intersect(rownames(ASD_expr[[1]]$datExpr), rownames(ASD_expr[[2]]$datExpr))
for(i in 2:n) genes = intersect(genes, rownames(ASD_expr[[i]]$datExpr))
for(i in 1:n) ASD_expr[[i]]$datExpr = ASD_expr[[i]]$datExpr[match(genes,rownames(ASD_expr[[i]]$datExpr)),]

ASD_datExpr = data.frame(row.names = genes)
ASD_datMeta = data.frame(Subject_ID=NA, Study=NA, Dx=NA, Age=NA, Sex=NA, Brain_Region=NA)

for (i in 1:n) {
  ASD_datExpr = cbind(ASD_datExpr, ASD_expr[[i]]$datExpr)
  ASD_datMeta = rbind(ASD_datMeta, ASD_expr[[i]]$datMeta[,c("Subject_ID", "Study", "Dx", 
                                                            "Age", "Sex", "Brain_Region")])
}

ASD_datMeta = ASD_datMeta[-1,]

keep = (ASD_datMeta$Dx %in% c("Normal", "ASD"))
ASD_datExpr = ASD_datExpr[,keep]
ASD_datMeta = ASD_datMeta[keep,]
ASD_datMeta$Dx = factor(ASD_datMeta$Dx, levels = c("Normal", "ASD"))
ASD_datMeta$Subject_ID = as.factor(ASD_datMeta$Subject_ID)
ASD_datMeta$Brain_Region = as.factor(ASD_datMeta$Brain_Region)
ASD_datMeta$Sex = as.factor(ASD_datMeta$Sex)
ASD_datMeta$Study = as.factor(ASD_datMeta$Study)
ASD_meta = matrix(NA, nrow=nrow(ASD_datExpr), ncol=3)
for(i in 1:nrow(ASD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(ASD_datExpr[i,])
  tryCatch({
    ASD_meta[i,] = summary(lme(expr~ Dx + Study ,data = ASD_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

ASD_meta=as.data.frame(ASD_meta)
colnames(ASD_meta) = c("logFC", "SE", "p.value")
rownames(ASD_meta) = genes
ASD_meta$fdr = p.adjust(ASD_meta$p.value, "fdr")
ASD_meta$symbol=ASD_expr[[1]]$datProbes$external_gene_name[match(genes, ASD_expr[[1]]$datProbes$ensembl_gene_id_version)]
ASD_meta$biotype <- ASD_expr[[1]]$datProbes$gene_biotype[match(genes, ASD_expr[[1]]$datProbes$ensembl_gene_id_version)]
sig.ASD_meta <-ASD_meta[ASD_meta$p.value<0.05,] 
write.csv(file="./results/ASD_metaanalysis.csv",ASD_meta)
write.csv(file="./results/ASD_sig_metaanalysis.csv",sig.ASD_meta)
#------------------------------------------------------------------------------------------------------------------------#

#--normalized PD data

#--It is a single dataset so we load sumstats data
PD.files = dir("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/", pattern = "_PD_")
load("./codes/pj283498_PD_normalized.Rdata")
PD_datExpr = datExpr; rm(datExpr)
PD_datMeta = datMeta; rm(datMeta)
PD_datProbes = datProbes; rm(datProbes)
genes= rownames(PD_datExpr)
PD_datMeta = PD_datMeta[,c("Subject_ID", "Study", "Dx", "PMI", 
                           "RIN", "Age", "Sex", "Brain_Region")]

PD_meta= read.csv(file ="results/pj283498_sumstats.csv", row.names = 1)
PD_meta = PD_meta %>% select(logFC,t, P.Value, adj.P.Val,gene,biotype)
colnames(PD_meta) <- c("logFC", "t", "p.value", "fdr", "symbol", "biotype")
sig.PD_meta <-PD_meta[PD_meta$p.value<0.05,] 
write.csv(file="./results/PD_metaanalysis.csv",PD_meta)
write.csv(file="./results/PD_sig_metaanalysis.csv",sig.PD_meta)



#----------------------------------------------------------------------------------------------------------------#

#--normalized PA data
PA.files = dir("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/", pattern = "_PA_")
n=length(PA.files)
PA_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/",PA.files[[i]],sep = ""))
  PA_expr[[i]]$datExpr = datExpr
  PA_expr[[i]]$datMeta = datMeta
  PA_expr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}


#----correct metadata
genes= intersect(rownames(PA_expr[[1]]$datExpr), rownames(PA_expr[[2]]$datExpr))
for(i in 1:n) PA_expr[[i]]$datExpr = PA_expr[[i]]$datExpr[match(genes,rownames(PA_expr[[i]]$datExpr)),]

PA_datExpr = data.frame(row.names = genes)
PA_datMeta = data.frame(Subject_ID=NA, Study=NA, Dx=NA, PMI=NA, RIN=NA, Age=NA, Sex=NA, Brain_Region=NA)

for (i in 1:n) {
  PA_datExpr = cbind(PA_datExpr, PA_expr[[i]]$datExpr)
  PA_datMeta = rbind(PA_datMeta, PA_expr[[i]]$datMeta[,c("Subject_ID", "Study", "Dx", "PMI", 
                                                         "RIN", "Age", "Sex", "Brain_Region")])
}

PA_datMeta = PA_datMeta[-1,]

keep = (PA_datMeta$Dx %in% c("Normal", "PA"))
PA_datExpr = PA_datExpr[,keep]
PA_datMeta = PA_datMeta[keep,]
PA_datMeta$Dx = factor(PA_datMeta$Dx, levels = c("Normal", "PA"))
PA_datMeta$Subject_ID = as.factor(PA_datMeta$Subject_ID)
PA_datMeta$Study = as.factor(PA_datMeta$Study)
PA_datMeta$Brain_Region = as.factor(PA_datMeta$Brain_Region)
PA_datMeta$Sex = as.factor(PA_datMeta$Sex)
PA_datMeta$PMI[is.na(PA_datMeta$PMI)]= mean(PA_datMeta$PMI, na.rm = T)

PA_meta = matrix(NA, nrow=nrow(PA_datExpr), ncol=3)
for(i in 1:nrow(PA_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(PA_datExpr[i,])
  tryCatch({
    PA_meta[i,] = summary(lme(expr~ Dx + Study, data = PA_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

PA_meta=as.data.frame(PA_meta)
colnames(PA_meta) = c("logFC", "SE", "p.value")
rownames(PA_meta) = genes
PA_meta$fdr = p.adjust(PA_meta$p.value, "fdr")
PA_meta$symbol=PA_expr[[1]]$datProbes$external_gene_name[match(genes, PA_expr[[1]]$datProbes$ensembl_gene_id_version)]
PA_meta$biotype <- PA_expr[[1]]$datProbes$gene_biotype[match(genes, PA_expr[[1]]$datProbes$ensembl_gene_id_version)]
sig.PA_meta <-PA_meta[PA_meta$p.value<0.05,] 
write.csv(file="./results/PA_metaanalysis.csv",PA_meta)
write.csv(file="./results/PA_sig_metaanalysis.csv",sig.PA_meta)

#-------------------------------------------------------------------------------------------------------------#
#--normalized PSP data
PSP.files = dir("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/", pattern = "_PSP_")
n=length(PSP.files)
PSP_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/codes/",PSP.files[[i]],sep = ""))
  PSP_expr[[i]]$datExpr = datExpr
  PSP_expr[[i]]$datMeta = datMeta
  PSP_expr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}


genes= intersect(rownames(PSP_expr[[1]]$datExpr), rownames(PSP_expr[[2]]$datExpr))
for(i in 1:n) PSP_expr[[i]]$datExpr = PSP_expr[[i]]$datExpr[match(genes,rownames(PSP_expr[[i]]$datExpr)),]

PSP_datExpr = data.frame(row.names = genes)
PSP_datMeta = data.frame(Subject_ID=NA, Study=NA, Dx=NA, PMI=NA, RIN=NA, Age=NA, Sex=NA, Brain_Region=NA)

for (i in 1:n) {
  PSP_datExpr = cbind(PSP_datExpr, PSP_expr[[i]]$datExpr)
  PSP_datMeta = rbind(PSP_datMeta, PSP_expr[[i]]$datMeta[,c("Subject_ID", "Study", "Dx", "PMI", 
                                                            "RIN", "Age", "Sex", "Brain_Region")])
}

PSP_datMeta = PSP_datMeta[-1,]

keep = (PSP_datMeta$Dx %in% c("Normal", "PSP"))
PSP_datExpr = PSP_datExpr[,keep]
PSP_datMeta = PSP_datMeta[keep,]
PSP_datMeta$Dx = factor(PSP_datMeta$Dx, levels = c("Normal", "PSP"))
PSP_datMeta$Subject_ID = as.factor(PSP_datMeta$Subject_ID)
PSP_datMeta$PMI[is.na(PSP_datMeta$PMI)] = mean(PSP_datMeta$PMI, na.rm = T)

PSP_meta = matrix(NA, nrow=nrow(PSP_datExpr), ncol=3)
for(i in 1:nrow(PSP_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(PSP_datExpr[i,])
  tryCatch({
    PSP_meta[i,] = summary(lme(expr~ Dx +Study,data = PSP_datMeta, random=~1|Subject_ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

PSP_meta=as.data.frame(PSP_meta)
colnames(PSP_meta) = c("logFC", "SE", "p.value")
rownames(PSP_meta) = genes
PSP_meta$fdr = p.adjust(PSP_meta$p.value, "fdr")
PSP_meta$symbol=PSP_expr[[1]]$datProbes$external_gene_name[match(genes, PSP_expr[[1]]$datProbes$ensembl_gene_id_version)]
PSP_meta$biotype <- PSP_expr[[1]]$datProbes$gene_biotype[match(genes, PSP_expr[[1]]$datProbes$ensembl_gene_id_version)]
sig.PSP_meta <-PSP_meta[PSP_meta$p.value<0.05,] 
write.csv(file="./results/PSP_metaanalysis.csv",PSP_meta)
write.csv(file="./results/PSP_sig_metaanalysis.csv",sig.PSP_meta)


#---------------------------------------------------------------------------------------------------------------------------#


multiExpr = vector(mode="list", length=8)
multiExpr[[1]]$datExpr = AD_datExpr;  multiExpr[[1]]$datMeta = AD_datMeta 
multiExpr[[2]]$datExpr = Scz_datExpr; multiExpr[[2]]$datMeta = Scz_datMeta 
multiExpr[[3]]$datExpr = BP_datExpr; multiExpr[[3]]$datMeta = BP_datMeta 
multiExpr[[4]]$datExpr = MDD_datExpr; multiExpr[[4]]$datMeta = MDD_datMeta 
multiExpr[[5]]$datExpr = ASD_datExpr; multiExpr[[5]]$datMeta = ASD_datMeta 
multiExpr[[6]]$datExpr = PD_datExpr; multiExpr[[6]]$datMeta = PD_datMeta
multiExpr[[7]]$datExpr = PSP_datExpr; multiExpr[[7]]$datMeta = PSP_datMeta 
multiExpr[[8]]$datExpr = PA_datExpr; multiExpr[[8]]$datMeta = PA_datMeta
names(multiExpr) = c("AD","Scz","BP","MDD","ASD","PD","PSP","PA")

##Save compiled expression and metadata for linear model analysis
lm_multiExpr = multiExpr
save(file = "./codes/st2_lm_multiExpr.Rdata", lm_multiExpr)

##Save compiled expression and metadata for permutation testing
genes = intersect(rownames(multiExpr[[1]]$datExpr), rownames(multiExpr[[2]]$datExpr))
for(i in 3:8) genes = intersect(genes, rownames(multiExpr[[i]]$datExpr))
for(i in 1:8) multiExpr[[i]]$datExpr = multiExpr[[i]]$datExpr[match(genes, rownames(multiExpr[[i]]$datExpr)),]
save(file="./codes/st2_Permutation_data.RData", multiExpr)



save.image(file = "./codes/st2_metaAnalysis.Rdata")
