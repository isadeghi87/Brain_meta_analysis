#step 2------combine datasets per Diagnosis
#compile_data.R
#------------------------------------------#
library(WGCNA); library(nlme); library(reshape); library(corrplot);
library(ggplot2); library(corrplot); library(biomaRt); library(cqn); library(limma)
library(readr); library(ggplot2); library(purrr); library(stringr); library(dplyr);
library(devtools); library(statmod); library(pheatmap);
library(RColorBrewer)
#-----------------------------------------------#
rm(list = ls())
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")


#### AD data #### 
AD.files = dir("../normalized_studies/codes/", pattern = "_AD_")
n=length(AD.files)
AD_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("../normalized_studies/codes/",AD.files[[i]],sep = ""))
  AD_expr[[i]]$datExpr = datExpr
  AD_expr[[i]]$datMeta = datMeta
  rm(datExpr,datMeta)
}


genes= intersect(rownames(AD_expr[[1]]$datExpr), rownames(AD_expr[[2]]$datExpr))
genes = intersect(genes, rownames(AD_expr[[3]]$datExpr))
for(i in 1:n) AD_expr[[i]]$datExpr = AD_expr[[i]]$datExpr[match(genes,rownames(AD_expr[[i]]$datExpr)),]

AD_datExpr = data.frame(row.names = genes)
AD_datMeta = data.frame(Sample=NA, Subject_ID=NA, Study=NA, Dx=NA, PMI=NA, RIN=NA, 
                        Age=NA, Sex=NA, Brain_Region=NA, Brain_Lobe=NA)

for (i in 1:n) {
  AD_datExpr = cbind(AD_datExpr, AD_expr[[i]]$datExpr)
  AD_datMeta = rbind(AD_datMeta, AD_expr[[i]]$datMeta[,c("Sample", 
                                                         "Subject_ID", "Study", "Dx", "PMI", 
                                                         "RIN", "Age", "Sex", "Brain_Region","Brain_Lobe")])
}

AD_datMeta = AD_datMeta[-1,]

keep = (AD_datMeta$Dx %in% c("CTL", "AD"))
AD_datExpr = AD_datExpr[,keep]
AD_datMeta = AD_datMeta[keep,]
AD_datMeta$Dx = factor(AD_datMeta$Dx, levels = c("CTL", "AD"))
AD_datMeta$Subject_ID = as.factor(AD_datMeta$Subject_ID)
AD_datMeta$Brain_Region = as.factor(AD_datMeta$Brain_Region)
AD_datMeta$PMI[is.na(AD_datMeta$PMI)] = mean(AD_datMeta$PMI, na.rm = T)
AD_datMeta$RIN[is.na(AD_datMeta$RIN)] = mean(AD_datMeta$RIN, na.rm = T)
AD_datMeta$Sex = as.factor((AD_datMeta$Sex))

# check duplicated samples
idx = !duplicated(AD_datMeta$Sample)
AD_datMeta = AD_datMeta[idx,]
AD_datExpr = AD_datExpr[,idx]
all(colnames(AD_datExpr)==rownames(AD_datMeta))

# check mds plot
pdf("./results/figures/AD_PCA2.pdf")
plotMDS(log2(AD_datExpr), 
        col = as.numeric(as.factor(AD_datMeta$Study)),
        pch = as.numeric(as.factor(AD_datMeta$Dx)))
dev.off()
#---------------------------------------------------------------------------------------------------------------------------------------#

#### PD data ####

#--It is a single dataset so we load sumstats data
PD.files = dir("../normalized_studies/codes/", pattern = "_PD_")
load(paste("../normalized_studies/codes/",PD.files,sep = ""))

PD_datExpr = datExpr; rm(datExpr)
PD_datMeta = datMeta; rm(datMeta)

PD_datMeta = PD_datMeta[,c("Sample","Subject_ID", "Study", "Dx", "PMI", 
                           "RIN", "Age", "Sex", "Brain_Region","Brain_Lobe")]
all(colnames(PD_datExpr)==rownames(PD_datMeta))

#----------------------------------------------------------------------------------------------------------------#

####  PA data ####

PA.files = dir("../normalized_studies/codes/", pattern = "_PA_")
n=length(PA.files)
PA_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("../normalized_studies/codes/",PA.files[[i]],sep = ""))
  PA_expr[[i]]$datExpr = datExpr
  PA_expr[[i]]$datMeta = datMeta
  rm(datExpr,datMeta)
}


#----correct metadata
genes= intersect(rownames(PA_expr[[1]]$datExpr), rownames(PA_expr[[2]]$datExpr))
for(i in 1:n) PA_expr[[i]]$datExpr = PA_expr[[i]]$datExpr[match(genes,rownames(PA_expr[[i]]$datExpr)),]

PA_datExpr = data.frame(row.names = genes)
PA_datMeta = data.frame(Sample =NA,
                        Subject_ID=NA, 
                        Study=NA,
                        Dx=NA, 
                        PMI=NA,
                        RIN=NA, 
                        Age=NA,
                        Sex=NA, 
                        Brain_Region=NA,
                        Brain_Lobe=NA)

for (i in 1:n) {
  PA_datExpr = cbind(PA_datExpr, PA_expr[[i]]$datExpr)
  PA_datMeta = rbind(PA_datMeta, PA_expr[[i]]$datMeta[,c("Sample","Subject_ID",
                                                         "Study", 
                                                         "Dx", 
                                                         "PMI", 
                                                         "RIN",
                                                         "Age",
                                                         "Sex", 
                                                         "Brain_Region",
                                                         "Brain_Lobe")])
}

PA_datMeta = PA_datMeta[-1,]

keep = (PA_datMeta$Dx %in% c("CTL", "PA"))
PA_datExpr = PA_datExpr[,keep]
PA_datMeta = PA_datMeta[keep,]
PA_datMeta$Dx = factor(PA_datMeta$Dx, levels = c("CTL", "PA"))
PA_datMeta$Subject_ID = as.factor(PA_datMeta$Subject_ID)
PA_datMeta$Study = as.factor(PA_datMeta$Study)
PA_datMeta$Brain_Region = as.factor(PA_datMeta$Brain_Region)
PA_datMeta$Sex = as.factor(PA_datMeta$Sex)
PA_datMeta$PMI[is.na(PA_datMeta$PMI)]= mean(PA_datMeta$PMI, na.rm = T)

# check duplicated samples
idx = !duplicated(PA_datMeta$Sample)
PA_datMeta = PA_datMeta[idx,]
PA_datExpr = PA_datExpr[,idx]
all(colnames(PA_datExpr)==rownames(PA_datMeta))

pdf("./results/figures/PA_PCA.pdf")
plotMDS(log2(PA_datExpr), 
        col = as.numeric(as.factor(PA_datMeta$Brain_Region)),
        pch = as.numeric(as.factor(PA_datMeta$Dx)))
dev.off()

#-------------------------------------------------------------------------------------------------------------#
#### normalized PSP data ####

PSP.files = dir("../normalized_studies/codes/", pattern = "_PSP_")
n=length(PSP.files)
PSP_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("../normalized_studies/codes/",PSP.files[[i]],sep = ""))
  PSP_expr[[i]]$datExpr = datExpr
  PSP_expr[[i]]$datMeta = datMeta
  rm(datExpr,datMeta)
}


genes= intersect(rownames(PSP_expr[[1]]$datExpr), rownames(PSP_expr[[2]]$datExpr))
for(i in 1:n) PSP_expr[[i]]$datExpr = PSP_expr[[i]]$datExpr[match(genes,rownames(PSP_expr[[i]]$datExpr)),]

PSP_datExpr = data.frame(row.names = genes)
PSP_datMeta = data.frame(Sample =NA, 
                         Subject_ID=NA,
                         Study=NA,
                         Dx=NA,
                         PMI=NA, 
                         RIN=NA,
                         Age=NA,
                         Sex=NA, 
                         Brain_Region=NA,
                         Brain_Lobe = NA)

for (i in 1:n) {
  PSP_datExpr = cbind(PSP_datExpr, PSP_expr[[i]]$datExpr)
  PSP_datMeta = rbind(PSP_datMeta, PSP_expr[[i]]$datMeta[,c("Sample","Subject_ID",
                                                            "Study",
                                                            "Dx",
                                                            "PMI", 
                                                            "RIN",
                                                            "Age",
                                                            "Sex",
                                                            "Brain_Region",
                                                            "Brain_Lobe")])
}

PSP_datMeta = PSP_datMeta[-1,]

keep = (PSP_datMeta$Dx %in% c("CTL", "PSP"))
PSP_datExpr = PSP_datExpr[,keep]
PSP_datMeta = PSP_datMeta[keep,]
PSP_datMeta$Dx = factor(PSP_datMeta$Dx, levels = c("CTL", "PSP"))
PSP_datMeta$Subject_ID = as.factor(PSP_datMeta$Subject_ID)
PSP_datMeta$PMI[is.na(PSP_datMeta$PMI)] = mean(PSP_datMeta$PMI, na.rm = T)

# check duplicated samples
idx = !duplicated(PSP_datMeta$Sample)
PSP_datMeta = PSP_datMeta[idx,]
PSP_datExpr = PSP_datExpr[,idx]
all(colnames(PSP_datExpr)==rownames(PSP_datMeta))

##### normalized Scz data ####

Scz.files = dir("../normalized_studies/codes/", pattern = "_Scz_")
n=length(Scz.files)
Scz_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("../normalized_studies/codes/",Scz.files[[i]],sep = ""))
  if ("datExpr_scz" %in% ls()) {
    Scz_expr[[i]]$datExpr = datExpr_scz
    Scz_expr[[i]]$datMeta = datMeta_scz
    rm(datExpr_scz,datMeta_scz)
  } else {
  Scz_expr[[i]]$datExpr = datExpr
  Scz_expr[[i]]$datMeta = datMeta
  rm(datExpr,datMeta)
  }
}

Scz_expr[[6]]$datMeta$Sample = rownames(Scz_expr[[6]]$datMeta)

genes= intersect(rownames(Scz_expr[[1]]$datExpr), rownames(Scz_expr[[2]]$datExpr))
for(i in 3:n) genes = intersect(genes, rownames(Scz_expr[[i]]$datExpr))
for(i in 1:n) Scz_expr[[i]]$datExpr = Scz_expr[[i]]$datExpr[match(genes,rownames(Scz_expr[[i]]$datExpr)),]

#---merge datexpr and datmeta
Scz_datExpr = data.frame(row.names = genes)
Scz_datMeta = data.frame(Sample = NA, 
                         Subject_ID=NA, 
                         Study=NA,
                         Dx=NA,
                         Age=NA,
                         Sex=NA,
                         Brain_Region=NA,
                         Brain_Lobe=NA)

for (i in 1:n) {
  Scz_datExpr = cbind(Scz_datExpr, Scz_expr[[i]]$datExpr)
  Scz_datMeta = rbind(Scz_datMeta, Scz_expr[[i]]$datMeta[,c("Sample",
                                                            "Subject_ID", 
                                                            "Study", 
                                                            "Dx",
                                                            "Age",
                                                            "Sex", 
                                                            "Brain_Region",
                                                            "Brain_Lobe")])
}

Scz_datMeta = Scz_datMeta[-1,]

keep = (Scz_datMeta$Dx %in% c("CTL", "Scz"))
Scz_datExpr = Scz_datExpr[,keep]
Scz_datMeta = Scz_datMeta[keep,]
Scz_datMeta$Dx = factor(Scz_datMeta$Dx, levels = c("CTL", "Scz"))
Scz_datMeta$Subject_ID = factor(Scz_datMeta$Subject_ID)
Scz_datMeta$Brain_Region = factor(Scz_datMeta$Brain_Region)

# check duplicated samples
idx = !duplicated(Scz_datMeta$Sample)
id = which(rownames(Scz_datMeta)!= Scz_datMeta$Sample)
Scz_datMeta = Scz_datMeta[idx,]
Scz_datExpr = Scz_datExpr[,idx]
all(colnames(Scz_datExpr)== rownames(Scz_datMeta))

#-------------------------------------------------------------------------------------------------------------------------------#
#### ASD data ####

ASD.files = dir("../normalized_studies/codes/", pattern = "_ASD_")

# exclude pj398545 study
ASD.files = ASD.files[!grepl("398545",ASD.files)]
n=length(ASD.files)
ASD_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("../normalized_studies/codes/",ASD.files[[i]],sep = ""))
  ASD_expr[[i]]$datExpr = datExpr
  ASD_expr[[i]]$datMeta = datMeta
  rm(datExpr,datMeta)
}

#---keep common genes only
genes = intersect(rownames(ASD_expr[[1]]$datExpr), rownames(ASD_expr[[2]]$datExpr))
for(i in 3:n) genes = intersect(genes, rownames(ASD_expr[[i]]$datExpr))
for(i in 1:n) ASD_expr[[i]]$datExpr = ASD_expr[[i]]$datExpr[match(genes,rownames(ASD_expr[[i]]$datExpr)),]


ASD_datExpr = data.frame(row.names = genes)
ASD_datMeta = data.frame(Sample=NA,
                         Subject_ID=NA, 
                         Study=NA,
                         Dx=NA, 
                         Age=NA,
                         Sex=NA, 
                         Brain_Region=NA,
                         Brain_Lobe=  NA)

for (i in 1:n) {
  ASD_datExpr = cbind(ASD_datExpr, ASD_expr[[i]]$datExpr)
  ASD_datMeta = rbind(ASD_datMeta, ASD_expr[[i]]$datMeta[,c("Sample",
                                                            "Subject_ID",
                                                            "Study",
                                                            "Dx", 
                                                            "Age",
                                                            "Sex", 
                                                            "Brain_Region",
                                                            "Brain_Lobe")])
}

ASD_datMeta = ASD_datMeta[-1,]

keep = (ASD_datMeta$Dx %in% c("CTL", "ASD"))
ASD_datExpr = ASD_datExpr[,keep]
ASD_datMeta = ASD_datMeta[keep,]
ASD_datMeta$Dx = factor(ASD_datMeta$Dx, levels = c("CTL", "ASD"))
ASD_datMeta$Subject_ID = as.factor(ASD_datMeta$Subject_ID)
ASD_datMeta$Brain_Region = as.factor(ASD_datMeta$Brain_Region)
ASD_datMeta$Sex = as.factor(ASD_datMeta$Sex)
ASD_datMeta$Study = as.factor(ASD_datMeta$Study)

# check duplicated samples
idx = !duplicated(ASD_datMeta$Sample)
ASD_datMeta = ASD_datMeta[idx,]
ASD_datExpr = ASD_datExpr[,idx]
all(colnames(ASD_datExpr)==rownames(ASD_datMeta))

plotMDS(log2(ASD_datExpr),
        col = as.numeric(as.factor(ASD_datMeta$Study)),
        pch = as.numeric(as.factor(ASD_datMeta$Dx)))


#### BP data ####

BP.files = dir("../normalized_studies/codes/", pattern = "_BP_")
n=length(BP.files)
BP_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("../normalized_studies/codes/",BP.files[[i]],sep = ""))
  if ("datExpr_bd" %in% ls()) {
    BP_expr[[i]]$datExpr = datExpr_bd
    BP_expr[[i]]$datMeta = datMeta_bd
    rm(datExpr_bd,datMeta_bd)
  } else {
  BP_expr[[i]]$datExpr = datExpr
  BP_expr[[i]]$datMeta = datMeta
  rm(datExpr,datMeta)
}
}

BP_expr[[4]]$datMeta$Sample = rownames(BP_expr[[4]]$datMeta)
genes= intersect(rownames(BP_expr[[1]]$datExpr), rownames(BP_expr[[2]]$datExpr))
for(i in 3:n) genes = intersect(genes, rownames(BP_expr[[i]]$datExpr))
for(i in 1:n) BP_expr[[i]]$datExpr = BP_expr[[i]]$datExpr[match(genes,rownames(BP_expr[[i]]$datExpr)),]

BP_datExpr = data.frame(row.names = genes)
BP_datMeta = data.frame(Sample=NA,
                        Subject_ID=NA,
                        Study=NA,
                        Dx=NA,
                        Age=NA,
                        Sex=NA,
                        Brain_Region=NA,
                        Brain_Lobe=NA)

for (i in 1:n) {
  BP_datExpr = cbind(BP_datExpr, BP_expr[[i]]$datExpr)
  BP_datMeta = rbind(BP_datMeta, BP_expr[[i]]$datMeta[,c("Sample",
                                                         "Subject_ID",
                                                         "Study",
                                                         "Dx",
                                                         "Age",
                                                         "Sex",
                                                         "Brain_Region",
                                                         "Brain_Lobe")])
}

BP_datMeta = BP_datMeta[-1,]

keep = (BP_datMeta$Dx %in% c("CTL", "BP"))
BP_datExpr = BP_datExpr[,keep]
BP_datMeta = BP_datMeta[keep,]
BP_datMeta$Dx = factor(BP_datMeta$Dx, levels = c("CTL", "BP"))
BP_datMeta$Subject_ID = as.factor(BP_datMeta$Subject_ID)

# check duplicated samples
idx = !duplicated(BP_datMeta$Sample)
BP_datMeta = BP_datMeta[idx,]
BP_datExpr = BP_datExpr[,idx]
all(colnames(BP_datExpr)==rownames(BP_datMeta))

plotMDS.default(log2(BP_datExpr), col = as.numeric(as.factor(BP_datMeta$Study)),
                pch = as.numeric(BP_datMeta$Dx))

#---------------------------------------------------------------------------------------------------------------------------#

####  MDD data  #### 

MDD.files = dir("../normalized_studies/codes/", pattern = "_MDD_")
n=length(MDD.files)
MDD_expr = vector(mode = "list", length = n)

for (i in 1:n) {
  load(paste("../normalized_studies/codes/",MDD.files[[i]],sep = ""))
  MDD_expr[[i]]$datExpr = datExpr
  MDD_expr[[i]]$datMeta = datMeta
  rm(datExpr,datMeta)
}


#---keep common genes only

genes= intersect(rownames(MDD_expr[[1]]$datExpr), rownames(MDD_expr[[2]]$datExpr))
for(i in 1:n) MDD_expr[[i]]$datExpr = MDD_expr[[i]]$datExpr[match(genes,rownames(MDD_expr[[i]]$datExpr)),]

MDD_datExpr = data.frame(row.names = genes)
MDD_datMeta = data.frame(Sample=NA,
                         Subject_ID=NA, 
                         Study=NA,
                         Dx=NA,
                         PMI=NA,
                         Age=NA,
                         Sex=NA, 
                         Brain_Region=NA,
                         Brain_Lobe = NA)

for (i in 1:n) {
  MDD_datExpr = cbind(MDD_datExpr, MDD_expr[[i]]$datExpr)
  MDD_datMeta = rbind(MDD_datMeta, MDD_expr[[i]]$datMeta[,c("Sample","Subject_ID",
                                                            "Study", 
                                                            "Dx",
                                                            "PMI", 
                                                            "Age",
                                                            "Sex",
                                                            "Brain_Region",
                                                            "Brain_Lobe")])
}

MDD_datMeta = MDD_datMeta[-1,]
keep = (MDD_datMeta$Dx %in% c("CTL", "MDD"))
MDD_datExpr = MDD_datExpr[,keep]
MDD_datMeta = MDD_datMeta[keep,]
MDD_datMeta$Dx = factor(MDD_datMeta$Dx,levels = c("CTL", "MDD"))
MDD_datMeta$Subject_ID = as.factor(MDD_datMeta$Subject_ID)
MDD_datMeta$Brain_Region = as.factor(MDD_datMeta$Brain_Region)
MDD_datMeta$Sex = as.factor(MDD_datMeta$Sex)
MDD_datMeta$Study = factor(MDD_datMeta$Study)

# check duplicated samples
idx = !duplicated(MDD_datMeta$Sample)
MDD_datMeta = MDD_datMeta[idx,]
MDD_datExpr = MDD_datExpr[,idx]

all(colnames(MDD_datExpr)==rownames(MDD_datMeta))


plotMDS.default(log2(MDD_datExpr), 
                col = as.numeric(as.factor(MDD_datMeta$Brain_Region)),
                pch = as.numeric(MDD_datMeta$Dx))
#--------------------------------------------------------------------------------------------------------------------#





## compile data into one object

multiExpr = vector(mode="list", length=8)
names(multiExpr) = c("AD","PD","PA","PSP","Scz","ASD","BP","MDD")

multiExpr[[1]]$datExpr = AD_datExpr;  multiExpr[[1]]$datMeta = AD_datMeta 
multiExpr[[2]]$datExpr = PD_datExpr; multiExpr[[2]]$datMeta = PD_datMeta 
multiExpr[[3]]$datExpr = PA_datExpr; multiExpr[[3]]$datMeta = PA_datMeta 
multiExpr[[4]]$datExpr = PSP_datExpr; multiExpr[[4]]$datMeta = PSP_datMeta 
multiExpr[[5]]$datExpr = Scz_datExpr; multiExpr[[5]]$datMeta = Scz_datMeta 
multiExpr[[6]]$datExpr = ASD_datExpr; multiExpr[[6]]$datMeta = ASD_datMeta
multiExpr[[7]]$datExpr = BP_datExpr; multiExpr[[7]]$datMeta = BP_datMeta 
multiExpr[[8]]$datExpr = MDD_datExpr; multiExpr[[8]]$datMeta = MDD_datMeta

# compile all in one data.frame

cov = c("Sample",
        "Subject_ID", 
        "Study",
        "Dx",
        "Age",
        "Sex", 
        "Brain_Region",
        "Brain_Lobe")

genes= intersect(rownames(multiExpr[[1]]$datExpr), rownames(multiExpr[[2]]$datExpr))
for(i in 3:8) genes = intersect(genes , rownames(multiExpr[[i]]$datExpr))

# union of samples
samples = colnames(AD_datExpr)
for(i in 2:length(multiExpr)) {
samples = union(samples,colnames(multiExpr[[i]]$datExpr))
  }

# dataframes
allExp = as.data.frame(matrix(NA,nrow = length(genes),
                              ncol = length(samples)))
colnames(allExp) = samples
rownames(allExp) = genes
allMeta = as.data.frame(matrix(NA,ncol = length(cov)))
colnames(allMeta) = cov
rownames(allMeta) = samples 

for(i in 1:length(multiExpr)) {
  meta = multiExpr[[i]]$datMeta
  exp  =  multiExpr[[i]]$datExpr
  exp  = exp[match(genes,rownames(exp)),]
  
  # expression
  # check for duplicated samples
  for (s in colnames(exp) ) {
    allExp[,s] = exp[,s] 
    if(!(s %in% allMeta$Sample)){
    allMeta = rbind(allMeta,meta[s,cov])
    }
  }
  # id = na.omit(match(colnames(allExp),colnames(exp)))
  # if(length(id)>0){ 
  #   exp = exp[,-id]
  #   } else {
  #     exp = exp 
  #   }
  # 
  # allExp = cbind(allExp, exp[match(genes,rownames(Exp)),])
  # 
  # rownames(meta) = NULL
  # meta = meta[match(colnames(exp),meta$Sample),]
  # allMeta = rbind(allMeta, dat[,cov],deparse.level = 0)
}

allMeta = allMeta[-1,]

# check dimnames
all(colnames(allExp) == rownames(allMeta))

# remove duplicated samples
idx = !duplicated(colnames(allExp))


##Save compiled expression and metadata for linear model analysis

save(file = "./codes/st2_lm_multiExpr.Rdata", multiExpr,allExp,allMeta)





