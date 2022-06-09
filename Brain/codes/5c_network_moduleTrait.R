#5c_network_moduleTrait.R

library(pSI); library(pSI.data)
library(WGCNA);library(ggplot2); 
library(reshape); library(nlme)
condition = TRUE

setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")
# load final network
load("./codes/finalNetwork.RData")

all_colors = unique(colors)
all_colors = all_colors[!grepl("grey",all_colors)]
all_genes = colors
names(all_genes) = rownames(datExpr)
datMeta$Age[is.na(datMeta$Age)] = mean(datMeta$Age, na.rm = T)


#Step 1 - Calculate module-trait P values, beta, and SEM
moduleTraitP = matrix(NA, nrow = length(all_colors), ncol=12)
colnames(moduleTraitP) = c("ANOVA", 
                           "AD", "ASD", "Scz", "BP", "MDD", "PD", "PSP", "PA",
                           "Sex",  "PMI", "Age")
rownames(moduleTraitP) = all_colors
moduleTraitB = moduleTraitSE = moduleTraitP

if(condition){
for (m in all_colors) {
  me_name = paste("ME", m, sep="")
  me = MEs$eigengenes[[me_name]]
  i = which(m == rownames(moduleTraitP))
  
  mixedmodel = lme(me ~ Dx +  Sex +  PMI, data = datMeta, random = ~1|Subject_ID)
  moduleTraitP[i,"ANOVA"] = anova(mixedmodel)["Dx","p-value"]
  mixedmodel = summary(mixedmodel)$tTable
  for(var in c("AD", "ASD", "Scz", "BP", "MDD","PD", "PSP", "PA", "Sex")) {
    moduleTraitP[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 5]
    moduleTraitB[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 1]
    moduleTraitSE[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 2]
  }
  for(cov in c("Sex", "PMI", "Age")) {
    mod = summary(lm(me ~ datMeta[,cov]))
    moduleTraitP[i,cov] = mod$coefficients[2, 4]
    moduleTraitB[i,cov] = mod$coefficients[2, 1]
    moduleTraitSE[i,cov] = mod$coefficients[2, 2]
  }
}
}

moduleTraitP.fdr = p.adjust(moduleTraitP, "fdr")
dim(moduleTraitP.fdr) = dim(moduleTraitP)
dimnames(moduleTraitP.fdr) = dimnames(moduleTraitP)

out = cbind(moduleTraitB, moduleTraitP, moduleTraitP.fdr)
colnames(out)[1:12] = paste("Beta.", colnames(out)[1:12], sep="")
colnames(out)[13:24] = paste("P.", colnames(out)[13:24], sep="")
colnames(out)[25:36] = paste("FDR.", colnames(out)[25:36], sep="")
out = out[,-1]
write.csv(file= "./results/tables/Table_moduleTraitB.csv",out)


col_idx = c(2:9)
#bar plot
bpdata = melt(moduleTraitB[, col_idx])
semdata = melt(moduleTraitSE[, col_idx])
pdata = melt(moduleTraitP.fdr[, col_idx])
bpdata$sem = semdata$value
bpdata$p = pdata$value
bpdata$p.symbol = ""
bpdata$p.symbol[bpdata$p<0.05] = "*"
bpdata$p.symbol[bpdata$p<0.01] = "**"
bpdata$p.symbol[bpdata$p<0.001] = "***"

colnames(bpdata) = c("Module", "Diagnosis", "beta", "SEM", "p", "p.symbol")
bpdata$Module = factor(bpdata$Module)


#----------------------------plot the differential expression for modules
byModule = ggplot(bpdata,
                        aes(x = Diagnosis,
                            y = beta,
                            fill = Module,
                            group = Module,
                            label = p.symbol))+
  facet_wrap(~ Module,
             scales = "free",ncol = 4) + 
  geom_bar(stat="identity", 
           position = position_dodge(),
           color="black") +   
  geom_errorbar(aes( ymin = (beta - SEM),
                     ymax=(beta + SEM)),
                position = position_dodge(.9),
                size = 0.3,
                width =.3) +
  scale_fill_manual(name="Diagnosis", values = levels(bpdata$Module)) + 
  labs(y="beta", 
       x="",
       title = "Differential expression of modules") + 
  geom_text(color = "red",
            size = 4, 
            aes(y = beta + sign(beta)*SEM + sign(beta)*.01), 
            position=position_dodge(.1))  + 
  scale_x_discrete() + 
  theme(legend.position = "none", 
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(angle=45),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(size = 13))

byModule
ggsave(plot = byModule,
       "./results/figures/ExpressionByModule.pdf", 
       height = 8, 
       width = 10,
       device = "pdf")

#SEM of modules
SEMbyModule = ggplot(bpdata,
                            aes(
                              x = Diagnosis,
                              y = as.numeric(SEM),
                              fill = Module,
                              group = Module,
                              label = p.symbol
                            )) + facet_wrap( ~ Module, ncol = 4, scales = "free_x") +
  geom_bar(stat = "identity",
           position = position_dodge(),
           color = "black") +
  theme_minimal() + 
  scale_fill_manual(name = "Diagnosis", values = levels(bpdata$Module)) +
  labs(y = "Std Error of Beta", x = "") +
  geom_text(color="red", size=4, aes(y=SEM *1.1), position=position_dodge(.9))  +
  scale_x_discrete() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 13)
  )


ggsave(plot = SEMbyModule,
       "./results/figures/ExpressionByModuleSEM.pdf", 
       height = 8, 
       width = 10,
       device = "pdf")


#----plot by disease
byDisease = ggplot(bpdata, 
                   aes(x = Module,
                       y = beta,
                       fill = Module,
                       group = Module,
                       label = p.symbol
                       )) + 
  facet_wrap( ~ Diagnosis, ncol = 1) +
  geom_bar(stat = "identity",
           position = position_dodge(),
           color = "black") +
  geom_errorbar(
    aes(ymin = (beta - SEM), ymax = (beta + SEM)),
    position = position_dodge(.9),
    size = 0.3,
    width = .3
  ) +
  #theme_minimal() + 
  scale_fill_manual(name = "Dx", values = levels(bpdata$Module)) +
  labs(y = "beta", 
       x = "",
       title = "Differential expression of modules") +
  geom_text(
    color = "red",
    size = 4,
    aes(y = beta + sign(beta) * SEM + sign(beta) * .01),
    position = position_dodge(.9)
  )  +
  scale_x_discrete() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # text = element_text(size = 14),
    strip.text.x = element_text(size = 16)
  )

byDisease 
saveRDS(byDisease, "./codes/Fig4d_ModuleExpression.rds")
  ggsave(plot = byDisease,
    filename = "./results/figures/ExpressionModbyDisease.pdf",
    width = 4,
    height = 10)


  #---------Cell type specific module expression
cellModule = vector(mode = "list", length = 4)
names(cellModule) = c("Oligodendrocytes", "Astrocytes","Neurons", "Microglia")
cellModule[["Oligodendrocytes"]] = bpdata[bpdata$Module =="brown",]
cellModule[["Astrocytes"]] = bpdata[bpdata$Module =="yellow",]
cellModule[["Neurons"]] = bpdata[bpdata$Module %in% c("magenta","turquoise"),]
cellModule[["Microglia"]] = bpdata[bpdata$Module =="purple",]

for ( i in names(cellModule)){
  cellModule[[i]]$cell = i
}

cellmod = do.call("rbind",cellModule)  
cellmod$cell = factor(cellmod$cell, levels = names(cellModule))
# file = data.frame(image = list.files("./results/figures/cell_schematics", 
#                                         full.names = T), 
#                       stringsAsFactors = F)
# file$name =  c("Astrocytes", "Neurons","Oligodendrocytes" )
# cellmod2$image = as.character(cellmod2$cell)
# cellmod2$name = as.character(cellmod2$cell)
# cellmod2$image= gsub("Oligodendrocytes",file$image[3], cellmod2$image)
# 
# cellmod2$image= gsub("Astrocytes",file$image[1], cellmod2$image)
# 
# cellmod2$image= gsub("Neurons",file$image[2], cellmod2$image)

p <-  ggplot(cellmod,
             aes(x = Diagnosis,
                 y = beta,
                 group = Module,
                 label = p.symbol))+
  geom_bar(stat="identity", 
           fill = cellmod$Module,
           position = position_dodge(),
           color="black") +   
  geom_errorbar(aes( ymin = (beta - SEM),
                     ymax=(beta + SEM)),
                position = position_dodge(.9),
                size = 0.3,
                width =.3) +
  labs(y="beta", 
       x="",
       title = "Brain cell types modules") + 
  geom_text(color = "red",
            size = 2, 
            aes(y = beta+ sign(beta)*SEM + sign(beta)*.005), 
            position=position_dodge(1))  + 
  scale_x_discrete() + 
  scale_fill_manual(values =cellmod$Module)+
  theme(legend.position = "none", 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle=45),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        text = element_text(size = 12)) +
  # theme_minimal()+
facet_wrap(~cell, ncol = 1)

p

# name= unique(as.character(cellmod2$cell))
# g = list()
# for(i in 1:length(name)){
#   dat = cellmod2[cellmod2$cell==name[i],]
#   img = readPNG(source = dat$image[i])
#   g[[i]] =  rasterGrob(img, interpolate=TRUE)
# 
#   PLOT = p +
#     annotation_custom(grob=g[[i]],
#                       xmin=i-.5,
#                       xmax=i+.5,
#                       ymin=-Inf,
#                       ymax=-0.03)
# }
# 
# PLOT + theme(axis.text.x = element_text(angle=60, vjust=-.001))
# 
ggsave(filename = "./results/figures/CellTypesModule.pdf", 
       height = 7, width = 3 )  


save.image("./codes/st5c_NetModTrait.Rdata")
