#5c_network_moduleTrait.R

library(pSI); library(pSI.data)
library(WGCNA);library(ggplot2); 
library(reshape); library(nlme);library(dplyr)
library(magrittr);library(ggthemes)
condition = TRUE

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")
# load final network
load("./codes/finalNetwork.RData")

all_colors = unique(mods)
all_colors = all_colors[!grepl("M0",all_colors)]
all_genes = colors
names(all_genes) = rownames(datExpr)
allMeta$Age[is.na(allMeta$Age)] = mean(allMeta$Age, na.rm = T)

#Step 1 - Calculate module-trait P values, beta, and SEM
moduleTraitP = matrix(NA, nrow = length(all_colors),
                      ncol=11)
colnames(moduleTraitP) = c("ANOVA", 
                           "AD","PD", "PSP", "PA", "Scz","ASD", "BP", "MDD",
                           "Sex", "Age")
rownames(moduleTraitP) = all_colors
moduleTraitB = moduleTraitSE = moduleTraitP
if(condition){
  for (m in all_colors) {
    me_name = gsub("M","ME",m)
    me = MEs$eigengenes[[me_name]]
    i = which(m == rownames(moduleTraitP))
    
    mixedmodel = lme(me ~ Dx + Sex + Age, data = allMeta, random = ~1|Subject_ID)
    moduleTraitP[i,"ANOVA"] = anova(mixedmodel)["Dx","p-value"]
    mixedmodel = summary(mixedmodel)$tTable
    for(var in c("AD", "ASD", "Scz", "BP", "MDD","PD", "PSP", "PA")) {
      moduleTraitP[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 5]
      moduleTraitB[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 1]
      moduleTraitSE[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 2]
    }
    for(cov in c("Sex", "Age")) {
      mod = summary(lm(me ~ allMeta[,cov]))
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
colnames(out)[1:11] = paste("Beta.", colnames(out)[1:11], sep="")
colnames(out)[12:23] = paste("P.", colnames(out)[12:23], sep="")
colnames(out)[24:33] = paste("FDR.", colnames(out)[24:33], sep="")
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
bpdata = bpdata %>%  arrange(.,Module)
bpdata$Module = as.character(bpdata$Module)

# plot the differential expression for disease
col_df = as.data.frame(cols)
col_df$module = unique(mods) %>% sort
bpdata$color = col_df$cols[match(bpdata$Module,col_df$module)]
bpdata$order = as.numeric(gsub("M","",bpdata$Module))
bpdata = bpdata[order(bpdata$order),]
bpdata$Module = factor(bpdata$Module,levels = unique(bpdata$Module))
bpdata$color = factor(bpdata$color,levels = unique(bpdata$color))

byModule = ggplot(bpdata,
                  aes(x = Diagnosis,
                      y = beta,
                      fill = Module,
                      # group = Module,
                      label = p.symbol))+
  facet_wrap(~ Module,
             scales = "free",ncol = 2) + 
  ylim(c(-0.0035,0.0035))+
  geom_bar(stat="identity", 
           position = position_dodge(),
           color="black") +   
  geom_errorbar(aes( ymin = (beta - SEM),
                     ymax=(beta + SEM)),
                position = position_dodge(.009),
                size = 0.3,
                width =.3) +
  scale_fill_manual(name="Diagnosis", values = cols[-1]) + 
  labs(y="beta", 
       x="") + 
  geom_text(color = "red",
            size = 4, 
            aes(y = beta + sign(beta)*SEM + sign(beta)*.0005), 
            position=position_dodge(0.1))  + 
   scale_x_discrete() + 
  theme_bw()+
  theme(legend.position = "none", 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle=45),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 13))

byModule
ggsave(plot = byModule,
       "./results/figures/network_analysis/ExpressionByModule.pdf", 
       height = 10, 
       width = 8,
       device = "pdf")

#----plot by disease
byDisease = ggplot(bpdata, 
                   aes(x = Module,
                       y = beta,
                       fill = Module,
                       # group = Module,
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
  scale_fill_manual(name = "Dx", values = cols) +
  labs(y = "beta", 
       x = "",
       title = "") +
  geom_text(
    color = "red",
    size = 4,
    aes(y = beta + sign(beta) * SEM + sign(beta) * .0005),
    position = position_dodge(.9)
  )  +
  scale_x_discrete() +
  ylim(-0.0035,0.0035)+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 14,angle = 90),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # text = element_text(size = 14),
    strip.text.x = element_text(size = 15)
  )

byDisease 
ggsave(plot = byDisease,
       filename = "./results/figures/network_analysis/ExpressionModbyDisease.pdf",
       width = 6,
       height = 10)


save.image("./codes/st5c_NetModTrait.Rdata")
