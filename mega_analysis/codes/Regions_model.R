##Region_DxCorr.R
##Here we calculate the correlation of brain regions with disease
## to see which region has more effect 
rm(list=ls())
library(dplyr)
library(pheatmap)
library(grid)
library(caret)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(Rtsne)
library(pROC)
library(ggrepel)
library(ggthemes)
library(WGCNA);

# library(biomaRt);
# library(ggpubr);
# library(ComplexHeatmap)
##
error = F; warnings = F; StringsAsFactors = F

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")
plot = TRUE

#######################
# Loading data
#######################
#load normalized data from each disease
files = dir("./codes/", pattern = "_metadata")
n=length(files)
multiExpr = vector(mode = "list", length = n)
names(multiExpr) = c("AD","ASD","BP","MDD","PA","PD","PSP","Scz")
for (i in seq_len(length(multiExpr))) {
  load(paste("./codes/",files[[i]],sep = ""))
  multiExpr[[i]]$datExpr = datExpr
  multiExpr[[i]]$datMeta = datMeta
  rm(datExpr,datMeta)
}
multiExpr[["PD"]]$datMeta$Brain_lobe = "Frontal"
load("./codes/st2_lm_multiExpr.Rdata")
rm(allMeta,allExp)

##rename normal
for (i in 1:length(multiExpr)) {
  multiExpr[[i]]$datMeta$Brain_Lobe = gsub("Basal_ganglia","Basal ganglia",multiExpr[[i]]$datMeta$Brain_Lobe)
  }

## load allexpr data to filter DE genes
allExpr = readRDS("./codes/regions_allExpr.rds")

##rename basal ganglia
allExpr$Region = gsub("Basal.ganglia", "Basal ganglia",allExpr$Region)

##options
enableWGCNAThreads()
allowWGCNAThreads()
condition =TRUE

# 
# ##parallel computation
require(doParallel)
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



# First we build a raw glm
res = data.frame(dx = NA,lobe = NA,accuracy = NA,
                 kappa = NA, sensitivity = NA, specificity = NA)

if(condition){
  for (dx  in names(multiExpr)) {
    ## filter the real data
    datMeta = multiExpr[[dx]]$datMeta
    datExpr = multiExpr[[dx]]$datExpr
    #
    for (lobe in unique(datMeta$Brain_Lobe)) {
      id = which(datMeta$Brain_Lobe == lobe)
      dat = datMeta[id,]
      exp = datExpr[,id]
      #look for significant genes first
      sig = allExpr[allExpr$Disease == dx & allExpr$Region == lobe,]
      sig = subset(sig, p.value < 0.05 & abs(logFC)>0.5)
      idx = unique(sig$gene)
      exp = exp[idx,]
      
  metadata = dat[,"Dx", drop = FALSE]
  data = cbind(metadata,t(exp))
# Partition data into train and test
set.seed(1234)
train.idx <- createDataPartition(y = data$Dx, p = 0.5, list = FALSE)[,1]
X_train <- data[train.idx, ]
X_test <- data[-train.idx, ]

# Check performance just on the partition (model fitted on train and tested over test)
# mod_fit <- train(Dx ~. ,  data = X_train, 
                 # method="rf")
# confusionMatrix(data = predict(mod_fit, X_test), reference = as.factor(X_test$Dx))

# Perform crossvalidation
ctrl1 <- trainControl(method = "cv",
                      number = 10,
                      # summaryFunction = twoClassSummary,
                      verboseIter = TRUE,
                      allowParallel = T,
                      savePredictions = TRUE,
                      classProbs = TRUE)
mod_fit_cv <- train(Dx ~ ., data = X_train, method = "rf",
                    trControl = ctrl1)

sens = confusionMatrix(data = predict(mod_fit_cv, X_test), reference = as.factor(X_test$Dx))


# -- ! CHECK THIS ! --
# selectedIndices <- mod_fit_cv$pred$mtry ==2
# roctest <- roc(mod_fit_cv$pred$obs[selectedIndices], mod_fit_cv$pred$Control[selectedIndices])
# plot(roctest)
accuracy = round(unname(sens$overall["Accuracy"]),2)
kappa = round(unname(sens$overall["Kappa"]),2)
# auc = round(roctest$auc,2)
sensitivity = round(unname(sens$byClass["Sensitivity"]),2)
specificity = round(unname(sens$byClass["Specificity"]),2)
res = rbind(res,c(dx,lobe,accuracy,kappa,sensitivity,specificity))

# mod[lobe,i] = accuracy
# auc.df[lobe,i] = auc
# kapp[lobe,i] = kappa

}
}
}

res = res[-1,]

# ########
####tSNE
# ########

col = c("blue","green","red","purple",
        "brown","magenta","lightgreen","black","grey")
names(col) = c(unique(allExpr$Disease),"CTL")

plots = vector("list",length = 56)
dx = c("AD","PD","PA","PSP","Scz","ASD","BP","MDD")
lobe = unique(allExpr$Region)
names(plots) = paste0(rep(dx,1,each=7),"-",rep(lobe,7))

for (dx  in names(multiExpr)) {
  ## filter the real data
  datMeta = multiExpr[[dx]]$datMeta
  datExpr = multiExpr[[dx]]$datExpr
  #
  for (lobe in unique(datMeta$Brain_Lobe)) {
    id = which(datMeta$Brain_Lobe== lobe)
    dat = datMeta[id,]
    exp = datExpr[,id]
    #look for significant genes first
    sig = allExpr[allExpr$Disease == dx & allExpr$Region == lobe,]
    sig = subset(sig, p.value < 0.05 & abs(logFC)>0.5)
    idx = unique(sig$gene)
    exp = exp[idx,]
    #tsne
      set.seed(12)
      fit.tsne = prcomp(t(exp), center = F, scale. = F)
      n = nrow(t(exp))
      fit.tsne <- Rtsne(fit.tsne$x[,1:5], 
                        dims = 2, theta = 0,check_duplicates = F,
                        max_iter = 3000,
                        perplexity =n/4)
      tsne.coordinates <- fit.tsne$Y[,1:2,drop = F]
      tsne.coordinates <- as.data.frame(tsne.coordinates)
      tsne.coordinates$Group <- dat$Dx
      #plot
      name = paste0(dx,"-",lobe)
      p1 = ggplot(tsne.coordinates, aes(x = V1, y = V2)) +
          geom_point(aes(color =  Group))+
        scale_color_manual(values = col)+
          theme_classic2()+
          labs(x = "tSNE1",
               y = "tSNE2",
               title = paste(dx,"-",lobe))+
          theme(text = element_text(size = 10),
                plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))

      plots[[name]] = p1
      
      }

    }
  
plots = plots[-which(sapply(plots, is.null))]
g1 = ggarrange(plotlist = plots,ncol = 4,nrow = 6,common.legend = F)
  ggsave(g1, file = "./results/figures/model/tsne-plots.pdf",
              width = 15 , height = 10, useDingbats = FALSE)

results = res
results$dx = factor(results$dx, levels = c("AD","PD","PA","PSP",
                                             "Scz","ASD","BP","MDD"))


for(i in 3:6){
  results[,i]= as.numeric(results[,i])
}

colnames(results)= gsub("lobe","region",colnames(results))

#plot

p2 = ggplot(results,aes(x = specificity, y = sensitivity, group = 1))+
  geom_point(aes(color = region), size = 3,show.legend = F)+
  geom_text_repel(aes(label = paste(region,"\n(",accuracy*100,"%)")),
            size = 2.5, seed = 12)+
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.25),limits = c(0,1.2))+
  facet_wrap(~dx, nrow = 2)+
  theme_clean(base_size = 14)+
  theme(text = element_text(size = 14))


ggsave(p2, filename = "./results/figures/model/regionsAccuracyPlots.pdf",
       width = 9, height = 5)

## MEAN ACCURACY FOR EACH DISEASE

mean.Acc = results %>% group_by(dx) %>% summarise(mean = mean(accuracy)) %>%
  arrange(desc(mean))

mean.Acc.reg = results %>% group_by(region) %>% summarise(mean = mean(accuracy)) %>%
  arrange(desc(mean))

## 
write.csv(results,file = "./results/tables/model_accuracy.csv")

##save the data
save.image("./codes/regions_model.Rdata")
