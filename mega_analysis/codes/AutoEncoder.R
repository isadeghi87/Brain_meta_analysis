## AutoEncoder.R

rm(list=ls()); options(stringsAsFactors=F)

library(ggplot2);library(sva); 
library(WGCNA);library(dplyr);
library(stringr);library(RColorBrewer)
library(ggpubr);library(edgeR)
library(Rmagic); library(Rtsne);
library(colorspace);library(sva);

### set working directory

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")

# options
enableWGCNAThreads()
allowWGCNAThreads()
condition = TRUE

### load normalized data from each disease

load("./codes/st2_lm_multiExpr.Rdata")
rm(multiExpr)

allMeta$Brain_Lobe = gsub("Basal_ganglia","Basal ganglia",
                          allMeta$Brain_Lobe)
allMeta$Brain_Lobe = as.factor(allMeta$Brain_Lobe)

var = c("Brain_Lobe","Brain_Region","Study","Sex")
for (i in var){allMeta[,i] = factor(allMeta[,i])}
allMeta$Dx = factor(allMeta$Dx, levels = c("CTL","AD","PD","PA","PSP","Scz","ASD","BP","MDD"))


#### filter genes ####

keep = rowSums(cpm(allExp)>0.5)>= 0.3 * ncol(allExp)
datExp = allExp[keep,]


# train and test data are the same here
## ideally should split into separate train and test sets
## to avoid overfitting
m = as.matrix(t(datExp))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
dat <- apply(m, 2, range01)
rownames(dat) <- rownames(m)
range(dat)

x_train <- x_test <- dat

## Load Keras (https://keras.rstudio.com/)
library(keras)
K <- keras::backend()

## Deep learning model
input_size <- ncol(x_train) ## 1000 genes
hidden_size <- 10 ## 10 dimensional hidden layer
code_size <- 2 ## 2 dimensional encoding
input <- layer_input(shape=c(input_size))
hidden_1 <- layer_dense(input, hidden_size) %>% 
  layer_activation_leaky_relu() %>%
  layer_dropout(rate=0.1)
code <- layer_dense(hidden_1, code_size) %>% 
  layer_activation_leaky_relu()
hidden_2 <- layer_dense(code, units=hidden_size) %>% 
  layer_activation_leaky_relu()
output <- layer_dense(hidden_2, units=input_size, activation="sigmoid")

## input and output should be the same
autoencoder <- keras_model(input, output)
## encoder from input to code space
encoder <- keras_model(input, code)

## Learn
autoencoder  %>% compile(optimizer='adam', 
                         loss='cosine_proximity',
                         metrics='mae')
autoencoder %>% fit(
  x_train, x_train, 
  shuffle=TRUE, 
  epochs=1000, 
  batch_size=100, 
  validation_data=list(x_test, x_test)
)

############### Plot

## predict code space using deep learning  model
x_test_encoded <- predict(encoder, 
                          x_test, 
                          batch_size=100)
emb2 <- x_test_encoded
rownames(emb2) <- rownames(x_test)
colnames(emb2) <- c('latent 1', 'latent 2')

## plot
plot(emb2, 
     col=rainbow(length(levels(group2)))[group2], 
     pch=16,
     main='Autoencoder: 2D code layer')


library(devtools)
install_github('rstudio/reticulate',force=T)
library(reticulate)
library(tensorflow)
install_tensorflow(version= "1.1.0")
install_github("rstudio/keras",force=T)
library(keras)
keras::install_keras()