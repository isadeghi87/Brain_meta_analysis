## Here we do linear discrimination analysis usign lda from MASS package
### set working directory

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")

## Load the library containing the linear discriminant analysis function
library(MASS)

## options
cols <- c("lightgrey","aquamarine","red","blue","chartreuse",
          "brown","blueviolet","yellow","orange")

# load top variable genes ranking
gene.ranks = read.csv("./results/tables/allExpr_gene_ranks.tsv",
                      sep = "\t", row.names = 1)

## Create a list of gene names sorted by decreasing variance
sorted.names <- rownames(gene.ranks)[order(gene.ranks$var, decreasing=TRUE)]

# load data
load("./codes/Combined.Rdata")

## Select the 100 top-ranking genes sorted by decreasing variance
top.variables <- 200
selected.genes <- sorted.names[1:top.variables]

## Train the classifier
sample.labels = as.character(allMeta$Dx)
names(sample.labels) = colnames(datExp.comb)

lda.classifier <- lda(t(datExp.comb[selected.genes,]),sample.labels,CV=FALSE) 

## Use the MASS:lda() function with the cross-validation option
lda.loo <- lda(t(datExp.comb[selected.genes,]),sample.labels,CV=TRUE) 

## Collect the LOO prediction result in a vector
loo.predicted.class <- as.vector(lda.loo$class)
print(loo.predicted.class)
table(loo.predicted.class)

pred.lda = predict(object = lda.classifier,
                   newdata = t(datExp.comb[selected.genes,]))

dataset = as.data.frame(pred.lda$x)
dataset = dataset %>% mutate(Dx = allMeta$Dx,
                             study = allMeta$Study,
                             region = allMeta$Brain_Lobe)

p = ggplot(dataset, aes(x = LD3, y = LD2,
                        color = Dx))+
  geom_point()+
  scale_color_manual(values = cols)

ggsave("./results/figures/pca/LDA_topVar_allExp.pdf",plot = p,
       width = 8,height = 5)

#3d plot
library(scatterplot3d)

pdf("./results/figures/pca/lda_3d.pdf", width = 10, height = 10)
par(mfrow=c(2,1))
with(dataset, {
  s3d <- scatterplot3d(
    x = LD1, 
    y = LD3, 
    z = LD2,
    color = cols[as.numeric(as.factor(allMeta$Dx))], 
    pch = 19, 
    # type = "h",
    lty.hplot = 1,
    angle = 45,
    scale.y = .7,
    main = "LDA-combined datasets",
    xlab = "LD1",
    ylab = "LD3",
    zlab = "LD2")
  #add legend
  legend(#location  
    "topright",
    inset=  0.05,
    # suppress legend box, shrink text 50%
    bty="n",
    cex=.7,
    title="Diagnosis",
    levels(allMeta$Dx),
    fill=cols)
    })
  
  with(dataset, {
    s3d <- scatterplot3d(
      x = LD1, 
      y = LD3, 
      z = LD2,
      color = cols[as.numeric(as.factor(allMeta$Dx))], 
      pch = 19, 
      # type = "h",
      lty.hplot = 1,
      angle = -45,
      scale.y = 0.7,
      main = "LDA-combined datasets",
      xlab = "LD1",
      ylab = "LD3",
      zlab = "LD2") 
    
  
  
})

dev.off()
## Build a contingency table of known versus predicted class
lda.loo.xtab <- table(sample.labels, loo.predicted  )
print(lda.loo.xtab)

## Display the contingency table as a heat map
image(lda.loo.xtab)
library(lattice)
levelplot(lda.loo.xtab)

## Compute the number of hits 
## Compute the hit rate
hits <- sample.labels == loo.predicted.class
errors <- sample.labels != loo.predicted.class
## (we need to omit NA values because LDA fails to assign a group to some objects).
(nb.hits <- sum(na.omit(hits))) ## this should give 2688
(nb.pred <- length(na.omit(hits))) ## This should give 4452
(hit.rate <- nb.hits / nb.pred ) ## This should give 0.57

