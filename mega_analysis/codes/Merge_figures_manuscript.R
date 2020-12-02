# Document: Merge_figrues_manuscript.R
# This is file contains scripts for merging figures for the manuscritp
# The data from different Rdata files are loaded here for each figures
# and we merge them by ggplot
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis")

#libraries
library(ggplot2);library(pheatmap); library(ComplexHeatmap)
library(RColorBrewer); library(ggpubr); library(circlize)
library(igraph);library(WGCNA); library(ggplotify)

fig2a = readRDS("./codes/Fig2a_scplot.rds")
fig2c =readRDS("./codes/fig2c.rds")
ggarrange(fig2a,fig2c,labels = c("A","C"))