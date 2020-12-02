## sample_size_deg_corr.R
## check correlation of DEGs numbers with sample size 
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/mega_analysis/")


# load DEGs table 
degTable = read.csv("./results/tables/Dis_specific_DEGs.csv",
                    row.names = 1)
degTable = as.data.frame(t(degTable[,1:8]))

# load sample size 
sampleSize = read.csv("./results/tables/sample_size.csv",row.names = 1)
sampleSize = as.numeric(t(sampleSize))

degTable$sample = sampleSize

# correlation test ####
cor_res = as.data.frame(matrix(NA,ncol = 4, nrow = 2))
colnames(cor_res) = colnames(degTable[,1:4])
rownames(cor_res) = c("rho","p.value")

n = 4 # number of comparison for each deg type
for (i in 1:n){
cor = cor.test(degTable[,5],
                   degTable[,i],
                   method = "spearman",
                   use = "pairwise.complete.obs")
 p = as.numeric(cor$p.value)
 rho = as.numeric(cor$estimate)
 cor_res["rho",i] = signif(rho,2)
 cor_res["p.value",i] = signif(p,2)
 
}

p = ggtexttable(cor_res,theme = ttheme(base_style = "lGreenWhite",
                                       base_size = 8))
ggsave(p,filename = "./results/figures/differential_expression/sample_deg_cor.pdf",
       width = 5 ,height = 2)
write.csv("./results/tables/sampleSize_deg_corr.csv",x = cor_res)
