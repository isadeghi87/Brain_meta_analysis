#projection_score.R

options(stringsAsFactors=FALSE)
library(ggplot2); library(mada); library(reshape); library(RColorBrewer)
library(ggrepel); library(ggpubr)
setwd("/nfs/users2/rg/isadeghi/data/combined/mega_analysis/")

###########
# OPTIONS #
###########


thresholds = c(0.1,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8)#thresholds
B= 10#Number of row permutations 
S = 3#Number of principal components to compute the projection score
  
##############
# BEGIN
##############

# Read the data

data = read.csv("./results/tables/datExpr.csv", header = T, row.names = 1)
colnames(data) = gsub("X","",colnames(data))

# Read the statistics
stats = read.table("./results/tables/GeneVariance.tsv", header = T)
# Order them to be in the same order as the input matrix row names
stats = stats[match(stats[,1], rownames(data)),]

# Normalize the variance by the maximum variance
variancen = stats[,2]/max(stats[,2], na.rm=T)
#print(head(variancen))
names(variancen) <- stats[,1]

# Function to compute the alpha_2 measure
alpha_2 = function(lambda, S) {
  sqrt(sum(lambda[1:S]^2)/sum(lambda^2))
}

# Initialize the vector of projection scores
proj_scores = array(numeric(0))

# Iterate over different variance thresholds
for (var_t in thresholds) {
  
  cat("Variance threshold: ", var_t, "\n")
  
  # Get the actual submatrix for this variance threshold
  m_t = data[which(variancen>var_t),]
  
  # Compute the pca on this submatrix
  pca1 = prcomp(t(m_t), center=FALSE, scale.=FALSE)
  
  # Obtain a matrix of lambdas (sdev) for each permutation
  # Rows are the components and columns are the iterations
  set.seed(123)
  lambda = replicate(B, prcomp(apply(m_t, 1, sample), center=FALSE, scale.=FALSE)$sdev) 
  # Count how many times the stdev for each component in the permutation is lower than the observed one
  lambda_counts = rowSums(apply(lambda, 2, function(x) pca1$sdev>=x))
  
  # If at least one observed stdev is not higher than 95% of the permutated lambdas 
  # the submatrix does not support S, and the projection score is assigned to NA
  if (sum(lambda_counts[1:S] < 0.05*10) != 0) {
    proj_score = NA
  } else {
    exp_l = mean(apply(lambda, 2, alpha_2, S=S))
    obs_l = alpha_2(pca1$sdev, S)
    proj_score = obs_l - exp_l
  }
  proj_scores = c(proj_scores, proj_score)
}


###############
# OUTPUT
###############

selected = names(which(variancen >= thresholds[which.max(proj_scores)]))

write.table(selected, file = "./results/tables/projectionScore.txt",
            quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
#subset variance data for selected genes
subdat = stats[match(selected, stats$gene),] %>% arrange(desc(var))
top = subdat[1:50,]

write.table(top, file = "./results/tables/projection_top50.txt",
            quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
# Plot the projection score as a function of the variance
df = data.frame(thresholds, proj_scores)
theme_set(theme_bw(base_size=18))

gp = ggplot(df, aes(x=as.numeric(thresholds), y=proj_scores)) + geom_point(size=3) + geom_line()
gp = gp + labs(x="Thresholds (fraction of max)")





