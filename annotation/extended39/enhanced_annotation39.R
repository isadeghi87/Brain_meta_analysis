###--enhanced_annotation39.R
## Here we prepare enhanced annotation file containing
## Gencode version 39 and long-read seq from PacBIO
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/annotation/extended39/")

#libraries ####

library(reshape); library(stringr); library(dplyr);library(edgeR);
library(readr);  library(tidyverse);library(stringr);library(sva);
library(WGCNA);library(biomaRt)

#options
options(stringsAsFactors = F)

## new annotation
new_annot = read.delim(file = gzfile("./merged_annotation39.gtf.gz"),header = F)

## remove spike-in from annotation file
new_annot = new_annot[!grepl("SIRVome_isoforms",new_annot$V1),]

##  transcript id
new_annot$transcript_id = str_extract(new_annot$V9,"ENST[0-9]++\\.[0-9]++")


### Gencode39 annotation
gencod = read.delim(file = gzfile("./gencode.v39.annotation.gtf.gz"),header = F,skip = 5)
## keep gene ID from gencode annotation for those with ENST ####
gencod$gene_id = str_extract(gencod$V9,"ENSG[0-9]++\\.[0-9]++")

##  transcript id
gencod$transcript_id = str_extract(gencod$V9,"ENST[0-9]++\\.[0-9]++")

# gene type
gencod$gene_type = str_split(gencod$V9,"gene_type ",simplify = T)[,2]
gencod$gene_type = str_split(gencod$gene_type,";",simplify = T)[,1]

# gene name
gencod$gene_name = str_split(gencod$V9,"gene_name ",simplify = T)[,2]
gencod$gene_name = str_split(gencod$gene_name,";",simplify = T)[,1]

saveRDS(gencod,file = "C:/Users/crgcomu/Desktop/Iman/Brain_meta/annotation/gencode.v39.rds")

## remove those without transcript id 
# g  = gencod[!is.na(gencod$transcript_id),]

new_annot$gene_id = gencod$gene_id[match(new_annot$transcript_id,gencod$transcript_id)]

id = !(is.na(new_annot$transcript_id))
new_annot$V9[id] = str_replace_all(new_annot$V9[id], pattern = "gene_id TM_[0-9]++",
                                   replacement = paste0("gene_id ",
                                                        new_annot$gene_id[id]))

new_annot$V9[id] = str_replace_all(new_annot$V9[id], pattern = "transcript_id TM_[0-9]++",
                                   replacement = paste0("transcript_id ", new_annot$transcript_id[id]))

# remove gene id and transcript id columns to export file
new_annot$gene_id = NULL
new_annot$transcript_id = NULL

## quote values of transcript id and gene id
new_annot$V9 = gsub("(\\S+);", "\"\\1\";", new_annot$V9)

write.table(new_annot,file = "./enhanced_annotation39.gtf",col.names = FALSE, row.names = FALSE, quote = FALSE,sep = "\t")
