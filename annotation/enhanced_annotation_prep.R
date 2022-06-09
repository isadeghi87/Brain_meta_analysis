###--enhanced_annotation_prep.R
## Here we prepare enhanced annotation file containing
## Gencode version 28 and long-read seq from PacBIO
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/annotation/")

#libraries ####

library(reshape); library(stringr); library(dplyr);library(edgeR);
library(readr);  library(tidyverse);library(stringr);library(sva);
library(WGCNA);library(biomaRt)

#options

options(stringsAsFactors = F)

## read bed file 
# biotype = read.delim("./Hv3_primary_targets.exons.reduced.gene_type.segments.gtf",header = F)
# biotype$biotype = gsub(".*.gene_type ","",biotype$V9)
# biotype$biotype = gsub(";","",biotype$biotype)

## new annotation
new_annot = read.delim(file = gzfile("./precap_merged_annotation.gtf.gz"),header = F)
## remove spike-in froma annotation file
new_annot = new_annot[!grepl("SIRVome_isoforms",new_annot$V1),]

##  transcript id
new_annot$transcript_id = str_extract(new_annot$V9,"ENST[0-9]++\\.[0-9]++")


### old Gencode28 annotation
old_annot = read.delim("./old_gencode28_annotation.gtf",header = F,skip = 6)
## keep gene ID from old annotation for those with ENST ####
old_annot$gene_id = str_extract(old_annot$V9,"ENSG[0-9]++\\.[0-9]++")

##  transcript id
old_annot$transcript_id = str_extract(old_annot$V9,"ENST[0-9]++\\.[0-9]++")

## remove those without transcript id 
old_annot  = old_annot[!is.na(old_annot$transcript_id),]

new_annot$gene_id = old_annot$gene_id[match(new_annot$transcript_id,old_annot$transcript_id)]


id = !(is.na(new_annot$transcript_id))
new_annot$V9[id] = str_replace_all(new_annot$V9[id], pattern = "gene_id TM_[0-9]++",
                                   replacement = paste0("gene_id ", new_annot$gene_id[id]))


# remove gene id and transcript id columns to export file
new_annot$gene_id = NULL
new_annot$transcript_id = NULL

## quote values of transcript id and gene id

new_annot$V9 = gsub("(\\S+);", "\"\\1\";", new_annot$V9)

# write.table(new_annot,file = "./new_enhanced_annotation.gtf",col.names = FALSE, row.names = FALSE, quote = FALSE,sep = "\t")
write.table(new_annot,file = "./precap_enhanced_annotation.gtf",col.names = FALSE, row.names = FALSE, quote = FALSE,sep = "\t")
