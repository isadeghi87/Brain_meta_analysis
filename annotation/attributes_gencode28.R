## preparing attributes of gencode28 annotation
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/annotation/")

annot = read.delim("./old_gencode28_annotation.gtf",header = F,skip = 6)
## gene id ####
annot$gene_id = str_extract(annot$V9,"ENSG[0-9]++\\.[0-9]++")
annot$gene_id = str_split(annot$gene_id,"\\.",simplify = T)[,1]

##  transcript id
annot$transcript_id = str_extract(annot$V9,"ENST[0-9]++\\.[0-9]++")

# gene type
annot$gene_type = str_split(annot$V9,"gene_type ",simplify = T)[,2]
annot$gene_type = str_split(annot$gene_type,";",simplify = T)[,1]

# gene name
annot$gene_name = str_split(annot$V9,"gene_name ",simplify = T)[,2]
annot$gene_name = str_split(annot$gene_name,";",simplify = T)[,1]

saveRDS(annot,"./gencode28.rds")
