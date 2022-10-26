# check , install and load dependencies
if( ! "tidyverse" %in% installed.packages()){
    install.packages("tidyverse")
}

if( ! "biomaRt" %in% installed.packages()){
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install("biomaRt")
}

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(biomaRt))

# Create and set working directory
dir.create("~/interaction/", showWarnings = FALSE)
setwd("~/interaction/")

protein <- "./9606.protein.links.full.v11.5.txt"
if( ! file.exists(protein)){
    download.file(url = "https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz",
                  destfile = protein)
    
}
save_protein <- "./copexpressed_proteins.RDS"


# Protein correlation
protein <- read.csv(file = protein,
                    header = TRUE,
                    sep = "")

protein$protein1 <- str_remove(protein$protein1, "9606.")
protein$protein2 <- str_remove(protein$protein2, "9606.")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_id <- getBM(filters = "ensembl_peptide_id", 
                attributes = c("ensembl_peptide_id","hgnc_symbol", "description"),
                values = ensp, mart = mart)

protein$protein1_gene <- "."
protein$protein2_gene <- "."

protein2 <- protein[which(protein$coexpression != 0),]

tmp <- unique(c(unique(protein2$protein1),unique(protein2$protein2)))
for (i in tmp) {
   j <- which(gene_id$ensembl_peptide_id == i)
   k <- which(protein2$protein1 == i)
   l <- which(protein2$protein2 == i)
   if (length(j) > 0) {
       if(length(k) > 0) {
           protein2$protein1_gene[k] <- gene_id$hgnc_symbol[j]   
       }
       if(length(l)>0){
           protein2$protein2_gene[l] <- gene_id$hgnc_symbol[j]   
       }
   }
}
saveRDS(object = protein2,file = save_protein)

protein <- readRDS(file = save_protein)

