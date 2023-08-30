setwd("~/JP_JL/gtex/")

# Get and load gtex
gtexPath <- "./gtex_gene_mean.gct.gz"
if ( ! file.exists(gtexPath)) {
    download.file(url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
                  destfile = gtexPath)
}

gtex <- read.delim(file = gtexPath,
                   header = T,
                   skip = 2,
                   sep = "\t")


# load RNA file
RNA <- read.delim("../pubmed_count/ALL/Complexity-RNA-sum-unit-0.7-PN.dat")


# Mutual filter
RNA <- RNA[which(RNA$RNA %in% gtex$Description),]

gtex <- gtex[which(gtex$Description %in% RNA$RNA),]
row.names(gtex) <- gtex$Description
gtex <- as.matrix(gtex[,3:ncol(gtex)])

breastEndogenousTissus <- c("Adipose...Subcutaneous",
                            "Adipose...Visceral..Omentum.",
                            "Breast...Mammary.Tissue")

rnaGtex  <- data.frame(RNA = RNA$RNA,
                       gtex = sapply(RNA$RNA,
                                     function(x){
                                         if("Breast...Mammary.Tissue" %in% colnames(gtex)[order(gtex[x,],decreasing = T)][1:5]){
                                             return("breast_high_expression")
                                         }else{
                                             return("breast_low_expression")
                                         }
                                     }))

nsample = 100
top_vs_bottom <- data.frame(RNA = c(head(rnaGtex$RNA,n = nsample),
                                    tail(rnaGtex$RNA,n = nsample)),
                            gtex = c(head(rnaGtex$gtex,n = nsample),
                                     tail(rnaGtex$gtex,n = nsample)),
                            type = c(rep("high_complexity",nsample),
                                     rep("low_complexity",nsample)))
chisq.test(table(top_vs_bottom$gtex,top_vs_bottom$type))
