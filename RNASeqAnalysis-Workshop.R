# Please write following lines to download necessary packages!
library("DESeq2")
library("stringr") #for splitting strings
library("ggplot2") #for nice plots
library("pheatmap") 
library("EnhancedVolcano") #for volcano plots
library("biomaRt") #for fetching gene annotations from Ensembl
library("gage") #for gene set enrichment analysis

# Set the correct working directory. Here this is for the HPC, but change as required for your system
setwd("~/RNA_Seq")

countsfile <- "counts.txt"
counts_table <- read.table(countsfile, header=TRUE, row.names=1)
head(counts_table) #check first few entries to see if command has worked as intended
dim(counts_table)

samples <- data.frame(row.names= colnames(counts_table))
#below splits column names of counts_table into two by severing at the _ within each column name
samples[ ,c("patient", "treatment")] <- str_split_fixed(colnames(counts_table), "_", 2) #str_split_fixed is from the stringr package

lowcov_threshold <- 3
lowcov <- (rowMaxs(as.matrix(counts_table)) <= lowcov_threshold)
table(lowcov)

filtered_counts <- counts_table[!lowcov, ]
dim(filtered_counts) #result should be the 'FALSE' bit of table lowcov, meaning removal of low coverage.

dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = samples,
                              design = ~ treatment)
dds <- DESeq(dds)

res<- results(dds, contrast = c("treatment", "treated", "control"))
head(res)
table(res$padj<0.1)
table(res$padj<0.05)
summary(res)

res_ordered = res[order(res$padj), ]

counts_normalised = counts(dds, normalized=TRUE)
counts_normalised[rownames(res_ordered)[1:5], ]

rld <- rlogTransformation(dds, blind = TRUE)



