

####################################
########    Libraries    ###########
####################################
## Downloading relevant tools to analyse RNA Seq data
library("DESeq2")
library("stringr")
library("ggplot2") #for nice plots
library("pheatmap") 
library("EnhancedVolcano") #for volcano plots
library("biomaRt") #for fetching gene annotations from Ensembl
library("gage") #for gene set enrichment analysis
library("RColorBrewer") #can make pretty colour

# Ensure in correct working directory, otherwise setwd(directorypathgoeshere)




####################################
########    Data Prep     ##########
####################################
countsfile <- "counts.txt" #saves this string into var name
counts_table <- read.table(countsfile, header=TRUE, row.names=1) #uses string defined above as the filename to be searched
head(counts_table)
dim(counts_table) #no. genes, no. patients(control + treated)

samples <- data.frame(row.names = colnames(counts_table)) 
samples[,c("patient", "treatment")] <- str_split_fixed(colnames(counts_table), "_", 2) #splits column names of the counts table into two at the '_'. from stringr package

lowcoverage_threshold <- 3
lowcoverage <- (rowMaxs(as.matrix(counts_table)) <= lowcoverage_threshold) #saves all entries with coverage lower than threshold into 'lowcoverage'
table(lowcoverage) #True/False corresponding result

filtered_countstable <- counts_table[!lowcoverage,] # bye bye low coverage data for the new set
dim(filtered_countstable)

dds <- DESeqDataSetFromMatrix(countData = filtered_countstable, #function for data loaded in as a matrix, other functions for HTSeq or txiimport
                              colData = samples,
                              design = ~ treatment) # Determines through which lens we view this data
dds <- DESeq(dds) # Runs DESeq
dds

res <- results(dds, contrast = c("treatment", "treated", "control"))
table(res$padj<0.05) # p-adjusted value less than 0.05
summary(res)

# Order genes by significance of differential expression
res_ordered = res[order(res$padj), ]
head(res_ordered, 5)

# Checks normalised counts for top 5 significant genes
counts_normalised = counts(dds, normalized=TRUE)
counts_normalised[rownames(res_ordered)[1:5], ]



#######################################
########  Log transformation  #########
#######################################
rld <- rlogTransformation(dds, blind = TRUE) # Regularised log trasnformation by DESeq2







####################################
########     PCA Plot     ##########
####################################
pcaData <- plotPCA(rld, intgroup = c("patient", "treatment"), ntop = 1000, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

(pca=ggplot(pcaData, aes(PC1, PC2, colour= patient, shape= treatment)) +
  geom_point(size=3) + #this is the part that generally decides type of plotting: bar/point/col
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() #combo of bar and coord_polar makes a pie :)
)
ggsave(pca, filename= "PCA_1000genesonly.png", width = 6, height = 4) #saves plot as png, in working directory







####################################
########      MA Plot     ##########
####################################
# MA plot visualises gene expression changes throughout, where x is mean counts (normalised), y is the LFC 
# Helps determine genes with large LFC
plotMA(res, # When using dds, because is metadata, the colours of the plot revert to default. so have used res here
       main = "MA Plot",
       xlab = "Mean of Normalised Counts",
       ylim = c(-8,8),
       colNonSig="purple",
       colSig="red", # Significant data in red
       colLine="gray")





dds_paired <- DESeqDataSetFromMatrix(countData = filtered_countstable,
                              colData = samples,
                              design = ~ patient+treatment)
dds_paired <- DESeq(dds_paired) # Runs DESeq2

res_paired <- results(dds_paired, contrast = c("treatment", "treated", "control"))
table(res_paired$padj<0.05) # p-adjusted value less than 0.05
summary(res_paired)





####################################
########      Gene ID     ##########
####################################
# Converting Ensembl number to a gene name, using biomaRt.
listEnsembl() # Requires internet connection
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      GRCh=37) #creates database for hg19 assembly

gene_label_conversion <- getBM(attributes = c("ensembl_gene_id",
                                              "external_gene_name",
                                              "description",
                                              "entrezgene_id"),
                               mart = ensembl) # Finds each of the specified values from the specified assembly
dim(gene_label_conversion)
head(gene_label_conversion)





####################################
########      Cleaning    ##########
####################################
res_paired_annotated <- merge(gene_label_conversion[!(duplicated(gene_label_conversion$ensembl_gene_id)), c(1,2,4)],
                              as.data.frame(res_paired),
                              by.x=1, by.y=0, all.y=T) # Merge two datasets, not including duplicated ids, get rid of column 3. Leave second set as is. Deal with inconsistencies in rows/columns

res_paired_annotated_ordered <- res_paired_annotated[order(res_paired_annotated$padj), ] #orders by p adujsted value

head(res_paired_annotated)
head(res_paired_annotated_ordered) # Heading both allows to visualise the change made by ordering. No longer 1,2,3,.., but by asc p-adj

output=res_paired_annotated_ordered
# Cleaning with sapply()
output$pvalue <- sapply(output$pvalue, format, digits=3)
output$padj <- sapply(output$padj, format, digits=3)
output$baseMean <- sapply(output$baseMean, round, digits=1) # Round to 1 decimal place
output$log2FoldChange <- sapply(output$log2FoldChange, round, digits=3)
output$stat <- sapply(output$stat, round, digits=3)
output$lfcSE <- NULL

head(output) #Compared to prev head. the number values are shortened to the decided decimal places e.g. 3.204553 becomes 3.205

write.table(output, file = "Kiki_DESeqFIRSTEVER.txt", sep = "\t", row.names = F, quote = F) # Save output as table in txt file, tab delimited




####################################
########      Volcano!    ##########
####################################
pdf(file = "RNASeq_volcano_plot.pdf", width= 12, height = 12)
(volc=EnhancedVolcano(res_paired_annotated,
                      lab= res_paired_annotated$external_gene_name,
                      x= 'log2FoldChange',
                      y= 'padj',
                      pCutoff= 0.01,
                      ylab= bquote(~-Log[10]~italic(P[adj])),
                      title= "Volcano Plot",
                      subtitle= "Treated vs Control",
                      titleLabSize =18,
                      subtitleLabSize= 14,
                      axisLabSize= 12,
                      captionLabSize= 9,
                      col = c("pink", "purple", "orange", "lightblue")
                      ))
dev.off() #important to be able to export





####################################
########      Heat Map    ##########
####################################

topgenes=res_paired_annotated_ordered[1:50, 1:2]
matr= assay(rld)[topgenes$ensembl_gene_id, ]
rownames(matr)=topgenes$external_gene_name
head(matr)

pdf(file = "RNASeq_heatmap.pdf", width= 12, height = 12)
(hm=pheatmap(matr,
             scale="row",
             annotation_col=samples[ ,"treatment",drop=F],
             fontsize=10,
             fontsize_row = 8,
             fontsize_col = 8,
             treeheight_row = 20,
             treeheight_col = 20))
dev.off() #important to be able to export

