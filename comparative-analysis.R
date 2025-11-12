library(VariantAnnotation)
library(vcfR)
library(ggplot2)

gatk_vcf <- readVcf("gatk_output.vcf", "hg38")
freebayes_vcf <- readVcf("freebayes_output.vcf", "hg38")

gatk_pos <- paste(seqnames(rowRanges(gatk_vcf)), start(rowRanges(gatk_vcf)))
freebayes_pos <- paste(seqnames(rowRanges(freebayes_vcf)), start(rowRanges(freebayes_vcf)))
overlap <- intersect(gatk_pos, freebayes_pos)
unique_calls <- union(gatk_pos, freebayes_pos)
concordance_rate <- length(overlap) / length(unique_calls)
concordance_rate

unique_gatk <- length(setdiff(gatk_pos, freebayes_pos))
unique_freebayes <- length(setdiff(freebayes_pos, gatk_pos))
bar_data <- data.frame(
  Category= c("GATK only", "FreeBayes only", "Shared"),
  Count= c(unique_gatk, unqiue_freebayes, overlap)
  )
ggplot(bar_data, aes(x= Category, y= Count, fill= Category)) +
geom_col() +
theme_minimal() +
labs(title= "Overlap of variants between GATK and FreeBayes", y= "Number of variants", x= "")

gatk_qual <- info(gatk_vcf)$QD
freebayes_qual <- info(freebayes_vcf)$QUAL
quals_df <- data.frame(
  Quality= c(gatk_qual, freebayes_qual),
  Tool= rep(c("GATK", "FreeBayes"), c(length(gatk_qual), length(freebayes_qual)))
  )
ggplot(quals_df, aes(x= Tool, y= Quality, fill= Tool)) +
geom_boxplot() +
theme_minimal() +
labs(title= "Variant quality distribution in GATK and FreeBayes", y= "Quality score", x= "")

gatk_types <- as.character(info(gatk_vcf)$TYPE)
freebayes_types <- as.character(info(freebayes_vcf)$TYPE)
types_df <- data.frame(
  Caller= rep(c("GATK", "FreeBayes"), c(length(gatk_types), length(freebayes_types))),
  Type= c(gatk_types, freebayes_types)
  )
ggplot(types_df, aes(x= Caller, fill= Type))+
geom_bar(position= "fill")+
theme_minimal()+
labs(title= "Proportion of SNVs vs indels in GATK and FreeBayes", y= "Proportion", x= "")

