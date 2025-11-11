library(vcfR)
vcf <- read.vcfR("recessive_hits.vcf")

# extract variants of interest amd filter
vcf_dataframe <- as.data.frame(vcf@fix)
vcf_dataframe$QUAL <- as.numeric(vcf_dataframe$QUAL)
filtered_var <- subset(vcf_dataframe, QUAL>30)

# compare number of variants before and after filtering
summary(vcf) 
summary(filtered_var)
