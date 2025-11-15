library("VariantAnnotation")
library("vcfR")
library("ggplot2")
library("GenomicRanges")
library("dplyr")

# setwd("~/Desktop/WES-data")

gatk_vcf <- readVcf("proband_gatk.vcf", genome = "hg19") 
filtered_gatkvcf <- subset(gatk_vcf, QUAL>30) #filtering
fb_vcf <- readVcf("proband_fb.vcf", genome = "hg19")
filtered_fbvcf <- subset(fb_vcf, QUAL>30) #filtering
gr_gatk <- rowRanges(filtered_gatkvcf)
gr_fb <- rowRanges(filtered_fbvcf)


# Finding unique and shared variants 

overlap <- (intersect(names(gr_gatk), names(gr_fb))) #finding shared variants of each
head(overlap)

unique_gatk <- setdiff(names(gr_gatk), names(gr_fb)) 
uniquegr_gatk <- gr_gatk[names(gr_gatk) %in% unique_gatk] #lists names of variants unique to gatk-hc calls
unique_fb <- setdiff(names(gr_fb), names(gr_gatk)) 
uniquegr_fb <- gr_fb[names(gr_fb) %in% unique_fb] #lists names of variants unique to freebayes calls
numunique_gatk <- length(uniquegr_gatk)
numunique_fb <- length(unique_fb)

overlap_df <- data.frame(
  Category= c("GATK only", "FreeBayes only", "Overlap"),
  Count= c(numunique_gatk, numunique_fb, length(overlap))
)
(chart1=ggplot(overlap_df, aes(x="", y=Count, fill=Category)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Overlap of Variants Discovered", y= "Number of Variants", x="")
)
ggsave(chart1, filename= "overlaps.pdf", width = 6, height = 6) #saves first chart as pdf






# Assessing types of variant found

get_variant_type <- function(vcf) {
  ref_alleles <- as.character(ref(vcf))
  alt_first <- sapply(alt(vcf), function(x) {
    if (length(x) == 0) return(NA_character_)
    as.character(x[1])
  })
  ref_len <- nchar(ref_alleles)
  alt_len <- nchar(alt_first)
  
  variant_type <- ifelse(is.na(alt_first), NA_character_,
                  ifelse(ref_len == 1 & alt_len == 1, "SNV",
                  ifelse(ref_len < alt_len, "Insertion",
                  ifelse(ref_len > alt_len, "Deletion", "MNP/Other"))))
  return(variant_type)
} #creating a function that find the variant type (SNV/indel/etc) for each in list

gatk_types <- get_variant_type(filtered_gatkvcf) #apply function
fb_types <- get_variant_type(filtered_fbvcf) #apply function
gatk_dftype <- data.frame(Caller="GATK", Type= gatk_types)
fb_dftype <- data.frame(Caller="FreeBayes", Type= fb_types)
types_df <- bind_rows(gatk_dftype, fb_dftype) %>%
  filter(!is.na(Type))

type_summary <- types_df %>%
  group_by(Caller, Type) %>%
  summarise(Count= n(), .groups = "drop") %>%
  group_by(Caller) %>%
  mutate(Proportion = Count / sum(Count)) #reformatting data frame for visualisation

cols <- c("SNV" = "#56B4E9",
          "Insertion" = "#009E73",
          "Deletion" = "#E69F00",
          "MNP/Other" = "#D55E00") #choosing colours!

(pie_gatk=ggplot(filter(type_summary, Caller== "GATK"),
                 aes(x="", y= Proportion, fill= Type)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols) +
    labs(title = "GATK variant types") +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
)
ggsave(pie_gatk, filename= "piechart_gatk_variants.pdf") #saves pie chart for gatk variants

(pie_fb=ggplot(filter(type_summary, Caller == "FreeBayes"),
                 aes(x="", y= Proportion, fill= Type)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols) +
    labs(title = "FreeBayes variant types") +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
)
ggsave(pie_fb, filename= "piechart_fb_variants.pdf") #saves pie chart for freebayes variants


type_summary %>%
  mutate(Percent = round (Proportion *100, 2)) %>%
  select(Caller, Type, Count, Percent)
type_summary
         





## attempt Ti/Tv ratio

get_titv <- function(vcf) {
  ref_alleles <- as.character(ref(vcf))
  alt_alleles <- sapply(alt(vcf), function(x) {
    if (length(x) == 0) return(NA_character_)
    as.character(x[1])
  }) 
  snv_mask <- nchar(ref_alleles) == 1 &nchar(alt_alleles) == 1
  ref_alleles <- ref_alleles[snv_mask]
  alt_alleles <- alt_alleles[snv_mask]
  transitions <- ((ref_alleles == "A" & alt_alleles == "G") |
                    (ref_alleles == "G" & alt_alleles == "A") |
                    (ref_alleles == "C" & alt_alleles == "T") |
                    (ref_alleles == "T" & alt_alleles == "C"))
  transversions <- !(transitions)
  ti_count <- sum(transitions, na.rm = TRUE)
  tv_count <- sum(transversions, na.rm = TRUE)
  ratio <- ifelse(tv_count == 0, NA, ti_count / tv_count)
  return(list(
    total_snvs = length(ref_alleles),
    transitions = ti_count,
    transversions = tv_count,
    TiTv_ratio = ratio
  ))
} #create function to find transition transversion ratio

gatk_titv <- get_titv(filtered_gatkvcf)
fb_titv <- get_titv(filtered_fbvcf)
titv_summary <- bind_rows(
  data.frame(Caller="GATK", gatk_titv),
  data.frame(Caller="FreeBayes", fb_titv))
titv_summary

(titv_plot=ggplot(titv_summary, aes(x=Caller, y=TiTv_ratio, fill=Caller))+
    geom_bar(stat="identity", width = 0.5) +
    geom_text(aes(label= round(TiTv_ratio, 2)), vjust = -0.5) +
    ylim(0, max(titv_summary$TiTv_ratio) * 1.2) +
    labs(title = "Transition/Transversion ratio per variant caller",
         y="Ti/Tv ratio", x=NULL) +
    theme_minimal()+
    theme(legend.position = "none")
)
ggsave(titv_plot, filename="titv_summary.pdf") #saves TiTv ratio summary as pdf
