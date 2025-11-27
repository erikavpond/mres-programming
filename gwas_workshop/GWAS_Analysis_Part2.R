###############################################
# LOAD REQUIRED LIBRARIES
###############################################

# data.table: fast loading of large GWAS results
library(data.table)

# qqman: standard Manhattan & QQ plots (simple)
library(qqman)

###############################################
# SET WORKING DIRECTORY AND LOAD PLINK RESULTS
###############################################

# Set the folder where your PLINK output files are located
setwd("~/Downloads/GWAS_Data/")

# Read the association results file produced by PLINK logistic regression
# fread() is faster than read.table()

df <- fread("wgas1_logistic_header.assoc.logistic", data.table = FALSE)

###############################################
# FILTER TO KEEP ONLY ADDITIVE MODEL RESULTS
###############################################

# PLINK produces several tests (ADD, DOM, REC).
# For GWAS we usually use the ADD (additive) model.
df_add <- subset(df, TEST == "ADD")

###############################################
# FIX P-VALUES (REMOVE NA AND ZERO)
###############################################

# Convert P-values to numeric
summary(df_add$P)
df_add$P <- as.numeric(df_add$P)

# Remove rows where P-value is NA
df_add_clean <- df_add[!is.na(df_add$P), ]

# Replace any P-values equal to 0 with the smallest possible number
# to avoid problems with -log10(0)
df_add_clean$P[df_add_clean$P == 0] <- .Machine$double.xmin

###############################################
# BASIC MANNHATTAN PLOT USING QQMAN
###############################################

# Basic Manhattan plot (no colours)
manhattan(df_add_clean,
          chr = "CHR",
          bp  = "BP",
          snp = "SNP",
          p   = "P",
          main = "GWAS Manhattan Plot (Logistic Regression)",
          genomewideline = -log10(5e-8),   # genome-wide significance
          suggestiveline = -log10(1e-5))   # suggestive threshold

# Manhattan plot with alternating colours
manhattan(df_add_clean,
          chr = "CHR",
          bp  = "BP",
          snp = "SNP",
          p   = "P",
          col = c("skyblue3", "orange3"),  # alternate colours per chromosome
          main = "GWAS Manhattan Plot (Logistic Regression)",
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-5))

# Manhattan plot with 4 custom colours
manhattan(df_add_clean,
          chr = "CHR",
          bp  = "BP",
          snp = "SNP",
          p   = "P",
          col = c("#4C72B0", "#55A868", "#C44E52", "#8172B2"),
          main = "GWAS Manhattan Plot (Coloured)",
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-5))

###############################################
# QQ PLOT TO ASSESS P-VALUE DISTRIBUTION
###############################################

# QQ plot checks for inflation or deflation of p-values
qq(df_add_clean$P,
   main = "QQ Plot of GWAS P-values")

# Calculate genomic inflation factor (lambda)
pvals <- df_add_clean$P
chisq <- qchisq(1 - pvals, df = 1)
lambda <- median(chisq) / qchisq(0.5, df = 1)
lambda

###############################################
# VIEW TOP SIGNIFICANT SNPs
###############################################

# Sort by P-value and view the most associated SNPs
top_hits <- df_add[order(df_add$P), ]
head(top_hits[, c("SNP", "CHR", "BP", "OR", "P")], 10)

###############################################
# CUSTOM GGPLOT2 MANNHATTAN PLOT
# (THIS VERSION IS NICE FOR FIGURES AND REPORTS)
###############################################

library(dplyr)
library(ggplot2)
library(ggrepel)

###############################################
# PREPARE DATA FOR GGPLOT2
###############################################

# Select relevant columns and compute -log10(p)
df <- df_add_clean %>%
  select(CHR, BP, SNP, P) %>%
  mutate(CHR = as.numeric(CHR),
         BP  = as.numeric(BP),
         P   = as.numeric(P),
         logp = -log10(P)) %>%
  filter(is.finite(logp) & logp > 0)

###############################################
# CREATE CUMULATIVE POSITIONS SO CHROMOSOMES LINE UP
###############################################

df <- df %>%
  arrange(CHR, BP) %>%         # sort by chromosome then by position
  group_by(CHR) %>%
  mutate(chr_len = max(BP)) %>%   # length of each chromosome
  ungroup() %>%
  mutate(tot = cumsum(lag(chr_len, default = 0)),   # cumulative offset
         BP_cum = BP + tot)                         # final genome position

# Create chromosome label positions (middle of each chromosome)
axis_df <- df %>%
  group_by(CHR) %>%
  summarise(center = mean(BP_cum))

###############################################
# IDENTIFY TOP 5 MOST SIGNIFICANT SNPs
###############################################

top_idx  <- order(df$P)[1:5]  # indices of lowest p-values
top_snps <- df$SNP[top_idx]   # SNP names
top_df   <- df[df$SNP %in% top_snps, ]  # subset for plotting

###############################################
# SAVE HIGH-QUALITY PDF OUTPUT
###############################################

pdf("manhattan_final.pdf", width = 30, height = 10)

###############################################
# FINAL GGPLOT2 MANNHATTAN PLOT WITH LABELLED TOP SNPS
###############################################

ggplot(df, aes(x = BP_cum, y = logp)) +

  # Plot the GWAS points
  geom_point(aes(color = as.factor(CHR)), alpha = 0.75, size = 1) +
  scale_color_manual(values = rep(c("#4C72B0", "#55A868"),
                                  length(unique(df$CHR)))) +

  # Add genome-wide significance lines
  geom_hline(yintercept = -log10(5e-8), colour = "red",
             linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = -log10(1e-5), colour = "grey40",
             linetype = "dotted", size = 0.5) +

  # Add highlighted top SNPs
  geom_point(data = top_df, aes(x = BP_cum, y = logp),
             color = "red", size = 2.3) +
  geom_text_repel(data = top_df,
                  aes(label = SNP),
                  color = "red", size = 3.2,
                  nudge_y = 0.4,
                  min.segment.length = 0,
                  segment.size = 0.3,
                  segment.color = "red",
                  box.padding = 0.3,
                  point.padding = 0.15) +

  # Chromosome numbers on x-axis
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +

  # Axis titles and overall title
  labs(title = "GWAS Manhattan Plot – Top SNPs Labeled (Final)",
       x = "Chromosome",
       y = expression(-log[10](p))) +

  # Improve plot appearance
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11, margin = margin(r = 12)),
    axis.title.y = element_text(size = 13, margin = margin(r = 18)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(20, 20, 20, 90)   # wide left margin
  )

dev.off()

##############################
##Combined Manhattan + QQ Plot
##############################
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

df <- df_add_clean %>%
  select(CHR, BP, SNP, P) %>%
  mutate(CHR = as.numeric(CHR),
         BP  = as.numeric(BP),
         P   = as.numeric(P),
         logp = -log10(P)) %>%
  filter(is.finite(logp) & logp > 0)

# Order SNPs by chromosome and genomic location
df <- df %>%
  arrange(CHR, BP) %>%
  group_by(CHR) %>%
  mutate(chr_len = max(BP)) %>%
  ungroup() %>%
  mutate(tot = cumsum(lag(chr_len, default = 0)),
         BP_cum = BP + tot)

# Chromosome label positions
axis_df <- df %>%
  group_by(CHR) %>%
  summarize(center = mean(BP_cum))

# Identify top SNPs to label
top_idx  <- order(df$P)[1:5]
top_snps <- df$SNP[top_idx]
top_df   <- df[df$SNP %in% top_snps, ]


#Build Manhattan plot
manhattan_plot <- ggplot(df, aes(x = BP_cum, y = logp)) +

  geom_point(aes(color = as.factor(CHR)), alpha = 0.75, size = 1) +
  scale_color_manual(values = rep(c("#4C72B0", "#55A868"), length(unique(df$CHR)))) +

  geom_hline(yintercept = -log10(5e-8), colour = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(1e-5), colour = "grey40", linetype = "dotted") +

  geom_point(data = top_df, color = "red", size = 2.5) +
  geom_text_repel(data = top_df,
                  aes(label = SNP),
                  color = "red", size = 3,
                  nudge_y = 0.3,
                  max.overlaps = Inf) +

  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  labs(title = "Manhattan Plot", x = "Chromosome", y = expression(-log[10](p))) +

  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(15,15,5,15)
  )


#Build QQ plot
# Expected and observed p-values
exp <- -log10(ppoints(length(df$P)))
obs <- -log10(sort(df$P))

lambda <- median(qchisq(1 - df$P, 1)) / qchisq(0.5, 1)

qq_plot <- ggplot(data.frame(exp, obs), aes(x = exp, y = obs)) +
  geom_point(size = 1, alpha = 0.6, color = "#4C72B0") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = paste0("QQ Plot (λ = ", round(lambda, 3), ")"),
       x = "Expected -log10(p)",
       y = "Observed -log10(p)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.margin = margin(5,15,15,15)
  )

combined_plot <- manhattan_plot / qq_plot +
  plot_annotation(title = "GWAS Manhattan + QQ Plot",
                  theme = theme(plot.title = element_text(size = 18, hjust = 0.5)))


ggsave("GWAS_Manhattan_QQ_Combined.pdf", combined_plot, width = 25, height = 10)
ggsave("GWAS_Manhattan_QQ_Combined.png", combined_plot, width = 14, height = 10, dpi = 300)


##########################################
# END OF THE SCRIPT
#ALL THE BEST
#DR MO ELASRAG
##########################################
