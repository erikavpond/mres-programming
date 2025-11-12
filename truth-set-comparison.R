## add to comparative-analysis.R script if needed

truthset <- readVcf("truth_set.vcf", "hg38")
truth_pos <- paste(seqnames(rowRanges(truthset)), ":", start(rowRanges(truthset)))

gatk_tp <- length(intersect(gatk_pos, truth_pos))
gatk_fp <- length(setdiff(gatk_pos, truth_pos))
gatk_fn <- length(setdiff(truth_pos, gatk_pos))

freebayes_tp <- length(intersect(freebayes_pos, truth_pos))
freebayes_fp <- length(setdiff(freebayes_pos, truth_pos))
freebayes_fn <- length(setdiff(truth_pos, freebayes_pos))

gatk_precision <- gatk_tp / (gatk_tp + gatk_fp)
gatk_recall <- gatk_tp / (gatk_tp + gatk_fn)
gatk_f1 <- 2 * (gatk_precison * gatk_recall) / (gatk_precision + gatk_recall)

freebayes_precision <- freebayes_tp / (freebayes_tp + freebayes_fp)
freebayes_recall <- freebayes_tp / (freebayes_tp + freebayes_fn)
freebayes_f1 <- 2 * (freebayes_precison * freebayes_recall) / (freebayes_precision + freebayes_recall)

## represent results as a table
