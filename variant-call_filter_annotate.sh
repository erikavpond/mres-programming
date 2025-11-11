## Using original trio of sorted BAM files. Ensure all files to start inc ref genome are in same directory
## As done in class


## CALLING using freebayes
freebayes -f hg19_ch8.fa Father_bwa.sorted.bam Mother_bwa.sorted.bam Proband_bwa.sorted.bam > Trio_raw.vcf 
# -f denotes file, so the reference file which is Chr8 on hg19, here.


## FILTERING using bcftools
bcftools filter -i 'QUAL>30 && MIN(FMT/DP)>10' Trio_raw.vcf -o Trio_filtered.vcf
grep -v '^#' trio.freebayes.vcf | awk '($10 ~ /1\/1/ && $11 ~ /0\/1/ && $12 ~ /0\/1/)' > recessive_hits.vcf
# specific to this case of osteopetrosis autosomal recessive. Find all in file (besides header) that meets following conditions: homozygous in column 10 AND heterozygous in row11, ...


## ANNOTATION using vep
vep -i Trio_filtered.vcf > Trio_annotated.txt --cache --symbol
--canonical --vcf
#
