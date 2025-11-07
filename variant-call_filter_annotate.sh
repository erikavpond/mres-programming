## Using original trio of sorted BAM files.


## CALLING
freebayes -f hg19_ch8.fa Father_bwa.sorted.bam Mother_bwa.sorted.bam Proband_bwa.sorted.bam > Trio_raw.vcf 
# -f denotes file, so the reference file which is Chr8 on hg19, here.



## FILTERING
bcftools filter -i 'QUAL>30 && MIN(FMT/DP)>10' Trio_raw.vcf -o Trio_filtered.vcf

grep -v '^#' trio.freebayes.vcf | awk '($10 ~ /1\/1/ && $11 ~ /0\/1/ && $12 ~ /0\/1/)' > recessive_hits.vcf
# specific to this case of osteopetrosis autosomal recessive



## ANNOTATION
vep -i Trio_filtered.vcf -o Trio_annotated.txt --cache --symbol
--canonical --vcf
