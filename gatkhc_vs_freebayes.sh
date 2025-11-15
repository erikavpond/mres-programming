## This script was used to perform variant calling on the same .bam file, to give two different .vcf files for further analysis


# FreeBayes
freebayes -f hg19_ch8.fa Proband_bwa.sorted.bam -o proband_fb.vcf

# GATK HaploTypeCaller
gatk HaplotypeCaller -R hg19_ch8 -I Proband_bwa.sorted.bam -O proband_gatk.vcf
