mkdir playingwithGATK
cd playingwithGATK
# ensure files needed are moved into new directory

module avail
module load (for GATK, BWA, samtools, picard)
module list

srun #open the node with 16g
module list

wget ##link to vcf files for known indels and snps in chr8

bwa index hg19_ch8.fa
samtools faidx hg19_ch8.fa

java -jar CreateSequenceDictionary.jar
  REFERENCE=hg19_ch8.fa
  OUTPUT=reference.dict
@RG\tID:Proband\tSM:Proband\tPL:Illumina\tLB:lib1\tPU:unit1

bwa mem -R "<read group info>" -p hg19_ch8.fa raw_reads.fq -o aligned_reads.sam

java -jar MarkDuplicates.jar
  -T RealignerTargetCreato
  -R hg19_ch8.fa
  -I dedupe_reads.bam
  -L 8
  -known ###vcf known indels##
  -o target_intervals.list

java -jar GenomeAnalysisTK.jar
  -T IndelRealigner
  -R hg19_ch8.fa
  -I dedupe_reads.bam
  -targetIntervals target_intervals.list
  -known ###vcf known indels##
  -o realigned_reads.bam

java -jar GenomeAnalysisTK.jar
  -T BaseRecalibrator
  -R hg19_ch8.fa
  -I realigned_reads.bam
  -L 8 
  -knownSites ###vcf known snvs##
  -knownSites ###vcf known indels##
  -o recal_data.grp
  -plots before_recal.pdf

java -jar GenomeAnalysisTK.jar
  -T BaseRecalibrator
  -R hg19_ch8.fa
  -I realigned_reads.bam
  -L 20
  -BQSR recal_data.grp
  -o post_recal_data.grp
  -plots after_recal.pdf

java -jar GenomeAnalysisTK.jar
  -T PrintReads
  -R hg19_ch8.fa
  -I realigned_reads.bam
  -L 8
  -BQSR recal_data.grp
  -o recal_reads.bam

java -jar GenomeAnalysisTK.jar
  -T ReduceReads
  -R hg19_ch8.fa
  -I recal_reads.bam
  -L 8
  -o reduced_reads.bam

java -jar GenomeAnalysisTK.jar
  -T HaplotypeCaller
  -R hg19_ch8.fa
  -I reduced_reads.bam
  -L 8
  --genotyping_mode DISCOVERY
  --output_mode EMIT_VARIANTS_ONLY
  --stand_emit_conf 10
  --stand_call_conf 30
  -o raw_variants.vcf

  ##next step is filtering variants
