#Made for command line use

#Sit in directory with all files needed please
#gPLINK.jar, pop.cov, wgas1.map, wgas.ped, 

head wgas1.map #.map includes base no, SNP id, zygosity, 
head wgas1.bed #.ped includes letter per SNP, so for first .map entry perhaps a C.

plink --file wgas1 --make-bed --out wgas1 #converts existing text files into binary format for faster processing 
ls -lts #sees for all files in directory that have been made, time sorted

plink --bfile wgas1 --geno 0.05 --make-bed --out wgas1_step1 #QC step that removes SNPs with too many missing genotypes (here, >5%), output file_step1

plink --bfile wgas1_step1 --mind 0.05 --make-bed --out wgas1_step2 #QC step 2: removes samples with tpp many missing genotypes, (here, >5% again)

plink --bfile wgas1_step2 --maf 0.01 --make-bed --out wgas1_step3 #QC step 3: removes ultra rare variants by only keeping SNPs that have AF>=1%

plink --bfile wgas1_step3 --hwe 1e-6 --make-bed --out wgas1_qc #\\\\\qc step4: remove SNPs that break HWE, with extreme genotype patterns

plink --bfile wgas1_qc --logistic --ci 0.95 --out wgas1_logistic #Runs basic GWAS, without covariates, with CI 95% for odds ratios

head wgas1_logistic.assoc.logistic #view first vits of results file

awk '{print $1, $2, int(rand()*2)}' wgas1_qc.fam > wgas1_cov.cov #creates simple population group covariate for each sample, may not ever have to do for us

tr '\t' ' ' < wgas1_cov.cov | awk '{$1=$1; print}' > covar_fixed.cov #cleans covariate file so separators are space, not tab delim

awk '{print NR, NF, $0}' covar_fixed.cov | head #optional check

echo -e "FID IID POP" | cat - covar_fixed.cov > covar_header.cov #adds header line to the covariate file

head covar_header.cov #check that header has been added

plink --bfile wgas1_qc \
      --covar covar_header.cov \
      --covar-name POP \
      --logistic \
      --ci 0.95 \
      --out wgas1_logistic_header
#Runs GWAS again, now adjusting for population covariates, POP means use POP column as covariate, logistic is for logistical casevcontrol regression





