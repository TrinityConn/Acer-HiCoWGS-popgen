#!/bin/bash
#SBATCH --job-name=plinkrecode
#SBATCH --account=ckenkel_26
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH -o plink2.out
#SBATCH -e plink2.error

source ~/.bashrc

conda activate sambcf

cd /scratch1/tlc_975 

bcftools query -l snps.5bpind.DPqual.AC3.nomiss.bi.vcf.gz > samplelist.txt

bcftools view --samples-file pop1.txt snps.5bpind.DPqual.AC3.nomiss.bi.vcf.gz > pop1.vcf

bcftools view --samples-file pop2.txt snps.5bpind.DPqual.AC3.nomiss.bi.vcf.gz > pop2.vcf 
bcftools view --samples-file pop3.txt snps.5bpind.DPqual.AC3.nomiss.bi.vcf.gz > pop3.vcf 
