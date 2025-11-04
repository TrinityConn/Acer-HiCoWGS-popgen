#!/bin/bash
#SBATCH --job-name=garlic
#SBATCH --account=ckenkel_26
#SBATCH --partition=epyc-64
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100gb
#SBATCH --time=40:00:00
#SBATCH -o garlic.out
#SBATCH -e garlic.error

#set up environment 
source ~/.bashrc
conda activate garlic 

cd /scratch1/tlc_975 

#define variables 
TPED="pop3.vcf_transformed.tped"
TFAM="pop3.vcf_transformed.tfam"
CENT="cent2.file"


#run garlic 
/home1/tlc_975/linux/garlic --tped $TPED --tfam $TFAM --centromere $CENT --error 0.001 --auto-winsize --winsize 10 --auto-overlap-frac --tped-missing . --gl-type GQ --resample 30 --out pop3
 
#end
