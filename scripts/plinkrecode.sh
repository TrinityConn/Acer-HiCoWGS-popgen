#!/bin/bash
#SBATCH --job-name=plinkrecode
#SBATCH --account=ckenkel_26
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --array=0-2
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH -o plink2.out
#SBATCH -e plink2.error

source ~/.bashrc

conda activate garlic 

cd /scratch1/tlc_975 
#list files
FILES=("pop1.vcf" "pop2.vcf" "pop3.vcf")
INPUT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

plink2 --vcf $INPUT_FILE  --recode transpose --chr OZ035972.1,OZ035966.1,OZ035967.1,OZ035968.1,OZ035969.1,OZ035971.1,OZ035973.1,OZ035975.1,OZ035976.1,OZ035977.1,OZ035978.1,OZ035979.1,OZ035974.1 --out ${INPUT_FILE}_transformed

echo "tped/tfam files created"
