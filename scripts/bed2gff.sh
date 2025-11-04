#!/bin/bash
#SBATCH --job-name=cm
#SBATCH --account=ckenkel_26
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100gb
#SBATCH --time=40:00:00
#SBATCH -o cm.out
#SBATCH -e cm.error

source ~/.bashrc
conda activate bed2gff
cd /scratch1/tlc_975


#create a gff3 file with TE annotations from existing files 

#TE annotations 

export BED="/project/ckenkel_26/RefSeqs/AcerGenome/Sanger_Inst_GCA_964034985.1/Acervicornis.repeatmask.bed"
export BED2="/scratch1/tlc_975/acerrepeat12.bed"

gt bed_to_gff3 $BED2
