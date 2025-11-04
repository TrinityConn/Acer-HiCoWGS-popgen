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


#set up environment 
source ~/.bashrc
conda activate quarTeTdependencies 

GENOME="/project/ckenkel_26/RefSeqs/AcerGenome/Sanger_Inst_GCA_964034985.1/GCA_964034985.1_jaAcrCerv1.1_genomic.fna"
GENE="/project/ckenkel_26/RefSeqs/AcerGenome/Sanger_Inst_GCA_964034985.1/Acervicornis.gff3"
TE="/scratch1/tlc_975/acerrepeat12.bed"

python3 /home1/tlc_975/quarTeT/quartet.py cm -i $GENOME --gene $GENE --TE $TE
