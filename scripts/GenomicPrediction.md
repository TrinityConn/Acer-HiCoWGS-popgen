---
layout: post
title: Genomic Prediction 
date: '2025 Feb 10'
categories: Analysis, Genome Assembly,
tags: [Bioinformatics, Genome, Assembly]
projects: Acropora cervicornis haplotype panel
---


# Genomic Prediction of *Acropora cervicornis* 

## In Collaboration with USC 


Genomic Prediction of ED50's of *Acropora cervicornis*

## Issue 
Complex traits like thermal resilience, particularly as measured by a proxy measurement like Fv/Fm, are unlikely to be one or two loci of large effect, particularly at the sample sizes being measured. Therefore, the next step is to build a genomic prediction model and measure a polygenic score to characterize how variation could be predictive of phenotype. 

## Basic Approach 
The idea is to use the GWAS summary file – the list of p-values of SNPs associated with the trait, to then predict the scores and trait of un-phenotyped individuals. Because we are building polygenic scores from the only individuals we are analyzing, 

## Files Needed 

### GWAS Summary File -- for summary-based analysis 

```
#DRTO Only 
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/DRTO_GWAS_rel/cbass_gwas_KING_ed50_norm.BLUP_just.DRTO_NO.WQ_retry.assoc.txt

#FLK Only
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/cbass_gwas_KING_raw.ed50_norm.BLUP.NO.WQ_nursery_NOT.DRTO.assoc.txt

#Merged
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/merged_GWAS_rel/cbass_gwas_franken.merge_post.BLUP_215.genets_after.quant.norm_W.POP.PCs.assoc.txt


```
### Bed/Bim/Fam File  -- for individual - level analysis 

```

#DRTO Only 
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/DRTO_GWAS_rel/gwas.LD97_just.DRTO_LoCo.imp_Dec.2025_111.genets

#FLK Only
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/FLKE_GWAS_rel/gwas.LD97_NOT.DRTO_LoCo.imp_Dec.2025_105.genets

#Merged
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/merged_GWAS_rel/gwas.LD97_all.sites_LoCo.imp_Dec.2025_215.genets_force.merge 

```
### Phenotype files 

```
#DRTO Only
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/DRTO_GWAS_rel/phenotype_main_just.DRTO_NO.WQ.txt
#FLK Only
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/FLKE_GWAS_rel/phenotype_main_NOT.DRTO_nursery_NO.WQ_raw.ed50.txt

#Merged Only 
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/merged_GWAS_rel/phenotype_main_all_force.merge.txt


```
### Covariates 

```
#DRTO only
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/DRTO_GWAS_rel/covars.pop.2PCs_just.DRTO_NO.WQ_retry.txt

#FLK Only
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/FLKE_GWAS_rel/covars.pop.2PCs_NOT.DRTO_nursery_NO.WQ_raw.ed50.txt 

#Merged 
/project2/ckenkel_26/Acer_WGS/GWAS.rel_for.TC/merged_GWAS_rel/covars.pop.2PCs_NOT.DRTO_nursery_NO.WQ_raw.ed50.txt
```


1. Calculate x different training sets so that every sample is in the training data set multiple times but in the test set
2. Build Genomic Prediction Model
3. Measure Model Accuracy relative to measured values 


### We will want to 
### Starting with BayesR 

```
#install bayesR
git clone https://github.com/syntheke/bayesR.git

#compile
gfortran -o bayesR -O2 -cpp RandomDistributions.f90 baymods.f90 bayesR.f90
gfortran -o bayesRv2 -O2 -cpp -Dblock -fopenmp RandomDistributions.f90 baymods.f90 bayesR.f90
ifort -o bayesR -O3 -fpp RandomDistributions.f90 baymods.f90 bayesR.f90
ifort -o bayesRv2 -O3 -fpp -Dblock -openmp -static RandomDistributions.f90 baymods.f90 bayesR.f90
```

```
#running bayesr

#first we need to incorporate the phenotype file into the fam file
conda activate plink
#
plink --bfile gwas.LD97_just.DRTO_LoCo.imp_Dec.2025_111.genets --pheno phenotype_main_just.DRTO_NO.WQ.txt --make-bed --out gwas.LD97_just.DRTO_LoCo.imp_Dec.2025_111.genets_sub_updated

```
### Running bayesR with 5-fold CV for genomic prediction 
#!/bin/bash

# 5-Fold Cross-Validation for BayesR


# Input arguments
BFILE=gwas.LD97_just.DRTO_LoCo.imp_Dec.2025_111.genets_sub_updated   # Path to PLINK bfile (without .bed/.bim/.fam extension) this file should include phenotype as column 6
OUTPUT=gwas.LD97_just.DRTO_LoCo.imp_Dec.2025_111.genets_sub_update_bayesr # Output prefix
NFOLDS=5           # Number of folds

# BayesR parameters
BURNIN=10000
ITERATIONS=50000
THIN=10

echo "Starting 5-Fold Cross-Validation for BayesR"
echo "Input: ${BFILE}"
echo "Output: ${OUTPUT}"

# Create output directory
mkdir -p ${OUTPUT}_cv_results

# Get list of individuals
plink --bfile ${BFILE} --write-snplist --out ${OUTPUT}_temp
awk '{print $1, $2}' ${BFILE}.fam > ${OUTPUT}_individuals.txt
N_IND=$(wc -l < ${OUTPUT}_individuals.txt)

echo "Total individuals: ${N_IND}"
FOLD_SIZE=$((N_IND / NFOLDS))

# Shuffle individuals and assign to folds
shuf ${OUTPUT}_individuals.txt > ${OUTPUT}_individuals_shuffled.txt

# Split into folds
split -l ${FOLD_SIZE} -d -a 1 ${OUTPUT}_individuals_shuffled.txt ${OUTPUT}_fold_

# Run BayesR for each fold
for FOLD in $(seq 0 $((NFOLDS-1))); do
    echo ""
    echo "========================================"
    echo "Processing Fold ${FOLD}"
    echo "========================================"
    
    # Test set for this fold
    TEST_FILE="${OUTPUT}_fold_${FOLD}"
    TRAIN_FILE="${OUTPUT}_train_fold_${FOLD}.txt"
    
    # Create training set (all folds except current)
    for OTHER_FOLD in $(seq 0 $((NFOLDS-1))); do
        if [ ${OTHER_FOLD} -ne ${FOLD} ]; then
            cat ${OUTPUT}_fold_${OTHER_FOLD} >> ${TRAIN_FILE}
        fi
    done
    
    echo "Training set: $(wc -l < ${TRAIN_FILE}) individuals"
    echo "Test set: $(wc -l < ${TEST_FILE}) individuals"
    
    # Create PLINK files for training set
    plink --bfile ${BFILE} \
          --keep ${TRAIN_FILE} \
          --make-bed \
          --out ${OUTPUT}_train_fold_${FOLD}
    
    # Extract training phenotypes
    awk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' ${TRAIN_FILE} ${PHENO} > ${OUTPUT}_train_pheno_fold_${FOLD}.txt
    
    # Run BayesR on training set
    echo "Running BayesR on training set..."
    
    # BayesR command - adjust path to your BayesR executable
    bayesR \
        -bfile ${OUTPUT}_train_fold_${FOLD} \
        -out ${OUTPUT}_cv_results/fold_${FOLD} \
        -n ${BURNIN} \
        -s ${ITERATIONS} \
        -thin ${THIN} \
        -seed $((12345 + FOLD))
    
    # Generate predictions for test set
    echo "Generating predictions for test set..."
    
    # Extract SNP effects from BayesR output
    # Format: SNP_ID effect
    tail -n +2 ${OUTPUT}_cv_results/fold_${FOLD}.SnpEffects > ${OUTPUT}_cv_results/fold_${FOLD}_effects.txt
    
    # Calculate genomic predictions for test set using PLINK --score
    plink --bfile ${BFILE} \
          --keep ${TEST_FILE} \
          --score ${OUTPUT}_cv_results/fold_${FOLD}_effects.txt 1 2 3 sum \
          --out ${OUTPUT}_cv_results/fold_${FOLD}_predictions
    
    # Extract test phenotypes
    awk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' ${TEST_FILE} ${PHENO} > ${OUTPUT}_test_pheno_fold_${FOLD}.txt
    
    echo "Fold ${FOLD} complete"
done

```






