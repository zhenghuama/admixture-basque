#!/bin/bash

########### PS2 - Global ancestry - ADMIXTURE ###############
# This script contains template code for running ADMIXTURE
# for CSE284 PS2
# Usage: ./ps2_run_admixture.sh

# Paths to data files you will need
VCF=./data_and_results/step1/dataset/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
SAMPLEINFO=./data_and_results/step1/dataset/igsr_samples.tsv

# Output prefix to use
PREFIX=step1_admixture
DATA_DIR="./data_and_results/step1"
INTERMEDIATE_DIR="${DATA_DIR}/intermediate"
mkdir -p "$INTERMEDIATE_DIR"


# Keep all preprocessing artifacts under step1/intermediate/
BASE_PREFIX="${INTERMEDIATE_DIR}/${PREFIX}"
PRUNED_PREFIX="${INTERMEDIATE_DIR}/${PREFIX}.pruned"
#############################################################

# Step 1: Extract samples from the CEU, PEL, GWD, ASW, and PUR populations.
# Recall these are in igsr_samples.tsv. 
# Not all samples in that file are included in our VCF.
# You can restrict to samples listed as "Genomes phase 3 release".
# Note plink will need a file with 2 columns (fam
# ily ID and sample ID).
# We do not have a family ID so you can use the sample ID as the family ID
#echo "Not implemented - extract samples"
KEEPFILE="${INTERMEDIATE_DIR}/${PREFIX}.keep"
awk -F'\t' '
  BEGIN{OFS="\t"}
  NR==1 {next}  # skip header
  ($4=="CEU" || $4=="PEL" || $4=="GWD" || $4=="ASW" || $4=="PUR") && ($9 ~ /phase 3 release/) {
    print $1, $1
  }
' "$SAMPLEINFO" > "$KEEPFILE"

# Step 2: Use plink to preprocess the dataset.
# (1) filter variants with MAF<1% (--maf 0.01)
# (2) keep only samples extracted in step 1 (--keep <sample list> --double-id)
# (3) Output a binary plink format (--make-bed)
# Use --out $PREFIX, which should output files ps2_admixture.bed, etc.
#echo "Not implemented - preprocess with plink"
plink --vcf "$VCF" \
  --keep "$KEEPFILE" \
  --double-id \
  --maf 0.01 \
  --geno 0.05 \
  --mind 0.05 \
  --make-bed \
  --out "$BASE_PREFIX"


# Step 3: Use plink to perform LD-pruning.
# You can use the plink --indep-pairwise option to get a list of 
# variants to include (plink.prune.in)
# Then use plink --extract plink.prune.in --make-bed --out $PREFIX.pruned
# to output a new dataset (ps2_admixture.pruned.bed, etc.) with only the
# pruned variants
# See example in https://dalexander.github.io/admixture/admixture-manual.pdf
#echo "Not implemented - LD pruning"
plink --bfile "$BASE_PREFIX" \
  --indep-pairwise 50 10 0.1 \
  --out "$BASE_PREFIX"

# Apply pruning list
plink --bfile "$BASE_PREFIX" \
  --extract "${BASE_PREFIX}.prune.in" \
  --make-bed \
  --out "$PRUNED_PREFIX"

# Step 4: Run admixture for values of K=3 or K=5
# Can try higher K but it takes a long time!
#for K in 3 5
#do
    #echo "Not implemented - ADMIXTURE for K=$K"
#    admixture "${PRUNED_PREFIX}.bed" "$K" > "${PRUNED_PREFIX}.${K}.0"
#done
