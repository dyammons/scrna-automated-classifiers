#!/usr/bin/env bash

### PBMC analysis
# Run classifier on "new" dataset
singularity run \
  -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer:/data \
  -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer/output/celltypist/PBMC:/opt/celltypist/data/models \
  --overlay ../software/overlay.img \
  --fakeroot \
  ../software/celltypist-latest.sif \
  celltypist --indata /data/internal_data/canine_PBMC_query.h5ad \
  --model k9_pbmc_celltypeL3.pkl --outdir /data/output/celltypist/PBMC
  
# Run classifier using pbmc data on OS query
singularity run \
  -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer:/data \
  -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer/output/celltypist/PBMC:/opt/celltypist/data/models \
  --overlay ../software/overlay.img \
  --fakeroot \
  ../software/celltypist-latest.sif \
  celltypist --indata /data/internal_data/canine_OS_query.h5ad \
  --model k9_pbmc_celltypeL3.pkl --outdir /data/output/celltypist/OS_PBMC
  
# Run classifier using pbmc data on OS query using prob_match
singularity run \
  -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer:/data \
  -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer/output/celltypist/PBMC:/opt/celltypist/data/models \
  --overlay ../software/overlay.img \
  --fakeroot \
  ../software/celltypist-latest.sif \
  celltypist --indata /data/internal_data/canine_OS_query.h5ad \
  --model k9_pbmc_celltypeL3.pkl --outdir /data/output/celltypist/OS_PBMC_prob_match --mode prob_match


### OS analysis
# Run classifier on the dataset used to train the classifier
singularity run \
  -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer:/data \
  -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer/output:/opt/celltypist/data/models \
  --overlay ../software/overlay.img \
  --fakeroot \
  ../software/celltypist-latest.sif \
  celltypist --indata /data/external_data/canine_naive_n6_annotated.h5ad \
  --model k9_OS_celltypeL3.pkl --outdir /data/output/celltypist/OS
  
# Run classifier on the half "new" dataset in which ~50% of the dataset was used
# to train the classifier and the other ~50% is new
singularity run \
  -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer:/data \
  -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer/output:/opt/celltypist/data/models \
  --overlay ../software/overlay.img \
  --fakeroot \
  ../software/celltypist-latest.sif \
  celltypist --indata /data/internal_data/canine_OS_query.h5ad \
  --model k9_OS_celltypeL3.pkl --outdir /data/output/celltypist/OS
  