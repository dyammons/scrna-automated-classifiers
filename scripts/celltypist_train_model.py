#!/usr/bin/env python3

### get the software installed - run the following in ../software
# singularity pull celltypist-latest.sif docker://quay.io/teichlab/celltypist:latest
# dd if=/dev/zero of=overlay.img bs=1M count=500 && mkfs.ext3 overlay.img
# Execute code with the following (bind mount path will vary): 
# singularity exec -B /pl/active/CSUClinHeme/users/dylan/proj04_label_transfer:/mnt --fakeroot --overlay ../software/overlay.img ../software/celltypist-latest.sif python3 celltypist_train_model.py

## Load modules
import scanpy as sc
import celltypist
## Train the OS tumor model
new_model = celltypist.train("/mnt/external_data/canine_naive_n6_annotated.h5ad", labels = "celltype.l3")
new_model.write('/mnt/output/celltypist/OS/k9_OS_celltypeL3.pkl')
## Train using healthy only annotated PBMC dataset
new_model = celltypist.train("/mnt/external_data/GSE225599_final_dataSet_H.h5ad", labels = "celltype.l3")
new_model.write('/mnt/output/celltypist/PBMC/k9_pbmc_celltypeL3.pkl')