#!/bin/bash
#
#Job Script
#
#SBATCH -J CM018_Loom_object_generation              # Job name
#SBATCH -N 1                   # Total number of nodes requested
#SBATCH -n 8                 # Total number of cpu requested
#SBATCH -t 24:00:00           # Run time (hh:mm:ss) - 24 hours
#SBATCH --mem=50000MB  #memory in MB i.e. 50GB
#SBATCH --mail-user=adr2189@cumc.columbia.edu 
#SBATCH --mail-type=ALL

source activate base

conda activate /users/adr2189/.conda/envs/velocyto-env

HDF5_USE_FILE_LOCKING='FALSE'

velocyto run -b /groups/mb928_gp/data_backup/220502_CM_rerelease/genomecenter.columbia.edu/CM018_cellranger_count_outs/filtered_feature_bc_matrix/barcodes.tsv.gz -o /groups/mb928_gp/adr2189/velocyto_visium/loom/ -m /groups/mb928_gp/adr2189/velocyto_visium/GRCh38_rmsk.gtf /groups/mb928_gp/data_backup/220502_CM_rerelease/genomecenter.columbia.edu/CM018_cellranger_count_outs/CM018_possorted_genome_bam.bam /groups/mb928_gp/adr2189/velocyto_visium/genes.gtf
