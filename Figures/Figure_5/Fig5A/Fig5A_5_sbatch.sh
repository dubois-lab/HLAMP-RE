#!/bin/bash
#
#SBATCH --job-name=Fig5a_20260914_CCND1_TF_BRCA
#SBATCH --output=logs/Fig5a_20260914_CCND1_TF_BRCA.out
#SBATCH --error=logs/Fig5a_20260914_CCND1_TF_BRCA.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=0-2:00:00
#SBATCH --mem=40G

# Set these for your system, then submit from this folder after creating logs/
IN_DIR=/path/to/input
OUT_DIR=/path/to/output
REPO_DIR=/path/to/HLAMP-RE
image=/path/to/jabba_gxg.sif
export IN_DIR OUT_DIR

apptainer exec --bind "$IN_DIR","$OUT_DIR","$REPO_DIR" "$image" bash -c "

Rscript Fig5A_5.r
"
