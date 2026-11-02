#!/bin/bash
#
#SBATCH --job-name=normHiC_BRCA_NSCLC_ESCC_EGFR
#SBATCH --output=logs/normHiC_20261004_BRCA_NSCLC_ESCC_EGFR.out
#SBATCH --error=logs/normHiC_20261004_BRCA_NSCLC_ESCC_EGFR.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=0-2:00
#SBATCH --mem=96G

# Set these for your system, then submit from this folder after creating logs/
IN_DIR=/path/to/input
OUT_DIR=/path/to/output
REPO_DIR=/path/to/HLAMP-RE
image=/path/to/jabba_gxg.sif
export IN_DIR OUT_DIR

apptainer exec --bind "$IN_DIR","$OUT_DIR","$REPO_DIR" "$image" bash -c "
Rscript FigS7F_normHiC_BRCA_NSCLC_ESCC_EGFR.R
"
