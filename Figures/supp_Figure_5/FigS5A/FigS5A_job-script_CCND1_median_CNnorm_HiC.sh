#!/bin/bash
#
#SBATCH --job-name=FigS5A_median_CNnorm_HiC
#SBATCH --output=logs/FigS5A_median_CNnorm_HiC_%j.out
#SBATCH --error=logs/FigS5A_median_CNnorm_HiC_%j.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-4:00
#SBATCH --mem=48G

# Set these for your system, then submit from this folder after creating logs/
IN_DIR=/path/to/input
OUT_DIR=/path/to/output
REPO_DIR=/path/to/HLAMP-RE
image=/path/to/jabba_gxg.sif
export IN_DIR OUT_DIR

# median CN-normalized HiC contact map of the 7 nf-core HiC CCND1 HLAMP tumors (3 BRCA, 4 LUSC), chr11:69-70.5Mb, with the gencode track;
# CN effect of the 20Mb background models; the 7 raw maps read in parallel
apptainer exec --bind "$IN_DIR","$OUT_DIR","$REPO_DIR" "$image" bash -c "
Rscript FigS5A_CCND1_median_CNnorm_HiC_BRCA_LUSC.R
"
