#!/bin/bash
#
#SBATCH --job-name=FigS7A_B_median_CNnorm_HiC
#SBATCH --output=logs/FigS7A_B_median_CNnorm_HiC_%a.out
#SBATCH --error=logs/FigS7A_B_median_CNnorm_HiC_%a.err
#SBATCH --array=0-1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-4:00
#SBATCH --mem=64G

# Set these for your system, then submit from this folder after creating logs/
IN_DIR=/path/to/input
OUT_DIR=/path/to/output
REPO_DIR=/path/to/HLAMP-RE
image=/path/to/jabba_gxg.sif
export IN_DIR OUT_DIR

# median CN-normalized HiC contact maps, BRCA and LUSC pages, with the gencode track; array index -> gene: 0 NSD3 (chr8:37.9-39.4Mb, nf-core HiC
# tumors), 1 EGFR (chr7:54.5-56Mb, linkprep tumors); CN effect of the 20Mb background models; the 8 maps of
# a gene read in parallel
genes=(NSD3 EGFR)
gene=${genes[$SLURM_ARRAY_TASK_ID]}
echo "gene $gene"
apptainer exec --bind "$IN_DIR","$OUT_DIR","$REPO_DIR" "$image" bash -c "
Rscript FigS7A_B_NSD3_EGFR_median_CNnorm_HiC_BRCA_LUSC_binfix.R $gene
"
