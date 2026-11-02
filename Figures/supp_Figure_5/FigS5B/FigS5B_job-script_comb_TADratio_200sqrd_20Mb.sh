#!/bin/bash
#
#SBATCH --job-name=TADratio_200sqrd_20Mb_comb
#SBATCH --output=logs/TADratio_200sqrd_20Mb_comb_%j.out
#SBATCH --error=logs/TADratio_200sqrd_20Mb_comb_%j.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-24:00
#SBATCH --mem=128G

# Set these for your system, then submit from this folder after creating logs/
IN_DIR=/path/to/input
OUT_DIR=/path/to/output
REPO_DIR=/path/to/HLAMP-RE
image=/path/to/jabba_gxg.sif
export IN_DIR OUT_DIR

# normals (partners <= 20Mb, two sets of 200 random windows each across chr7/chr8/chr11; the per-bin tables
# do not depend on the tumor background model) + TAD vs random window plots and permutation test for one gene; needs the
# per-tumor TSVs from FigS5B_job-script_tumors_TADratio_200sqrd_20Mb.sh
# usage: sbatch FigS5B_job-script_comb_TADratio_200sqrd_20Mb.sh <gene> <tumors, comma separated>
gene=$1
tumors=$2
echo "gene $gene; tumors $tumors"

apptainer exec --bind "$IN_DIR","$OUT_DIR","$REPO_DIR" "$image" bash -c "
Rscript FigS5B_lungTAD_inout_ratio_normals20Mb_vs_tumors_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_200sqrd.R $gene $OUT_DIR/TADratio/$gene/ $tumors
"
