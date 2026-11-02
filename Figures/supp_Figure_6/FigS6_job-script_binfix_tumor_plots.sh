#!/bin/bash
#
#SBATCH --job-name=binfix_tumor_plots
#SBATCH --output=logs/tumor_plots_%A_%a.out
#SBATCH --error=logs/tumor_plots_%A_%a.err
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-3:00
#SBATCH --mem=64G

# Set these for your system, then submit from this folder after creating logs/
IN_DIR=/path/to/input
OUT_DIR=/path/to/output
REPO_DIR=/path/to/HLAMP-RE
image=/path/to/jabba_gxg.sif
export IN_DIR OUT_DIR

# per-tumor plots of the NSD3 tumors with 0-based contact maps, with the refitted logistic 200 x 200 models
# array index -> tumor (LUSC): 0 1IDEVD, 1 NODTWC (map NODWTC), 2 K0HCJ3
samples=(1IDEVD NODTWC K0HCJ3)
map_ids=(1IDEVD NODWTC K0HCJ3)
sid=${samples[$SLURM_ARRAY_TASK_ID]}
map_id=${map_ids[$SLURM_ARRAY_TASK_ID]}
hic=$IN_DIR/hic/nfcore_results
model=$IN_DIR/CNmodels_binfix/NSD3/20261023_${sid}_chr7_8_11_20Mb_HiCcov_logistic_df8_x200_d200.rds
script=FigS6_NSD3_EGFR_tumor_plots_hiccovlogisticnorm_DMR_20Mb_binfix.R
outdir=$OUT_DIR/binfix/NSD3/
echo "Plots ($script): $sid chr8 lung $hic (map $map_id); gene NSD3; model $model; output in $outdir"
apptainer exec --bind "$IN_DIR","$OUT_DIR","$REPO_DIR" "$image" bash -c "
Rscript $script $sid chr8 lung $hic NSD3 $outdir $model $map_id
"
