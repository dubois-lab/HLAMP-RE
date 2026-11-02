#!/bin/bash
#
#SBATCH --job-name=FigS6_NSD3_EGFR_tumor_plots
#SBATCH --output=logs/FigS6_tumor_plots_%A_%a.out
#SBATCH --error=logs/FigS6_tumor_plots_%A_%a.err
#SBATCH --array=0-16
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

# per-tumor plots of the NSD3 and EGFR HLAMP tumors with DMR tracks and the CN-normalized v4C saved as rds
# (array 0-8 NSD3, 9-16 EGFR)
samples=(C12 C132 C17 Cb5 1IDEVD 1E2AVA NODTWC K0HCJ3 C168 C109 C15 C168 C208 1HGXBV 1WNRWL AJZ5JO XO7JTV)
genes=(NSD3 NSD3 NSD3 NSD3 NSD3 NSD3 NSD3 NSD3 NSD3 EGFR EGFR EGFR EGFR EGFR EGFR EGFR EGFR)
msets=(NSD3 NSD3 CCND1 NSD3 NSD3 NSD3 NSD3 NSD3 EGFR EGFR EGFR EGFR EGFR EGFR EGFR EGFR EGFR)
ampchrs=(chr8 chr8 chr8 chr8 chr8 chr8 chr8 chr8 chr11 chr11 chr11 chr11 chr11 chr11 chr11 chr11 chr11)
tissues=(breast breast breast breast lung lung lung lung breast breast breast breast breast lung lung lung lung)
i=$SLURM_ARRAY_TASK_ID
sid=${samples[$i]}
gene=${genes[$i]}
map_id=$sid
case $sid in
    C12|C132|C17|Cb5|1IDEVD|1E2AVA|K0HCJ3) hic=$IN_DIR/hic/nfcore_results ;;
    NODTWC) hic=$IN_DIR/hic/nfcore_results; map_id=NODWTC ;;
    *) hic=$IN_DIR/hic/linkprep/${sid}_1.hic ;;
esac
model=$IN_DIR/CNmodels/${msets[$i]}/20261023_${sid}_chr7_8_11_20Mb_HiCcov_logistic_df8_x200_d200.rds
outdir=$OUT_DIR/${gene}/
echo "Plots: $sid ${ampchrs[$i]} ${tissues[$i]} $hic (map $map_id); gene $gene; model $model; output in $outdir"

apptainer exec --bind "$IN_DIR","$OUT_DIR","$REPO_DIR" "$image" bash -c "
Rscript FigS6_NSD3_EGFR_tumor_plots_hiccovlogisticnorm_DMR_20Mb.R $sid ${ampchrs[$i]} ${tissues[$i]} $hic $gene $outdir $model $map_id
"
