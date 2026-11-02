#!/bin/bash
#
#SBATCH --job-name=FigS4_CCND1_tumor_plots
#SBATCH --output=logs/FigS4_CCND1_tumor_plots_%a.out
#SBATCH --error=logs/FigS4_CCND1_tumor_plots_%a.err
#SBATCH --array=0-8
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

# per-tumor plots of the CCND1 HLAMP tumors with DMR and significance tracks (array 0-8:
# C169, C17, C3 (BRCA), H3D44S, KDAMIF, Y57U8R, F57J9F (LUSC), C15 (BRCA, linkprep), 1WNRWL (LUSC, linkprep))
samples=(C169 C17 C3 H3D44S KDAMIF Y57U8R F57J9F C15 1WNRWL)
msets=(CCND1 CCND1 CCND1 CCND1 CCND1 CCND1 CCND1 EGFR EGFR)
ampchrs=(chr11 chr8 chr11 chr11 chr11 chr11 chr11 chr11 chr11)
tissues=(breast breast breast lung lung lung lung breast lung)
i=$SLURM_ARRAY_TASK_ID
sid=${samples[$i]}
case $sid in
    C169|C17) hic=$IN_DIR/hic/nfcore_results ;;
    H3D44S) hic=$IN_DIR/hic/H3D44S_results ;;
    C3|KDAMIF|Y57U8R|F57J9F) hic=$IN_DIR/hic/arima_results ;;
    *) hic=$IN_DIR/hic/linkprep/${sid}_1.hic ;;
esac
model=$IN_DIR/CNmodels/${msets[$i]}/20261023_${sid}_chr7_8_11_20Mb_HiCcov_logistic_df8_x200_d200.rds
outdir=$OUT_DIR/
echo "Plots: $sid ${ampchrs[$i]} ${tissues[$i]} $hic; model $model; output in $outdir"

apptainer exec --bind "$IN_DIR","$OUT_DIR","$REPO_DIR" "$image" bash -c "
Rscript FigS4_CCND1_tumor_plots_hiccovlogisticnorm_DMR_sig_20Mb.R $sid ${ampchrs[$i]} ${tissues[$i]} $hic CCND1 $outdir $model $sid
"
