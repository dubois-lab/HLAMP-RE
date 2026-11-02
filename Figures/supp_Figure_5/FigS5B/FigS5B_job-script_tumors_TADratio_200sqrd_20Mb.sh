#!/bin/bash
#
#SBATCH --job-name=TADratio_200sqrd_20Mb_tumors
#SBATCH --output=logs/TADratio_200sqrd_20Mb_tumor_%a.out
#SBATCH --error=logs/TADratio_200sqrd_20Mb_tumor_%a.err
#SBATCH --array=0-25
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-16:00
#SBATCH --mem=128G

# Set these for your system, then submit from this folder after creating logs/
IN_DIR=/path/to/input
OUT_DIR=/path/to/output
REPO_DIR=/path/to/HLAMP-RE
image=/path/to/jabba_gxg.sif
export IN_DIR OUT_DIR

# lung TAD vs two random window sets, partners up to 20Mb, counts / CN effect of the 200 x 200 cell background model (pairs up to 20Mb)
# array index -> tumor and gene arm: 0-8 CCND1 (C169, C17, C3, H3D44S, KDAMIF, Y57U8R, F57J9F, C15, 1WNRWL) -> CCND1/;
#   9-17 NSD3 (C12, C132, C17, Cb5, 1IDEVD, 1E2AVA, NODTWC, K0HCJ3, C168) -> NSD3/; 18-25 EGFR (C109, C15, C168, C208, 1HGXBV, 1WNRWL, AJZ5JO,
#   XO7JTV) -> EGFR/
# background model: fitted once per tumor (CCND1 folder for C169-F57J9F incl. C17, NSD3 folder for the other NSD3 tumors, EGFR folder for the
# linkprep tumors); JaBbA: ampChr chr11 -> CCND1 folder, chr8 -> FGFR1 folder; HiC: raw 5kb maps of the nf-core runs or linkprep .hic files
# (NODTWC map named NODWTC)
samples=(C169 C17 C3 H3D44S KDAMIF Y57U8R F57J9F C15 1WNRWL C12 C132 C17 Cb5 1IDEVD 1E2AVA NODTWC K0HCJ3 C168 C109 C15 C168 C208 1HGXBV 1WNRWL AJZ5JO XO7JTV)
genes=(CCND1 CCND1 CCND1 CCND1 CCND1 CCND1 CCND1 CCND1 CCND1 NSD3 NSD3 NSD3 NSD3 NSD3 NSD3 NSD3 NSD3 NSD3 EGFR EGFR EGFR EGFR EGFR EGFR EGFR EGFR)
msets=(CCND1 CCND1 CCND1 CCND1 CCND1 CCND1 CCND1 EGFR EGFR NSD3 NSD3 CCND1 NSD3 NSD3 NSD3 NSD3 NSD3 EGFR EGFR EGFR EGFR EGFR EGFR EGFR EGFR EGFR)
ampchrs=(chr11 chr8 chr11 chr11 chr11 chr11 chr11 chr11 chr11 chr8 chr8 chr8 chr8 chr8 chr8 chr8 chr8 chr11 chr11 chr11 chr11 chr11 chr11 chr11 chr11 chr11)
tissues=(breast breast breast lung lung lung lung breast lung breast breast breast breast lung lung lung lung breast breast breast breast breast lung lung lung lung)
i=$SLURM_ARRAY_TASK_ID
sid=${samples[$i]}
gene=${genes[$i]}
map_id=$sid
case $sid in
    C169|C17|C12|C132|Cb5|1IDEVD|1E2AVA|K0HCJ3) hic=$IN_DIR/hic/nfcore_results ;;
    NODTWC) hic=$IN_DIR/hic/nfcore_results; map_id=NODWTC ;;
    H3D44S) hic=$IN_DIR/hic/H3D44S_results ;;
    C3|KDAMIF|Y57U8R|F57J9F) hic=$IN_DIR/hic/arima_results ;;
    *) hic=$IN_DIR/hic/linkprep/${sid}_1.hic ;;
esac
model=$IN_DIR/CNmodels/${msets[$i]}/20261023_${sid}_chr7_8_11_20Mb_HiCcov_logistic_df8_x200_d200.rds
outdir=$OUT_DIR/TADratio/${gene}/
echo "Processing sample: $sid ${ampchrs[$i]} ${tissues[$i]} $hic (map $map_id); gene $gene; model $model; output in $outdir"

apptainer exec --bind "$IN_DIR","$OUT_DIR","$REPO_DIR" "$image" bash -c "
Rscript FigS5B_lungTAD_inout_ratio_tumor_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_20Mb_200sqrd.R $sid ${ampchrs[$i]} ${tissues[$i]} $hic $gene $outdir $model $map_id
"
