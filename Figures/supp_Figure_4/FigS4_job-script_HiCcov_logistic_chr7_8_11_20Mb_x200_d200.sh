#!/bin/bash
#
#SBATCH --job-name=HiCcov_logistic_20Mb_x200_d200
#SBATCH --output=logs/HiCcov_logistic_20Mb_x200_d200_%a.out
#SBATCH --error=logs/HiCcov_logistic_20Mb_x200_d200_%a.err
#SBATCH --array=0-22
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-16:00
#SBATCH --mem=192G

# Set these for your system, then submit from this folder after creating logs/
IN_DIR=/path/to/input
OUT_DIR=/path/to/output
REPO_DIR=/path/to/HLAMP-RE
image=/path/to/jabba_gxg.sif
export IN_DIR OUT_DIR

# joint chr7/chr8/chr11 logistic HiC coverage CN model, pairs up to 20Mb, distance decay spline df 8, fit set from 200 x 200 cells
# array index -> tumor and gene set:
#   0-6 CCND1 (C169, C17, C3, H3D44S, KDAMIF, Y57U8R, F57J9F; raw 5kb maps of the nf-core runs) -> CCND1/
#   7-14 NSD3 (BRCA: C12, C132, C17, Cb5; LUSC: 1IDEVD, 1E2AVA, NODTWC, K0HCJ3; raw 5kb maps of the nf-core run, NODTWC map named NODWTC) -> NSD3/
#   15-22 EGFR (BRCA: C109, C15, C168, C208; LUSC: 1HGXBV, 1WNRWL, AJZ5JO, XO7JTV; linkprep .hic files) -> EGFR/
# the gene only sets the TAD marked in the plots
samples=(C169 C17 C3 H3D44S KDAMIF Y57U8R F57J9F C12 C132 C17 Cb5 1IDEVD 1E2AVA NODTWC K0HCJ3 C109 C15 C168 C208 1HGXBV 1WNRWL AJZ5JO XO7JTV)
map_ids=(C169 C17 C3 H3D44S KDAMIF Y57U8R F57J9F C12 C132 C17 Cb5 1IDEVD 1E2AVA NODWTC K0HCJ3 C109 C15 C168 C208 1HGXBV 1WNRWL AJZ5JO XO7JTV)
ccnd1_dirs=($IN_DIR/hic/nfcore_results $IN_DIR/hic/nfcore_results
    $IN_DIR/hic/arima_results $IN_DIR/hic/H3D44S_results
    $IN_DIR/hic/arima_results $IN_DIR/hic/arima_results $IN_DIR/hic/arima_results)
i=$SLURM_ARRAY_TASK_ID
sid=${samples[$i]}
map_id=${map_ids[$i]}
if [ $i -lt 7 ]; then
    gene=CCND1; hic=${ccnd1_dirs[$i]}
elif [ $i -lt 15 ]; then
    gene=NSD3; hic=$IN_DIR/hic/nfcore_results
else
    gene=EGFR; hic=$IN_DIR/hic/linkprep/${sid}_1.hic
fi
echo "Processing sample: $sid $hic (map $map_id); gene $gene; df 8; pairs up to 20Mb; 200 x 200 cells"

apptainer exec --bind "$IN_DIR","$OUT_DIR","$REPO_DIR" "$image" bash -c "
Rscript FigS4_chr7_8_11_20Mb_HiCcov_logistic_CNmodel_x200_d200_tumor_gene_hic.R $sid $hic 8 $OUT_DIR/CNmodels/$gene/ $gene $map_id
"
