#!/bin/bash
#
#SBATCH --job-name=linkprep_hic
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=0-12:00:00
#SBATCH --mem=96G

# Builds the .hic contact map of one sample from the pairs file written by 01_align_and_make_pairs.sh.
# Usage: 02_pairs_to_hic.sh <sample>
# Set OUT_DIR, CHROM_SIZES and JUICER_JAR in the environment or edit the defaults here

set -euo pipefail

sample=$1

OUT_DIR="${OUT_DIR:-output}"
CHROM_SIZES="${CHROM_SIZES:-hs1.genome}"
JUICER_JAR="${JUICER_JAR:-juicer_tools_1.22.01.jar}"
threads="${SLURM_CPUS_PER_TASK:-24}"

echo "Processing sample: $sample"

java -Xmx80g -Djava.awt.headless=true -jar "$JUICER_JAR" pre --threads "$threads" \
    "$OUT_DIR/${sample}mapped.pairs" "$OUT_DIR/${sample}_1.hic" "$CHROM_SIZES"
