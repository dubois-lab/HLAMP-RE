#!/bin/bash
#
#SBATCH --job-name=linkprep_pairs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96
#SBATCH --time=1-12:00:00
#SBATCH --mem=384G

# Aligns one LinkPrep Hi-C library and writes its pairs file, duplicate statistics and BAM file.
# Usage: 01_align_and_make_pairs.sh <sample> <fastq_1> <fastq_2>
# Set OUT_DIR, TMP_DIR, REF_FASTA and CHROM_SIZES in the environment or edit the defaults here

set -euo pipefail

sample=$1
fastq_1=$2
fastq_2=$3

OUT_DIR="${OUT_DIR:-output}"
TMP_DIR="${TMP_DIR:-$OUT_DIR/tmp}/$sample"
REF_FASTA="${REF_FASTA:-hs1.fa}"       # indexed with bwa index
CHROM_SIZES="${CHROM_SIZES:-hs1.genome}"  # chromosome name and length, tab-separated
threads="${SLURM_CPUS_PER_TASK:-96}"
half=$((threads / 2))

mkdir -p "$OUT_DIR" "$TMP_DIR"
echo "Processing sample: $sample"

bwa mem -5SP -T0 -t"$threads" "$REF_FASTA" <(zcat "$fastq_1") <(zcat "$fastq_2") | \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
    --nproc-in "$half" --nproc-out "$half" --chroms-path "$CHROM_SIZES" | \
pairtools sort --tmpdir="$TMP_DIR" --nproc "$threads" | \
pairtools dedup --nproc-in "$half" --nproc-out "$half" --mark-dups --output-stats "$OUT_DIR/${sample}stats.txt" | \
pairtools split --nproc-in "$half" --nproc-out "$half" --output-pairs "$OUT_DIR/${sample}mapped.pairs" --output-sam - | \
samtools view -bS -@"$threads" | \
samtools sort -@"$threads" -o "$OUT_DIR/${sample}mapped.PT.bam"

samtools index "$OUT_DIR/${sample}mapped.PT.bam"
