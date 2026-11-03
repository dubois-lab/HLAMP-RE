#!/bin/bash
#SBATCH --job-name=sarek_controlfreec
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=2
#SBATCH --output=logs/sarek_controlfreec.out
#SBATCH --error=logs/sarek_controlfreec.err

# Set these for your system, then submit from this folder after creating logs/
# Needs Nextflow and Apptainer
SAMPLESHEET=/path/to/samplesheet_sarek.csv # see samplesheet_sarek_blank.csv
OUT_DIR=/path/to/output/sarek_controlfreec
WORK_DIR=/path/to/scratch/sarek_controlfreec
CLUSTER_CONFIG=/path/to/cluster.config
CONTAINER_CACHE=/path/to/singularity
export NXF_APPTAINER_CACHEDIR="$CONTAINER_CACHE"
export NXF_SINGULARITY_CACHEDIR="$CONTAINER_CACHE"

mkdir -p "$OUT_DIR" "$WORK_DIR"
cd "$WORK_DIR"

# Starts from aligned BAM files
nextflow run nf-core/sarek -profile apptainer --input "$SAMPLESHEET" \
    -c "$CLUSTER_CONFIG" --step variant_calling \
    --outdir "$OUT_DIR" --genome GATK.GRCh38 --tools controlfreec --cf_window 10000 \
    --save_output_as_bam
