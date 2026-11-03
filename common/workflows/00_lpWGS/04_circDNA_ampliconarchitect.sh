#!/bin/bash
#SBATCH --job-name=circDNA_ampliconarchitect
#SBATCH --time=2-00:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=2
#SBATCH --output=logs/circDNA_ampliconarchitect.out
#SBATCH --error=logs/circDNA_ampliconarchitect.err

# Set these for your system, then submit from this folder after creating logs/
# Needs Nextflow, Singularity, the AmpliconArchitect data repository and a MOSEK license
SAMPLESHEET=/path/to/samplesheet_circdna.csv # see samplesheet_circdna_blank.csv
OUT_DIR=/path/to/output/circDNA
WORK_DIR=/path/to/scratch/circDNA
CLUSTER_CONFIG=/path/to/cluster.config
AA_DATA_REPO=/path/to/aa_data_repo
MOSEK_LICENSE_DIR=/path/to/mosek_license_dir
CONTAINER_CACHE=/path/to/singularity
export NXF_VER=23.10.1
export NXF_APPTAINER_CACHEDIR="$CONTAINER_CACHE"
export NXF_SINGULARITY_CACHEDIR="$CONTAINER_CACHE"

mkdir -p "$OUT_DIR" "$WORK_DIR"
cd "$WORK_DIR"

# Skips QC and duplicate marking, which already ran in sarek
nextflow run nf-core/circdna --input "$SAMPLESHEET" \
    --outdir "$OUT_DIR" --genome GRCh38 \
    -profile singularity --circle_identifier ampliconarchitect --input_format 'BAM' --bam_sorted true --save_markduplicates_bam false \
    -c "$CLUSTER_CONFIG" \
    --aa_data_repo "$AA_DATA_REPO" --mosek_license_dir "$MOSEK_LICENSE_DIR" \
    --reference_build GRCh38 -resume
