#!/bin/bash
#SBATCH --job-name=ichorCNA_snakemake
#SBATCH --time=4-10:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=4
#SBATCH --output=logs/ichorCNA_snakemake.out
#SBATCH --error=logs/ichorCNA_snakemake.err

# Set these for your system, then submit from this folder after creating logs/
# Needs Snakemake with the SLURM executor plugin, and Singularity
WORK_DIR=/path/to/lpWGS
SNAKEFILE=/path/to/Snakefile
WORKFLOW_PROFILE=/path/to/workflow_profile
SINGULARITY_DIR=/path/to/singularity_images

cd "$WORK_DIR"
mkdir -p results

snakemake --executor slurm -j unlimited --workflow-profile "$WORKFLOW_PROFILE" --rerun-incomplete \
    -s "$SNAKEFILE" --use-singularity --singularity-prefix "$SINGULARITY_DIR" --latency-wait 3
