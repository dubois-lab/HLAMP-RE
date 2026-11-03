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
CONFIG=/path/to/config.yaml # see config_blank.yaml
WORKFLOW_PROFILE=/path/to/workflow_profile # SLURM settings for your cluster
SINGULARITY_DIR=/path/to/singularity_images

SCRIPT_DIR="$PWD"
mkdir -p "$WORK_DIR/results"
cp 01_ichorCNA_runscript.R "$WORK_DIR/"
cd "$WORK_DIR"

snakemake --executor slurm -j unlimited --workflow-profile "$WORKFLOW_PROFILE" --rerun-incomplete \
    -s "$SCRIPT_DIR/Snakefile" --configfile "$CONFIG" --use-singularity --singularity-prefix "$SINGULARITY_DIR" --latency-wait 3
