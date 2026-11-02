#!/bin/bash
#
##SBATCH --job-name=nf-core_hic
#SBATCH --output=logs/output.out
#SBATCH --error=logs/error.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8

#export TOWER_ACCESS_TOKEN=
#export NXF_APPTAINER_CACHEDIR=
#export NXF_SINGULARITY_CACHEDIR=

# mkdir ./results

nextflow run nf-core/hic -profile singularity -with-tower -params-file params.yaml -resume
