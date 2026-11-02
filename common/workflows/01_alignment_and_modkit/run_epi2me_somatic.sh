#!/bin/bash

# obtain the workflow
nextflow run epi2me-labs/wf-somatic-variation --help


input_file="samplesheet_blank.csv"

while IFS=, read -r col1 col2
do
# Skip the first row (header)
    if $header; then
        header=false
        continue
    fi

    echo "$col1"
    mkdir -p somaticWF/${col1}
    cd somaticWF/${col1}
    mkdir results
    mkdir ${col1}_cache
cat > job-script_${col1}.sh <<EOF
#!/bin/bash
#SBATCH --job-name=epi2me-somaticWF-${col1}
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8
#SBATCH --output=epi2me-somaticWF-${col1}.out
#SBATCH --error=epi2me-somaticWF-${col1}.err

#export TOWER_ACCESS_TOKEN
#export NXF_APPTAINER_CACHEDIR
#export NXF_SINGULARITY_CACHEDIR
#export NXF_ASSETS=somaticWF/${col1}/${col1}_cache

nextflow run epi2me-labs/wf-somatic-variation --sample_name ${col1} \
	--mod \
	--ref ../../data/refGenome/USCS_T2T-CHM13v1.1_genomic.fna \
	--bam_tumor ${col2} \
	--out_dir somaticWF/${col1}/results \
	--annotation false --tumor_min_coverage 4 \
	--include_all_ctgs true \
	--modkit_threads 32 --ubam_map_threads 64 --ubam_sort_threads 24 --ubam_bam2fq_threads 8 \
	-with-tower -profile singularity \
	-resume
EOF
sbatch job-script_${col1}.sh
done < "$input_file"
