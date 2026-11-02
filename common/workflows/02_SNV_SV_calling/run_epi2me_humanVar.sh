#!/bin/bash

# obtain the workflow
nextflow run epi2me-labs/wf-human-variation --help


input_file="samplesheet_blank.csv"

while IFS=, read -r col1 col2
do
# Skip the first row (header)
    if $header; then
        header=false
        continue
    fi

    echo "$col1"
    mkdir -p humanVar/${col1}
    cd humanVar/${col1}
    mkdir results
    mkdir ${col1}_cache
cat > job-script_${col1}.sh <<EOF
#!/bin/bash
#SBATCH --job-name=epi2me-humanVar-${col1}
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8
#SBATCH --output=epi2me-humanVar-${col1}.out
#SBATCH --error=epi2me-humanVar-${col1}.err

#export TOWER_ACCESS_TOKEN
#export NXF_APPTAINER_CACHEDIR
#export NXF_SINGULARITY_CACHEDIR
#export NXF_ASSETS=humanVar/${col1}/${col1}_cache

nextflow run epi2me-labs/wf-human-variation --sample_name ${col1} \
	--snp --phased \
	--sv \
	--bam ${col2} \
	--ref ../../data/refGenome/USCS_T2T-CHM13v1.1_genomic.fna \
	--out_dir humanVar/${col1}/results \
	--annotation false \
	--include_all_ctgs true \
	--threads 96 \
	-with-tower -profile singularity \
	-resume
EOF
sbatch job-script_${col1}.sh
done < "$input_file"
