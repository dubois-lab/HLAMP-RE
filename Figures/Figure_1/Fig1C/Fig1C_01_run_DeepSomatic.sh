#!/bin/bash

input_file="./samplesheet_blank.csv"
REF="../../../common/data/refGenome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
TUMOR="/path/to/input"
OUT="/path/to/out"
image="/path/to/deepsomatic_1.10.0.sif"

while IFS= read -r sid; do
# Skip the first row (header)
    if $header; then
        header=false
        continue
    fi

    echo "$sid"
    mkdir -p ${OUT}/${sid}

    apptainer exec -B /path/local/:/path/mount/ \
                ${image} \
        	run_deepsomatic \
        	--model_type=ONT_TUMOR_ONLY \
        	--ref=${REF} \
        	--reads_tumor=${TUMOR}/${sid}.haplotagged.bam \
        	--output_vcf=${OUT}/${sid}/${sid}_tumor_only.vcf.gz \
        	--sample_name_tumor=${sid} \
        	--num_shards=32 \
        	--logging_dir=${OUT}/${sid}/logs \
        	--intermediate_results_dir=${OUT}/${sid}/intermediate_results \
        	--use_default_pon_filtering=true
    
done < "$input_file"
