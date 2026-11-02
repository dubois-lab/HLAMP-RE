#!/bin/bash

input_file="./samplesheet_blank.csv"
REF="../../../common/data/refGenome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
dir="./"

mkdir -p ${dir}/input ${dir}/output/vep

while IFS= read -r sid; do
# Skip the first row (header)
    if $header; then
        header=false
        continue
    fi

    echo "$sid"
    #filter
    bcftools view -f PASS -Oz -o ${dir}/input/${sid}_tumor_only.PASS.vcf.gz /path/to/out/${sid}_tumor_only.vcf.gz
    
    vep -i ${dir}/input/${sid}_tumor_only.PASS.vcf.gz \
     -o ${dir}/output/vep/${sid}_tumor_only.vep.vcf.gz \
     --vcf \
     --compress_output bgzip \
     --cache --offline --dir_cache ${dir}/input/vep_cache \
     --fasta ${REF} --assembly GRCh38 \
     --symbol --canonical --protein --biotype --hgvs \
     --fork 4

done < "$input_file"
