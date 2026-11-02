#!/bin/bash

input_file="samplesheet_normals_blank.csv"
touch normals_vcf_list.txt

while IFS=, read -r col1 col2
do
# Skip the first row (header)
    if $header; then
        header=false
        continue
    fi

    echo "$col1"
    ls ${col2} >> normals_vcf_list.txt
done < "$input_file"

# merge normal germline snv .vcf files
xargs bcftools merge -m all -Oz -o hets_merged_normals.vcf.gz < vcf_list.txt

# filter for hets
bcftools view -i 'COUNT(GT="het") >= 2' hets_merged_normals.vcf.gz -Oz -o PoN_hets.vcf.gz
gunzip -k PoN_hets.vcf.gz
