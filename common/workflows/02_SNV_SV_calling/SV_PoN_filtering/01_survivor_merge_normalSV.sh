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

# merge SV normals
SURVIVOR merge normals_vcf_list.txt 5 1 1 1 0 500 sniffles_survivor_merged_normals.unsorted.vcf

# sort vcf
bcftools sort -O z -o sniffles_survivor_merged_normals.sorted.vcf.gz sniffles_survivor_merged_normals.unsorted.vcf
tabix -p vcf sniffles_survivor_merged_normals.sorted.vcf.gz
gunzip -k sniffles_survivor_merged_normals.sorted.vcf.gz
