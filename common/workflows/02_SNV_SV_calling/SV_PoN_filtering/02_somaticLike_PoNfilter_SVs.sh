#!/bin/bash

input_file="samplesheet_tumor_blank.csv"
logfile="logfile.log"
sv_dir="../humanVar"
pon_dir="."

mkdir output

while IFS=, read -r col1 col2
do
# Skip the first row (header)
    if $header; then
        header=false
        continue
    fi
    
  echo "Starting ${col1}.." | tee -a "$logfile"
  gunzip -k ${sv_dir}/${col1}/results/${col1}.wf_sv.vcf.gz 2>&1 | tee -a "$logfile"
  total_variants=$(grep -v '^#' "${sv_dir}/${col1}/results/${col1}.wf_sv.vcf" | wc -l)
  echo "Total variants = $total_variants" | tee -a "$logfile"
  ls ${pon_dir}/sniffles_survivor_merged_normals.sorted.vcf > output/tmp_list.txt
  ls ${sv_dir}/${col1}/results/${col1}.wf_sv.vcf >> output/tmp_list.txt 

  SURVIVOR merge output/tmp_list.txt 5 1 1 1 0 500 output/${col1}.tmp_merged.vcf 2>&1 | tee -a "$logfile"

  bcftools query -i 'INFO/SUPP_VEC="01"' -f '%ID\n' output/${col1}.tmp_merged.vcf > output/${col1}.somaticIDs.txt 2>&1 | tee -a "$logfile"
  bcftools view -i "ID=@output/${col1}.somaticIDs.txt" "${sv_dir}/${col1}/results/${col1}.wf_sv.vcf" -O z -o "output/${col1}.sniffles_sv.PoNfiltered.vcf.gz" 2>&1 | tee -a "$logfile"
  tabix -p vcf "output/${col1}.sniffles_sv.PoNfiltered.vcf.gz" 2>&1 | tee -a "$logfile"
  somatic_variants=$(zgrep -v '^#' "output/${col1}.sniffles_sv.PoNfiltered.vcf.gz" | wc -l)
  filtered_out=$((total_variants - somatic_variants))
  echo "Number of SV's filtered out = $filtered_out" | tee -a "$logfile"

  # CleanUp
  rm ${sv_dir}/${col1}/results/${col1}.wf_sv.vcf 2>&1 | tee -a "$logfile" # remove uncompressed .vcf
  rm output/${col1}.tmp_merged.vcf output/${col1}.somaticIDs.txt output/tmp_list.txt

done < "$input_file"
