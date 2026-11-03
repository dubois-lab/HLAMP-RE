#!/bin/bash

# Builds the peak overlap tables read by 01_peak_enrichment_fisher.R: one per row group in the sample sheet, and one across all peak files.
# Usage: bash 00_peak_overlaps.sh <samplesheet.csv> <output file for all samples>
# Sample sheet columns: output_file (shared by the samples of one tumor type, as named in gene_settings.R), peak_file.
# Use the .narrowPeak file of each sample, sorted by chromosome and start. Set OUT_DIR in the environment or edit the default here

set -euo pipefail

samplesheet=$1
all_samples_output=$2

OUT_DIR="${OUT_DIR:-output}"
mkdir -p "$OUT_DIR"

tail -n +2 "$samplesheet" | cut -d, -f1 | awk '!seen[$0]++' | while IFS= read -r output_file; do
  peak_files=()
  while IFS= read -r peak_file; do
    peak_files+=("$peak_file")
  done < <(awk -F, -v output_file="$output_file" 'NR > 1 && $1 == output_file {print $2}' "$samplesheet")
  bedtools multiinter -i "${peak_files[@]}" > "$OUT_DIR/$output_file"
done

all_peak_files=()
while IFS= read -r peak_file; do
  all_peak_files+=("$peak_file")
done < <(tail -n +2 "$samplesheet" | cut -d, -f2)
bedtools multiinter -i "${all_peak_files[@]}" > "$OUT_DIR/$all_samples_output"
