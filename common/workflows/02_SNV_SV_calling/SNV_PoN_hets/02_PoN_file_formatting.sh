#!/bin/bash

pon="./PoN_hets.vcf"

awk '
BEGIN { OFS="\t"; print "seqnames","start","end","alt_count","ref_count" }
!/^#/ {
  alt_sum = 0; ref_sum = 0;
  for (i = 10; i <= NF; i++) {
    n = split($i, f, ":");
    if (n >= 4) {
      split(f[4], a, ",");
      if (a[1] ~ /^[0-9]+$/ && a[2] ~ /^[0-9]+$/) {
        ref_sum += a[1];
        alt_sum += a[2];
      }
    }
  }
  print $1, $2, $2, alt_sum, ref_sum;
}' ${pon} > PoN_hets_counts.tsv
