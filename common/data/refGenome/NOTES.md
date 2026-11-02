Main assembly: T2T-CHM13 v1.1

Used by `common/workflows/06_consensus_amplicon/HLAMP_walks_combined_hs1`:

- `hg19ToHs1.over.chain` — UCSC liftOver chain hg19 → hs1 (T2T-CHM13 v2.0),
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHs1.over.chain.gz (decompressed)
- `hs1.chrom.sizes` — UCSC hs1 chromosome sizes (incl. chrY)
- `cgc_T2T_gr_consAmpWalks.rds` — Cancer Gene Census GRanges on T2T used to place the
  gene windows (CCND1, NSD3, EGFR) in that workflow; gene coordinates differ slightly
  from `cgc_T2T_gr_UCSC.rds`
