# H3K27ac and ATAC peak statistics

These scripts test which H3K27ac and ATAC peaks in the consensus amplified regions of a gene are enriched in one tumor type, and which of those peaks overlap differentially methylated regions (DMRs). Input file names and cutoffs for each gene are in `gene_settings.R`. To analyze another gene, add a block there.

Fill in `samplesheet_peaks_blank.csv` with one row per peak file, then build the overlap tables for one gene, with the output name that `gene_settings.R` gives for all samples:

```bash
bash 00_peak_overlaps.sh samplesheet_peaks_blank.csv all_NSD3_withCLs_merged.bed
```

Copy the tables from `OUT_DIR` to `IN_DIR`, then run from this folder:

```bash
Rscript 01_peak_enrichment_fisher.R CCND1
Rscript 01_peak_enrichment_fisher.R NSD3
Rscript 01_peak_enrichment_fisher.R EGFR
Rscript 02_peak_dmr_overlap.R NSD3
Rscript 02_peak_dmr_overlap.R EGFR
```

`01_peak_enrichment_fisher.R` writes the `*_linDiffs_stats.csv` tables read by `02_peak_dmr_overlap.R`, so copy them from `OUT_DIR` to `IN_DIR` before running the second script. `N_CORES` sets the number of cores (default 10).

## Inputs in `IN_DIR`

- Peak overlap tables for each tumor type and for all samples together, written by `00_peak_overlaps.sh` with `bedtools multiinter` from the narrowPeak file of each sample, not from MACS summits. The columns are chromosome, start, end, number of samples with a peak, sample list, then one presence column per sample. Coordinates are hg38.
- Consensus amplification tracks: rds files holding a gTrack whose first data element is a GRanges in hs1 with an `AmpFreq` column. For NSD3 and EGFR, the `_fisher` and `_overlap` tracks are separate files, the versions used by each script.
- DMR table from `common/workflows/07_differential_methylation`, in hs1, with columns `chr`, `start`, `end` and `diff.Methy` (the first DSS group, LUSC, minus the second, BRCA).
- `hs1ToHg38.over.chain`, which `01_peak_enrichment_fisher.R` downloads from UCSC if it is missing.

## Test details

The group sizes in the contingency table come from the tables rather than a sample list: `n_tt` is the highest sample count among the tested peaks of the tumor type, and the size of the other group is the highest count across all samples minus `n_tt`. These equal the number of samples if at least one interval is shared by all of them. Peaks with no interval sharing more than 5 bp in the all-sample table get p = 1 and an empty `compttype_count`.

# BEDTools (https://github.com/arq5x/bedtools2)

Quinlan, A.R. and Hall, I.M. *BEDTools: a flexible suite of utilities for comparing genomic features*. Bioinformatics, 2010. 26(6): p. 841-842.

# qvalue (https://bioconductor.org/packages/qvalue)

Storey, J.D. and Tibshirani, R. *Statistical significance for genomewide studies*. PNAS, 2003. 100(16): p. 9440-9445.
