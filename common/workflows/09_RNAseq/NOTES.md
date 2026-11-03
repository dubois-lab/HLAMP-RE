# Differential expression by copy number

`01_DESeq2_copy_number_groups.R` compares samples within one tumor type by the copy number status of one gene. Run it from this folder:

```bash
# Tables read by Fig. S1D
for gene in CCND1 EGFR FGFR1; do
  for tumor in BRCA NSCLC; do
    Rscript 01_DESeq2_copy_number_groups.R "$gene" "$tumor" AMP
  done
done

# CCND1 gained or amplified, against normal copy number
Rscript 01_DESeq2_copy_number_groups.R CCND1 BRCA GAIN
Rscript 01_DESeq2_copy_number_groups.R CCND1 NSCLC GAIN
```

The script needs:

- `salmon.merged.gene_counts*.tsv` in `IN_DIR`, one per batch, as written by nf-core/rnaseq
- `samplesheet_blank.csv` filled in with one row per sample: the column name in the count tables, `group` (tumor type), `subgroup`, and the status of each gene (`AMP`, `GAIN`, or empty for normal copy number)

# DESeq2 (https://bioconductor.org/packages/DESeq2)

Love, M.I. et al. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology, 2014. 15(12): p. 550.

# Salmon (https://github.com/COMBINE-lab/salmon)

Patro, R. et al. *Salmon provides fast and bias-aware quantification of transcript expression*. Nature Methods, 2017. 14(4): p. 417-419.

# nf-core/rnaseq (https://nf-co.re/rnaseq)

Ewels, P.A. et al. *The nf-core framework for community-curated bioinformatics pipelines*. Nature Biotechnology, 2020. 38(3): p. 276-278.
