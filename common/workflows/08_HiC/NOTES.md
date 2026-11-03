# nf-core/hic pipeline (https://github.com/nf-core/hic.git)
Servant, N., et al., nf-core/hic: nf-core/hic v2.1.0. Zenodo, 2023. doi:10.5281/zenodo.2669512

# HiCDCPlus (https://github.com/mervesa/HiCDCPlus.git)
Sahin, M., et al., *HiC-DC+ enables systematic 3D interaction calls and differential analysis for Hi-C and HiChIP*. Nat Commun, 2021. 12(1): p. 3366.

# Virtual 4C tracks

`05_interactions_table.R`, `06_virtual_4C_tracks.R` and `07_combine_virtual_4C_pvalues.R` turn the HiCDCPlus results of `02_HiCDCPlus_signInteractions.R` into the virtual 4C tracks at a gene's TSS. Run them from this folder:

```bash
Rscript 05_interactions_table.R C3 chr7 chr8 chr11
Rscript 06_virtual_4C_tracks.R C3 CCND1 NSD3
Rscript 07_combine_virtual_4C_pvalues.R CCND1 BRCA C3 C169 C17
```

Run the first two scripts for each sample, and the third for each gene and group. Besides the results of `02_HiCDCPlus_signInteractions.R`, `IN_DIR` needs:

- `normalized_abs.bed`: the 5 kb bins of the contact matrix in hs1, with chromosome, start (0-based) and end in the first three columns
- `<gene>_TSSgr_T2T.rds`: a GRanges with the TSS of the gene in hs1

The figure scripts read the per-sample tracks from `individuals/` under `IN_DIR`.

# gGnome (https://github.com/mskilab-org/gGnome)

Hadi, K. et al. *Distinct classes of complex structural variation uncovered across thousands of cancer genome graphs*. Cell, 2020. 183(1): p. 197-210.

# poolr (https://github.com/ozancinar/poolR)

Cinar, O. and Viechtbauer, W. *The poolr package for combining independent and dependent p values*. Journal of Statistical Software, 2022. 101(1): p. 1-42.
