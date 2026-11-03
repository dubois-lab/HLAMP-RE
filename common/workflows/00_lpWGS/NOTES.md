# Low-pass whole-genome sequencing

| Script | Step |
|---|---|
| `01_ichorCNA_snakemake.sh` | copy number segments per sample with ichorCNA |
| `02_sarek_controlfreec.sh` | copy number with Control-FREEC, from aligned BAM files |
| `03_ichorCNA_amplified_oncogenes.R` | amplified oncogenes in the ichorCNA segments, with IGV segment files per gene |
| `04_circDNA_ampliconarchitect.sh` | amplicon reconstruction with AmpliconArchitect, from aligned BAM files |

`03_ichorCNA_amplified_oncogenes.R` reads these files from `IN_DIR`:

- `ichorCNA/results/<sample>/<sample>.seg` and `<sample>.seg.txt`
- `20240621_ichorCNA_IGV.seg` and `20240621_ichorCNA.seg.rds`, the segments of an earlier ichorCNA run

# ichorCNA (https://github.com/broadinstitute/ichorCNA)

Adalsteinsson, V.A. et al. *Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors*. Nature Communications, 2017. 8(1): p. 1324.

# nf-core/sarek (https://nf-co.re/sarek)

Garcia, M. et al. *Sarek: A portable workflow for whole-genome sequencing analysis of germline and somatic variants*. F1000Research, 2020. 9: p. 63.

# Control-FREEC (https://github.com/BoevaLab/FREEC)

Boeva, V. et al. *Control-FREEC: a tool for assessing copy number and allelic content using next-generation sequencing data*. Bioinformatics, 2012. 28(3): p. 423-425.

# nf-core/circdna (https://nf-co.re/circdna)

https://doi.org/10.5281/zenodo.6685250

# AmpliconArchitect (https://github.com/AmpliconSuite/AmpliconArchitect)

Deshpande, V. et al. *Exploring the landscape of focal amplifications in cancer using AmpliconArchitect*. Nature Communications, 2019. 10(1): p. 392.
