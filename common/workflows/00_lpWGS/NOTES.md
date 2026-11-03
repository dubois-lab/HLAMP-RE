# Low-pass whole-genome sequencing

| Script | Step |
|---|---|
| `01_ichorCNA_snakemake.sh` | copy number segments per sample with ichorCNA, using `Snakefile` and `01_ichorCNA_runscript.R` |
| `02_sarek_controlfreec.sh` | copy number with Control-FREEC, from aligned BAM files |
| `03_ichorCNA_amplified_oncogenes.R` | amplified oncogenes in the ichorCNA segments, with IGV segment files per gene |
| `04_circDNA_ampliconarchitect.sh` | amplicon reconstruction with AmpliconArchitect, from aligned BAM files |
| `05_ichorCNA_gene_alterations.R` | copy number alteration calls per sample and gene, written as `lpWGS_genomic.tsv` for Fig. 1B and Fig. S1A |

`03_ichorCNA_amplified_oncogenes.R` reads these files from `IN_DIR`:

- `results/ichorCNA/<sample>/<sample>.seg` and `<sample>.seg.txt`
- `20240621_ichorCNA_IGV.seg` and `20240621_ichorCNA.seg.rds`, the segments of an earlier ichorCNA run

`05_ichorCNA_gene_alterations.R` reads `20240717_ichorCNA_cgc.seg.rds` from `IN_DIR`, which `03_ichorCNA_amplified_oncogenes.R` writes to `OUT_DIR`.

# ichorCNA (https://github.com/broadinstitute/ichorCNA)

Adalsteinsson, V.A. et al. *Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors*. Nature Communications, 2017. 8(1): p. 1324.

The `Snakefile` runs ichorCNA in `docker://frankdubois/ichorcna_optparse:v1`, a derivative of the `fredhutch/ichorcna` image.

# nf-core/sarek (https://nf-co.re/sarek)

Garcia, M. et al. *Sarek: A portable workflow for whole-genome sequencing analysis of germline and somatic variants*. F1000Research, 2020. 9: p. 63.

# Control-FREEC (https://github.com/BoevaLab/FREEC)

Boeva, V. et al. *Control-FREEC: a tool for assessing copy number and allelic content using next-generation sequencing data*. Bioinformatics, 2012. 28(3): p. 423-425.

# nf-core/circdna (https://nf-co.re/circdna)

https://doi.org/10.5281/zenodo.6685250

# AmpliconArchitect (https://github.com/AmpliconSuite/AmpliconArchitect)

Deshpande, V. et al. *Exploring the landscape of focal amplifications in cancer using AmpliconArchitect*. Nature Communications, 2019. 10(1): p. 392.
