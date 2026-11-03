# HLAMP-RE

Analysis code and source data for the HLAMP-RE manuscript. The code covers
long-read cancer genome reconstruction, consensus amplicons, differential
methylation and Hi-C around *CCND1*, *NSD3* and *EGFR*.

> **Status:** this repository is still being prepared for publication. Folder
> and file names may change.

## Layout

```
Figures/            one folder per figure, with plotting scripts and source data
common/data/        reference genome files (see common/data/refGenome/NOTES.md)
common/workflows/   upstream analysis pipeline, in run order
```

## Running the scripts

Run each script from its own folder. Files included in this repository are read
through relative paths. All other inputs are read from `IN_DIR`, and outputs are
written to `OUT_DIR`:

```bash
export IN_DIR=/path/to/input
export OUT_DIR=/path/to/output
Rscript Fig3A_diffmethyl_CCND1_T2T_vis.R
```

If these are not set, the scripts use `input/` and `output/` in the current
folder. Each script names the files it needs. Most sit directly in `IN_DIR`; these
are kept in subfolders:

Some inputs are outputs of earlier steps, so copy them to `IN_DIR` before running
the next step. For example, the Fig. S4 job writes the models in `CNmodels`, and
the Fig. S5B job writes the tables in `TADratio`.

Scripts named `*job-script*` or `*sbatch*` run on a SLURM cluster with Apptainer.
Set the paths at the top of each one, create a `logs/` folder, and submit it from
its own folder. They need a container with R and the packages the R script loads.

Patient-level and third-party input data are not part of this repository.

## Figures

### Main figures

| Figure | Folder | Scripts |
|---|---|---|
| Fig. 1 | `Figures/Figure_1/Fig1B` | [`Fig1B.R`](Figures/Figure_1/Fig1B/Fig1B.R) |
| Fig. 1 | `Figures/Figure_1/Fig1C` | [`Fig1C_01_run_DeepSomatic.sh`](Figures/Figure_1/Fig1C/Fig1C_01_run_DeepSomatic.sh)<br>[`Fig1C_02_vep_annotation.sh`](Figures/Figure_1/Fig1C/Fig1C_02_vep_annotation.sh)<br>[`Fig1C_03_make_cbioportal_files.py`](Figures/Figure_1/Fig1C/Fig1C_03_make_cbioportal_files.py)<br>[`Fig1C_jabba_cbio.R`](Figures/Figure_1/Fig1C/Fig1C_jabba_cbio.R) |
| Fig. 2 | `Figures/Figure_2/Fig2A` | [`Fig2A.R`](Figures/Figure_2/Fig2A/Fig2A.R) |
| Fig. 2 | `Figures/Figure_2/Fig2B` | [`Fig2B.R`](Figures/Figure_2/Fig2B/Fig2B.R) |
| Fig. 2 | `Figures/Figure_2/Fig2C_D` | [`Fig2C_D_NSD3_EGFR_genomeWide.R`](Figures/Figure_2/Fig2C_D/Fig2C_D_NSD3_EGFR_genomeWide.R) |
| Fig. 2 | `Figures/Figure_2/Fig2E` | [`Fig2E_BRCA_NSCLC_CCND1_LeftVsRight_FisherStats_cons80.R`](Figures/Figure_2/Fig2E/Fig2E_BRCA_NSCLC_CCND1_LeftVsRight_FisherStats_cons80.R) |
| Fig. 3 | `Figures/Figure_3/Fig3A` | [`Fig3A_diffmethyl_CCND1_T2T_vis.R`](Figures/Figure_3/Fig3A/Fig3A_diffmethyl_CCND1_T2T_vis.R) |
| Fig. 3 | `Figures/Figure_3/Fig3B` | [`Fig3B_CCND1_vis.R`](Figures/Figure_3/Fig3B/Fig3B_CCND1_vis.R) |
| Fig. 4 | `Figures/Figure_4` | [`Fig4A_vis_hicdcdiff_results_covCovariate.R`](Figures/Figure_4/Fig4A_vis_hicdcdiff_results_covCovariate.R)<br>[`Fig4B_TCGA_C8_A12M_HiChIP_BRCA_CCND1_loc.R`](Figures/Figure_4/Fig4B_TCGA_C8_A12M_HiChIP_BRCA_CCND1_loc.R)<br>[`Fig4C_TCGA-77-A5G6_HiChIP_LUSC_CCND1_loc.R`](Figures/Figure_4/Fig4C_TCGA-77-A5G6_HiChIP_LUSC_CCND1_loc.R) |
| Fig. 5 | `Figures/Figure_5/Fig5A` | [`Fig5A_5.r`](Figures/Figure_5/Fig5A/Fig5A_5.r)<br>[`Fig5A_5_sbatch.sh`](Figures/Figure_5/Fig5A/Fig5A_5_sbatch.sh) |
| Fig. 5 | `Figures/Figure_5/Fig5B` | [`Fig5B_TFs_LUSC_CCND1.R`](Figures/Figure_5/Fig5B/Fig5B_TFs_LUSC_CCND1.R) |
| Fig. 5 | `Figures/Figure_5/Fig5C` | [`Fig5C_CCND1_E1_MCF7_Reporter_boxplot_log2_paired_ttest_ymax0.R`](Figures/Figure_5/Fig5C/Fig5C_CCND1_E1_MCF7_Reporter_boxplot_log2_paired_ttest_ymax0.R) |
| Fig. 5 | `Figures/Figure_5/Fig5D` | [`Fig5D_CCND1_E1_EPLC_boxplot_log2_paired_ttest_ymax0.R`](Figures/Figure_5/Fig5D/Fig5D_CCND1_E1_EPLC_boxplot_log2_paired_ttest_ymax0.R) |
| Fig. 5 | `Figures/Figure_5/Fig5E` | [`Fig5E_ZR751_CCND1_CRISPRi_3rep_merged_curve_no_g1.R`](Figures/Figure_5/Fig5E/Fig5E_ZR751_CCND1_CRISPRi_3rep_merged_curve_no_g1.R) |
| Fig. 5 | `Figures/Figure_5/Fig5F` | [`Fig5F_ZR751_AUC_LMM_adjP.R`](Figures/Figure_5/Fig5F/Fig5F_ZR751_AUC_LMM_adjP.R) |
| Fig. 5 | `Figures/Figure_5/Fig5G` | [`Fig5G_EPLC272H_CCND1_CRISPRi_3rep_merged_curve_no_g1.R`](Figures/Figure_5/Fig5G/Fig5G_EPLC272H_CCND1_CRISPRi_3rep_merged_curve_no_g1.R) |
| Fig. 5 | `Figures/Figure_5/Fig5H` | [`Fig5H_EPLC_AUC_LMM_adjp.R`](Figures/Figure_5/Fig5H/Fig5H_EPLC_AUC_LMM_adjp.R) |

### Supplementary figures

| Figure | Folder | Scripts |
|---|---|---|
| Fig. S1 | `Figures/supp_Figure_1/FigS1B` | [`FigS1B_hmfpcawg_cDNA_prev.R`](Figures/supp_Figure_1/FigS1B/FigS1B_hmfpcawg_cDNA_prev.R) |
| Fig. S1 | `Figures/supp_Figure_1/FigS1C` | [`FigS1C.R`](Figures/supp_Figure_1/FigS1C/FigS1C.R) |
| Fig. S1 | `Figures/supp_Figure_1/FigS1D` | [`FigS1D_DESeq2_ampVSnonamp_CCND1_EGFR_FGFR1NSD3_BRCA_LUSC.R`](Figures/supp_Figure_1/FigS1D/FigS1D_DESeq2_ampVSnonamp_CCND1_EGFR_FGFR1NSD3_BRCA_LUSC.R) |
| Fig. S2 | `Figures/supp_Figure_2` | [`FigS2A_B_CCND1_NSD3_EGFR_genomeWide.R`](Figures/supp_Figure_2/FigS2A_B_CCND1_NSD3_EGFR_genomeWide.R)<br>[`FigS2C_D_NSD3_EGFR_vis_cons80.R`](Figures/supp_Figure_2/FigS2C_D_NSD3_EGFR_vis_cons80.R) |
| Fig. S3 | `Figures/supp_Figure_3` | [`FigS3C_CCND1_hg38_lindiff_vis_with_peaks_and_medianwigs.R`](Figures/supp_Figure_3/FigS3C_CCND1_hg38_lindiff_vis_with_peaks_and_medianwigs.R)<br>[`FigS3D_diffmethyl_NSD3_T2T_vis.R`](Figures/supp_Figure_3/FigS3D_diffmethyl_NSD3_T2T_vis.R)<br>[`FigS3E_H3K27acATAC_NSD3_BRCAvsNSCLC_with_peaks_and_medianwigs.R`](Figures/supp_Figure_3/FigS3E_H3K27acATAC_NSD3_BRCAvsNSCLC_with_peaks_and_medianwigs.R)<br>[`FigS3F_diffmethyl_EGFR_T2T_vis.R`](Figures/supp_Figure_3/FigS3F_diffmethyl_EGFR_T2T_vis.R)<br>[`FigS3G_H3K27acATAC_EGFR_BRCAvsNSCLC_with_peaks_and_medianwigs.R`](Figures/supp_Figure_3/FigS3G_H3K27acATAC_EGFR_BRCAvsNSCLC_with_peaks_and_medianwigs.R) |
| Fig. S3 | `Figures/supp_Figure_3/FigS3A` | [`FigS3A_01_plot_HCC1954_data.R`](Figures/supp_Figure_3/FigS3A/FigS3A_01_plot_HCC1954_data.R)<br>[`FigS3A_02_ROC_meth_H3K27acChIP_covLevels.R`](Figures/supp_Figure_3/FigS3A/FigS3A_02_ROC_meth_H3K27acChIP_covLevels.R)<br>[`FigS3A_03_plot_ROCanalysis.R`](Figures/supp_Figure_3/FigS3A/FigS3A_03_plot_ROCanalysis.R) |
| Fig. S3 | `Figures/supp_Figure_3/FigS3B` | [`FigS3B_01_create_meth_matrix.R`](Figures/supp_Figure_3/FigS3B/FigS3B_01_create_meth_matrix.R)<br>[`FigS3B_02_UMAP.R`](Figures/supp_Figure_3/FigS3B/FigS3B_02_UMAP.R) |
| Fig. S3 | `Figures/supp_Figure_3/atac_H3K27ac_stats` | [`00_peak_overlaps.sh`](Figures/supp_Figure_3/atac_H3K27ac_stats/00_peak_overlaps.sh)<br>[`01_peak_enrichment_fisher.R`](Figures/supp_Figure_3/atac_H3K27ac_stats/01_peak_enrichment_fisher.R)<br>[`02_peak_dmr_overlap.R`](Figures/supp_Figure_3/atac_H3K27ac_stats/02_peak_dmr_overlap.R)<br>[`gene_settings.R`](Figures/supp_Figure_3/atac_H3K27ac_stats/gene_settings.R) |
| Fig. S4 | `Figures/supp_Figure_4` | [`FigS4_CCND1_tumor_plots_hiccovlogisticnorm_DMR_sig_20Mb.R`](Figures/supp_Figure_4/FigS4_CCND1_tumor_plots_hiccovlogisticnorm_DMR_sig_20Mb.R)<br>[`FigS4_chr7_8_11_20Mb_HiCcov_logistic_CNmodel_x200_d200_tumor_gene_hic.R`](Figures/supp_Figure_4/FigS4_chr7_8_11_20Mb_HiCcov_logistic_CNmodel_x200_d200_tumor_gene_hic.R)<br>[`FigS4_job-script_CCND1_tumor_plots.sh`](Figures/supp_Figure_4/FigS4_job-script_CCND1_tumor_plots.sh)<br>[`FigS4_job-script_HiCcov_logistic_chr7_8_11_20Mb_x200_d200.sh`](Figures/supp_Figure_4/FigS4_job-script_HiCcov_logistic_chr7_8_11_20Mb_x200_d200.sh) |
| Fig. S5 | `Figures/supp_Figure_5/FigS5A` | [`FigS5A_CCND1_median_CNnorm_HiC_BRCA_LUSC.R`](Figures/supp_Figure_5/FigS5A/FigS5A_CCND1_median_CNnorm_HiC_BRCA_LUSC.R)<br>[`FigS5A_job-script_CCND1_median_CNnorm_HiC.sh`](Figures/supp_Figure_5/FigS5A/FigS5A_job-script_CCND1_median_CNnorm_HiC.sh) |
| Fig. S5 | `Figures/supp_Figure_5/FigS5B` | [`FigS5B_job-script_comb_TADratio_200sqrd_20Mb.sh`](Figures/supp_Figure_5/FigS5B/FigS5B_job-script_comb_TADratio_200sqrd_20Mb.sh)<br>[`FigS5B_job-script_tumors_TADratio_200sqrd_20Mb.sh`](Figures/supp_Figure_5/FigS5B/FigS5B_job-script_tumors_TADratio_200sqrd_20Mb.sh)<br>[`FigS5B_lungTAD_inout_ratio_normals20Mb_vs_tumors_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_200sqrd.R`](Figures/supp_Figure_5/FigS5B/FigS5B_lungTAD_inout_ratio_normals20Mb_vs_tumors_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_200sqrd.R)<br>[`FigS5B_lungTAD_inout_ratio_tumor_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_20Mb_200sqrd.R`](Figures/supp_Figure_5/FigS5B/FigS5B_lungTAD_inout_ratio_tumor_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_20Mb_200sqrd.R) |
| Fig. S5 | `Figures/supp_Figure_5/FigS5C` | [`FigS5C_job-script_normHiC_combBRCA_NSCLC_chr11_69M_70M.sh`](Figures/supp_Figure_5/FigS5C/FigS5C_job-script_normHiC_combBRCA_NSCLC_chr11_69M_70M.sh)<br>[`FigS5C_normHiC_BRCA_NSCLC_chr11_69M_70M.R`](Figures/supp_Figure_5/FigS5C/FigS5C_normHiC_BRCA_NSCLC_chr11_69M_70M.R) |
| Fig. S6 | `Figures/supp_Figure_6` | [`FigS6_NSD3_EGFR_tumor_plots_hiccovlogisticnorm_DMR_20Mb.R`](Figures/supp_Figure_6/FigS6_NSD3_EGFR_tumor_plots_hiccovlogisticnorm_DMR_20Mb.R)<br>[`FigS6_NSD3_EGFR_tumor_plots_hiccovlogisticnorm_DMR_20Mb_binfix.R`](Figures/supp_Figure_6/FigS6_NSD3_EGFR_tumor_plots_hiccovlogisticnorm_DMR_20Mb_binfix.R)<br>[`FigS6_job-script_NSD3_EGFR_tumor_plots.sh`](Figures/supp_Figure_6/FigS6_job-script_NSD3_EGFR_tumor_plots.sh)<br>[`FigS6_job-script_binfix_tumor_plots.sh`](Figures/supp_Figure_6/FigS6_job-script_binfix_tumor_plots.sh) |
| Fig. S7 | `Figures/supp_Figure_7` | [`FigS7A_B_NSD3_EGFR_median_CNnorm_HiC_BRCA_LUSC_binfix.R`](Figures/supp_Figure_7/FigS7A_B_NSD3_EGFR_median_CNnorm_HiC_BRCA_LUSC_binfix.R)<br>[`FigS7A_B_job-script_median_CNnorm_HiC.sh`](Figures/supp_Figure_7/FigS7A_B_job-script_median_CNnorm_HiC.sh)<br>[`FigS7C_job-script_binfix_TADratio_logistic.sh`](Figures/supp_Figure_7/FigS7C_job-script_binfix_TADratio_logistic.sh)<br>[`FigS7C_lungTAD_inout_ratio_tumor_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_20Mb_200sqrd_binfix.R`](Figures/supp_Figure_7/FigS7C_lungTAD_inout_ratio_tumor_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_20Mb_200sqrd_binfix.R)<br>[`FigS7E_job-script_normHiC_BRCA_NSCLC_ESCC_NSD3.sh`](Figures/supp_Figure_7/FigS7E_job-script_normHiC_BRCA_NSCLC_ESCC_NSD3.sh)<br>[`FigS7E_normHiC_BRCA_NSCLC_ESCC_NSD3.R`](Figures/supp_Figure_7/FigS7E_normHiC_BRCA_NSCLC_ESCC_NSD3.R)<br>[`FigS7F_job-script_normHiC_BRCA_NSCLC_ESCC_EGFR.sh`](Figures/supp_Figure_7/FigS7F_job-script_normHiC_BRCA_NSCLC_ESCC_EGFR.sh)<br>[`FigS7F_normHiC_BRCA_NSCLC_ESCC_EGFR.R`](Figures/supp_Figure_7/FigS7F_normHiC_BRCA_NSCLC_ESCC_EGFR.R) |
| Fig. S7 | `Figures/supp_Figure_7/FigS7D` | [`FigS7D_job-script_comb_TADratio_200sqrd_20Mb.sh`](Figures/supp_Figure_7/FigS7D/FigS7D_job-script_comb_TADratio_200sqrd_20Mb.sh)<br>[`FigS7D_job-script_tumors_TADratio_200sqrd_20Mb.sh`](Figures/supp_Figure_7/FigS7D/FigS7D_job-script_tumors_TADratio_200sqrd_20Mb.sh)<br>[`FigS7D_lungTAD_inout_ratio_normals20Mb_vs_tumors_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_200sqrd.R`](Figures/supp_Figure_7/FigS7D/FigS7D_lungTAD_inout_ratio_normals20Mb_vs_tumors_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_200sqrd.R)<br>[`FigS7D_lungTAD_inout_ratio_tumor_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_20Mb_200sqrd.R`](Figures/supp_Figure_7/FigS7D/FigS7D_lungTAD_inout_ratio_tumor_hiccovlogisticnorm_CNgt9_randwin_chr7_8_11_20Mb_200sqrd.R) |
| Fig. S8 | `Figures/supp_Figure_8` | [`FigS8A_vis_hicdcdiff_results_covCovariate_NSD3_counts_sigbars.R`](Figures/supp_Figure_8/FigS8A_vis_hicdcdiff_results_covCovariate_NSD3_counts_sigbars.R)<br>[`FigS8B_vis_hicdcdiff_results_covCovariate_EGFR_counts_sigbars.R`](Figures/supp_Figure_8/FigS8B_vis_hicdcdiff_results_covCovariate_EGFR_counts_sigbars.R) |
| Fig. S9 | `Figures/supp_Figure_9/FigS9C` | [`FigS9C_MFM223_CCND1_CRISPRi_3rep_merged_curve.R`](Figures/supp_Figure_9/FigS9C/FigS9C_MFM223_CCND1_CRISPRi_3rep_merged_curve.R) |
| Fig. S9 | `Figures/supp_Figure_9/FigS9D` | [`FigS9D_MFM223AUC_LMM_adjp.R`](Figures/supp_Figure_9/FigS9D/FigS9D_MFM223AUC_LMM_adjp.R) |

## Workflows

| Step | Folder | Tools |
|---|---|---|
| 0 | [`00_lpWGS`](common/workflows/00_lpWGS) | ichorCNA, nf-core/sarek with Control-FREEC, nf-core/circdna with AmpliconArchitect |
| 1 | [`01_alignment_and_modkit`](common/workflows/01_alignment_and_modkit) | epi2me wf-somatic-variation |
| 2 | [`02_SNV_SV_calling`](common/workflows/02_SNV_SV_calling) | epi2me wf-human-variation, panel-of-normals filtering |
| 3 | [`03_CNV_calling`](common/workflows/03_CNV_calling) | CNVkit |
| 4 | [`04_JaBbA_genome_reconstruction`](common/workflows/04_JaBbA_genome_reconstruction) | JaBbA |
| 5 | [`05_JaBbA_genome_walks`](common/workflows/05_JaBbA_genome_walks) | JaBbA walks |
| 6 | [`06_consensus_amplicon`](common/workflows/06_consensus_amplicon/HLAMP_walks_combined_hs1) | consensus amplicons and Fisher tests ([method notes](common/workflows/06_consensus_amplicon/HLAMP_walks_combined_hs1/NOTES.md)) |
| 7 | [`07_differential_methylation`](common/workflows/07_differential_methylation) | DSS |
| 8 | [`08_HiC`](common/workflows/08_HiC) | nf-core/hic, HiCDCPlus, virtual 4C tracks; [`linkprep`](common/workflows/08_HiC/linkprep): BWA, pairtools, Juicer tools |
| 9 | [`09_RNAseq`](common/workflows/09_RNAseq) | nf-core/rnaseq, DESeq2 |

Tool citations are in the `NOTES.md` file in each folder. Sample sheets are
provided as blank templates (`samplesheet_*_blank.csv`). Patient-level inputs
to these workflows are not part of this repository.

## License

The code in this repository is licensed under the
[GNU General Public License v3.0](LICENSE). Source data and third-party reference
files are not covered by this license and remain under their own terms.

Copyright (C) 2026 HLAMP-RE contributors

## Contributors

Code contributed by Petrov, YingXu11, tisop99 and FrankDubois.
