# Input files and cutoffs for each amplified gene, read by the two scripts in this folder.
# File names are relative to IN_DIR. To analyze another gene, add a block with the same fields.

# Tumor types tested, and the amplification group whose consensus regions each one is restricted to
tumor_types <- c(LUAD = 'NSCLC', 'BRCA_HR-' = 'BRCA', 'BRCA_HR+' = 'BRCA', LUSC = 'NSCLC', HNSC = 'NSCLC', UCEC = 'BRCA')

gene_settings <- list(
  CCND1 = list(
    fisher_amp_tracks = c(NSCLC = 'Fig.2C_NSCLC_CCND1_AmpFreq_track.rds',
                          BRCA = 'Fig.2C_BRCA_CCND1_AmpFreq_track.rds'),
    peak_files = c(LUAD = 'LUAD_CCND1_withCLs_merged.bed',
                   'BRCA_HR-' = 'BRCA_HRneg_CCND1_withCLs_merged.bed',
                   'BRCA_HR+' = 'BRCA_HRpos_CCND1_withCLs_merged.bed',
                   LUSC = 'LUSC_CCND1_withCLs_merged.bed',
                   HNSC = 'HNSC_CCND1_merged.bed',
                   UCEC = 'UCEC_CCND1_merged.bed'),
    all_peaks_file = 'all_CCND1_withCLs_merged.bed',
    min_samples = c(LUAD = 5, 'BRCA_HR-' = 5, 'BRCA_HR+' = 5, LUSC = 3, HNSC = 3, UCEC = 3),
    # CCND1 has no settings for 02_peak_dmr_overlap.R
    stats_file = 'CCND1_H3K27ac&ATAC_linDiffs_stats.csv'
  ),
  NSD3 = list(
    # 01_peak_enrichment_fisher.R
    fisher_amp_tracks = c(NSCLC = 'Fig_S3_NSCLC_NSD3_AmpFreq_track_fisher.rds',
                          BRCA = 'Fig_S3_BRCA_NSD3_AmpFreq_track_fisher.rds'),
    peak_files = c(LUAD = 'LUAD_NSD3_withCLs_merged.bed',
                   'BRCA_HR-' = 'BRCA_HRneg_NSD3_withCLs_merged.bed',
                   'BRCA_HR+' = 'BRCA_HRpos_NSD3_withCLs_merged.bed',
                   LUSC = 'LUSC_NSD3_withCLs_merged.bed',
                   HNSC = 'HNSC_NSD3_merged.bed',
                   UCEC = 'UCEC_NSD3_merged.bed'),
    all_peaks_file = 'all_NSD3_withCLs_merged.bed',
    # a peak is tested if at least this many samples of the tumor type have it
    min_samples = c(LUAD = 5, 'BRCA_HR-' = 5, 'BRCA_HR+' = 5, LUSC = 3, HNSC = 3, UCEC = 3),
    # written by 01_peak_enrichment_fisher.R, read by 02_peak_dmr_overlap.R
    stats_file = 'NSD3_H3K27ac_ATAC_linDiffs_stats.csv',
    # 02_peak_dmr_overlap.R
    overlap_amp_tracks = c(NSCLC = 'Fig_S3_NSCLC_NSD3_AmpFreq_track_overlap.rds',
                           BRCA = 'Fig_S3_BRCA_NSD3_AmpFreq_track_overlap.rds'),
    dmr_file = 'NSD3_TM_dmrs.csv',
    min_diff_methyl = 0.3,
    q_cutoff = 0.1,
    preferential_amp_file = 'BRCAvsNSCLC_NSD3_chr8_FisherStats_CN10_cons80_gr.rds'
  ),
  EGFR = list(
    fisher_amp_tracks = c(NSCLC = 'Fig_S3_NSCLC_EGFR_AmpFreq_track_fisher.rds',
                          BRCA = 'Fig_S3_BRCA_EGFR_AmpFreq_track_fisher.rds'),
    peak_files = c(LUAD = 'LUAD_EGFR_withCLs_merged.bed',
                   'BRCA_HR-' = 'BRCA_HRneg_EGFR_withCLs_merged.bed',
                   'BRCA_HR+' = 'BRCA_HRpos_EGFR_withCLs_merged.bed',
                   LUSC = 'LUSC_EGFR_withCLs_merged.bed',
                   HNSC = 'HNSC_EGFR_merged.bed',
                   UCEC = 'UCEC_EGFR_merged.bed'),
    all_peaks_file = 'all_EGFR_withCLs_merged.bed',
    min_samples = c(LUAD = 5, 'BRCA_HR-' = 5, 'BRCA_HR+' = 4, LUSC = 3, HNSC = 3, UCEC = 3),
    stats_file = 'EGFR_H3K27ac_ATAC_linDiffs_stats.csv',
    overlap_amp_tracks = c(NSCLC = 'Fig_S3_NSCLC_EGFR_AmpFreq_track_overlap.rds',
                           BRCA = 'Fig_S3_BRCA_EGFR_AmpFreq_track_overlap.rds'),
    dmr_file = 'EGFR_TM_dmrs_d_delt_0.4.csv',
    min_diff_methyl = 0.4,
    q_cutoff = 0.15,
    preferential_amp_file = NULL
  )
)
