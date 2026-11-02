# CCND1 TCGA-77-A5G6 HiChIP (LUSC) locus visualisation
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Load packages tolerantly. A transitive dependency (biomaRt, pulled in via gGnome/khtools)
## runs an Ensembl connectivity check at load time that can fail. biomaRt is not used, so suppress the
## load-time messages/warnings and ignore any such load-time error.
suppressWarnings(suppressMessages({
  try(library(parallel), silent = TRUE)
  try(library(GxG), silent = TRUE)
  try(library(rtracklayer), silent = TRUE)
}))
CGC <- readRDS('../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
gencode_hg38 = readRDS('../supp_Figure_3/FigS3A/source_data/gencode_hg38.rds')

# hs1 -> hg38 liftOver chain (CCND1 amp / lindiff / methylation / v4C tracks are in hs1/T2T coordinates)
hs1ToHg38_chain <- import.chain(file.path(IN_DIR, 'hs1ToHg38.over.chain'))

xG = 'CCND1'
xT = 'TCGA-77-A5G6'
Oncgr <- CGC %Q%  (Gene_Symbol == xG)
window_oi = Oncgr + 12e5

# CCND1 NSCLC amplification frequency + left-vs-right and BRCA-vs-NSCLC significance tracks (hs1 -> hg38)
NSCLC_xGtad_per_amp_gt <- readRDS(file.path(IN_DIR, '20260907_Fig.2C_NSCLC_CCND1_AmpFreq_track.rds'))
NSCLC_xGtad_per_amp_gt_gr <- unlist(liftOver(NSCLC_xGtad_per_amp_gt@data[[1]], hs1ToHg38_chain))
NSCLC_xGtad_per_amp_gt <- gTrack(data = NSCLC_xGtad_per_amp_gt_gr, y.field = "AmpFreq", name = "%amp", col = "blue")
NSCLC_xGtad_per_amp_gt$bars <- F


BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- readRDS(file.path(IN_DIR, '20260907c_Fig.2C_BRCAvsNSCLC_CCND1_chr11_sig_track_cons80.rds'))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr <- unlist(liftOver(BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt@data[[1]], hs1ToHg38_chain))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- gTrack(data = BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr, name = "nlgp", col = "blue", height = 0.5)

# read in HiChIP data from https://gdc.cancer.gov/about-data/publications/TCGA-HiChIP-2024 (TCGA-77-A5G6 sample)
TCGA_xT_HiChiP = straw(file.path(IN_DIR, 'LUSC-5AC6ED60-B4C0-42C1-BAB1-A173859818AE-X019-S04-B1-T1_H3K27ac.allValidPairs.hic'),  gr = window_oi, res = 5e3)

# TCGA-77-A5G6 H3K27ac 1D signal (normalized) bigWig
TCGA_xT_HiChIP_1D_signal_bw_norm_counts_gt <- gTrack(data = file.path(IN_DIR, 'LUSC-5AC6ED60-B4C0-42C1-BAB1-A173859818AE-X019-S04-B1-T1_H3K27ac_1D-signal.norm.bw'), bars = T)

# TCGA-77-A5G6 WGS ASCAT copy-number segments
TCGA_xT_seg_dt <- fread(file.path(IN_DIR, 'afbf285e-0560-4e1d-9c6e-9741368aa9f5.wgs.ASCAT.copy_number_variation.seg.txt'))
TCGA_xT_seg_gt <- gTrack(data = dt2gr(TCGA_xT_seg_dt), y.field = 'Copy_Number')

# make 4C track from the HiChIP matrix at the CCND1 TSS
# (hg38; 'chromosome_name' column is non-'chr' -> matches the straw HiChIP gr directly)
CCND1_TSS <- fread(file.path(IN_DIR, 'CCND1_TSS_hg38.csv'))
colnames(CCND1_TSS)[4] = 'seqnames'
colnames(CCND1_TSS)[8] = 'start'
CCND1_TSS[, end := start + 1]
CCND1_TSS[, strand := '*']
CCND1_TSSgr_hg38 <- dt2gr(CCND1_TSS)
TCGA_xT_HiChiP_4C <- TCGA_xT_HiChiP[TCGA_xT_HiChiP$gr %^% CCND1_TSSgr_hg38, drop = T]
TCGA_xT_HiChiP_4C_gt <- gTrack(data = TCGA_xT_HiChiP_4C, y.field = 'value', bars = T, col = 'red', y1 = 30)

# our tumor methylation data: DMRs, NSCLC direction (hs1 -> hg38)
xD = 0.4
dmrs_d4 <- readRDS(file.path(IN_DIR, '20260512_CCND1_TM_dmrs_d.rds'))
colnames(dmrs_d4)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d4_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs_d4)), hs1ToHg38_chain))
dmrs_d4 <- as.data.table(dmrs_d4_gr)
dmrs_d4_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -xD])))
dmrs_d5_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -0.5])))
dmrs_d6_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -0.6])))


# significant lineage-differential H3K27ac & ATAC peaks (LUSC)
ATACfisher <- fread(file.path(IN_DIR, '20250614_CCND1_H3K27ac&ATAC_linDiffs_stats.csv'))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = "[:]"))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[3]), by = peakID]
fishres_A2gr <- dt2gr(ATACfisher)
fishres_A_sigs_LUSC <- fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "LUSC")
fishres_sigs_LUSC_gt <- gTrack(data = reduce(fishres_A_sigs_LUSC))

# CCND1 LUSC virtual-4C combined (Fisher) p-values (hs1 -> hg38)
comb_Pvals_v4C_CCND1_LUSC_gr <- readRDS(file.path(IN_DIR, '20260702_comb_Pvals_v4C_CCND1_LUSC_gr.rds'))
comb_Pvals_v4C_CCND1_LUSC_gr_hg38 <- grl.unlist(liftOver(comb_Pvals_v4C_CCND1_LUSC_gr, hs1ToHg38_chain))
comb_Pvals_v4C_CCND1_LUSC_gt <- gTrack(comb_Pvals_v4C_CCND1_LUSC_gr_hg38, y.field = 'nlg_comb_q', bars = T, col = 'blue', sep.draw =F)

# plot window (window_oiz)
# stripped to non-'chr' seqnames to match the straw HiChIP matrix
window_oi <- gr.nochr(GRanges("chr11:69000000-70000000"))
window_oiz <- gr.nochr(GRanges("chr11:69000000-70000000") - 2e5)

pdf(file = paste0(OUT_DIR, '/Fig.4C_20260908_', xT, '_HiChIP_LUSC_CCND1.pdf')
    ,paper = 'a4', height = 70)
  plot(c(gencode_hg38,
         fishres_sigs_LUSC_gt,
         dmrs_d4_NSCLCgt,
         NSCLC_xGtad_per_amp_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
         comb_Pvals_v4C_CCND1_LUSC_gt,
         TCGA_xT_seg_gt,
         TCGA_xT_HiChIP_1D_signal_bw_norm_counts_gt,
         TCGA_xT_HiChiP_4C_gt,
         TCGA_xT_HiChiP$gtrack(name = paste0(xT, "_HiChiP"), height = 90, clim = c(0, 30))),
       window_oi, cex.label = 0.1)
  plot(c(gencode_hg38,
         fishres_sigs_LUSC_gt,
         dmrs_d4_NSCLCgt,dmrs_d5_NSCLCgt,dmrs_d6_NSCLCgt,
         NSCLC_xGtad_per_amp_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
         comb_Pvals_v4C_CCND1_LUSC_gt,
         TCGA_xT_seg_gt,
         TCGA_xT_HiChIP_1D_signal_bw_norm_counts_gt,
         TCGA_xT_HiChiP_4C_gt,
         TCGA_xT_HiChiP$gtrack(name = paste0(xT, "_HiChiP"), height = 90, clim = c(0, 30))),
       window_oiz, cex.label = 0.1)
  plot(c(gencode_hg38,
         fishres_sigs_LUSC_gt,
         dmrs_d4_NSCLCgt,dmrs_d5_NSCLCgt,dmrs_d6_NSCLCgt,
         NSCLC_xGtad_per_amp_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
         comb_Pvals_v4C_CCND1_LUSC_gt,
         TCGA_xT_seg_gt,
         TCGA_xT_HiChIP_1D_signal_bw_norm_counts_gt,
         TCGA_xT_HiChiP_4C_gt,
         TCGA_xT_HiChiP$gtrack(name = paste0(xT, "_HiChiP"), height = 90, clim = c(0, 30))),
       window_oi, cex.label = 0.1)
dev.off()
