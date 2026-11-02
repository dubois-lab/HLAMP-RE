# Fig S9B: NSD3 locus visualisation with TF tracks
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
  try(library(data.table), silent = TRUE)
}))

CGC <- readRDS('../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
gencode_hg38 <- readRDS('../supp_Figure_3/FigS3A/source_data/gencode_hg38.rds')

# hs1 -> hg38 liftOver chain (NSD3 amp / diffmethyl / v4C tracks are in hs1/T2T coordinates)
hs1ToHg38_chain <- import.chain(file.path(IN_DIR, 'hs1ToHg38.over.chain'))

xG <- 'NSD3'
xT <- 'TCGA-77-A5G6'
Oncgr <- gr.chr(CGC %Q% (Gene_Symbol == xG))
window_oi <- Oncgr + 12e5

# amp/sig tracks (T2T/hs1 -> hg38)

# BRCA vs NSCLC NSD3 differential amp-frequency significance track (chr8)
BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt <- readRDS(paste0(IN_DIR, "/Fig_S3_20260624_BRCAvsNSCLC_NSD3_chr8_sig_track_cons80.rds"))
BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt_gr <- unlist(liftOver(BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt@data[[1]], hs1ToHg38_chain))
BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt <- gTrack(data = BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt_gr, name = "nlgp", col = "black", height = 0.5)

# NSCLC NSD3 left-vs-right (amp boundary) significance track
NSCLC_NSD3_LvsR_sig_LR_gt <- readRDS(paste0(IN_DIR, "/Fig_S3_20260624_NSCLC_NSD3_LvsR_sig_LR_cons80.rds"))
NSCLC_NSD3_LvsR_sig_LR_gt_gr <- unlist(liftOver(NSCLC_NSD3_LvsR_sig_LR_gt@data[[1]], hs1ToHg38_chain))
NSCLC_NSD3_LvsR_sig_LR_gt <- gTrack(data = NSCLC_NSD3_LvsR_sig_LR_gt_gr, name = "nlgp", col = "blue", y0 = 0, y1 = 1, height = 0.5)

# NSCLC NSD3 amp-frequency track
NSCLC_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, "/Fig_S3_20260624_NSCLC_NSD3_AmpFreq_track.rds"))
NSCLC_xG_amps_grc <- NSCLC_xGtad_per_amp_gt@data[[1]]
NSCLC_xGtad_amps_hg38gr <- unlist(liftOver(NSCLC_xG_amps_grc, hs1ToHg38_chain))
NSCLC_xGtad_per_amp_gt <- gTrack(data = NSCLC_xGtad_amps_hg38gr, y.field = "AmpFreq", name = "%amp", col = "blue", y0 = 0, y1 = 1)
NSCLC_xGtad_per_amp_gt$bars <- F

# Fisher & diffmethyl tracks
fishres_A2 <- fread(file.path(IN_DIR, "20250619_NSD3_H3K27ac_ATAC_linDiffs_stats.csv"))
fishres_A2[, seqnames := unlist(strsplit(x = peakID, split = "[:]"))[1], by = peakID]
fishres_A2[, start := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[2]), by = peakID]
fishres_A2[, end := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[3]), by = peakID]
fishres_A2gr <- dt2gr(fishres_A2)
fishres_A_sigs_BRCA_HR <- fishres_A2gr %Q% (q_val < 0.1 & tstdttype == "BRCA_HR+")
fishres_sigs_BRCA_HR_gt <- gTrack(data = reduce(fishres_A_sigs_BRCA_HR))
fishres_sigs_LUSCgt <- gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.1 & tstdttype == "LUSC")))
fishres_sigs_LUSC_gt <- fishres_sigs_LUSCgt

# diffmethyl
xD <- 0.3
dmrs_d <- fread(paste0(IN_DIR, "/20260615_NSD3_TM_dmrs_d.csv"))
colnames(dmrs_d)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs_d)), hs1ToHg38_chain))
dmrs_d <- as.data.table(dmrs_d_gr)
dmrs_d_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d[diff.Methy < -xD])))
dmrs_d_BRCAgt <- gTrack(dt2gr(dmrs_d[diff.Methy > xD]))
dmrs <- fread(paste0(IN_DIR, "/20260615_NSD3_TM_dmrs.csv"))
colnames(dmrs)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs)), hs1ToHg38_chain))
dmrs <- gr2dt(dmrs_gr)
dmrs3_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > 0.3]))
dmrs3_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -0.3]))

# NSD3 LUSC virtual-4C combined (Fisher) p-values (hs1 -> hg38)
comb_Pvals_v4C_NSD3_LUSC_gr <- readRDS(file.path(IN_DIR, '20260702_comb_Pvals_v4C_NSD3_LUSC_gr.rds'))
comb_Pvals_v4C_NSD3_LUSC_gr_hg38 <- grl.unlist(liftOver(comb_Pvals_v4C_NSD3_LUSC_gr, hs1ToHg38_chain))
comb_Pvals_v4C_NSD3_LUSC_gt <- gTrack(comb_Pvals_v4C_NSD3_LUSC_gr_hg38, y.field = 'nlg_comb_q', bars = TRUE, col = 'blue', sep.draw = FALSE)

# TF FE tracks (preprocessed RDS files from bdgcmp)
HCC2814_SOX2_FE <- readRDS(file.path(IN_DIR, 'HCC2814_SOX2_FE_gr.rds'))
HCC2814_SOX2_FE_gt <- gTrack(data = HCC2814_SOX2_FE, y.field = 'score', bars = TRUE, y1 = 28, name = "HCC2814_SOX2_FE")

HCC95_SOX2_FE <- readRDS(file.path(IN_DIR, 'HCC95_SOX2_FE_gr.rds'))
HCC95_SOX2_FE_gt <- gTrack(data = HCC95_SOX2_FE, y.field = 'score', bars = TRUE, y1 = 14, name = "HCC95_SOX2_FE")

HCC95_p63_FE <- readRDS(file.path(IN_DIR, 'HCC95_p63_FE_gr.rds'))
HCC95_p63_FE_gt <- gTrack(data = HCC95_p63_FE, y.field = 'score', bars = TRUE, y1 = 14, name = "HCC95_p63_FE")

# TF narrowPeak tracks
macs_dir <- file.path(IN_DIR, '20250630_CL_chipseqTF_macs')

HCC2814_SOX2_peaks <- fread(file.path(macs_dir, 'HCC2814_SOX2_peaks.narrowPeak'))
colnames(HCC2814_SOX2_peaks)[1:3] <- c("chr", "start", "end")
HCC2814_SOX2_peaks_gt <- gTrack(dt2gr(HCC2814_SOX2_peaks), height = 0.5, name = "HCC2814_SOX2_peaks")

HCC95_SOX2_peaks <- fread(file.path(macs_dir, 'HCC95_SOX2_peaks.narrowPeak'))
colnames(HCC95_SOX2_peaks)[1:3] <- c("chr", "start", "end")
HCC95_SOX2_peaks_gt <- gTrack(dt2gr(HCC95_SOX2_peaks), height = 0.5, name = "HCC95_SOX2_peaks")

HCC95_p63_peaks <- fread(file.path(macs_dir, 'HCC95_p63_peaks.narrowPeak'))
colnames(HCC95_p63_peaks)[1:3] <- c("chr", "start", "end")
HCC95_p63_peaks_gt <- gTrack(dt2gr(HCC95_p63_peaks), height = 0.5, name = "HCC95_p63_peaks")

# plot window (window_oi, window_oiz) for NSD3
window_oi <- Oncgr + 2e5
window_oiz <- Oncgr + 5e5

reporters_dt <- data.table(rbind(cbind('chr8', 38261122, 38261673), 
                                 cbind('chr8', 38562585, 38562946)))
colnames(reporters_dt) <- c('chr', 'start', 'end')
reporters_dt$start <- as.numeric(reporters_dt$start)
reporters_dt$end <- as.numeric(reporters_dt$end)
deadreporters_gt <- gTrack(dt2gr(reporters_dt[1,]))
active_reporters_gt <- gTrack(dt2gr(reporters_dt[2,]))


pdf(file = file.path(OUT_DIR, 'FigS9B_LUSC_NSD3.pdf'), paper = 'a4', height = 70)
  plot(c(gencode_hg38, deadreporters_gt, active_reporters_gt,
         fishres_sigs_LUSC_gt,
         dmrs_d_NSCLCgt,
         comb_Pvals_v4C_NSD3_LUSC_gt,
         HCC2814_SOX2_peaks_gt,
         HCC2814_SOX2_FE_gt,
         HCC95_p63_peaks_gt,
         HCC95_p63_FE_gt,
         BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt,
         NSCLC_xGtad_per_amp_gt),
       window_oi, cex.label = 0.1)
  plot(c(gencode_hg38, deadreporters_gt, active_reporters_gt,
         fishres_sigs_LUSC_gt,
         dmrs_d_NSCLCgt,
         comb_Pvals_v4C_NSD3_LUSC_gt,
         HCC2814_SOX2_peaks_gt,
         HCC2814_SOX2_FE_gt,
         HCC95_p63_peaks_gt,
         HCC95_p63_FE_gt,
         BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt,
         NSCLC_xGtad_per_amp_gt),
       window_oiz, cex.label = 0.1)
dev.off()
