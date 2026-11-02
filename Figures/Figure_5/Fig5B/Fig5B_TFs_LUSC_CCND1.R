# CCND1 TCGA-77-A5G6 HiChIP (LUSC) locus visualisation with TF tracks
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

CGC <- readRDS('../../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
gencode_hg38 <- readRDS('../../supp_Figure_3/FigS3A/source_data/gencode_hg38.rds')

# hs1 -> hg38 liftOver chain (CCND1 amp / lindiff / methylation / v4C tracks are in hs1/T2T coordinates)
hs1ToHg38_chain <- import.chain(file.path(IN_DIR, 'hs1ToHg38.over.chain'))

xG <- 'CCND1'
Oncgr <- CGC %Q% (Gene_Symbol == xG)
window_oi <- Oncgr + 12e5

# CCND1 NSCLC amplification frequency + BRCA-vs-NSCLC significance tracks (hs1 -> hg38)
NSCLC_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, "/20260907_Fig.2C_NSCLC_CCND1_AmpFreq_track.rds"))
NSCLC_xG_amps_grc <- NSCLC_xGtad_per_amp_gt@data[[1]]
NSCLC_xGtad_amps_hg38gr <- unlist(liftOver(NSCLC_xG_amps_grc, hs1ToHg38_chain))
NSCLC_xGtad_per_amp_gt <- gTrack(data = NSCLC_xGtad_amps_hg38gr, y.field = "AmpFreq", name = "%amp", col = "blue", y0 = 0, y1 = 1)
NSCLC_xGtad_per_amp_gt$bars <- FALSE

BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- readRDS(paste0(IN_DIR, "/20260907c_Fig.2C_BRCAvsNSCLC_CCND1_chr11_sig_track_cons80.rds"))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr <- unlist(liftOver(BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt@data[[1]], hs1ToHg38_chain))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- gTrack(data = BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr, name = "nlgp", col = "blue", height = 0.5)

xD <- 0.4
dmrs_d4 <- fread(paste0(IN_DIR, "/20260512_CCND1_TM_dmrs_d.csv"))
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
comb_Pvals_v4C_CCND1_LUSC_gt <- gTrack(comb_Pvals_v4C_CCND1_LUSC_gr_hg38, y.field = 'nlg_comb_q', bars = TRUE, col = 'blue', sep.draw = FALSE)

# TF FE tracks (preprocessed RDS files from bdgcmp)
HCC2814_SOX2_FE <- readRDS(file.path(IN_DIR, 'HCC2814_SOX2_FE_gr.rds'))
HCC2814_SOX2_FE_gt <- gTrack(data = HCC2814_SOX2_FE, y.field = 'score', bars = TRUE, y1 = 28, name = "HCC2814_SOX2_FE")

HCC95_SOX2_FE <- readRDS(file.path(IN_DIR, 'HCC95_SOX2_FE_gr.rds'))
HCC95_SOX2_FE_gt <- gTrack(data = HCC95_SOX2_FE, y.field = 'score', bars = TRUE, y1 = 28, name = "HCC95_SOX2_FE")

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

# plot window (window_oiz)
# stripped to non-'chr' seqnames to match the straw HiChIP matrix
window_oi <- gr.nochr(GRanges("chr11:69000000-70000000"))
window_oiz <- gr.nochr(GRanges("chr11:69000000-70000000") - 2e5)

pdf(file = file.path(OUT_DIR, 'Fig5B_LUSC_CCND1.pdf'), paper = 'a4', height = 70)
  plot(c(gencode_hg38,
         fishres_sigs_LUSC_gt,
         dmrs_d4_NSCLCgt,dmrs_d5_NSCLCgt,dmrs_d6_NSCLCgt,
         comb_Pvals_v4C_CCND1_LUSC_gt,
         HCC2814_SOX2_peaks_gt,
         HCC2814_SOX2_FE_gt,
         HCC95_p63_peaks_gt,
         HCC95_p63_FE_gt,
         BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
         NSCLC_xGtad_per_amp_gt),
       window_oi, cex.label = 0.1)
  plot(c(gencode_hg38,
         fishres_sigs_LUSC_gt,
         dmrs_d4_NSCLCgt,dmrs_d5_NSCLCgt,dmrs_d6_NSCLCgt,
         comb_Pvals_v4C_CCND1_LUSC_gt,
         HCC2814_SOX2_peaks_gt,
         HCC2814_SOX2_FE_gt,
         HCC95_p63_peaks_gt,
         HCC95_p63_FE_gt,
         BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
         NSCLC_xGtad_per_amp_gt),
       window_oiz, cex.label = 0.1)
dev.off()
