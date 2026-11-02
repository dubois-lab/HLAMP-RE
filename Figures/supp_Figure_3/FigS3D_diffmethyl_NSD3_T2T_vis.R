# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(gGnome)
library(gTrack)
library(rtracklayer)

# NSD3 version of the Fig. 3A script (red/blue scheme, median + BRCA-LUSC median difference tracks, shade-coded DMR tracks)

# same formatting as passed to the gTrack() calls below; applied to the tracks loaded from rds and to all methylation tracks
fmt_gt <- function(gt) {
    gt$draw.paths = FALSE
    gt$sep.draw = FALSE
    gt$sep.bg.col = 'white'
    gt$sep.lwd = 0
    gt$stack.gap = 0
    gt$y.grid.col = "white"
    gt
}
# recolor tracks loaded from rds
col_gt <- function(gt, col) {
    gt$col = col
    gt$border = col
    gt
}

# BRCA tracks in shades of red, LUSC/NSCLC tracks in shades of blue (light -> dark); BRCA-LUSC median difference stays black
BRCA_cols = c(light = '#FCBBA1', mid = '#FB6A4A', dark = '#CB181D', darkest = '#67000D')
LUSC_cols = c(light = '#C6DBEF', mid = '#6BAED6', dark = '#2171B5', darkest = '#08306B')

CGC_gr_T2T <- readRDS(file.path(IN_DIR, 'CGC_gr_T2T.rds'))
gencode_T2T_gt <- fmt_gt(readRDS(file.path(IN_DIR, 'gencodeTrack_t2t.rds')))

conversion_chain <- import.chain(file.path(IN_DIR, 'hg38ToHs1.over.chain'))

xG = 'NSD3'
Oncgr <- CGC_gr_T2T %Q%  (Gene_Symbol == xG)
out_dir = paste0(OUT_DIR, '/')
BRCA_smpl_col = BRCA_cols[['mid']] # individual sample methylation tracks; the BRCA/LUSC median tracks use the darkest shade
LUSC_smpl_col = LUSC_cols[['mid']]

# individual sample methylation track, same CpG filter and formatting as in the CCND1 script (Fig. 3A)
smpl_gt <- function(bm, col) {
    fmt_gt(gTrack(dt2gr(bm[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = col, border = col, stack.gap = 0))
}

# # amp dists ####
BRCA_NSD3_CNMatm <- readRDS(file.path(IN_DIR, '20260910_BRCA_NSD3_CNMatm_combined_hs1.rds'))
NSCLC_NSD3_CNMatm <- readRDS(file.path(IN_DIR, '20260910_NSCLC_LUSC_NSD3_CNMatm_combined_hs1.rds'))
BRCA_xGtad_per_amp_gt <- fmt_gt(gTrack(dt2gr(BRCA_NSD3_CNMatm), y.field = 'AmpPercent', bars = F, y0=0, y1=1, col = BRCA_cols[['dark']], border = BRCA_cols[['dark']]))
NSCLC_xGtad_per_amp_gt <- fmt_gt(gTrack(dt2gr(NSCLC_NSD3_CNMatm), y.field = 'AmpPercent', bars = F, y0=0, y1=1, col = LUSC_cols[['dark']], border = LUSC_cols[['dark']]))
window_oi <- reduce(dt2gr(rbind(BRCA_NSD3_CNMatm[AmpPercent > 0.8, 1:10], NSCLC_NSD3_CNMatm[AmpPercent > 0.8, 1:10])))
window_oi2 = window_oi + 1e5
window_oi3 = window_oi - 1e5

BRCAvsNSCLC_NSD3_fisher_gr <- readRDS(file.path(IN_DIR, '20260910_BRCAvsNSCLC_NSD3_chr8_FisherStats_CN10_cons80_gr.rds'))
BRCAvsNSCLC_NSD3_sig <- BRCAvsNSCLC_NSD3_fisher_gr[BRCAvsNSCLC_NSD3_fisher_gr$Fisher_pval < 0.05]
mcols(BRCAvsNSCLC_NSD3_sig)$tissue <- ifelse(mcols(BRCAvsNSCLC_NSD3_sig)$BRCA_AmpFreq > mcols(BRCAvsNSCLC_NSD3_sig)$NSCLC_AmpFreq, "BRCA", "NSCLC")
colormap_BRCAvsNSCLC_NSD3 <- list(tissue = c(NSCLC = LUSC_cols[['dark']], BRCA = BRCA_cols[['dark']]))
BRCAvsNSCLC_NSD3_sig_track_cons80_gt <- fmt_gt(gTrack(BRCAvsNSCLC_NSD3_sig, gr.colorfield = "color", colormaps = colormap_BRCAvsNSCLC_NSD3, height = 0.5, stack.gap = 0))

# no BRCA LeftVsRight NSD3 data available -> only the NSCLC LvsR track (placed next to the NSCLC amp track, mirroring BRCA LvsR in the CCND1 figure)
NSCLC_NSD3_LvsR_fisher <- fread(file.path(IN_DIR, '20260910_NSCLC_NSD3_LeftVsRight_FisherStats_CN10_cons80.txt'))
NSCLC_NSD3_LvsR_fisher_gr <- dt2gr(NSCLC_NSD3_LvsR_fisher[, ':='(seqnames = fifelse(right_AmpFreq > left_AmpFreq, right_chr, left_chr),
                                             start = fifelse(right_AmpFreq > left_AmpFreq, right_start, left_start),
                                             end = fifelse(right_AmpFreq > left_AmpFreq, right_end, left_end)
)])
NSCLC_NSD3_LvsR_sig <- NSCLC_NSD3_LvsR_fisher_gr[NSCLC_NSD3_LvsR_fisher_gr$Fisher_pval < 0.05]
NSCLC_NSD3_LvsR_sig_LR_cons80_gt <- fmt_gt(gTrack(NSCLC_NSD3_LvsR_sig, col = LUSC_cols[['dark']], border = LUSC_cols[['dark']], height = 0.5, stack.gap = 0))

# median bedmethyl, median window only ####
# median window: NSD3 +-1Mb (as for CCND1), extended to cover the widest plot window (window_oi2)
# per sample: same CpG filter as the sample tracks, restricted to the median window;
# percent_modified is recomputed from neglog_mth (= -log10(percent_modified/100)), which every bedmethyl rds has
# seqnames style differs between bedmethyl sources ('chr8' or '8'): match both, and write win_chr ('chr8') for all samples so CpGs line up in median_bm()
win_chr = as.character(seqnames(Oncgr))[1]
win_chr_alias = unique(c(win_chr, sub('^chr', '', win_chr), paste0('chr', sub('^chr', '', win_chr))))
win_start = min(start(Oncgr + 10e5), start(window_oi2))
win_end = max(end(Oncgr + 10e5), end(window_oi2))
win_bm <- function(bm) bm[modified_base_code == 'm' & Nvalid_cov > 5 & seqnames %in% win_chr_alias & start >= win_start & end <= win_end,
                          .(seqnames = win_chr, start, end, Nvalid_cov, percent_modified = 100 * 10^(-neglog_mth))]
# median per CpG over the samples covering it (n_samples); neglog_mth recomputed from the median percent_modified
median_bm <- function(win_bms) {
    n_cpg <- sapply(win_bms, nrow)
    print(n_cpg)
    if (any(n_cpg == 0)) stop('no CpGs in the median window for: ', paste(names(n_cpg)[n_cpg == 0], collapse = ', '))
    med <- rbindlist(win_bms)[, .(Nvalid_cov = median(as.numeric(Nvalid_cov)), percent_modified = median(percent_modified), n_samples = .N),
                              by = .(seqnames, start, end)]
    med[, ':='(modified_base_code = 'm', neglog_mth = -log10(percent_modified/100))]
    setorder(med, seqnames, start)
    # CpGs per number of contributing samples; most CpGs should be covered by (nearly) all samples
    print(med[, .N, keyby = n_samples])
    med
}
BRCA_win_bm <- list()
LUSC_win_bm <- list()

#BRCA ####
#Cb5####
Cb5.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'Cb5wf_mods.ungrouped.bedmethyl.doradoSup.rds'))
Cb5.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(Cb5.wf_mods.ungrouped.bedmethyl, BRCA_smpl_col)
BRCA_win_bm$Cb5 <- win_bm(Cb5.wf_mods.ungrouped.bedmethyl)
rm(Cb5.wf_mods.ungrouped.bedmethyl)
#C17####
C17.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C17.wf_mods.ungrouped.bedmethyl.rds'))
C17.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(C17.wf_mods.ungrouped.bedmethyl, BRCA_smpl_col)
BRCA_win_bm$C17 <- win_bm(C17.wf_mods.ungrouped.bedmethyl)
rm(C17.wf_mods.ungrouped.bedmethyl)
#Cb13####
Cb13.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'Cb13wf_mods.ungrouped.bedmethyl.doradoSup.rds'))
Cb13.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(Cb13.wf_mods.ungrouped.bedmethyl, BRCA_smpl_col)
BRCA_win_bm$Cb13 <- win_bm(Cb13.wf_mods.ungrouped.bedmethyl)
rm(Cb13.wf_mods.ungrouped.bedmethyl)
#C132####
C132.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C132wf_mods.ungrouped.bedmethyl.doradoSup.rds'))
C132.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(C132.wf_mods.ungrouped.bedmethyl, BRCA_smpl_col)
BRCA_win_bm$C132 <- win_bm(C132.wf_mods.ungrouped.bedmethyl)
rm(C132.wf_mods.ungrouped.bedmethyl)
#C95####
C95.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C95wf_mods.ungrouped.bedmethyl.doradoSup.rds'))
C95.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(C95.wf_mods.ungrouped.bedmethyl, BRCA_smpl_col)
BRCA_win_bm$C95 <- win_bm(C95.wf_mods.ungrouped.bedmethyl)
rm(C95.wf_mods.ungrouped.bedmethyl)
#C30####
C30_bedmethyl_ungrouped <- readRDS(file.path(IN_DIR, 'C30wf_mods.ungrouped.bedmethyl.doradoSup.rds'))
C30.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(C30_bedmethyl_ungrouped, BRCA_smpl_col)
BRCA_win_bm$C30 <- win_bm(C30_bedmethyl_ungrouped)
rm(C30_bedmethyl_ungrouped)
#C12####
C12wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C12wf_mods.ungrouped.bedmethyl.doradoSup.rds'))
C12wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(C12wf_mods.ungrouped.bedmethyl, BRCA_smpl_col)
BRCA_win_bm$C12 <- win_bm(C12wf_mods.ungrouped.bedmethyl)
rm(C12wf_mods.ungrouped.bedmethyl)
#C81####
C81_bedmethyl_ungrouped <- readRDS(file.path(IN_DIR, 'C81wf_mods.ungrouped.bedmethyl.regularDoradoHac.rds'))
C81.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(C81_bedmethyl_ungrouped, BRCA_smpl_col)
BRCA_win_bm$C81 <- win_bm(C81_bedmethyl_ungrouped)
rm(C81_bedmethyl_ungrouped)
#BRCA median####
print('read in all BRCA')
saveRDS(BRCA_win_bm, paste0(out_dir, '20261019_BRCA_NSD3_win_bm.rds'))
BRCA_median_bedmethyl <- median_bm(BRCA_win_bm)
rm(BRCA_win_bm)
saveRDS(BRCA_median_bedmethyl, paste0(out_dir, '20261019_BRCA_NSD3_window_oi_median_bedmethyl.rds'))
BRCA_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(BRCA_median_bedmethyl), y.field = 'neglog_mth', name = 'BRCA_median_-log5mC_%',
    bars = T, y1= 2.5, col = BRCA_cols[['darkest']], border = BRCA_cols[['darkest']], stack.gap = 0, height =20))
print('created tracks for all BRCA')
#NSCLC ####
#1E2AVA####
H1E2AVA.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, '1E2AVAwf_mods.ungrouped.bedmethyl.doradoSup.rds'))
H1E2AVA.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(H1E2AVA.wf_mods.ungrouped.bedmethyl, LUSC_smpl_col)
LUSC_win_bm$H1E2AVA <- win_bm(H1E2AVA.wf_mods.ungrouped.bedmethyl)
rm(H1E2AVA.wf_mods.ungrouped.bedmethyl)
#1IDEVD####
H1IDEVD.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, '1IDEVDwf_mods.ungrouped.bedmethyl.doradoSup.rds'))
H1IDEVD.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(H1IDEVD.wf_mods.ungrouped.bedmethyl, LUSC_smpl_col)
LUSC_win_bm$H1IDEVD <- win_bm(H1IDEVD.wf_mods.ungrouped.bedmethyl)
rm(H1IDEVD.wf_mods.ungrouped.bedmethyl)
#K0HCJ3####
HK0HCJ3.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'K0HCJ3wf_mods.ungrouped.bedmethyl.doradoSup.rds'))
HK0HCJ3.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(HK0HCJ3.wf_mods.ungrouped.bedmethyl, LUSC_smpl_col)
LUSC_win_bm$HK0HCJ3 <- win_bm(HK0HCJ3.wf_mods.ungrouped.bedmethyl)
rm(HK0HCJ3.wf_mods.ungrouped.bedmethyl)
#1AEWDG####
H1AEWDG.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, '1AEWDGwf_mods.ungrouped.bedmethyl.regularDoradoHac.rds'))
H1AEWDG.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(H1AEWDG.wf_mods.ungrouped.bedmethyl, LUSC_smpl_col)
LUSC_win_bm$H1AEWDG <- win_bm(H1AEWDG.wf_mods.ungrouped.bedmethyl)
rm(H1AEWDG.wf_mods.ungrouped.bedmethyl)
#PR84W0####
HPR84W0.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'PR84W0wf_mods.ungrouped.bedmethyl.doradoSup.rds'))
HPR84W0.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(HPR84W0.wf_mods.ungrouped.bedmethyl, LUSC_smpl_col)
LUSC_win_bm$HPR84W0 <- win_bm(HPR84W0.wf_mods.ungrouped.bedmethyl)
rm(HPR84W0.wf_mods.ungrouped.bedmethyl)
#H1MZ9S4####
H1MZ9S4wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, '1MZ9S4wf_mods.ungrouped.bedmethyl.regularDoradoHac.rds'))
H1MZ9S4wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(H1MZ9S4wf_mods.ungrouped.bedmethyl, LUSC_smpl_col)
LUSC_win_bm$H1MZ9S4 <- win_bm(H1MZ9S4wf_mods.ungrouped.bedmethyl)
rm(H1MZ9S4wf_mods.ungrouped.bedmethyl)
#H12093N####
H12093N_bedmethyl_ungrouped <- readRDS(file.path(IN_DIR, '12093Nwf_mods.ungrouped.bedmethyl.doradoSup.rds'))
H12093N_bedmethyl_ungrouped_logmCe3_gt <- smpl_gt(H12093N_bedmethyl_ungrouped, LUSC_smpl_col)
LUSC_win_bm$H12093N <- win_bm(H12093N_bedmethyl_ungrouped)
rm(H12093N_bedmethyl_ungrouped)
#HNODTWC####
HNODTWC_bedmethyl_ungrouped <- readRDS(file.path(IN_DIR, 'NODTWCwf_mods.ungrouped.bedmethyl.doradoSup.rds'))
HNODTWC_bedmethyl_ungrouped_logmCe3_gt <- smpl_gt(HNODTWC_bedmethyl_ungrouped, LUSC_smpl_col)
LUSC_win_bm$HNODTWC <- win_bm(HNODTWC_bedmethyl_ungrouped)
rm(HNODTWC_bedmethyl_ungrouped)
#LUSC median####
print('read in all LUSC')
saveRDS(LUSC_win_bm, paste0(out_dir, '20261019_LUSC_NSD3_win_bm.rds'))
LUSC_median_bedmethyl <- median_bm(LUSC_win_bm)
rm(LUSC_win_bm)
saveRDS(LUSC_median_bedmethyl, paste0(out_dir, '20261019_LUSC_NSD3_window_oi_median_bedmethyl.rds'))
LUSC_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(LUSC_median_bedmethyl), y.field = 'neglog_mth', name = 'LUSC_median_-log5mC_%', bars = T, y1= 2.5,
    col = LUSC_cols[['darkest']], border = LUSC_cols[['darkest']], stack.gap = 0, height =20))
print('created tracks for all LUSC')
#BRCA - LUSC median####
# per CpG covered by both median tracks: BRCA minus LUSC median value (diff_neglog_mth = plotted track value; diff_percent_modified kept for reference);
# neglog_mth is capped at nl_cap before subtracting, so CpGs with median 0% (neglog_mth = Inf) give finite differences
nl_cap = 2
BRCA_minus_LUSC_median_bedmethyl <- as.data.table(merge(
    BRCA_median_bedmethyl[, .(seqnames, start, end, BRCA_percent_modified = percent_modified, BRCA_neglog_mth = neglog_mth, BRCA_n_samples = n_samples)],
    LUSC_median_bedmethyl[, .(seqnames, start, end, LUSC_percent_modified = percent_modified, LUSC_neglog_mth = neglog_mth, LUSC_n_samples = n_samples)],
    by = c('seqnames', 'start', 'end'))) # merge() returns a plain data.frame here, so := below fails without as.data.table()
BRCA_minus_LUSC_median_bedmethyl[, ':='(diff_percent_modified = BRCA_percent_modified - LUSC_percent_modified,
                                        diff_neglog_mth = pmin(BRCA_neglog_mth, nl_cap) - pmin(LUSC_neglog_mth, nl_cap))]
setorder(BRCA_minus_LUSC_median_bedmethyl, seqnames, start)
saveRDS(BRCA_minus_LUSC_median_bedmethyl, paste0(out_dir, '20261019_BRCA_minus_LUSC_NSD3_window_oi_median_bedmethyl.rds'))
logBRCA_minus_LUSC_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(BRCA_minus_LUSC_median_bedmethyl), y.field = 'diff_neglog_mth', name = 'log_BRCA-LUSC_median_-log5mC_%',
    bars = T, y0 = -nl_cap, y1 = nl_cap, col = 'black', border = 'black', stack.gap = 0, height =20))
BRCA_minus_LUSC_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(BRCA_minus_LUSC_median_bedmethyl), y.field = 'diff_percent_modified', name = 'BRCA-LUSC_median_-log5mC_%',
    bars = T, y0 = -70, y1 = 70, col = 'black', border = 'black', stack.gap = 0, height =20))


#functional val####
GERCARres <- as.data.table(read.csv2(file.path(IN_DIR, 'GERCARres.csv')))
colnames(GERCARres)
GERCARres_gr <- dt2gr(GERCARres[chromosome == unique(as.character(seqnames(Oncgr)))])
GERCARres_gr_t2tgr <- grl.unlist(liftOver(gr.chr(GERCARres_gr), conversion_chain))
GERCARres_prolifd21_gt <- gTrack(GERCARres_gr_t2tgr, y.field = 'logFC_prolif_21', bars = T, col = BRCA_cols[['mid']], border = BRCA_cols[['mid']], stack.gap = 0)
GERCARres_gr_sig <- dt2gr(gr2dt(GERCARres_gr_t2tgr)[FDR_prolif_21 < 0.1 ])
GERCARres_prolifd21_q0.1_gt <- gTrack(GERCARres_gr_sig, col = BRCA_cols[['dark']], border = BRCA_cols[['dark']], stack.gap = 0)

SID_MCF7_E2 <- as.data.table(read.csv2(file.path(IN_DIR, 'cd-23-1157_supplementary_table_s3_suppst3_MCF7_E2.csv')))
SID_MCF7_E2_gr <- dt2gr(SID_MCF7_E2)
SID_MCF7_E2_prolifd21_gt <- gTrack(SID_MCF7_E2_gr, y.field = 'log2FC_day21_vs_day07', bars = TRUE, col = BRCA_cols[['mid']], border = BRCA_cols[['mid']], stack.gap = 0)
SID_MCF7_E2_sig <- SID_MCF7_E2[FDR_day21_vs_day07 < 0.1]
SID_MCF7_E2_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_MCF7_E2_sig), col = BRCA_cols[['dark']], border = BRCA_cols[['dark']], stack.gap = 0)

SID_T47D_E2 <- as.data.table(read.csv2(file.path(IN_DIR, 'cd-23-1157_supplementary_table_s4_suppst4_T47D_E2.csv')))
SID_T47D_E2_gr <- dt2gr(SID_T47D_E2)
SID_T47D_E2_prolifd21_gt <- gTrack(SID_T47D_E2_gr, y.field = 'log2FC_day21_vs_lib', bars = TRUE, col = BRCA_cols[['mid']], border = BRCA_cols[['mid']], stack.gap = 0)
SID_T47D_E2_sig <- SID_T47D_E2[FDR_day21_vs_lib < 0.1]
SID_T47D_E2_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_T47D_E2_sig), col = BRCA_cols[['dark']], border = BRCA_cols[['dark']], stack.gap = 0)

#####stats#####

fishres_A2 <- fread(file.path(IN_DIR, '20250619_NSD3_H3K27ac_ATAC_linDiffs_stats.csv'))
fishres_A2[, seqnames := unlist(strsplit(x = peakID, split = '[:]'))[1], by = peakID]
fishres_A2[, start := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[2]), by = peakID]
fishres_A2[, end := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[3]), by = peakID]
fishres_A2_gr<- dt2gr(fishres_A2)
ATACfisher_t2tgr <- grl.unlist(liftOver(gr.chr(fishres_A2_gr), conversion_chain))

ATACfisher_HRpBC_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.1 & tstdttype == 'BRCA_HR+')), col = BRCA_cols[['dark']], border = BRCA_cols[['dark']], stack.gap = 0)
ATACfisher_LUSC_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.1 & tstdttype == 'LUSC')), col = LUSC_cols[['dark']], border = LUSC_cols[['dark']], stack.gap = 0)
xD = 0.4

dmls <- fread(paste0(IN_DIR, '/20260615_NSD3_TM_dmls.csv'))
dmls[, end := pos +1, by =pos]
dmls_NSCLC_gt <- gTrack(reduce(dt2gr(dmls[diff < -1*(xD)])), col = LUSC_cols[['dark']], border = LUSC_cols[['dark']], stack.gap = 0)#, ygap = 0.1, height = 0.7)
dmls_BRCA_gt <- gTrack(reduce(dt2gr(dmls[diff > xD])), col = BRCA_cols[['dark']], border = BRCA_cols[['dark']], stack.gap = 0)#, ygap = 0.1, height = 0.7)
rm(dmls)
dmrs <- fread(paste0(IN_DIR, '/20260615_NSD3_TM_dmrs.csv'))
colnames(dmrs)[6:7] = c("meanMethy_NSCLC", "meanMethy_BRCA")
# DMR cutoff tracks shaded light (0.3) -> dark (0.5)
dmrs4_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -xD]), col = LUSC_cols[['mid']], border = LUSC_cols[['mid']], stack.gap = 0)
dmrs4_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > xD]), col = BRCA_cols[['mid']], border = BRCA_cols[['mid']], stack.gap = 0)
dmrs3_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -0.3]), col = LUSC_cols[['light']], border = LUSC_cols[['light']], stack.gap = 0)
dmrs3_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > 0.3]), col = BRCA_cols[['light']], border = BRCA_cols[['light']], stack.gap = 0)
dmrs5_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -0.5]), col = LUSC_cols[['darkest']], border = LUSC_cols[['darkest']], stack.gap = 0)
dmrs5_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > 0.5]), col = BRCA_cols[['darkest']], border = BRCA_cols[['darkest']], stack.gap = 0)
# one shade-coded track per tissue instead of the three dmrs3/4/5 cutoff tracks (same approach as dmrs_d below):
# each DMR is shaded by the highest |diff.Methy| cutoff (0.3/0.4/0.5) it passes, same shades as the individual cutoff tracks
dmrs[, methyl_level := ifelse(abs(diff.Methy) > 0.5, "d5", ifelse(abs(diff.Methy) > 0.4, "d4", "d3"))]
colormap_dmrs_BRCA <- list(methyl_level = c(d3 = BRCA_cols[['light']], d4 = BRCA_cols[['mid']], d5 = BRCA_cols[['darkest']]))
colormap_dmrs_NSCLC <- list(methyl_level = c(d3 = LUSC_cols[['light']], d4 = LUSC_cols[['mid']], d5 = LUSC_cols[['darkest']]))
dmrs_NSCLC_col_gt <- gTrack(dt2gr(dmrs[diff.Methy < -0.3]), gr.colorfield = "color", colormaps = colormap_dmrs_NSCLC, stack.gap = 0)
dmrs_BRCA_col_gt <- gTrack(dt2gr(dmrs[diff.Methy > 0.3]), gr.colorfield = "color", colormaps = colormap_dmrs_BRCA, stack.gap = 0)
rm(dmrs)

dmrs_d4 <- fread(paste0(IN_DIR, '/20260615_NSD3_TM_dmrs_d.csv'))
colnames(dmrs_d4)[6:7] = c("meanMethy_NSCLC", "meanMethy_BRCA")
# one shade-coded track per tissue (reds for BRCA, blues for NSCLC) instead of the three dmrs_d4/d5/d6 cutoff tracks:
# each DMR is shaded by the highest |diff.Methy| cutoff (0.4/0.5/0.6) it passes;
# categorical colormap as for the tissue-colored Fisher track (a continuous per-row color column does not render in gTrack)
dmrs_d4[, methyl_level := ifelse(abs(diff.Methy) > 0.6, "d6", ifelse(abs(diff.Methy) > 0.5, "d5", "d4"))]
colormap_dmrs_d_BRCA <- list(methyl_level = c(d4 = BRCA_cols[['light']], d5 = BRCA_cols[['dark']], d6 = BRCA_cols[['darkest']]))
colormap_dmrs_d_NSCLC <- list(methyl_level = c(d4 = LUSC_cols[['light']], d5 = LUSC_cols[['dark']], d6 = LUSC_cols[['darkest']]))
dmrs_d_NSCLC_col_gt <- gTrack(dt2gr(dmrs_d4[diff.Methy < -xD]), gr.colorfield = "color", colormaps = colormap_dmrs_d_NSCLC, stack.gap = 0)
dmrs_d_BRCA_col_gt <- gTrack(dt2gr(dmrs_d4[diff.Methy > xD]), gr.colorfield = "color", colormaps = colormap_dmrs_d_BRCA, stack.gap = 0)



##### plot #####
# track order as in the CCND1 figure: NSCLC amp/stats, NSCLC samples, LUSC median, NSCLC DMRs, BRCA-LUSC median difference, BRCA DMRs, BRCA median, BRCA samples, BRCA amp
NSCLC_smpl_gts <- c(H1E2AVA.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            H1IDEVD.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HK0HCJ3.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HNODTWC_bedmethyl_ungrouped_logmCe3_gt,
            HPR84W0.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            H1MZ9S4wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            H1AEWDG.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            H12093N_bedmethyl_ungrouped_logmCe3_gt)
BRCA_smpl_gts <- c(Cb5.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C17.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            Cb13.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C132.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C95.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C30.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C12wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C81.wf_mods.ungrouped.bedmethyl_logmCe3_gt)

# DMRs: shade-coded dmrs tracks (|diff.Methy| > 0.3 light / 0.4 mid / 0.5 darkest) instead of the shade-coded dmrs_d tracks
gts_logdiff <- c(gencode_T2T_gt,
            NSCLC_xGtad_per_amp_gt,
            NSCLC_NSD3_LvsR_sig_LR_cons80_gt,
            BRCAvsNSCLC_NSD3_sig_track_cons80_gt,
            NSCLC_smpl_gts,
            LUSC_median_bedmethyl_gt,
            dmrs_NSCLC_col_gt,
            logBRCA_minus_LUSC_median_bedmethyl_gt,
            dmrs_BRCA_col_gt,
            BRCA_median_bedmethyl_gt,
            BRCA_smpl_gts,
            BRCA_xGtad_per_amp_gt)
gts_pctdiff <- c(gencode_T2T_gt,
            NSCLC_xGtad_per_amp_gt,
            NSCLC_NSD3_LvsR_sig_LR_cons80_gt,
            BRCAvsNSCLC_NSD3_sig_track_cons80_gt,
            NSCLC_smpl_gts,
            LUSC_median_bedmethyl_gt,
            dmrs_NSCLC_col_gt,
            BRCA_minus_LUSC_median_bedmethyl_gt,
            dmrs_BRCA_col_gt,
            BRCA_median_bedmethyl_gt,
            BRCA_smpl_gts,
            BRCA_xGtad_per_amp_gt)
gts_logdiff_noamp <- c(gencode_T2T_gt,
            NSCLC_NSD3_LvsR_sig_LR_cons80_gt,
            BRCAvsNSCLC_NSD3_sig_track_cons80_gt,
            NSCLC_smpl_gts,
            LUSC_median_bedmethyl_gt,
            dmrs_NSCLC_col_gt,
            logBRCA_minus_LUSC_median_bedmethyl_gt,
            dmrs_BRCA_col_gt,
            BRCA_median_bedmethyl_gt,
            BRCA_smpl_gts)

pdf(file = paste0(out_dir, 'Fig.S3D_20261020_redblue_dmrs345col_diffmethyl_BRCA_vsNSCLC_NSD3amps.pdf')
        ,paper = 'a4', height = 64)
    # pages 1-3 as in the CCND1 figure, at the NSD3 amp-core window (AmpPercent > 0.8)
    plot(gts_logdiff, window_oi, cex.label = 0.1)
    plot(gts_pctdiff, window_oi, cex.label = 0.1)
    plot(gts_logdiff_noamp, window_oi, cex.label = 0.1)
    # wider / narrower windows
    plot(gts_logdiff, window_oi2, cex.label = 0.1)
    plot(gts_logdiff, window_oi3, cex.label = 0.1)
dev.off()
