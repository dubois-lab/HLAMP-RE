# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(gGnome)
library(gTrack)
library(rtracklayer)

# EGFR version of the Fig. 3A script (red/blue scheme, median + BRCA-LUSC median difference tracks, shade-coded DMR tracks)

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

xG = 'EGFR'
Oncgr <- CGC_gr_T2T %Q%  (Gene_Symbol == xG)
out_dir = paste0(OUT_DIR, '/')
BRCA_smpl_col = BRCA_cols[['mid']] # individual sample methylation tracks; the BRCA/LUSC median tracks use the darkest shade
LUSC_smpl_col = LUSC_cols[['mid']]

# individual sample methylation track, same CpG filter and formatting as in the CCND1 script (Fig. 3A)
smpl_gt <- function(bm, col) {
    fmt_gt(gTrack(dt2gr(bm[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = col, border = col, stack.gap = 0))
}

# # amp dists ####
BRCA_EGFR_CNMatm <- readRDS(file.path(IN_DIR, '20260910_BRCA_EGFR_CNMatm_combined_hs1.rds'))
NSCLC_EGFR_CNMatm <- readRDS(file.path(IN_DIR, '20260910_NSCLC_LUSC_EGFR_CNMatm_combined_hs1.rds'))
BRCA_xGtad_per_amp_gt <- fmt_gt(gTrack(dt2gr(BRCA_EGFR_CNMatm), y.field = 'AmpPercent', bars = F, y0=0, y1=1, col = BRCA_cols[['dark']], border = BRCA_cols[['dark']]))
NSCLC_xGtad_per_amp_gt <- fmt_gt(gTrack(dt2gr(NSCLC_EGFR_CNMatm), y.field = 'AmpPercent', bars = F, y0=0, y1=1, col = LUSC_cols[['dark']], border = LUSC_cols[['dark']]))
window_oi <- reduce(dt2gr(rbind(BRCA_EGFR_CNMatm[AmpPercent > 0.8, 1:10], NSCLC_EGFR_CNMatm[AmpPercent > 0.8, 1:10])))
window_oi2 = window_oi + 1e5
window_oi3 = window_oi - 1e5

# BRCA vs NSCLC EGFR differential amp-frequency significance track (chr7)
BRCAvsNSCLC_EGFR_fisher_gr <- readRDS(file.path(IN_DIR, '20260910_BRCAvsNSCLC_EGFR_chr7_FisherStats_CN10_cons80_gr.rds'))
BRCAvsNSCLC_EGFR_sig <- BRCAvsNSCLC_EGFR_fisher_gr[BRCAvsNSCLC_EGFR_fisher_gr$Fisher_pval < 0.05]
mcols(BRCAvsNSCLC_EGFR_sig)$tissue <- ifelse(mcols(BRCAvsNSCLC_EGFR_sig)$BRCA_AmpFreq > mcols(BRCAvsNSCLC_EGFR_sig)$NSCLC_AmpFreq, "BRCA", "NSCLC")
colormap_BRCAvsNSCLC_EGFR <- list(tissue = c(NSCLC = LUSC_cols[['dark']], BRCA = BRCA_cols[['dark']]))
BRCAvsNSCLC_EGFR_chr7_sig_track_cons80_gt <- fmt_gt(gTrack(BRCAvsNSCLC_EGFR_sig, gr.colorfield = "tissue", colormaps = colormap_BRCAvsNSCLC_EGFR, height = 0.5, stack.gap = 0))

# no EGFR LeftVsRight Fisher stats available (BRCA or NSCLC) -> no LvsR track here

# median bedmethyl, median window only ####
# median window: EGFR +-1Mb (as for CCND1), extended to cover the widest plot window (window_oi2)
# per sample: same CpG filter as the sample tracks, restricted to the median window;
# percent_modified is recomputed from neglog_mth (= -log10(percent_modified/100)), which every bedmethyl rds has
win_chr = as.character(seqnames(Oncgr))[1]
win_start = min(start(Oncgr + 10e5), start(window_oi2))
win_end = max(end(Oncgr + 10e5), end(window_oi2))
win_bm <- function(bm) bm[modified_base_code == 'm' & Nvalid_cov > 5 & seqnames == win_chr & start >= win_start & end <= win_end,
                          .(seqnames, start, end, Nvalid_cov, percent_modified = 100 * 10^(-neglog_mth))]
# median per CpG over the samples covering it (n_samples); neglog_mth recomputed from the median percent_modified
median_bm <- function(win_bms) {
    med <- rbindlist(win_bms)[, .(Nvalid_cov = median(as.numeric(Nvalid_cov)), percent_modified = median(percent_modified), n_samples = .N),
                              by = .(seqnames, start, end)]
    med[, ':='(modified_base_code = 'm', neglog_mth = -log10(percent_modified/100))]
    setorder(med, seqnames, start)
    med
}
BRCA_win_bm <- list()
LUSC_win_bm <- list()

#BRCA ####
#C109####
C109.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C109.wf_mods.ungrouped.bedmethyl.rds'))
C109.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(C109.wf_mods.ungrouped.bedmethyl, BRCA_smpl_col)
BRCA_win_bm$C109 <- win_bm(C109.wf_mods.ungrouped.bedmethyl)
rm(C109.wf_mods.ungrouped.bedmethyl)
#C208####
C208.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C208.wf_mods.ungrouped.bedmethyl.rds'))
C208.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(C208.wf_mods.ungrouped.bedmethyl, BRCA_smpl_col)
BRCA_win_bm$C208 <- win_bm(C208.wf_mods.ungrouped.bedmethyl)
rm(C208.wf_mods.ungrouped.bedmethyl)
#C168####
C168.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C168.wf_mods.ungrouped.bedmethyl.rds'))
C168.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(C168.wf_mods.ungrouped.bedmethyl, BRCA_smpl_col)
BRCA_win_bm$C168 <- win_bm(C168.wf_mods.ungrouped.bedmethyl)
rm(C168.wf_mods.ungrouped.bedmethyl)
#C15####
C15.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C15.wf_mods.ungrouped.bedmethyl.rds'))
C15.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(C15.wf_mods.ungrouped.bedmethyl, BRCA_smpl_col)
BRCA_win_bm$C15 <- win_bm(C15.wf_mods.ungrouped.bedmethyl)
rm(C15.wf_mods.ungrouped.bedmethyl)
#BRCA median####
BRCA_median_bedmethyl <- median_bm(BRCA_win_bm)
rm(BRCA_win_bm)
saveRDS(BRCA_median_bedmethyl, paste0(out_dir, '20261017_BRCA_EGFR_window_oi_median_bedmethyl.rds'))
BRCA_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(BRCA_median_bedmethyl), y.field = 'neglog_mth', name = 'BRCA_median_-log5mC_%',
    bars = T, y1= 2.5, col = BRCA_cols[['darkest']], border = BRCA_cols[['darkest']], stack.gap = 0, height =20))

#NSCLC ####
#1WNRWL####
H1WNRWL.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, '1WNRWL.wf_mods.ungrouped.bedmethyl.rds'))
H1WNRWL.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(H1WNRWL.wf_mods.ungrouped.bedmethyl, LUSC_smpl_col)
LUSC_win_bm$H1WNRWL <- win_bm(H1WNRWL.wf_mods.ungrouped.bedmethyl)
rm(H1WNRWL.wf_mods.ungrouped.bedmethyl)
#1HGXBV####
H1HGXBV.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, '1HGXBV.wf_mods.ungrouped.bedmethyl.rds'))
H1HGXBV.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(H1HGXBV.wf_mods.ungrouped.bedmethyl, LUSC_smpl_col)
LUSC_win_bm$H1HGXBV <- win_bm(H1HGXBV.wf_mods.ungrouped.bedmethyl)
rm(H1HGXBV.wf_mods.ungrouped.bedmethyl)
#XO7JTV####
HXO7JTV.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'XO7JTV.wf_mods.ungrouped.bedmethyl.rds'))
HXO7JTV.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(HXO7JTV.wf_mods.ungrouped.bedmethyl, LUSC_smpl_col)
LUSC_win_bm$HXO7JTV <- win_bm(HXO7JTV.wf_mods.ungrouped.bedmethyl)
rm(HXO7JTV.wf_mods.ungrouped.bedmethyl)
#AJZ5JO####
HAJZ5JO.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'AJZ5JO.wf_mods.ungrouped.bedmethyl.rds'))
HAJZ5JO.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- smpl_gt(HAJZ5JO.wf_mods.ungrouped.bedmethyl, LUSC_smpl_col)
LUSC_win_bm$HAJZ5JO <- win_bm(HAJZ5JO.wf_mods.ungrouped.bedmethyl)
rm(HAJZ5JO.wf_mods.ungrouped.bedmethyl)
#LUSC median####
LUSC_median_bedmethyl <- median_bm(LUSC_win_bm)
rm(LUSC_win_bm)
saveRDS(LUSC_median_bedmethyl, paste0(out_dir, '20261017_LUSC_EGFR_window_oi_median_bedmethyl.rds'))
LUSC_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(LUSC_median_bedmethyl), y.field = 'neglog_mth', name = 'LUSC_median_-log5mC_%', bars = T, y1= 2.5,
    col = LUSC_cols[['darkest']], border = LUSC_cols[['darkest']], stack.gap = 0, height =20))

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
saveRDS(BRCA_minus_LUSC_median_bedmethyl, paste0(out_dir, '20261017_BRCA_minus_LUSC_EGFR_window_oi_median_bedmethyl.rds'))
logBRCA_minus_LUSC_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(BRCA_minus_LUSC_median_bedmethyl), y.field = 'diff_neglog_mth', name = 'log_BRCA-LUSC_median_-log5mC_%',
    bars = T, y0 = -nl_cap, y1 = nl_cap, col = 'black', border = 'black', stack.gap = 0, height =20))
BRCA_minus_LUSC_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(BRCA_minus_LUSC_median_bedmethyl), y.field = 'diff_percent_modified', name = 'BRCA-LUSC_median_-log5mC_%',
    bars = T, y0 = -70, y1 = 70, col = 'black', border = 'black', stack.gap = 0, height =20))


#####stats#####

ATACfisher <- fread(file.path(IN_DIR, '20260623q_EGFR_H3K27ac_ATAC_linDiffs_stats.csv'))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = '[:]'))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[3]), by = peakID]
ATACfisher_gr <- dt2gr(ATACfisher)
ATACfisher_t2tgr <- grl.unlist(liftOver(gr.chr(ATACfisher_gr), conversion_chain))

ATACfisher_HRnBC_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.15 & tstdttype == 'BRCA_HR-')), col = BRCA_cols[['dark']], border = BRCA_cols[['dark']], stack.gap = 0)
ATACfisher_LUSC_gt <- gTrack(data = (ATACfisher_t2tgr %Q% (q_val < 0.15 & tstdttype == 'LUSC')), col = LUSC_cols[['dark']], border = LUSC_cols[['dark']], stack.gap = 0)
xD = 0.4

dmls <- fread(paste0(IN_DIR, '/20251003_EGFR_TM_dmls_delt_',xD ,'.csv'))
dmls[, end := pos +1, by =pos]
dmls_NSCLC_gt <- gTrack(reduce(dt2gr(dmls[diff < -1*(xD)])), col = LUSC_cols[['dark']], border = LUSC_cols[['dark']], stack.gap = 0)#, ygap = 0.1, height = 0.7)
dmls_BRCA_gt <- gTrack(reduce(dt2gr(dmls[diff > xD])), col = BRCA_cols[['dark']], border = BRCA_cols[['dark']], stack.gap = 0)#, ygap = 0.1, height = 0.7)
rm(dmls)
dmrs <- fread(paste0(IN_DIR, '/20251003_EGFR_TM_dmrs_delt_',xD ,'.csv'))
colnames(dmrs)[6:7] = c("meanMethy_NSCLC", "meanMethy_BRCA")
# DMR cutoff tracks shaded light (0.3) -> dark (0.5)
dmrs4_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -xD]), col = LUSC_cols[['mid']], border = LUSC_cols[['mid']], stack.gap = 0)
dmrs4_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > xD]), col = BRCA_cols[['mid']], border = BRCA_cols[['mid']], stack.gap = 0)
dmrs3_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -0.3]), col = LUSC_cols[['light']], border = LUSC_cols[['light']], stack.gap = 0)
dmrs3_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > 0.3]), col = BRCA_cols[['light']], border = BRCA_cols[['light']], stack.gap = 0)
dmrs5_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -0.5]), col = LUSC_cols[['darkest']], border = LUSC_cols[['darkest']], stack.gap = 0)
dmrs5_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > 0.5]), col = BRCA_cols[['darkest']], border = BRCA_cols[['darkest']], stack.gap = 0)
rm(dmrs)

dmrs_d4 <- fread(paste0(IN_DIR, '/20251003_EGFR_TM_dmrs_d_delt_',xD ,'.csv'))
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
NSCLC_smpl_gts <- c(H1WNRWL.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            H1HGXBV.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HXO7JTV.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HAJZ5JO.wf_mods.ungrouped.bedmethyl_logmCe3_gt)
BRCA_smpl_gts <- c(C109.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C208.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C168.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C15.wf_mods.ungrouped.bedmethyl_logmCe3_gt)

gts_logdiff <- c(gencode_T2T_gt,
            NSCLC_xGtad_per_amp_gt,
            BRCAvsNSCLC_EGFR_chr7_sig_track_cons80_gt,
            NSCLC_smpl_gts,
            LUSC_median_bedmethyl_gt,
            dmrs_d_NSCLC_col_gt,
            logBRCA_minus_LUSC_median_bedmethyl_gt,
            dmrs_d_BRCA_col_gt,
            BRCA_median_bedmethyl_gt,
            BRCA_smpl_gts,
            BRCA_xGtad_per_amp_gt)
gts_pctdiff <- c(gencode_T2T_gt,
            NSCLC_xGtad_per_amp_gt,
            BRCAvsNSCLC_EGFR_chr7_sig_track_cons80_gt,
            NSCLC_smpl_gts,
            LUSC_median_bedmethyl_gt,
            dmrs_d_NSCLC_col_gt,
            BRCA_minus_LUSC_median_bedmethyl_gt,
            dmrs_d_BRCA_col_gt,
            BRCA_median_bedmethyl_gt,
            BRCA_smpl_gts,
            BRCA_xGtad_per_amp_gt)
gts_logdiff_noamp <- c(gencode_T2T_gt,
            BRCAvsNSCLC_EGFR_chr7_sig_track_cons80_gt,
            NSCLC_smpl_gts,
            LUSC_median_bedmethyl_gt,
            dmrs_d_NSCLC_col_gt,
            logBRCA_minus_LUSC_median_bedmethyl_gt,
            dmrs_d_BRCA_col_gt,
            BRCA_median_bedmethyl_gt,
            BRCA_smpl_gts)

pdf(file = paste0(out_dir, 'Fig.S3F_20261017_redblue_diffmethyl_BRCA_vsNSCLC_EGFRamps.pdf')
        ,paper = 'a4', height = 64)
    # pages 1-3 as in the CCND1 figure, at the EGFR amp-core window (AmpPercent > 0.8)
    plot(gts_logdiff, window_oi, cex.label = 0.1)
    plot(gts_pctdiff, window_oi, cex.label = 0.1)
    plot(gts_logdiff_noamp, window_oi, cex.label = 0.1)
    # wider / narrower windows
    plot(gts_logdiff, window_oi2, cex.label = 0.1)
    plot(gts_logdiff, window_oi3, cex.label = 0.1)
dev.off()
