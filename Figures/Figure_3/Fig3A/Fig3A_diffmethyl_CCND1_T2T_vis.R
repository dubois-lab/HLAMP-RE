# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(gGnome)
library(gTrack)
library(rtracklayer)

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

xG = 'CCND1'
Oncgr <- CGC_gr_T2T %Q%  (Gene_Symbol == xG)
window_oi = reduce(Oncgr + 10e5)
out_dir = paste0(OUT_DIR, '/')
BRCA_smpl_col = BRCA_cols[['mid']] # individual sample methylation tracks; the BRCA/LUSC median tracks use the darkest shade
LUSC_smpl_col = LUSC_cols[['mid']]

# median bedmethyl, window_oi only ####
# seqnames style differs between bedmethyl sources ('chr11' or '11')
# per sample: same CpG filter as the sample tracks, restricted to window_oi, matching 'chr11' and '11' and writing 'chr11' for all samples so CpGs line up in median_bm();
# percent_modified is recomputed from neglog_mth (= -log10(percent_modified/100)), which every bedmethyl rds has
win_chr = as.character(seqnames(window_oi))[1]
win_chr_alias = unique(c(win_chr, sub('^chr', '', win_chr), paste0('chr', sub('^chr', '', win_chr))))
win_start = min(start(window_oi))
win_end = max(end(window_oi))
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

# # amp dists ####
NSCLC_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, '/20260907_Fig.2C_NSCLC_CCND1_AmpFreq_track.rds'))
NSCLC_xGtad_per_amp_gt$bars =F
NSCLC_xGtad_per_amp_gt <- col_gt(fmt_gt(NSCLC_xGtad_per_amp_gt), LUSC_cols[['dark']])
BRCA_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, '/20260907_Fig.2C_BRCA_CCND1_AmpFreq_track.rds'))
BRCA_xGtad_per_amp_gt$bars =F
BRCA_xGtad_per_amp_gt <- col_gt(fmt_gt(BRCA_xGtad_per_amp_gt), BRCA_cols[['dark']])
BRCA_CCND1_LvsR_sig_LR_gt <- col_gt(fmt_gt(readRDS(paste0(IN_DIR, '/20260907c_Fig.2C_BRCA_CCND1_LvsR_sig_LR_cons80.rds'))), BRCA_cols[['dark']])
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- col_gt(fmt_gt(readRDS(paste0(IN_DIR, '/20260907c_Fig.2C_BRCAvsNSCLC_CCND1_chr11_sig_track_cons80.rds'))), LUSC_cols[['dark']])
NSCLC_CCND1_LvsR_sig_LR_cons80_gt <- col_gt(fmt_gt(readRDS(paste0(IN_DIR, '/20260907c_Fig.2C_NSCLC_CCND1_LvsR_sig_LR_cons80.rds'))), LUSC_cols[['dark']])


#C169####
C169_jabba <- gG(jabba = file.path(IN_DIR, 'jabba/breast/C169/jabba.rds'))
C169.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C169.wf_mods.ungrouped.bedmethyl.rds'))
C169.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(C169.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = BRCA_smpl_col, border = BRCA_smpl_col, stack.gap = 0))
BRCA_win_bm$C169 <- win_bm(C169.wf_mods.ungrouped.bedmethyl)
rm(C169.wf_mods.ungrouped.bedmethyl)
#C17####
C17_jabba <- gG(jabba = file.path(IN_DIR, 'jabba/breast/C17/jabba.rds'))
C17.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C17.wf_mods.ungrouped.bedmethyl.rds'))
C17.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(C17.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = BRCA_smpl_col, border = BRCA_smpl_col, stack.gap = 0))
BRCA_win_bm$C17 <- win_bm(C17.wf_mods.ungrouped.bedmethyl)
rm(C17.wf_mods.ungrouped.bedmethyl)
#C3####
C3_jabba <- gG(jabba = file.path(IN_DIR, 'jabba/breast/C3/jabba.rds'))
C3.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C3.wf_mods.ungrouped.bedmethyl.rds'))
C3.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(C3.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = BRCA_smpl_col, border = BRCA_smpl_col, stack.gap = 0))
BRCA_win_bm$C3 <- win_bm(C3.wf_mods.ungrouped.bedmethyl)
rm(C3.wf_mods.ungrouped.bedmethyl)
#C97####
C97.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C97wf_mods.ungrouped.bedmethyl.regularDoradoHac.rds'))
C97.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(C97.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = BRCA_smpl_col, border = BRCA_smpl_col, stack.gap = 0))
BRCA_win_bm$C97 <- win_bm(C97.wf_mods.ungrouped.bedmethyl)
rm(C97.wf_mods.ungrouped.bedmethyl)
#C96####
C96_jabba <- gG(jabba = file.path(IN_DIR, 'jabba/breast/C96/jabba.rds'))
C96.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C96.wf_mods.ungrouped.bedmethyl.rds'))
C96.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(C96.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = BRCA_smpl_col, border = BRCA_smpl_col, stack.gap = 0))
BRCA_win_bm$C96 <- win_bm(C96.wf_mods.ungrouped.bedmethyl)
rm(C96.wf_mods.ungrouped.bedmethyl)
#C106####
C106.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'C106wf_mods.ungrouped.bedmethyl.regularDoradoHac.rds'))
C106.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(C106.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = BRCA_smpl_col, border = BRCA_smpl_col, stack.gap = 0))
BRCA_win_bm$C106 <- win_bm(C106.wf_mods.ungrouped.bedmethyl)
rm(C106.wf_mods.ungrouped.bedmethyl)
#Cb32## CN=10, no RNA
Cb32.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'Cb32wf_mods.ungrouped.bedmethyl.regularDoradoHac.rds'))
Cb32_gt <- fmt_gt(gTrack(dt2gr(Cb32.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field='neglog_mth', bars=T, y0=0,y1=3, col = BRCA_smpl_col, border = BRCA_smpl_col, stack.gap = 0))
BRCA_win_bm$Cb32 <- win_bm(Cb32.wf_mods.ungrouped.bedmethyl)
rm(Cb32.wf_mods.ungrouped.bedmethyl)
#Cb11## CN=10, no RNA
Cb11.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'Cb11wf_mods.ungrouped.bedmethyl.regularDoradoHac.rds'))
Cb11_gt <- fmt_gt(gTrack(dt2gr(Cb11.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field='neglog_mth', bars=T, y0=0,y1=3, col = BRCA_smpl_col, border = BRCA_smpl_col, stack.gap = 0))
BRCA_win_bm$Cb11 <- win_bm(Cb11.wf_mods.ungrouped.bedmethyl)
rm(Cb11.wf_mods.ungrouped.bedmethyl)
#BRCA median####
BRCA_median_bedmethyl <- median_bm(BRCA_win_bm)
rm(BRCA_win_bm)
saveRDS(BRCA_median_bedmethyl, paste0(out_dir, '20261019_BRCA_CCND1_window_oi_median_bedmethyl.rds'))
BRCA_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(BRCA_median_bedmethyl), y.field = 'neglog_mth', name = 'BRCA_median_-log5mC_%', 
    bars = T, y1= 2.5, col = BRCA_cols[['darkest']], border = BRCA_cols[['darkest']], stack.gap = 0, height =20))


#NSCLC ####
#1WNRWL####
H1WNRWL_jabba <- gG(jabba = file.path(IN_DIR, 'jabba/lung/1WNRWL/jabba.rds'))
H1WNRWL.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, '1WNRWL.wf_mods.ungrouped.bedmethyl.rds'))
H1WNRWL.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(H1WNRWL.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = LUSC_smpl_col, border = LUSC_smpl_col, stack.gap = 0))
LUSC_win_bm$H1WNRWL <- win_bm(H1WNRWL.wf_mods.ungrouped.bedmethyl)
rm(H1WNRWL.wf_mods.ungrouped.bedmethyl)
#Y57U8R####
HY57U8R_jabba <- gG(jabba = file.path(IN_DIR, 'jabba/lung/Y57U8R/jabba.rds'))
HY57U8R.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'Y57U8R.wf_mods.ungrouped.bedmethyl.rds'))
HY57U8R.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(HY57U8R.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = LUSC_smpl_col, border = LUSC_smpl_col, stack.gap = 0))
LUSC_win_bm$HY57U8R <- win_bm(HY57U8R.wf_mods.ungrouped.bedmethyl)
rm(HY57U8R.wf_mods.ungrouped.bedmethyl)
#998A7Z####
H998A7Z_jabba <- gG(jabba = file.path(IN_DIR, 'jabba/lung/998A7Z/jabba.rds'))
H998A7Z.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, '998A7Z.wf_mods.ungrouped.bedmethyl.rds'))
H998A7Z.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(H998A7Z.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = LUSC_smpl_col, border = LUSC_smpl_col, stack.gap = 0))
LUSC_win_bm$H998A7Z <- win_bm(H998A7Z.wf_mods.ungrouped.bedmethyl)
rm(H998A7Z.wf_mods.ungrouped.bedmethyl)
#KDAMIF####
HKDAMIF_jabba <- gG(jabba = file.path(IN_DIR, 'jabba/lung/KDAMIF/jabba.rds'))
HKDAMIF.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'KDAMIF.wf_mods.ungrouped.bedmethyl.rds'))
HKDAMIF.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(HKDAMIF.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = LUSC_smpl_col, border = LUSC_smpl_col, stack.gap = 0))
LUSC_win_bm$HKDAMIF <- win_bm(HKDAMIF.wf_mods.ungrouped.bedmethyl)
rm(HKDAMIF.wf_mods.ungrouped.bedmethyl)
#H3D44S####
HH3D44S_jabba <- gG(jabba = file.path(IN_DIR, 'jabba/lung/H3D44S/jabba.rds'))
HH3D44S.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'H3D44S.wf_mods.ungrouped.bedmethyl.rds'))
HH3D44S.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(HH3D44S.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = LUSC_smpl_col, border = LUSC_smpl_col, stack.gap = 0))
LUSC_win_bm$HH3D44S <- win_bm(HH3D44S.wf_mods.ungrouped.bedmethyl)
rm(HH3D44S.wf_mods.ungrouped.bedmethyl)
#F57J9F####
HF57J9F_jabba <- gG(jabba = file.path(IN_DIR, 'jabba/lung/F57J9F/jabba.rds'))
HF57J9F.wf_mods.ungrouped.bedmethyl <- readRDS(file.path(IN_DIR, 'F57J9Fwf_mods.ungrouped.bedmethyl.regularDoradoHac.rds'))
HF57J9F.wf_mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(HF57J9F.wf_mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3, col = LUSC_smpl_col, border = LUSC_smpl_col, stack.gap = 0))
LUSC_win_bm$HF57J9F <- win_bm(HF57J9F.wf_mods.ungrouped.bedmethyl)
rm(HF57J9F.wf_mods.ungrouped.bedmethyl)
#LUSC median####
LUSC_median_bedmethyl <- median_bm(LUSC_win_bm)
rm(LUSC_win_bm)
saveRDS(LUSC_median_bedmethyl, paste0(out_dir, '20261019_LUSC_CCND1_window_oi_median_bedmethyl.rds'))
LUSC_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(LUSC_median_bedmethyl), y.field = 'neglog_mth', name = 'LUSC_median_-log5mC_%', bars = T, y1= 2.5, 
    col = LUSC_cols[['darkest']], border = LUSC_cols[['darkest']], stack.gap = 0, height =20))

#BRCA - LUSC median####
# per CpG covered by both median tracks: BRCA minus LUSC median value (diff_neglog_mth = plotted track value; diff_percent_modified kept for reference);
# neglog_mth is capped at the median tracks' y1 (3) before subtracting, so CpGs with median 0% (neglog_mth = Inf) give finite differences
nl_cap = 2
BRCA_minus_LUSC_median_bedmethyl <- as.data.table(merge(
    BRCA_median_bedmethyl[, .(seqnames, start, end, BRCA_percent_modified = percent_modified, BRCA_neglog_mth = neglog_mth, BRCA_n_samples = n_samples)],
    LUSC_median_bedmethyl[, .(seqnames, start, end, LUSC_percent_modified = percent_modified, LUSC_neglog_mth = neglog_mth, LUSC_n_samples = n_samples)],
    by = c('seqnames', 'start', 'end'))) # merge() returns a plain data.frame here, so := below fails without as.data.table()
BRCA_minus_LUSC_median_bedmethyl[, ':='(diff_percent_modified = BRCA_percent_modified - LUSC_percent_modified,
                                        diff_neglog_mth = pmin(BRCA_neglog_mth, nl_cap) - pmin(LUSC_neglog_mth, nl_cap))]
setorder(BRCA_minus_LUSC_median_bedmethyl, seqnames, start)
# saveRDS(BRCA_minus_LUSC_median_bedmethyl, paste0(out_dir, '20261014_BRCA_minus_LUSC_CCND1_window_oi_median_bedmethyl.rds'))
logBRCA_minus_LUSC_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(BRCA_minus_LUSC_median_bedmethyl), y.field = 'diff_neglog_mth', name = 'log_BRCA-LUSC_median_-log5mC_%', 
    bars = T, y0 = -nl_cap, y1 = nl_cap, col = 'black', border = 'black', stack.gap = 0, height =20))
# saveRDS(BRCA_minus_LUSC_median_bedmethyl_gt, paste0(out_dir, '20261014_BRCA_minus_LUSC_CCND1_window_oi_median_bedmethyl_track.rds'))
BRCA_minus_LUSC_median_bedmethyl_gt <- fmt_gt(gTrack(dt2gr(BRCA_minus_LUSC_median_bedmethyl), y.field = 'diff_percent_modified', name = 'BRCA-LUSC_median_-log5mC_%', 
    bars = T, y0 = -70, y1 = 70, col = 'black', border = 'black', stack.gap = 0, height =20))


# #functional val####
GERCARres <- as.data.table(read.csv2(file.path(IN_DIR, 'GERCARres.csv')))
colnames(GERCARres)
GERCARres_gr <- dt2gr(GERCARres[chromosome == unique(as.character(seqnames(Oncgr)))])
GERCARres_gr_t2tgr <- grl.unlist(liftOver(gr.chr(GERCARres_gr), conversion_chain))
GERCARres_prolifd21_gt <- gTrack(GERCARres_gr_t2tgr, y.field = 'logFC_prolif_21', bars = T, col = BRCA_cols[['mid']], border = BRCA_cols[['mid']], stack.gap = 0)
GERCARres_gr_sig <- dt2gr(gr2dt(GERCARres_gr_t2tgr)[FDR_prolif_21 < 0.1 ])
GERCARres_prolifd21_q0.1_gt <- gTrack(GERCARres_gr_sig, col = BRCA_cols[['dark']], border = BRCA_cols[['dark']], stack.gap = 0)
#####stats#####

ATACfisher <- fread(file.path(IN_DIR, '20260614_CCND1_H3K27ac&ATAC_linDiffs_stats.csv'))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = '[:]'))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[3]), by = peakID]
ATACfisher_gr <- dt2gr(ATACfisher)
ATACfisher_t2tgr <- grl.unlist(liftOver(gr.chr(ATACfisher_gr), conversion_chain))

ATACfisher_HRpBC_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.05 & tstdttype == 'BRCA_HR+')), col = BRCA_cols[['dark']], border = BRCA_cols[['dark']], stack.gap = 0)
ATACfisher_LUSC_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.05 & tstdttype == 'LUSC')), col = LUSC_cols[['dark']], border = LUSC_cols[['dark']], stack.gap = 0)
xD = 0.4

dmls <- fread(paste0(IN_DIR, '/20260512_CCND1_TM_dmls.csv'))
dmls[, end := pos +1, by =pos]
dmls_NSCLC_gt <- gTrack(reduce(dt2gr(dmls[diff < -1*(xD)])), col = LUSC_cols[['dark']], border = LUSC_cols[['dark']], stack.gap = 0)#, ygap = 0.1, height = 0.7)
dmls_BRCA_gt <- gTrack(reduce(dt2gr(dmls[diff > xD])), col = BRCA_cols[['dark']], border = BRCA_cols[['dark']], stack.gap = 0)#, ygap = 0.1, height = 0.7)
rm(dmls)
dmrs <- fread(paste0(IN_DIR, '/20260512_CCND1_TM_dmrs.csv'))
colnames(dmrs)[6:7] = c("meanMethy_NSCLC", "meanMethy_BRCA")
# DMR cutoff tracks shaded light (0.3) -> dark (0.5)
dmrs4_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -xD]), col = LUSC_cols[['mid']], border = LUSC_cols[['mid']], stack.gap = 0)
dmrs4_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > xD]), col = BRCA_cols[['mid']], border = BRCA_cols[['mid']], stack.gap = 0)
dmrs3_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -0.3]), col = LUSC_cols[['light']], border = LUSC_cols[['light']], stack.gap = 0)
dmrs3_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > 0.3]), col = BRCA_cols[['light']], border = BRCA_cols[['light']], stack.gap = 0)
dmrs5_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -0.5]), col = LUSC_cols[['darkest']], border = LUSC_cols[['darkest']], stack.gap = 0)
dmrs5_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > 0.5]), col = BRCA_cols[['darkest']], border = BRCA_cols[['darkest']], stack.gap = 0)


dmrs_d4 <- fread(paste0(IN_DIR, '/20260512_CCND1_TM_dmrs_d.csv'))
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
window_oi = Oncgr + 10e5
window_oi2 = Oncgr + 5e5
window_oi3 = Oncgr + 2e5

pdf(file = paste0(out_dir, 'Fig.3A_20261019_redblue_diffmethyl_BRCA_vsNSCLC_CCND1amps_.pdf')
        ,paper = 'a4', height = 64)
    plot(c(gencode_T2T_gt,
            NSCLC_xGtad_per_amp_gt,
            BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
            H1WNRWL.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HY57U8R.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            H998A7Z.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HKDAMIF.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HH3D44S.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HF57J9F.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            LUSC_median_bedmethyl_gt,
            dmrs_d_NSCLC_col_gt,
            logBRCA_minus_LUSC_median_bedmethyl_gt,
            dmrs_d_BRCA_col_gt,
            BRCA_median_bedmethyl_gt,
            Cb11_gt,
            Cb32_gt,
            C169.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C17.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C97.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C3.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C96.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C106.wf_mods.ungrouped.bedmethyl_logmCe3_gt,  
            BRCA_CCND1_LvsR_sig_LR_gt, BRCA_xGtad_per_amp_gt  
    ), window_oi2, cex.label = 0.1)
    plot(c(gencode_T2T_gt,
            NSCLC_xGtad_per_amp_gt,
            BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
            H1WNRWL.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HY57U8R.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            H998A7Z.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HKDAMIF.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HH3D44S.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HF57J9F.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            LUSC_median_bedmethyl_gt,
            dmrs_d_NSCLC_col_gt,
            BRCA_minus_LUSC_median_bedmethyl_gt,
            dmrs_d_BRCA_col_gt,
            BRCA_median_bedmethyl_gt,
            Cb11_gt,
            Cb32_gt,
            C169.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C17.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C97.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C3.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C96.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C106.wf_mods.ungrouped.bedmethyl_logmCe3_gt,  
            BRCA_CCND1_LvsR_sig_LR_gt, BRCA_xGtad_per_amp_gt  
    ), window_oi2, cex.label = 0.1)
    plot(c(gencode_T2T_gt,
            BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
            H1WNRWL.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HY57U8R.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            H998A7Z.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HKDAMIF.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HH3D44S.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            HF57J9F.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            LUSC_median_bedmethyl_gt,
            dmrs_d_NSCLC_col_gt,
            logBRCA_minus_LUSC_median_bedmethyl_gt,
            dmrs_d_BRCA_col_gt,
            BRCA_median_bedmethyl_gt,
            Cb11_gt,
            Cb32_gt,
            C169.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C17.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C97.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C3.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C96.wf_mods.ungrouped.bedmethyl_logmCe3_gt,
            C106.wf_mods.ungrouped.bedmethyl_logmCe3_gt,  
            BRCA_CCND1_LvsR_sig_LR_gt
    ), window_oi2, cex.label = 0.1)
dev.off()


