# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(GxG)
library(gGnome)
library(rtracklayer)
conversion_chain <- import.chain(file.path(IN_DIR, 'hg38ToHs1.over.chain'))

outdir <- paste0(OUT_DIR, '/')

gencode_T2T_gt <- readRDS(file.path(IN_DIR, 'gencodeTrack_t2t.rds'))
CGC_gr_T2T <- readRDS('../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')
xG = 'EGFR'
Oncgr <- CGC_gr_T2T %Q%  (Gene_Symbol == xG)
window_oi = reduce(Oncgr + 10e5)
binsize =5e3

xT_bins <- fread(file.path(IN_DIR, 'normalized_abs.bed'))
colnames(xT_bins) <- c('chr', 'start', 'end')
xT_bins[, start := start +1, by = start]

# EGFR 10%
xG_hicdcdiff <- fread(paste0(IN_DIR, '/EGFR_hicdcdiff_breastVSlung_fitTypeLocal_10prctSig_LinkPipediff_resbreastoverlung_chr7.txt'))
xG_hicdcdiff[, nlg_padj := -log10(padj)]
xG_hicdcdiff[log2FoldChange > 0, nlg_padj_sgn := nlg_padj]
xG_hicdcdiff[log2FoldChange < 0, nlg_padj_sgn := -nlg_padj]
xG_hicdcdiff[, nlg_p := -log10(pvalue)]
xG_hicdcdiff[log2FoldChange > 0, nlg_p_sgn := nlg_p]
xG_hicdcdiff[log2FoldChange < 0, nlg_p_sgn := -nlg_p]
xG_hicdcdiff[, i := (startI/ binsize) + 1, by = startI]
xG_hicdcdiff[, j := (startJ/ binsize) + 1, by = startJ]

gmat_EGFR_10prct_diff_nlgq <- gM(gr=dt2gr(xT_bins[chr== 'chr7']),dat=xG_hicdcdiff[, .(i, j, value=nlg_padj_sgn)],field="nlg_padj_sgn",fill=0,agg.fun=min)
gmat_EGFR_10prct_diff_l2fc <- gM(gr=dt2gr(xT_bins[chr== 'chr7']),dat=xG_hicdcdiff[, .(i, j, value=log2FoldChange)],field="log2FoldChange",fill=0,agg.fun=min)

EGFR_TSSgr <- readRDS(file.path(IN_DIR, 'EGFR_TSS_T2T.rds'))

xG_hicdcdiff_4C_nlgq <- gmat_EGFR_10prct_diff_nlgq[gmat_EGFR_10prct_diff_nlgq$gr %^% EGFR_TSSgr, drop =T]
# saveRDS(xG_hicdcdiff_4C_nlgq, file = paste0(outdir,'EGFR_hicdcdiff_4C_nlgq.rds'))
xG_hicdcdiff_4C_nlgq_BRCA <- xG_hicdcdiff_4C_nlgq %Q% (!is.na(value) & value > 1)
xG_hicdcdiff_4C_nlgq_BRCA_gt <- gTrack(data = xG_hicdcdiff_4C_nlgq_BRCA, name = '4C_hicdcdiff_BRCA', height = 0.5, col = 'red')
xG_hicdcdiff_4C_nlgq_LUSC<- xG_hicdcdiff_4C_nlgq %Q% (!is.na(value) & value < -1)
xG_hicdcdiff_4C_nlgq_LUSC_gt <- gTrack(data = xG_hicdcdiff_4C_nlgq_LUSC, name = '4C_hicdcdiff_LUSC', height = 0.5, col = 'blue')

xG_hicdcdiff_4C_l2fc <- gmat_EGFR_10prct_diff_l2fc[gmat_EGFR_10prct_diff_l2fc$gr %^% EGFR_TSSgr, drop =T]
# saveRDS(xG_hicdcdiff_4C_l2fc, file = paste0(outdir,'EGFR_hicdcdiff_4C_l2fc.rds'))
xG_hicdcdiff_4C_l2fc_gt <- gTrack(data = xG_hicdcdiff_4C_l2fc, y.field = 'value', bars = T, name = '4C_hicdcdiff_BRCA', y0 = -2, y1 = 2)

BRCA_EGFR_CNMatm <- readRDS(file.path(IN_DIR, "20260910_BRCA_EGFR_CNMatm_combined_hs1.rds"))
NSCLC_EGFR_CNMatm <- readRDS(file.path(IN_DIR, "20260910_NSCLC_LUSC_EGFR_CNMatm_combined_hs1.rds"))
BRCA_xGtad_per_amp_gt <- gTrack(dt2gr(BRCA_EGFR_CNMatm), y.field = 'AmpPercent', col = 'brown3', bars = F, y1=1, y0=0)
NSCLC_xGtad_per_amp_gt <- gTrack(dt2gr(NSCLC_EGFR_CNMatm), y.field = 'AmpPercent', col = 'royalblue3', bars = F, y1=1, y0=0)

# BRCA vs NSCLC EGFR differential amp-frequency significance track (chr7)
# Note: No NSCLC EGFR left-vs-right Fisher stats available for EGFR
BRCAvsNSCLC_EGFR_fisher_gr <- readRDS(file.path(IN_DIR, "20260910_BRCAvsNSCLC_EGFR_chr7_FisherStats_CN10_cons80_gr.rds"))
BRCAvsNSCLC_EGFR_sig <- BRCAvsNSCLC_EGFR_fisher_gr[BRCAvsNSCLC_EGFR_fisher_gr$Fisher_pval < 0.05]
mcols(BRCAvsNSCLC_EGFR_sig)$tissue <- ifelse(mcols(BRCAvsNSCLC_EGFR_sig)$BRCA_AmpFreq > mcols(BRCAvsNSCLC_EGFR_sig)$NSCLC_AmpFreq, "BRCA", "NSCLC")
colormap_BRCAvsNSCLC_EGFR <- list(tissue = c(NSCLC = "royalblue3", BRCA = "brown3"))
BRCAvsNSCLC_EGFR_chr7_sig_track_cons80_gt <- gTrack(BRCAvsNSCLC_EGFR_sig, gr.colorfield = "color", colormaps = colormap_BRCAvsNSCLC_EGFR, height = 0.5)

# #functional val####
GERCARres <- as.data.table(read.csv2(file.path(IN_DIR, 'GERCARres.csv')))
GERCARres_gr <- dt2gr(GERCARres[chromosome == unique(as.character(seqnames(Oncgr)))])
GERCARres_gr_t2tgr <- grl.unlist(liftOver(gr.chr(GERCARres_gr), conversion_chain))
GERCARres_prolifd21_gt <- gTrack(GERCARres_gr_t2tgr, y.field = 'logFC_prolif_21', bars = T)
GERCARres_gr_sig <- dt2gr(gr2dt(GERCARres_gr_t2tgr)[FDR_prolif_21 < 0.1 ])
GERCARres_prolifd21_q0.1_gt <- gTrack(GERCARres_gr_sig)

SID_MCF7_E2 <- as.data.table(read.csv2(file.path(IN_DIR, 'cd-23-1157_supplementary_table_s3_suppst3_MCF7_E2.csv')))
SID_MCF7_E2_gr <- dt2gr(SID_MCF7_E2)
# add liftover to t2t
SID_MCF7_E2_gr_t2tgr <- grl.unlist(liftOver(gr.chr(SID_MCF7_E2_gr), conversion_chain))
SID_MCF7_E2_prolifd21_gt <- gTrack(SID_MCF7_E2_gr_t2tgr, y.field = 'log2FC_day21_vs_day07', bars = TRUE)
SID_MCF7_E2_sig <- SID_MCF7_E2[FDR_day21_vs_day07 < 0.1]
SID_MCF7_E2_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_MCF7_E2_sig))

SID_T47D_E2 <- as.data.table(read.csv2(file.path(IN_DIR, 'cd-23-1157_supplementary_table_s4_suppst4_T47D_E2.csv')))
SID_T47D_E2_gr <- dt2gr(SID_T47D_E2)
# add liftover to t2t
SID_T47D_E2_gr_t2tgr <- grl.unlist(liftOver(gr.chr(SID_T47D_E2_gr), conversion_chain))
SID_T47D_E2_prolifd21_gt <- gTrack(SID_T47D_E2_gr_t2tgr, y.field = 'log2FC_day21_vs_lib', bars = TRUE)
SID_T47D_E2_sig <- SID_T47D_E2[FDR_day21_vs_lib < 0.1]
SID_T47D_E2_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_T47D_E2_sig))


xD <- 0.3
dmrs_d <- fread(paste0(IN_DIR, "/20251003_EGFR_TM_dmrs_d_delt_0.4.csv"))
colnames(dmrs_d)[6:7] = c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d_NSCLCgt <- gTrack(dt2gr(dmrs_d[diff.Methy < -xD]))
dmrs_d_NSCLCgt$name = 'DMRsd0.3_NSCLC'
dmrs_d_BRCAgt <- gTrack(dt2gr(dmrs_d[diff.Methy > xD]))
dmrs_d_BRCAgt$name = 'DMRsd0.3_BRCA'

v4C_EGFR_NSCLC_fm_qVal <- readRDS(file.path(IN_DIR, '20260702_comb_Pvals_v4C_EGFR_NSCLC_gr.rds'))
v4C_EGFR_NSCLC_fm_qVal_gt <- gTrack(v4C_EGFR_NSCLC_fm_qVal, y.field = 'nlg_comb_q', name = 'fm_nlg10q', col = 'blue', y0=0, bars=TRUE, y1=10)
v4C_EGFR_BRCA_fm_qVal <- readRDS(file.path(IN_DIR, '20260702_comb_Pvals_v4C_EGFR_BRCA_gr.rds'))
v4C_EGFR_BRCA_fm_qVal_gt <- gTrack(v4C_EGFR_BRCA_fm_qVal, y.field = 'nlg_comb_q', name = 'fm_nlg10q', col = 'red', y0=0, bars=TRUE, y1=10)

ATACfisher <- fread(paste0(IN_DIR, "/20260623q_", xG, "_H3K27ac_ATAC_linDiffs_stats.csv"))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = '[:]'))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[3]), by = peakID]
ATACfisher_gr <- dt2gr(ATACfisher)
ATACfisher_t2tgr <- grl.unlist(liftOver(gr.chr(ATACfisher_gr), conversion_chain))
ATACfisher_LUSC_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.1 & tstdttype == 'LUSC')))
ATACfisher_BRCA_HR_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.1 & tstdttype == 'BRCA_HR+')))


window_oiz = window_oi -8e5
window_oiz2 = window_oi - 6e5
print('plotting')

datafield <- 'counts'

read_v4C <- function(sampleid, datafield) {
    readRDS(paste0(IN_DIR, '/individuals/', sampleid,
                   '_covCovariate_v4C_', xG, '_', datafield, '.rds'))
}

# little significance bar drawn under each 4C track:
# grey where nlg_p > 2, black where nlg_q > 1 (black wins where both hold)
v4C_sig_gt <- function(sampleid) {
    sig_p <- read_v4C(sampleid, 'nlg_p') %Q% (!is.na(value) & value > 2)
    sig_q <- read_v4C(sampleid, 'nlg_q') %Q% (!is.na(value) & value > 1)
    mcols(sig_p) <- DataFrame(sig = rep('nlg_p', length(sig_p)))
    mcols(sig_q) <- DataFrame(sig = rep('nlg_q', length(sig_q)))
    sig_gr <- sort(c(sig_p[!(sig_p %^% sig_q)], sig_q))
    gTrack(sig_gr, gr.colorfield = 'sig',
           colormaps = list(sig = c(nlg_p = 'grey60', nlg_q = 'black')),
           name = paste0(sampleid, '_sig'), height = 0.5)
}

H1HGXBV_v4C_gr <- read_v4C('1HGXBV', datafield) # CN19 # hg38 CN15 # circwalk
H1WNRWL_v4C_gr <- read_v4C('1WNRWL', datafield) # CN15 bampass # CN14 #
AJZ5JO_v4C_gr <- read_v4C('AJZ5JO', datafield) # CN49 bampass # CN39 # CN40 in hg38
XO7JTV_v4C_gr <- read_v4C('XO7JTV', datafield) # CN26 bampass # CN23 # CN20 in hg38 # circwalk

C109_v4C_gr <- read_v4C('C109', datafield) # CN70+ segmentation off, # CN160 in hg38
C15_v4C_gr <- read_v4C('C15', datafield) # CN35 # CN40 in hg38
C168_v4C_gr <- read_v4C('C168', datafield) # CN140 # CN150 in hg38 # circwalk
C208_v4C_gr <- read_v4C('C208', datafield) # CN23 # CN29 in hg38 # circwalk

y1_val <- 60

H1HGXBV_v4C <- gTrack(H1HGXBV_v4C_gr, y.field = 'value', bars = T, y0 = 0, y1 = y1_val, name = 'H1HGXBV', col = 'blue')
H1WNRWL_v4C <- gTrack(H1WNRWL_v4C_gr, y.field = 'value', bars = T, y0 = 0, y1 = y1_val, name = 'H1WNRWL', col = 'blue')
AJZ5JO_v4C <- gTrack(AJZ5JO_v4C_gr, y.field = 'value', bars = T, y0 = 0, y1 = 80, name = 'HAJZ5JO', col = 'blue')
XO7JTV_v4C <- gTrack(XO7JTV_v4C_gr, y.field = 'value', bars = T, y0 = 0, y1 = y1_val, name = 'HXO7JTV', col = 'blue')

C109_v4C <- gTrack(C109_v4C_gr, y.field = 'value', bars = T, y0 = 0, y1 = 560, name = 'C109', col = 'red')
C15_v4C <- gTrack(C15_v4C_gr, y.field = 'value', bars = T, y0 = 0, y1 = 20, name = 'C15', col = 'red')
C168_v4C <- gTrack(C168_v4C_gr, y.field = 'value', bars = T, y0 = 0, y1 = 180, name = 'C168', col = 'red')
C208_v4C <- gTrack(C208_v4C_gr, y.field = 'value', bars = T, y0 = 0, y1 = 180, name = 'C208', col = 'red')

H1HGXBV_sig <- v4C_sig_gt('1HGXBV')
H1WNRWL_sig <- v4C_sig_gt('1WNRWL')
AJZ5JO_sig <- v4C_sig_gt('AJZ5JO')
XO7JTV_sig <- v4C_sig_gt('XO7JTV')

C109_sig <- v4C_sig_gt('C109')
C15_sig <- v4C_sig_gt('C15')
C168_sig <- v4C_sig_gt('C168')
C208_sig <- v4C_sig_gt('C208')

l2fc_gt <- gmat_EGFR_10prct_diff_l2fc$gtrack(name=paste0("log2(FC)"),height=120,clim=c(-2,2), colormap = c('blue', 'white', 'red'))
nlgq_gt <- gmat_EGFR_10prct_diff_nlgq$gtrack(name=paste0("signed -log10(qvalue)"),height=120,clim=c(-2,2), colormap = c('blue', 'white', 'red'))

NSCLC_tracks <- c(gencode_T2T_gt, ATACfisher_LUSC_gt, dmrs_d_NSCLCgt, v4C_EGFR_NSCLC_fm_qVal_gt,
    H1HGXBV_sig, H1HGXBV_v4C, H1WNRWL_sig, H1WNRWL_v4C, AJZ5JO_sig, AJZ5JO_v4C, XO7JTV_sig, XO7JTV_v4C,
    NSCLC_xGtad_per_amp_gt,
    BRCAvsNSCLC_EGFR_chr7_sig_track_cons80_gt)
functional_tracks <- c(GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt,
    SID_T47D_E2_sig_prolifd21_q0.1_gt, SID_T47D_E2_prolifd21_gt,
    SID_MCF7_E2_sig_prolifd21_q0.1_gt, SID_MCF7_E2_prolifd21_gt)
BRCA_tracks <- c(ATACfisher_BRCA_HR_gt, dmrs_d_BRCAgt, v4C_EGFR_BRCA_fm_qVal_gt,
    C109_sig, C109_v4C, C15_sig, C15_v4C, C168_sig, C168_v4C, C208_sig, C208_v4C,
    BRCA_xGtad_per_amp_gt,
    xG_hicdcdiff_4C_nlgq_LUSC_gt, xG_hicdcdiff_4C_l2fc_gt, xG_hicdcdiff_4C_nlgq_BRCA_gt)

pdf(file =paste0(outdir, 'FigS8B_vis_hicdcdiff_', datafield,'_sigbars_EGFR_diff.pdf'),paper = 'a4', height = 48)
    plot(c(NSCLC_tracks, BRCA_tracks, l2fc_gt), window_oiz2, cex.label = 0.1)
    plot(c(NSCLC_tracks, functional_tracks, BRCA_tracks, l2fc_gt), window_oiz, cex.label = 0.1)
    plot(c(NSCLC_tracks, functional_tracks, BRCA_tracks, l2fc_gt), window_oiz2, cex.label = 0.1)
    plot(c(NSCLC_tracks, functional_tracks, BRCA_tracks, nlgq_gt), window_oiz, cex.label = 0.1)
    plot(c(NSCLC_tracks, functional_tracks, BRCA_tracks, nlgq_gt), window_oiz2, cex.label = 0.1)
dev.off()
