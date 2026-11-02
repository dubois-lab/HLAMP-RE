# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(GxG)
library(gGnome)
library(rtracklayer)
conversion_chain <- import.chain(file.path(IN_DIR, 'hg38ToHs1.over.chain'))

gencode_T2T_gt <- readRDS(file.path(IN_DIR, 'gencodeTrack_t2t.rds'))
CGC_gr_T2T <- readRDS('../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')
xG = 'NSD3'
Oncgr <- CGC_gr_T2T %Q%  (Gene_Symbol == xG)
binsize =5e3

xT_bins <- fread(file.path(IN_DIR, 'normalized_abs.bed'))
colnames(xT_bins) <- c('chr', 'start', 'end')
xT_bins[, start := start +1, by = start]

# NSD3 10%
xG_hicdcdiff <- fread(paste0(IN_DIR, '/FGFR1_hicdcdiff_breastVSlung_fitTypeLocal_10prctSigdiff_resbreastoverlung_chr8.txt'))
xG_hicdcdiff[, nlg_padj := -log10(padj)]
xG_hicdcdiff[log2FoldChange > 0, nlg_padj_sgn := nlg_padj]
xG_hicdcdiff[log2FoldChange < 0, nlg_padj_sgn := -nlg_padj]
xG_hicdcdiff[, nlg_p := -log10(pvalue)]
xG_hicdcdiff[log2FoldChange > 0, nlg_p_sgn := nlg_p]
xG_hicdcdiff[log2FoldChange < 0, nlg_p_sgn := -nlg_p]
xG_hicdcdiff[, i := (startI/ binsize) + 1, by = startI]
xG_hicdcdiff[, j := (startJ/ binsize) + 1, by = startJ]

gmat_NSD3_10prct_diff_nlgq <- gM(gr=dt2gr(xT_bins[chr== 'chr8']),dat=xG_hicdcdiff[, .(i, j, value=nlg_padj_sgn)],field="nlg_padj_sgn",fill=0,agg.fun=min)
gmat_NSD3_10prct_diff_l2fc <- gM(gr=dt2gr(xT_bins[chr== 'chr8']),dat=xG_hicdcdiff[, .(i, j, value=log2FoldChange)],field="log2FoldChange",fill=0,agg.fun=min)

NSD3_TSSgr <- readRDS(file.path(IN_DIR, 'NSD3_TSSgr_T2T.rds'))

xG_hicdcdiff_4C_nlgq <- gmat_NSD3_10prct_diff_nlgq[gmat_NSD3_10prct_diff_nlgq$gr %^% NSD3_TSSgr, drop =T]
xG_hicdcdiff_4C_nlgq_BRCA <- xG_hicdcdiff_4C_nlgq %Q% (!is.na(value) & value > 1)
xG_hicdcdiff_4C_nlgq_BRCA_gt <- gTrack(data = xG_hicdcdiff_4C_nlgq_BRCA, name = '4C_hicdcdiff_BRCA', height = 0.5, col = "red")
xG_hicdcdiff_4C_nlgq_LUSC<- xG_hicdcdiff_4C_nlgq %Q% (!is.na(value) & value < -1)
xG_hicdcdiff_4C_nlgq_LUSC_gt <- gTrack(data = xG_hicdcdiff_4C_nlgq_LUSC, name = '4C_hicdcdiff_LUSC', height = 0.5, col = "blue")

xG_hicdcdiff_4C_l2fc <- gmat_NSD3_10prct_diff_l2fc[gmat_NSD3_10prct_diff_l2fc$gr %^% NSD3_TSSgr, drop =T]
xG_hicdcdiff_4C_l2fc_gt <- gTrack(data = xG_hicdcdiff_4C_l2fc, y.field = 'value', bars = T, name = '4C_hicdcdiff_BRCA', y0 = -2, y1 = 2)

BRCA_NSD3_CNMatm <- readRDS(file.path(IN_DIR, "20260910_BRCA_NSD3_CNMatm_combined_hs1.rds"))
NSCLC_NSD3_CNMatm <- readRDS(file.path(IN_DIR, "20260910_NSCLC_LUSC_NSD3_CNMatm_combined_hs1.rds"))
BRCA_xGtad_per_amp_gt <- gTrack(dt2gr(BRCA_NSD3_CNMatm), y.field = 'AmpPercent', col = 'brown3', bars = F, y1=1, y0=0)
NSCLC_xGtad_per_amp_gt <- gTrack(dt2gr(NSCLC_NSD3_CNMatm), y.field = 'AmpPercent', col = 'royalblue3', bars = F, y1=1, y0=0)
window_oi <- reduce(dt2gr(rbind(BRCA_NSD3_CNMatm[AmpPercent > 0.8, 1:10], NSCLC_NSD3_CNMatm[AmpPercent > 0.8, 1:10])))

# BRCA vs NSCLC NSD3 differential amp-frequency significance track (chr8)
BRCAvsNSCLC_NSD3_fisher_gr <- readRDS(file.path(IN_DIR, "20260910_BRCAvsNSCLC_NSD3_chr8_FisherStats_CN10_cons80_gr.rds"))
BRCAvsNSCLC_NSD3_sig <- BRCAvsNSCLC_NSD3_fisher_gr[BRCAvsNSCLC_NSD3_fisher_gr$Fisher_pval < 0.05]
mcols(BRCAvsNSCLC_NSD3_sig)$tissue <- ifelse(mcols(BRCAvsNSCLC_NSD3_sig)$BRCA_AmpFreq > mcols(BRCAvsNSCLC_NSD3_sig)$NSCLC_AmpFreq, "BRCA", "NSCLC")
colormap_BRCAvsNSCLC_NSD3 <- list(tissue = c(NSCLC = "royalblue3", BRCA = "brown3"))
BRCAvsNSCLC_xG_sig_track_cons80_gt <- gTrack(BRCAvsNSCLC_NSD3_sig, gr.colorfield = "color", colormaps = colormap_BRCAvsNSCLC_NSD3, height = 0.5)

# NSCLC NSD3 left-vs-right (amp boundary) significance track
NSCLC_NSD3_LvsR_fisher <- fread(file.path(IN_DIR, "20260910_NSCLC_NSD3_LeftVsRight_FisherStats_CN10_cons80.txt"))
NSCLC_NSD3_LvsR_fisher_gr <- dt2gr(NSCLC_NSD3_LvsR_fisher[, ':='(seqnames = fifelse(right_AmpFreq > left_AmpFreq, right_chr, left_chr),
                                             start = fifelse(right_AmpFreq > left_AmpFreq, right_start, left_start),
                                             end = fifelse(right_AmpFreq > left_AmpFreq, right_end, left_end)
)])
NSCLC_NSD3_LvsR_sig <- NSCLC_NSD3_LvsR_fisher_gr[NSCLC_NSD3_LvsR_fisher_gr$Fisher_pval < 0.05]
NSCLC_xG_LvsR_sig_LR_gt <- gTrack(NSCLC_NSD3_LvsR_sig, col = "royalblue3", height = 0.5)

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
dmrs_d <- fread(paste0(IN_DIR, "/20260615_", xG, "_TM_dmrs.csv"))
colnames(dmrs_d)[6:7] = c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d_NSCLCgt <- gTrack(dt2gr(dmrs_d[diff.Methy < -xD]))
dmrs_d_NSCLCgt$name = 'DMRsd0.3_NSCLC'
dmrs_d_BRCAgt <- gTrack(dt2gr(dmrs_d[diff.Methy > xD]))
dmrs_d_BRCAgt$name = 'DMRsd0.3_BRCA'

v4C_NSD3_LUSC_fm_qVal <- readRDS(file.path(IN_DIR, '20260702_comb_Pvals_v4C_NSD3_LUSC_gr.rds'))
v4C_NSD3_LUSC_fm_qVal_gt <- gTrack(v4C_NSD3_LUSC_fm_qVal, y.field = 'nlg_comb_q', name = 'fm_nlg10q', col = 'blue', y0=0, bars=TRUE, y1=5)
v4C_NSD3_breast_fm_qVal <- readRDS(file.path(IN_DIR, '20260702_comb_Pvals_v4C_NSD3_BRCA_gr.rds'))
v4C_NSD3_breast_fm_qVal_gt <- gTrack(v4C_NSD3_breast_fm_qVal, y.field = 'nlg_comb_q', name = 'fm_nlg10q', col = 'red', y0=0, bars=TRUE, y1=5)

ATACfisher <- fread(paste0(IN_DIR, "/20260619_", xG, "_H3K27ac_ATAC_linDiffs_stats.csv"))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = '[:]'))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[3]), by = peakID]
ATACfisher_gr <- dt2gr(ATACfisher)
ATACfisher_t2tgr <- grl.unlist(liftOver(gr.chr(ATACfisher_gr), conversion_chain))
ATACfisher_LUSC_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.1 & tstdttype == 'LUSC')))
ATACfisher_BRCA_HR_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.1 & tstdttype == 'BRCA_HR+')))


window_oi2 = window_oi + 1e5
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

H1IDEVD_v4C_gr <- read_v4C('1IDEVD', datafield)
H1E2AVA_v4C_gr <- read_v4C('1E2AVA', datafield)
NODTWC_v4C_gr <- read_v4C('NODTWC', datafield)
K0HCJ3_v4C_gr <- read_v4C('K0HCJ3', datafield)

C12_v4C_gr <- read_v4C('C12', datafield)
C132_v4C_gr <- read_v4C('C132', datafield)
C17_v4C_gr <- read_v4C('C17', datafield)
Cb5_v4C_gr <- read_v4C('Cb5', datafield)

H1IDEVD_v4C <- gTrack(H1IDEVD_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'H1IDEVD', col = 'blue')#, y1 = 1.5)
H1E2AVA_v4C <- gTrack(H1E2AVA_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'H1E2AVA', col = 'blue')#, y1 = 1.5)
NODTWC_v4C <- gTrack(NODTWC_v4C_gr, y.field = 'value', bars = T, y0 = 0,name = 'NODTWC', col = 'blue')#, y1 = 1.5)
K0HCJ3_v4C <- gTrack(K0HCJ3_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'K0HCJ3', col = 'blue')#, y1 = 1.5)

C12_v4C <- gTrack(C12_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C12', col = 'red')#, y1 = 1.5)
C132_v4C <- gTrack(C132_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C132', col = 'red')#, y1 = 1.5)
C17_v4C <- gTrack(C17_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C17', col = 'red')#, y1 = 1.5)
Cb5_v4C <- gTrack(Cb5_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'Cb5', col = 'red')#, y1 = 1.5)

H1IDEVD_sig <- v4C_sig_gt('1IDEVD')
H1E2AVA_sig <- v4C_sig_gt('1E2AVA')
NODTWC_sig <- v4C_sig_gt('NODTWC')
K0HCJ3_sig <- v4C_sig_gt('K0HCJ3')

C12_sig <- v4C_sig_gt('C12')
C132_sig <- v4C_sig_gt('C132')
C17_sig <- v4C_sig_gt('C17')
Cb5_sig <- v4C_sig_gt('Cb5')

C12_v4C$y1 = 30
C132_v4C$y1 = 100
C17_v4C$y1 = 100
Cb5_v4C$y1 = 100
H1IDEVD_v4C$y1 = 100
H1E2AVA_v4C$y1 = 100
NODTWC_v4C$y1 = 100
K0HCJ3_v4C$y1 = 100

pdf(file =paste0(OUT_DIR, '/FigS8A_vis_hicdcdiff_', datafield,'_sigbars_NSD3_diff.pdf'),paper = 'a4', height = 48)
    plot(c(gencode_T2T_gt, ATACfisher_LUSC_gt, dmrs_d_NSCLCgt, v4C_NSD3_LUSC_fm_qVal_gt,
        H1IDEVD_sig, H1IDEVD_v4C, H1E2AVA_sig, H1E2AVA_v4C, NODTWC_sig, NODTWC_v4C, K0HCJ3_sig, K0HCJ3_v4C, 
        NSCLC_xGtad_per_amp_gt, NSCLC_xG_LvsR_sig_LR_gt,
        BRCAvsNSCLC_xG_sig_track_cons80_gt,
        ATACfisher_BRCA_HR_gt, dmrs_d_BRCAgt, v4C_NSD3_breast_fm_qVal_gt,
        C12_sig, C12_v4C, C17_sig, C17_v4C, C132_sig, C132_v4C, Cb5_sig, Cb5_v4C, 
        BRCA_xGtad_per_amp_gt,
        xG_hicdcdiff_4C_nlgq_LUSC_gt, xG_hicdcdiff_4C_l2fc_gt, xG_hicdcdiff_4C_nlgq_BRCA_gt,
        gmat_NSD3_10prct_diff_l2fc$gtrack(name=paste0("log2(FC)"),height=120,clim=c(-2,2), colormap = c('blue', 'white', 'red'))),
        window_oi, cex.label = 0.1)
    # zoomed out
    plot(c(gencode_T2T_gt, ATACfisher_LUSC_gt, dmrs_d_NSCLCgt, v4C_NSD3_LUSC_fm_qVal_gt,
        H1IDEVD_sig, H1IDEVD_v4C, H1E2AVA_sig, H1E2AVA_v4C, NODTWC_sig, NODTWC_v4C, K0HCJ3_sig, K0HCJ3_v4C, 
        NSCLC_xGtad_per_amp_gt, NSCLC_xG_LvsR_sig_LR_gt,
        BRCAvsNSCLC_xG_sig_track_cons80_gt,
        ATACfisher_BRCA_HR_gt, dmrs_d_BRCAgt, v4C_NSD3_breast_fm_qVal_gt,
        C12_sig, C12_v4C, C17_sig, C17_v4C, C132_sig, C132_v4C, Cb5_sig, Cb5_v4C, 
        BRCA_xGtad_per_amp_gt,
        xG_hicdcdiff_4C_nlgq_LUSC_gt, xG_hicdcdiff_4C_l2fc_gt, xG_hicdcdiff_4C_nlgq_BRCA_gt,
        gmat_NSD3_10prct_diff_l2fc$gtrack(name=paste0("log2(FC)"),height=120,clim=c(-2,2), colormap = c('blue', 'white', 'red'))),
        window_oi2, cex.label = 0.1)
dev.off()
