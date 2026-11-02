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
xG = 'CCND1'
Oncgr <- CGC_gr_T2T %Q%  (Gene_Symbol == xG)
window_oi = GRanges("chr11:69000000-70000000")
binsize =5e3

xT_bins <- fread(file.path(IN_DIR, 'normalized_abs.bed'))
colnames(xT_bins) <- c('chr', 'start', 'end')
xT_bins[, start := start +1, by = start]

# CCND1 10%
xG_hicdcdiff <- fread(paste0(IN_DIR, '/CCND1_hicdcdiff_breastVSlung_fitTypeLocal_10prctSigdiff_resbreastoverlung_chr11.txt'))
xG_hicdcdiff[, nlg_padj := -log10(padj)]
xG_hicdcdiff[log2FoldChange > 0, nlg_padj_sgn := nlg_padj]
xG_hicdcdiff[log2FoldChange < 0, nlg_padj_sgn := -nlg_padj]
xG_hicdcdiff[, nlg_p := -log10(pvalue)]
xG_hicdcdiff[log2FoldChange > 0, nlg_p_sgn := nlg_p]
xG_hicdcdiff[log2FoldChange < 0, nlg_p_sgn := -nlg_p]
xG_hicdcdiff[, i := (startI/ binsize) + 1, by = startI]
xG_hicdcdiff[, j := (startJ/ binsize) + 1, by = startJ]

gmat_CCND1_10prct_diff_nlgq <- gM(gr=dt2gr(xT_bins[chr== 'chr11']),dat=xG_hicdcdiff[, .(i, j, value=nlg_padj_sgn)],field="nlg_padj_sgn",fill=0,agg.fun=min)
gmat_CCND1_10prct_diff_l2fc <- gM(gr=dt2gr(xT_bins[chr== 'chr11']),dat=xG_hicdcdiff[, .(i, j, value=log2FoldChange)],field="log2FoldChange",fill=0,agg.fun=min)

CCND1_TSSgr <- readRDS(file.path(IN_DIR, 'CCND1_TSSgr_T2T.rds'))

xG_hicdcdiff_4C_nlgq <- gmat_CCND1_10prct_diff_nlgq[gmat_CCND1_10prct_diff_nlgq$gr %^% CCND1_TSSgr, drop =T]
xG_hicdcdiff_4C_nlgq_BRCA <- xG_hicdcdiff_4C_nlgq %Q% (!is.na(value) & value > 1.3)
xG_hicdcdiff_4C_nlgq_BRCA_gt <- gTrack(data = xG_hicdcdiff_4C_nlgq_BRCA, name = '4C_hicdcdiff_BRCA', height = 0.5, color = "red")
xG_hicdcdiff_4C_nlgq_LUSC<- xG_hicdcdiff_4C_nlgq %Q% (!is.na(value) & value < -0.7)
xG_hicdcdiff_4C_nlgq_LUSC_gt <- gTrack(data = xG_hicdcdiff_4C_nlgq_LUSC, name = '4C_hicdcdiff_LUSC', height = 0.5, color = "blue")

xG_hicdcdiff_4C_l2fc <- gmat_CCND1_10prct_diff_l2fc[gmat_CCND1_10prct_diff_l2fc$gr %^% CCND1_TSSgr, drop =T]
xG_hicdcdiff_4C_l2fc_gt <- gTrack(data = xG_hicdcdiff_4C_l2fc, y.field = 'value', bars = T, name = '4C_hicdcdiff_BRCA', y0 = -1, y1 = 2)

NSCLC_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, '/20260907_Fig.2C_NSCLC_CCND1_AmpFreq_track.rds'))
NSCLC_xGtad_per_amp_gt$bars =F
BRCA_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, '/20260907_Fig.2C_BRCA_CCND1_AmpFreq_track.rds'))
BRCA_xGtad_per_amp_gt$bars =F
BRCA_CCND1_LvsR_sig_LR_gt <- readRDS(paste0(IN_DIR, '/20260907c_Fig.2C_BRCA_CCND1_LvsR_sig_LR_cons80.rds'))
BRCA_CCND1_LvsR_sig_LR_gt <- gTrack(data = BRCA_CCND1_LvsR_sig_LR_gt@data[[1]], col = "red", height = 0.5)
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- readRDS(paste0(IN_DIR, '/20260907c_Fig.2C_BRCAvsNSCLC_CCND1_chr11_sig_track_cons80.rds'))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- gTrack(data = BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt@data[[1]], col = "blue", height = 0.5)


# #functional val####
GERCARres <- as.data.table(read.csv2(file.path(IN_DIR, 'GERCARres.csv')))
GERCARres_gr <- dt2gr(GERCARres[chromosome == unique(as.character(seqnames(Oncgr)))])
GERCARres_gr_t2tgr <- grl.unlist(liftOver(gr.chr(GERCARres_gr), conversion_chain))
GERCARres_prolifd21_gt <- gTrack(GERCARres_gr_t2tgr, y.field = 'logFC_prolif_21', bars = T)
GERCARres_gr_sig <- dt2gr(gr2dt(GERCARres_gr_t2tgr)[FDR_prolif_21 < 0.1 ])
GERCARres_prolifd21_q0.1_gt <- gTrack(GERCARres_gr_sig)

dmrs_d4 <- fread(paste0(IN_DIR, '/20260512_CCND1_TM_dmrs_d.csv'))
colnames(dmrs_d4)[6:7] = c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d4_NSCLCgt <- gTrack(dt2gr(dmrs_d4[diff.Methy < -0.4]), ygap =0.5)
dmrs_d4_NSCLCgt$name = 'DMRsd0.4_NSCLC'
dmrs_d4_BRCAgt <- gTrack(dt2gr(dmrs_d4[diff.Methy > 0.4]), ygap =0.5)
dmrs_d4_BRCAgt$name = 'DMRsd0.4_BRCA'
dmrs_d5_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -0.5])), ygap =0.5)
dmrs_d5_BRCAgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy > 0.5])), ygap =0.5)
dmrs_d6_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -0.6])), ygap =0.5)
dmrs_d6_BRCAgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy > 0.6])), ygap =0.5)


v4C_CCND1_LUSC_fm_qVal <- readRDS(file.path(IN_DIR, '20260702_comb_Pvals_v4C_CCND1_LUSC_gr.rds'))
v4C_CCND1_LUSC_fm_qVal_gt <- gTrack(v4C_CCND1_LUSC_fm_qVal, y.field = 'nlg_comb_q', name = 'fm_nlg10q', col = 'blue', y0=0, bars=TRUE, y1=15)
v4C_CCND1_breast_fm_qVal <- readRDS(file.path(IN_DIR, '20260702_comb_Pvals_v4C_CCND1_BRCA_gr.rds'))
v4C_CCND1_breast_fm_qVal_gt <- gTrack(v4C_CCND1_breast_fm_qVal, y.field = 'nlg_comb_q', name = 'fm_nlg10q', col = 'red', y0=0, bars=TRUE, y1=15)

ATACfisher <- fread(file.path(IN_DIR, '20260614_CCND1_H3K27ac&ATAC_linDiffs_stats.csv'))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = '[:]'))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[3]), by = peakID]
ATACfisher_gr <- dt2gr(ATACfisher)
ATACfisher_t2tgr <- grl.unlist(liftOver(gr.chr(ATACfisher_gr), conversion_chain))
ATACfisher_LUSC_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.05 & tstdttype == 'LUSC')), ygap = 0.5)
ATACfisher_BRCA_HR_gt <- gTrack(data = reduce(ATACfisher_t2tgr %Q% (q_val < 0.05 & tstdttype == 'BRCA_HR+')), ygap = 0.5)


window_oiz = window_oi -2e5
window_oiz2 = window_oi - 4e5
print('plotting')

lapply(c('counts', 'nlg_q'), function(datafield) {

    H3D44S_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','H3D44S','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    KDAMIF_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','KDAMIF','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    Y57U8R_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','Y57U8R','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    F57J9F_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','F57J9F','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    
    C3_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','C3','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    C169_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','C169','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    C17_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','C17','_covCovariate_v4C_', xG, '_',datafield,'.rds'))

    H3D44S_v4C <- gTrack(H3D44S_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'H3D44S', col = 'blue')#, y1 = 1.5)
    KDAMIF_v4C <- gTrack(KDAMIF_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'KDAMIF', col = 'blue')#, y1 = 1.5)
    Y57U8R_v4C <- gTrack(Y57U8R_v4C_gr, y.field = 'value', bars = T, y0 = 0,name = 'Y57U8R', col = 'blue')#, y1 = 1.5)
    F57J9F_v4C <- gTrack(F57J9F_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'F57J9F', col = 'blue')#, y1 = 1.5)
    
    C3_v4C <- gTrack(C3_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C3', col = 'red')#, y1 = 1.5)
    C169_v4C <- gTrack(C169_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C169', col = 'red')#, y1 = 1.5)
    C17_v4C <- gTrack(C17_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C17', col = 'red')#, y1 = 1.5)

        if(datafield == 'nlg_q') {
        C3_v4C$y1 = 3
        C169_v4C$y1 = 3
        C17_v4C$y1 = 3
        H3D44S_v4C$y1 = 3
        KDAMIF_v4C$y1 = 3
        Y57U8R_v4C$y1 = 3
        F57J9F_v4C$y1 = 3
    }
    pdf(file =paste0(OUT_DIR, '/Fig4A_vis_hicdcdiff_', datafield,'_CCND1_diff_gercar.pdf'),paper = 'a4', height = 48)
        plot(c(gencode_T2T_gt, ATACfisher_LUSC_gt, dmrs_d4_NSCLCgt, dmrs_d5_NSCLCgt, dmrs_d6_NSCLCgt, v4C_CCND1_LUSC_fm_qVal_gt,
            H3D44S_v4C, KDAMIF_v4C, Y57U8R_v4C, F57J9F_v4C, NSCLC_xGtad_per_amp_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
            GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt, dmrs_d4_BRCAgt, dmrs_d5_BRCAgt, dmrs_d6_BRCAgt, v4C_CCND1_breast_fm_qVal_gt,
            C3_v4C, C17_v4C, C169_v4C, BRCA_xGtad_per_amp_gt,BRCA_CCND1_LvsR_sig_LR_gt,
            xG_hicdcdiff_4C_nlgq_LUSC_gt, xG_hicdcdiff_4C_l2fc_gt, xG_hicdcdiff_4C_nlgq_BRCA_gt,
            gmat_CCND1_10prct_diff_l2fc$gtrack(name=paste0("log2(FC)"),height=120,clim=c(-2,2), colormap = c('blue', 'white', 'red'))),window_oi, cex.label = 0.1)
        plot(c(gencode_T2T_gt, ATACfisher_LUSC_gt, dmrs_d4_NSCLCgt, dmrs_d5_NSCLCgt, dmrs_d6_NSCLCgt, v4C_CCND1_LUSC_fm_qVal_gt,
            H3D44S_v4C, KDAMIF_v4C, Y57U8R_v4C, F57J9F_v4C, NSCLC_xGtad_per_amp_gt,BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
            GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt, dmrs_d4_BRCAgt, dmrs_d5_BRCAgt, dmrs_d6_BRCAgt, v4C_CCND1_breast_fm_qVal_gt,
            C3_v4C, C17_v4C, C169_v4C, BRCA_xGtad_per_amp_gt,BRCA_CCND1_LvsR_sig_LR_gt,
            xG_hicdcdiff_4C_nlgq_LUSC_gt, xG_hicdcdiff_4C_l2fc_gt,xG_hicdcdiff_4C_nlgq_BRCA_gt,
            gmat_CCND1_10prct_diff_l2fc$gtrack(name=paste0("log2(FC)"),height=120,clim=c(-2,2), colormap = c('blue', 'white', 'red'))),window_oiz, cex.label = 0.1)
    dev.off()
})

