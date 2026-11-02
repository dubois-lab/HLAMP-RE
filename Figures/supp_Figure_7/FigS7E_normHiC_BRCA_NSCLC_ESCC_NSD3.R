# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(GxG)
library(gGnome)
library(rtracklayer)
hs1ToHg38_chain <- import.chain(file.path(IN_DIR, 'hs1ToHg38.over.chain')) # get the reverse

gencode_hg38_gt <- readRDS('../supp_Figure_3/FigS3A/source_data/gencode_hg38.rds')
CGC_gr_hg38 <- readRDS('../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
xG = 'NSD3'
Oncgr <- CGC_gr_hg38 %Q%  (Gene_Symbol == xG)
window_oi = gr.chr(Oncgr) + 10e5
binsize =5e3

# white y.grid for all data tracks
fmt_gt <- function(gt) {
    gt$y.grid.col= "white"
    gt
}

BRCA_NSD3_CNMatm <- readRDS(file.path(IN_DIR, "20260910_BRCA_NSD3_CNMatm_combined_hs1.rds"))
NSCLC_NSD3_CNMatm <- readRDS(file.path(IN_DIR, "20260910_NSCLC_LUSC_NSD3_CNMatm_combined_hs1.rds"))
BRCA_NSD3_CNMatm_hg38gr <- unlist(liftOver(gr.chr(dt2gr(BRCA_NSD3_CNMatm)), hs1ToHg38_chain))
NSCLC_NSD3_CNMatm_hg38gr <- unlist(liftOver(gr.chr(dt2gr(NSCLC_NSD3_CNMatm)), hs1ToHg38_chain))
BRCA_xGtad_per_amp_gt <- fmt_gt(gTrack(BRCA_NSD3_CNMatm_hg38gr, y.field = "AmpPercent", col = "brown3", bars = F, y0 = 0, y1 = 100))
NSCLC_xGtad_per_amp_gt <- fmt_gt(gTrack(NSCLC_NSD3_CNMatm_hg38gr, y.field = "AmpPercent", col = "royalblue3", bars = F, y0 = 0, y1 =100))

# BRCA vs NSCLC NSD3 differential amp-frequency significance track (chr8)
BRCAvsNSCLC_NSD3_fisher_gr <- readRDS(file.path(IN_DIR, "20260910_BRCAvsNSCLC_NSD3_chr8_FisherStats_CN10_cons80_gr.rds"))
BRCAvsNSCLC_NSD3_sig <- BRCAvsNSCLC_NSD3_fisher_gr[BRCAvsNSCLC_NSD3_fisher_gr$Fisher_pval < 0.05]
mcols(BRCAvsNSCLC_NSD3_sig)$tissue <- ifelse(mcols(BRCAvsNSCLC_NSD3_sig)$BRCA_AmpFreq > mcols(BRCAvsNSCLC_NSD3_sig)$NSCLC_AmpFreq, "BRCA", "NSCLC")
colormap_BRCAvsNSCLC_NSD3 <- list(tissue = c(NSCLC = "royalblue3", BRCA = "brown3"))
BRCAvsNSCLC_NSD3_sig_hg38gr <- unlist(liftOver(gr.chr(BRCAvsNSCLC_NSD3_sig), hs1ToHg38_chain))
BRCAvsNSCLC_xG_chr11_sig_track_cons80_gt <- gTrack(BRCAvsNSCLC_NSD3_sig_hg38gr, gr.colorfield = "color", colormaps = colormap_BRCAvsNSCLC_NSD3, height = 0.5)

# NSCLC NSD3 left-vs-right (amp boundary) significance track
NSCLC_NSD3_LvsR_fisher <- fread(file.path(IN_DIR, "20260910_NSCLC_NSD3_LeftVsRight_FisherStats_CN10_cons80.txt"))
NSCLC_NSD3_LvsR_fisher_gr <- dt2gr(NSCLC_NSD3_LvsR_fisher[, ':='(seqnames = fifelse(right_AmpFreq > left_AmpFreq, right_chr, left_chr),
                                             start = fifelse(right_AmpFreq > left_AmpFreq, right_start, left_start),
                                             end = fifelse(right_AmpFreq > left_AmpFreq, right_end, left_end)
)])
NSCLC_NSD3_LvsR_sig <- NSCLC_NSD3_LvsR_fisher_gr[NSCLC_NSD3_LvsR_fisher_gr$Fisher_pval < 0.05]
NSCLC_NSD3_LvsR_sig_hg38gr <- unlist(liftOver(gr.chr(NSCLC_NSD3_LvsR_sig), hs1ToHg38_chain))
NSCLC_xG_LvsR_sig_LR_gt <- gTrack(NSCLC_NSD3_LvsR_sig_hg38gr, col = "royalblue3", height = 0.5)

# #functional val####
SID_MCF7_E2 <- as.data.table(read.csv2(file.path(IN_DIR, 'cd-23-1157_supplementary_table_s3_suppst3_MCF7_E2.csv')))
SID_MCF7_E2_gr <- dt2gr(SID_MCF7_E2)
SID_MCF7_E2_prolifd21_gt <- gTrack(SID_MCF7_E2_gr, y.field = 'log2FC_day21_vs_day07', bars = TRUE)
SID_MCF7_E2_sig <- SID_MCF7_E2[FDR_day21_vs_day07 < 0.1]
SID_MCF7_E2_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_MCF7_E2_sig))

SID_T47D_E2 <- as.data.table(read.csv2(file.path(IN_DIR, 'cd-23-1157_supplementary_table_s4_suppst4_T47D_E2.csv')))
SID_T47D_E2_gr <- dt2gr(SID_T47D_E2)
SID_T47D_E2_prolifd21_gt <- gTrack(SID_T47D_E2_gr, y.field = 'log2FC_day21_vs_lib', bars = TRUE)
SID_T47D_E2_sig <- SID_T47D_E2[FDR_day21_vs_lib < 0.1]
SID_T47D_E2_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_T47D_E2_sig))

# sigdiff ####
ATACfisher <- fread(paste0(IN_DIR, "/20260619_", xG, "_H3K27ac_ATAC_linDiffs_stats.csv"))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = "[:]"))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[3]), by = peakID]
fishres_A2gr <- dt2gr(ATACfisher)
# fishres_A_sigs_LUAD <- fishres_A2gr %Q% c(q_val < 0.05 & tstdttype == "LUAD")
# fishres_sigs_LUADgt <- gTrack(data = reduce(fishres_A_sigs_LUAD))
fishres_A_sigs_BRCA_HR <- fishres_A2gr %Q% (q_val < 0.1 & tstdttype == "BRCA_HR+")
fishres_sigs_BRCA_HR_gt <- gTrack(data = reduce(fishres_A_sigs_BRCA_HR))
# fishres_A_sigs_BRCA_HRN <- fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "BRCA_HR-")
# fishres_sigs_BRCA_HRNgt <- gTrack(data = reduce(fishres_A_sigs_BRCA_HRN))
ATACfisher_LUSC_gt <- gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.1 & tstdttype == "LUSC")))
print('Fisher test readin done')
xD <- 0.3
dmrs_d <- fread(paste0(IN_DIR, "/20260615_", xG, "_TM_dmrs.csv"))
colnames(dmrs_d)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs_d)), hs1ToHg38_chain))
dmrs_d <- as.data.table(dmrs_d_gr)
# each DMR is shaded by the highest |diff.Methy| cutoff (xD / xD+0.1 / xD+0.2) it passes;
# categorical colormap as for the tissue-colored Fisher track (a continuous per-row grey() color column does not render in gTrack)
dmrs_d[, methyl_level := ifelse(abs(diff.Methy) > xD + 0.2, "dhi", ifelse(abs(diff.Methy) > xD + 0.1, "dmid", "dlo"))]
colormap_dmrs_d_grey <- list(methyl_level = c(dlo = grey(0.7), dmid = grey(0.4), dhi = grey(0.1)))
dmrs_d_NSCLCgt <- gTrack(dt2gr(dmrs_d[diff.Methy < -xD]), gr.colorfield = "color", colormaps = colormap_dmrs_d_grey)
dmrs_d_BRCAgt <- gTrack(dt2gr(dmrs_d[diff.Methy > xD]), gr.colorfield = "color", colormaps = colormap_dmrs_d_grey)

xG_TSS <- fread(paste0(IN_DIR, '/', xG, '_TSS_hg38.csv'))
colnames(xG_TSS)[4] = 'seqnames'
colnames(xG_TSS)[8] = 'start'
xG_TSS[, end := start +1]
xG_TSSgr <- dt2gr(xG_TSS)

print('reading in breast HiC data')
xT_map = straw(file.path(IN_DIR, 'ENCFF924WYB_BRCA_hg38_norm.hic'),  gr = window_oi, res = 5e3) #53f
gr_map_dt <- gr2dt(xT_map$gr)
# gr_map_dt
gr_map_dt[, seqnames := gsub(pattern = paste0(':', start(window_oi), ':', end(window_oi)), replacement = '', x = seqnames), by = seqnames]
# gr_map_dt
xT_map$gr <- dt2gr(gr_map_dt)
BRCA_hg38_norm_HiC_4C <- xT_map[xT_map$gr %^% gr.chr(xG_TSSgr), drop =T]
BRCA_hg38_norm_HiC_v4C_gt <- fmt_gt(gTrack(data = BRCA_hg38_norm_HiC_4C, y.field = 'value', bars = T, name = '4C_BRCA', y0 = 0, y1 = 15))
print('reading in lung HiC data done')
xT_map_Lu = straw(file.path(IN_DIR, 'ENCFF094KKX_lung16f_hg38.hic'),  gr = window_oi, res = 5e3) #53f
gr_map_LU_dt <- gr2dt(xT_map_Lu$gr)
# gr_map_dt
gr_map_LU_dt[, seqnames := gsub(pattern = paste0(':', start(window_oi), ':', end(window_oi)), replacement = '', x = seqnames), by = seqnames]
# gr_map_dt
xT_map_Lu$gr <- dt2gr(gr_map_LU_dt)
lung16f_hg38_norm_HiC_4C <- xT_map_Lu[xT_map_Lu$gr %^% gr.chr(xG_TSSgr), drop =T]
lung16f_hg38_norm_HiC_v4C_gt <- fmt_gt(gTrack(data = lung16f_hg38_norm_HiC_4C, y.field = 'value', bars = T, name = '4C_lung16f', y0 = 0, y1 = 20))

print('reading in esophageal squamous HiC data')
xT_map_sqEsoph = straw(file.path(IN_DIR, 'ENCFF224HKR_squamousEpEsoph_hg38.hic'),  gr = window_oi, res = 5e3)
gr_map_sqEsoph_dt <- gr2dt(xT_map_sqEsoph$gr)
# gr_map_dt
gr_map_sqEsoph_dt[, seqnames := gsub(pattern = paste0(':', start(window_oi), ':', end(window_oi)), replacement = '', x = seqnames), by = seqnames]
# gr_map_dt
xT_map_sqEsoph$gr <- dt2gr(gr_map_sqEsoph_dt)
squamousEpEsoph_hg38_norm_HiC_4C <- xT_map_sqEsoph[xT_map_sqEsoph$gr %^% gr.chr(xG_TSSgr), drop =T]
squamousEpEsoph_hg38_norm_HiC_v4C_gt <- fmt_gt(gTrack(data = squamousEpEsoph_hg38_norm_HiC_4C, y.field = 'value', bars = T, name = '4C_squamousEpEsoph', y0 = 0, y1 = 20))
print('reading in HiC data done')

print('reading in breast ATAC wig data')
BRCA_atacFE_median_gt <- fmt_gt(gTrack(file.path(IN_DIR, 'BRCA_median.bw'), y.field = "score", name = "BRCA_atacFE_median", bars = T, y0 = 0, y1 = 10))
print('reading in breast ATAC bed data')
BRCA51f_atac_bed <- fread(file.path(IN_DIR, 'ENCFF544RPU_BRCA51f_hg38_atac_narrowPeaks.bed'))
colnames(BRCA51f_atac_bed)[1:3] <- c("seqnames", "start", "end")
BRCA51f_atac_bed$seqnames <- gsub(pattern = "chr", replacement = "", BRCA51f_atac_bed$seqnames)
BRCA51f_atac_bed_gr <- dt2gr(BRCA51f_atac_bed[V9 > -log10(0.05)])
BRCA51f_atac_bed_gt <- gTrack(data = reduce(BRCA51f_atac_bed_gr), name = 'BRCA51f_atac_peaks', height = 0.5)

BRCA53f_atac_bed <- fread(file.path(IN_DIR, 'ENCFF970MTO_BRCA53f_hg38_atac_narrowPeaks.bed'))
colnames(BRCA53f_atac_bed)[1:3] <- c("seqnames", "start", "end")
BRCA53f_atac_bed$seqnames <- gsub(pattern = "chr", replacement = "", BRCA53f_atac_bed$seqnames)
BRCA53f_atac_bed_gr <- dt2gr(BRCA53f_atac_bed[V9 > -log10(0.05)])
BRCA53f_atac_bed_gt <- gTrack(data = reduce(BRCA53f_atac_bed_gr), name = 'BRCA53f_atac_peaks', height = 0.5)

BRCA54f_atac_bed <- fread(file.path(IN_DIR, 'ENCFF738IZU_BRCA54f_hg38_atac_narrowPeaks.bed'))
colnames(BRCA54f_atac_bed)[1:3] <- c("seqnames", "start", "end")
BRCA54f_atac_bed$seqnames <- gsub(pattern = "chr", replacement = "", BRCA54f_atac_bed$seqnames)
BRCA54f_atac_bed_gr <- dt2gr(BRCA54f_atac_bed[V9 > -log10(0.05)])
BRCA54f_atac_bed_gt <- gTrack(data = reduce(BRCA54f_atac_bed_gr), name = 'BRCA54f_atac_peaks', height = 0.5)

BRCA51f_hg38_H3K4me1_p_gt = gTrack(file.path(IN_DIR, 'ENCFF901FSR_BRCA51f_hg38_H3K4me1_p.bigWig'), name = 'BRCA51f_hg38_H3K4me1_p',
    y.field = 'score', bars = T, y0 = 0, y1 = 10)

print('reading in lung ATAC wig data')
lung_atacFE_median_gt <- fmt_gt(gTrack(file.path(IN_DIR, 'lung_median.bw'), y.field = "score", name = "lung_atacFE_median", bars = T, , y0 = 0, y1 = 10))
lung16f_hg38_H3K4me1_p_gt = gTrack(file.path(IN_DIR, 'ENCFF235XOI_lung16f_hg38_H3K4me1_p.bigWig'), name = 'lung16f_hg38_H3K4me1_p',
    y.field = 'score', bars = T, y0 = 0, y1 = 10)
print('reading in lung ATAC bed data')
lung16f_atac_bed <- fread(file.path(IN_DIR, 'ENCFF138ISY_lung16f_hg38_ATAC_narrowPeaks.bed'))
colnames(lung16f_atac_bed)[1:3] <- c("seqnames", "start", "end")
lung16f_atac_bed$seqnames <- gsub(pattern = "chr", replacement = "", lung16f_atac_bed$seqnames)
lung16f_atac_bed_gr <- dt2gr(lung16f_atac_bed[V9 > -log10(0.05)])
lung16f_atac_bed_gt <- gTrack(data = reduce(lung16f_atac_bed_gr), name = 'lung16f_atac_peaks', height = 0.5)

lung65f_atac_bed <- fread(file.path(IN_DIR, 'ENCFF778HMH_lung65f_hg38_ATAC_narrowPeaks.bed'))
colnames(lung65f_atac_bed)[1:3] <- c("seqnames", "start", "end")
lung65f_atac_bed$seqnames <- gsub(pattern = "chr", replacement = "", lung65f_atac_bed$seqnames)
lung65f_atac_bed_gr <- dt2gr(lung65f_atac_bed[V9 > -log10(0.05)])
lung65f_atac_bed_gt <- gTrack(data = reduce(lung65f_atac_bed_gr), name = 'lung65f_atac_peaks', height = 0.5)

lung60m_atac_bed <- fread(file.path(IN_DIR, 'ENCFF047NIZ_lung60m_hg38_ATAC_narrowPeaks.bed'))
colnames(lung60m_atac_bed)[1:3] <- c("seqnames", "start", "end")
lung60m_atac_bed$seqnames <- gsub(pattern = "chr", replacement = "", lung60m_atac_bed$seqnames)
lung60m_atac_bed_gr <- dt2gr(lung60m_atac_bed[V9 > -log10(0.05)])
lung60m_atac_bed_gt <- gTrack(data = reduce(lung60m_atac_bed_gr), name = 'lung60m_atac_peaks', height = 0.5)

print('reading in squamous esophagus H3K27ac wig data')
squamousEsoph_median_gt <- fmt_gt(gTrack(file.path(IN_DIR, 'squamousEsoph_median.bw'), y.field = "score", name = "squamousEsoph_median", bars = T, y0 = 0, y1 = 10))
print('reading in squamous esophagus H3K27ac bed data')
squamousEsoph51f_H3K27ac_bed <- fread(file.path(IN_DIR, 'ENCFF673PCP_squamousEsoph51f_hg38_H3K27ac_FE.bed'))
colnames(squamousEsoph51f_H3K27ac_bed)[1:3] <- c("seqnames", "start", "end")
squamousEsoph51f_H3K27ac_bed$seqnames <- gsub(pattern = "chr", replacement = "", squamousEsoph51f_H3K27ac_bed$seqnames)
squamousEsoph51f_H3K27ac_bed_gr <- dt2gr(squamousEsoph51f_H3K27ac_bed[V9 > -log10(0.05)])
squamousEsoph51f_H3K27ac_bed_gt <- gTrack(data = reduce(squamousEsoph51f_H3K27ac_bed_gr), name = 'squamousEsoph51f_H3K27ac_peaks', height = 0.5)

squamousEsoph53f_H3K27ac_bed <- fread(file.path(IN_DIR, 'ENCFF057CRD_squamousEsoph53f_hg38_H3K27ac_FE.bed'))
colnames(squamousEsoph53f_H3K27ac_bed)[1:3] <- c("seqnames", "start", "end")
squamousEsoph53f_H3K27ac_bed$seqnames <- gsub(pattern = "chr", replacement = "", squamousEsoph53f_H3K27ac_bed$seqnames)
squamousEsoph53f_H3K27ac_bed_gr <- dt2gr(squamousEsoph53f_H3K27ac_bed[V9 > -log10(0.05)])
squamousEsoph53f_H3K27ac_bed_gt <- gTrack(data = reduce(squamousEsoph53f_H3K27ac_bed_gr), name = 'squamousEsoph53f_H3K27ac_peaks', height = 0.5)

squamousEsoph54m_H3K27ac_bed <- fread(file.path(IN_DIR, 'ENCFF253VVF_squamousEsoph54m_hg38_H3K27ac_FE.bed'))
colnames(squamousEsoph54m_H3K27ac_bed)[1:3] <- c("seqnames", "start", "end")
squamousEsoph54m_H3K27ac_bed$seqnames <- gsub(pattern = "chr", replacement = "", squamousEsoph54m_H3K27ac_bed$seqnames)
squamousEsoph54m_H3K27ac_bed_gr <- dt2gr(squamousEsoph54m_H3K27ac_bed[V9 > -log10(0.05)])
squamousEsoph54m_H3K27ac_bed_gt <- gTrack(data = reduce(squamousEsoph54m_H3K27ac_bed_gr), name = 'squamousEsoph54m_H3K27ac_peaks', height = 0.5)

print('reading in v4C data')
datafield='counts'
    C12_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','C12','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    C12_v4C_gr <- unlist(liftOver(gr.chr(C12_v4C_gr), hs1ToHg38_chain))
    C132_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','C132','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    C132_v4C_gr <- unlist(liftOver(gr.chr(C132_v4C_gr), hs1ToHg38_chain))
    C17_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','C17','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    C17_v4C_gr <- unlist(liftOver(gr.chr(C17_v4C_gr), hs1ToHg38_chain))
    Cb5_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','Cb5','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    Cb5_v4C_gr <- unlist(liftOver(gr.chr(Cb5_v4C_gr), hs1ToHg38_chain))
    C12_v4C <- fmt_gt(gTrack(C12_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C12', col = 'red', y1 = 30))
    C132_v4C <- fmt_gt(gTrack(C132_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C132', col = 'red', y1 = 100))
    C17_v4C <- fmt_gt(gTrack(C17_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C17', col = 'red', y1 = 100))
    Cb5_v4C <- fmt_gt(gTrack(Cb5_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'Cb5', col = 'red', y1 = 100))
    H1IDEVD_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','1IDEVD','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    H1IDEVD_v4C_gr <- unlist(liftOver(gr.chr(H1IDEVD_v4C_gr), hs1ToHg38_chain))
    H1E2AVA_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','1E2AVA','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    H1E2AVA_v4C_gr <- unlist(liftOver(gr.chr(H1E2AVA_v4C_gr), hs1ToHg38_chain))
    NODTWC_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','NODTWC','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    NODTWC_v4C_gr <- unlist(liftOver(gr.chr(NODTWC_v4C_gr), hs1ToHg38_chain))
    K0HCJ3_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','K0HCJ3','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    K0HCJ3_v4C_gr <- unlist(liftOver(gr.chr(K0HCJ3_v4C_gr), hs1ToHg38_chain))
    H1IDEVD_v4C <- fmt_gt(gTrack(H1IDEVD_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'H1IDEVD', col = 'blue', y1 = 100))
    H1E2AVA_v4C <- fmt_gt(gTrack(H1E2AVA_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'H1E2AVA', col = 'blue', y1 = 100))
    NODTWC_v4C <- fmt_gt(gTrack(NODTWC_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'NODTWC', col = 'blue', y1 = 100))
    K0HCJ3_v4C <- fmt_gt(gTrack(K0HCJ3_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'K0HCJ3', col = 'blue', y1 = 100))


window_oiz = window_oi -6e5
print('plotting')
pdf(file =paste0(OUT_DIR, '/FigS7E_normHiC_BRCA_NSCLC_ESCC_', xG, '.pdf'),paper = 'a4', height = 36)
        plot(c(gencode_hg38_gt, ATACfisher_LUSC_gt, NSCLC_xG_LvsR_sig_LR_gt, NSCLC_xGtad_per_amp_gt,
            dmrs_d_NSCLCgt,
            H1IDEVD_v4C, H1E2AVA_v4C, NODTWC_v4C, K0HCJ3_v4C,
            squamousEsoph51f_H3K27ac_bed_gt, squamousEsoph53f_H3K27ac_bed_gt, squamousEsoph54m_H3K27ac_bed_gt, squamousEsoph_median_gt,
            squamousEpEsoph_hg38_norm_HiC_v4C_gt,
            xT_map_sqEsoph$gtrack(name=paste0("contacts"),height=210,clim=c(0,20)),
            lung16f_atac_bed_gt, lung65f_atac_bed_gt, lung60m_atac_bed_gt, lung_atacFE_median_gt,
            lung16f_hg38_norm_HiC_v4C_gt,
            xT_map_Lu$gtrack(name=paste0("contacts"),height=210,clim=c(0,20)),
            fishres_sigs_BRCA_HR_gt,
            BRCAvsNSCLC_xG_chr11_sig_track_cons80_gt, BRCA_xGtad_per_amp_gt,
            dmrs_d_BRCAgt,
            C12_v4C, C132_v4C, C17_v4C, Cb5_v4C,
            BRCA51f_atac_bed_gt, BRCA53f_atac_bed_gt, BRCA54f_atac_bed_gt, BRCA_atacFE_median_gt,
            BRCA_hg38_norm_HiC_v4C_gt,
            xT_map$gtrack(name=paste0("contacts"),height=210,clim=c(0,15))),
            window_oiz, cex.label = 0.1)
        plot(c(gencode_hg38_gt, ATACfisher_LUSC_gt, NSCLC_xG_LvsR_sig_LR_gt, NSCLC_xGtad_per_amp_gt,
            dmrs_d_NSCLCgt,
            H1IDEVD_v4C, H1E2AVA_v4C, NODTWC_v4C, K0HCJ3_v4C,
            squamousEsoph51f_H3K27ac_bed_gt, squamousEsoph53f_H3K27ac_bed_gt, squamousEsoph54m_H3K27ac_bed_gt, squamousEsoph_median_gt,
            squamousEpEsoph_hg38_norm_HiC_v4C_gt,
            xT_map_sqEsoph$gtrack(name=paste0("contacts"),height=210,clim=c(0,20))),
            window_oiz, cex.label = 0.1)
        plot(c(gencode_hg38_gt, ATACfisher_LUSC_gt, NSCLC_xG_LvsR_sig_LR_gt, NSCLC_xGtad_per_amp_gt,
            dmrs_d_NSCLCgt,
            H1IDEVD_v4C, H1E2AVA_v4C, NODTWC_v4C, K0HCJ3_v4C,
            lung16f_atac_bed_gt, lung65f_atac_bed_gt, lung60m_atac_bed_gt, lung_atacFE_median_gt,
            lung16f_hg38_norm_HiC_v4C_gt,
            xT_map_Lu$gtrack(name=paste0("contacts"),height=210,clim=c(0,20))),
            window_oiz, cex.label = 0.1)
        plot(c(gencode_hg38_gt,             fishres_sigs_BRCA_HR_gt,
            BRCAvsNSCLC_xG_chr11_sig_track_cons80_gt, BRCA_xGtad_per_amp_gt,
            dmrs_d_BRCAgt,
            C12_v4C, C132_v4C, C17_v4C, Cb5_v4C,
            BRCA51f_atac_bed_gt, BRCA53f_atac_bed_gt, BRCA54f_atac_bed_gt, BRCA_atacFE_median_gt,
            BRCA_hg38_norm_HiC_v4C_gt,
            xT_map$gtrack(name=paste0("contacts"),height=210,clim=c(0,15))),
            window_oiz, cex.label = 0.1)
dev.off()


