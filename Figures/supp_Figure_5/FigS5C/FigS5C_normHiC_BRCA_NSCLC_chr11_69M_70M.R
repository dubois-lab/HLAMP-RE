# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(GxG)
library(gGnome)
library(rtracklayer)
hs1ToHg38_chain <- import.chain(file.path(IN_DIR, 'hs1ToHg38.over.chain')) # get the reverse

gencode_hg38_gt <- readRDS('../../supp_Figure_3/FigS3A/source_data/gencode_hg38.rds')
CGC_gr_hg38 <- readRDS('../../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
xG = 'CCND1'
Oncgr <- CGC_gr_hg38 %Q%  (Gene_Symbol == xG)
window_oi = GRanges("chr11:69000000-70000000")
binsize =5e3

# white y.grid for all data tracks
fmt_gt <- function(gt) {
    gt$y.grid.col= "white"
    gt
}


BRCA_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, '/20260907_Fig.2C_BRCA_CCND1_AmpFreq_track.rds'))
BRCA_xGtad_per_amp_gt_gr <- unlist(liftOver(BRCA_xGtad_per_amp_gt@data[[1]], hs1ToHg38_chain))
BRCA_xGtad_per_amp_gt <- gTrack(data = BRCA_xGtad_per_amp_gt_gr, y.field = "AmpFreq", name = "%amp", col = "red")
BRCA_xGtad_per_amp_gt <- fmt_gt(BRCA_xGtad_per_amp_gt)
BRCA_xGtad_per_amp_gt$bars <- F

BRCA_CCND1_LvsR_sig_LR_gt <- readRDS(paste0(IN_DIR, '/20260907c_Fig.2C_BRCA_CCND1_LvsR_sig_LR_cons80.rds'))
BRCA_CCND1_LvsR_sig_LR_gt_gr <- unlist(liftOver(BRCA_CCND1_LvsR_sig_LR_gt@data[[1]], hs1ToHg38_chain))
BRCA_CCND1_LvsR_sig_LR_gt <- gTrack(data = BRCA_CCND1_LvsR_sig_LR_gt_gr, name = "nlgp", col = "red", height = 0.5)

BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- readRDS(paste0(IN_DIR, '/20260907c_Fig.2C_BRCAvsNSCLC_CCND1_chr11_sig_track_cons80.rds'))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr <- unlist(liftOver(BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt@data[[1]], hs1ToHg38_chain))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- gTrack(data = BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr, name = "nlgp", col = "blue", height = 0.5)

NSCLC_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, '/20260907_Fig.2C_NSCLC_CCND1_AmpFreq_track.rds'))
NSCLC_xG_amps_grc <- NSCLC_xGtad_per_amp_gt@data[[1]]
NSCLC_xGtad_amps_hg38gr <- unlist(liftOver(NSCLC_xG_amps_grc, hs1ToHg38_chain))
NSCLC_xGtad_per_amp_gt <- gTrack(data = NSCLC_xGtad_amps_hg38gr, y.field = "AmpFreq", name = "%amp", col = "blue", y0 = 0, y1 = 1)
NSCLC_xGtad_per_amp_gt <- fmt_gt(NSCLC_xGtad_per_amp_gt)
NSCLC_xGtad_per_amp_gt$bars <- F

# #functional val####
GERCARres <- as.data.table(read.csv2(file.path(IN_DIR, 'GERCARres.csv')))
colnames(GERCARres)
GERCARres_gr <- dt2gr(GERCARres[chromosome == "chr11"])
GERCARres_prolifd21_gt <- gTrack(GERCARres_gr, y.field = "logFC_prolif_21", bars = T)
GERCARres_gr_sig <- dt2gr(GERCARres[FDR_prolif_21 < 0.1])
GERCARres_prolifd21_q0.1_gt <- gTrack(GERCARres_gr_sig)

# sigdiff ####
ATACfisher <- fread(file.path(IN_DIR, "20260614_CCND1_H3K27ac&ATAC_linDiffs_stats.csv"))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = "[:]"))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[3]), by = peakID]
fishres_A2gr <- dt2gr(ATACfisher)
fishres_A_sigs_BRCA_HR <- fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "BRCA_HR+")
fishres_sigs_BRCA_HR_gt <- gTrack(data = reduce(fishres_A_sigs_BRCA_HR))
ATACfisher_LUSC_gt <- gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "LUSC")))
print('Fisher test readin done')
xD <- 0.4
dmrs_d4 <- fread(paste0(IN_DIR, "/20260512_CCND1_TM_dmrs_d.csv"))
colnames(dmrs_d4)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d4_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs_d4)), hs1ToHg38_chain))
dmrs_d4 <- as.data.table(dmrs_d4_gr)
# each DMR is shaded by the highest |diff.Methy| cutoff (0.4/0.5/0.6) it passes;
# categorical colormap as for the tissue-colored Fisher track (a continuous per-row grey() color column does not render in gTrack)
dmrs_d4[, methyl_level := ifelse(abs(diff.Methy) > 0.6, "d6", ifelse(abs(diff.Methy) > 0.5, "d5", "d4"))]
colormap_dmrs_d_grey <- list(methyl_level = c(d4 = grey(0.7), d5 = grey(0.4), d6 = grey(0.1)))
dmrs_d_NSCLC_grey_gt <- gTrack(dt2gr(dmrs_d4[diff.Methy < -xD]), gr.colorfield = "color", colormaps = colormap_dmrs_d_grey)
dmrs_d_BRCA_grey_gt <- gTrack(dt2gr(dmrs_d4[diff.Methy > xD]), gr.colorfield = "color", colormaps = colormap_dmrs_d_grey)

CCND1_TSS <- fread(file.path(IN_DIR, 'CCND1_TSS_hg38.csv'))
colnames(CCND1_TSS)[4] = 'seqnames'
colnames(CCND1_TSS)[8] = 'start'
CCND1_TSS[, end := start +1]
CCND1_TSSgr <- dt2gr(CCND1_TSS)

print('reading in breast HiC data')
xT_map = straw(file.path(IN_DIR, 'ENCFF306HOC.hic'),  gr = window_oi, res = 5e3) #53f
gr_map_dt <- gr2dt(xT_map$gr)
# gr_map_dt
gr_map_dt[, seqnames := gsub(pattern = ':69000000:70000000', replacement = '', x = seqnames), by = seqnames]
# gr_map_dt
xT_map$gr <- dt2gr(gr_map_dt)
BRCA_hg38_norm_HiC_4C <- xT_map[xT_map$gr %^% gr.chr(CCND1_TSSgr), drop =T]
BRCA_hg38_norm_HiC_v4C_gt <- fmt_gt(gTrack(data = BRCA_hg38_norm_HiC_4C, y.field = 'value', bars = T, name = '4C_BRCA', y0 = 0, y1 = 30))

print('reading in lung HiC data done')
xT_map_Lu = straw(file.path(IN_DIR, 'ENCFF094KKX_lung16f_hg38.hic'),  gr = window_oi, res = 5e3) #53f
gr_map_LU_dt <- gr2dt(xT_map_Lu$gr)
# gr_map_dt
gr_map_LU_dt[, seqnames := gsub(pattern = ':69000000:70000000', replacement = '', x = seqnames), by = seqnames]
# gr_map_dt
xT_map_Lu$gr <- dt2gr(gr_map_LU_dt)
lung16f_hg38_norm_HiC_4C <- xT_map_Lu[xT_map_Lu$gr %^% gr.chr(CCND1_TSSgr), drop =T]
lung16f_hg38_norm_HiC_v4C_gt <- fmt_gt(gTrack(data = lung16f_hg38_norm_HiC_4C, y.field = 'value', bars = T, name = '4C_lung16f', y0 = 0, y1 = 40))
Lung_norm_16f_TAD_1 <- fread(file.path(IN_DIR, 'Lung_norm_16f_TAD_1_ENCFF581SYY.bedpe'))
colnames(Lung_norm_16f_TAD_1)[1:3] = c('chr', 'start', 'end')
Lung_norm_16f_TAD_1_gt <- gTrack(data = dt2gr(Lung_norm_16f_TAD_1) %&% gr.chr(CGC_gr_hg38), name = "Lung_norm_16f_TAD_1", col = "grey", height = 0.5)

print('reading in esophageal squamous HiC data')
xT_map_sqEsoph = straw(file.path(IN_DIR, 'ENCFF224HKR_squamousEpEsoph_hg38.hic'),  gr = window_oi, res = 5e3)
gr_map_sqEsoph_dt <- gr2dt(xT_map_sqEsoph$gr)
# gr_map_dt
gr_map_sqEsoph_dt[, seqnames := gsub(pattern = ':69000000:70000000', replacement = '', x = seqnames), by = seqnames]
# gr_map_dt
xT_map_sqEsoph$gr <- dt2gr(gr_map_sqEsoph_dt)
squamousEpEsoph_hg38_norm_HiC_4C <- xT_map_sqEsoph[xT_map_sqEsoph$gr %^% gr.chr(CCND1_TSSgr), drop =T]
squamousEpEsoph_hg38_norm_HiC_v4C_gt <- fmt_gt(gTrack(data = squamousEpEsoph_hg38_norm_HiC_4C, y.field = 'value', bars = T, name = '4C_squamousEpEsoph', y0 = 0, y1 = 40))
Esophageal_norm_53f_TAD2<- fread(file.path(IN_DIR, 'Esophageal_norm_53f_TAD_2_ENCFF735CNP.bedpe'))
colnames(Esophageal_norm_53f_TAD2)[1:3] = c('chr', 'start', 'end')
Esophageal_norm_53f_TAD2_gt <- gTrack(data = dt2gr(Esophageal_norm_53f_TAD2) %&% gr.chr(CGC_gr_hg38), name = "Esophageal_norm_53f_TAD2", col = "grey", height = 0.5)

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


print('reading in lung ATAC wig data')
lung_atacFE_median_gt <- fmt_gt(gTrack(file.path(IN_DIR, 'lung_median.bw'), y.field = "score", name = "lung_atacFE_median", bars = T, y0 = 0, y1 = 10))
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
    C3_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','C3','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    C3_v4C_gr <- unlist(liftOver(gr.chr(C3_v4C_gr), hs1ToHg38_chain))
    C169_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','C169','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    C169_v4C_gr <- unlist(liftOver(gr.chr(C169_v4C_gr), hs1ToHg38_chain))
    C17_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','C17','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    C17_v4C_gr <- unlist(liftOver(gr.chr(C17_v4C_gr), hs1ToHg38_chain))
    C3_v4C <- fmt_gt(gTrack(C3_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C3', col = 'red', y1 = 80))
    C169_v4C <- fmt_gt(gTrack(C169_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C169', col = 'red', y1 = 80))
    C17_v4C <- fmt_gt(gTrack(C17_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'C17', col = 'red', y1 = 80))
    H3D44S_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','H3D44S','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    H3D44S_v4C_gr <- unlist(liftOver(gr.chr(H3D44S_v4C_gr), hs1ToHg38_chain))
    KDAMIF_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','KDAMIF','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    KDAMIF_v4C_gr <- unlist(liftOver(gr.chr(KDAMIF_v4C_gr), hs1ToHg38_chain))
    Y57U8R_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','Y57U8R','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    Y57U8R_v4C_gr <- unlist(liftOver(gr.chr(Y57U8R_v4C_gr), hs1ToHg38_chain))
    F57J9F_v4C_gr <- readRDS(paste0(IN_DIR, '/individuals/','F57J9F','_covCovariate_v4C_', xG, '_',datafield,'.rds'))
    F57J9F_v4C_gr <- unlist(liftOver(gr.chr(F57J9F_v4C_gr), hs1ToHg38_chain))
    H3D44S_v4C <- fmt_gt(gTrack(H3D44S_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'H3D44S', col = 'blue', y1 = 80))
    KDAMIF_v4C <- fmt_gt(gTrack(KDAMIF_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'KDAMIF', col = 'blue', y1 = 80))
    Y57U8R_v4C <- fmt_gt(gTrack(Y57U8R_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'Y57U8R', col = 'blue', y1 = 80))
    F57J9F_v4C <- fmt_gt(gTrack(F57J9F_v4C_gr, y.field = 'value', bars = T, y0 = 0, name = 'F57J9F', col = 'blue', y1 = 80))

window_oiz = window_oi -2e5
window_oiz2 = window_oi - 4e5
print('plotting')
pdf(file =paste0(OUT_DIR, '/FigS5C.pdf'),paper = 'a4', height = 36)
        plot(c(gencode_hg38_gt, ATACfisher_LUSC_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt, NSCLC_xGtad_per_amp_gt,
            dmrs_d_NSCLC_grey_gt,
            H3D44S_v4C, KDAMIF_v4C, Y57U8R_v4C, F57J9F_v4C,
            squamousEsoph51f_H3K27ac_bed_gt, squamousEsoph53f_H3K27ac_bed_gt, squamousEsoph54m_H3K27ac_bed_gt, squamousEsoph_median_gt,
            squamousEpEsoph_hg38_norm_HiC_v4C_gt,
            Esophageal_norm_53f_TAD2_gt,
            xT_map_sqEsoph$gtrack(name=paste0("contacts"),height=210,clim=c(0,20)),
            lung16f_atac_bed_gt, lung65f_atac_bed_gt, lung60m_atac_bed_gt, lung_atacFE_median_gt,
            lung16f_hg38_norm_HiC_v4C_gt,
            Lung_norm_16f_TAD_1_gt,
            xT_map_Lu$gtrack(name=paste0("contacts"),height=210,clim=c(0,20)),
            fishres_sigs_BRCA_HR_gt,BRCA_CCND1_LvsR_sig_LR_gt, BRCA_xGtad_per_amp_gt,
            dmrs_d_BRCA_grey_gt, 
            C3_v4C, C17_v4C, C169_v4C,
            BRCA51f_atac_bed_gt, BRCA53f_atac_bed_gt, BRCA54f_atac_bed_gt, BRCA_atacFE_median_gt,
            BRCA_hg38_norm_HiC_v4C_gt,
            xT_map$gtrack(name=paste0("contacts"),height=210,clim=c(0,20))),
            window_oi, cex.label = 0.1)
        plot(c(gencode_hg38_gt, 
            fishres_sigs_BRCA_HR_gt,BRCA_CCND1_LvsR_sig_LR_gt, BRCA_xGtad_per_amp_gt,
            dmrs_d_BRCA_grey_gt, 
            C3_v4C, C17_v4C, C169_v4C,
            BRCA51f_atac_bed_gt, BRCA53f_atac_bed_gt, BRCA54f_atac_bed_gt, BRCA_atacFE_median_gt,
            BRCA_hg38_norm_HiC_v4C_gt,
            xT_map$gtrack(name=paste0("contacts"),height=210,clim=c(0,20))),
            window_oi, cex.label = 0.1)
        plot(c(gencode_hg38_gt, 
            ATACfisher_LUSC_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt, NSCLC_xGtad_per_amp_gt,
            dmrs_d_NSCLC_grey_gt,
            H3D44S_v4C, KDAMIF_v4C, Y57U8R_v4C, F57J9F_v4C,
            squamousEsoph51f_H3K27ac_bed_gt, squamousEsoph53f_H3K27ac_bed_gt, squamousEsoph54m_H3K27ac_bed_gt, squamousEsoph_median_gt,
            squamousEpEsoph_hg38_norm_HiC_v4C_gt,
            Esophageal_norm_53f_TAD2_gt,
            xT_map_sqEsoph$gtrack(name=paste0("contacts"),height=210,clim=c(0,20))),
            window_oi, cex.label = 0.1)
        plot(c(gencode_hg38_gt, 
            ATACfisher_LUSC_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt, NSCLC_xGtad_per_amp_gt,
            dmrs_d_NSCLC_grey_gt,
            H3D44S_v4C, KDAMIF_v4C, Y57U8R_v4C, F57J9F_v4C,
            lung16f_atac_bed_gt, lung65f_atac_bed_gt, lung60m_atac_bed_gt, lung_atacFE_median_gt,
            lung16f_hg38_norm_HiC_v4C_gt,
            Lung_norm_16f_TAD_1_gt,
            xT_map_Lu$gtrack(name=paste0("contacts"),height=210,clim=c(0,20))),
            window_oi, cex.label = 0.1)
        plot(c(gencode_hg38_gt, ATACfisher_LUSC_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt, NSCLC_xGtad_per_amp_gt,
            dmrs_d_NSCLC_grey_gt,
            H3D44S_v4C, KDAMIF_v4C, Y57U8R_v4C, F57J9F_v4C,
            squamousEsoph51f_H3K27ac_bed_gt, squamousEsoph53f_H3K27ac_bed_gt, squamousEsoph54m_H3K27ac_bed_gt, squamousEsoph_median_gt,
            squamousEpEsoph_hg38_norm_HiC_v4C_gt,
            xT_map_sqEsoph$gtrack(name=paste0("contacts"),height=240,clim=c(0,20)),
            lung16f_atac_bed_gt, lung65f_atac_bed_gt, lung60m_atac_bed_gt, lung_atacFE_median_gt,
            lung16f_hg38_norm_HiC_v4C_gt,
            Lung_norm_16f_TAD_1_gt,
            xT_map_Lu$gtrack(name=paste0("contacts"),height=240,clim=c(0,20)),
            fishres_sigs_BRCA_HR_gt,BRCA_CCND1_LvsR_sig_LR_gt, BRCA_xGtad_per_amp_gt,
            GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt, dmrs_d_BRCA_grey_gt, 
            C3_v4C, C17_v4C, C169_v4C,
            BRCA51f_atac_bed_gt, BRCA53f_atac_bed_gt, BRCA54f_atac_bed_gt, BRCA_atacFE_median_gt,
            BRCA_hg38_norm_HiC_v4C_gt,
            xT_map$gtrack(name=paste0("contacts"),height=240,clim=c(0,20))),
            window_oiz, cex.label = 0.1)
        plot(c(gencode_hg38_gt, ATACfisher_LUSC_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt, NSCLC_xGtad_per_amp_gt,
            dmrs_d_NSCLC_grey_gt,
            H3D44S_v4C, KDAMIF_v4C, Y57U8R_v4C, F57J9F_v4C,
            squamousEsoph51f_H3K27ac_bed_gt, squamousEsoph53f_H3K27ac_bed_gt, squamousEsoph54m_H3K27ac_bed_gt, squamousEsoph_median_gt,
            squamousEpEsoph_hg38_norm_HiC_v4C_gt,
            xT_map_sqEsoph$gtrack(name=paste0("contacts"),height=210,clim=c(0,20)),
            lung16f_atac_bed_gt, lung65f_atac_bed_gt, lung60m_atac_bed_gt, lung_atacFE_median_gt,
            lung16f_hg38_norm_HiC_v4C_gt,
            xT_map_Lu$gtrack(name=paste0("contacts"),height=210,clim=c(0,20)),        
            fishres_sigs_BRCA_HR_gt,BRCA_CCND1_LvsR_sig_LR_gt, BRCA_xGtad_per_amp_gt,
            GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt, dmrs_d_BRCA_grey_gt, 
            C3_v4C, C17_v4C, C169_v4C,
            BRCA51f_atac_bed_gt, BRCA53f_atac_bed_gt, BRCA54f_atac_bed_gt, BRCA_atacFE_median_gt,
            BRCA_hg38_norm_HiC_v4C_gt,
            xT_map$gtrack(name=paste0("contacts"),height=210,clim=c(0,20))),
            window_oiz2, cex.label = 0.1)
dev.off()
