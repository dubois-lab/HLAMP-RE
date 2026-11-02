# nolint start: line_length_linter.
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(gGnome)
library(gTrack)
library(rtracklayer)
# library(BSgenome.Hsapiens.UCSC.hg38)

# white y.grid for all data tracks
fmt_gt <- function(gt) {
    # gt$draw.paths = FALSE
    # gt$sep.draw = FALSE
    # gt$sep.bg.col = "white"
    # gt$sep.lwd = 0
    # gt$stack.gap = 0
    # gt$stack.gap = 0
    gt$y.grid.col= "white"
    gt
}

# hg38
CGC <- readRDS("../../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds")
gencode_hg38 <- readRDS("../../supp_Figure_3/FigS3A/source_data/gencode_hg38.rds")
conversion_chain <- import.chain(file.path(IN_DIR, "hg38ToHs1.over.chain")) # hg38ToHs1 chain

# Load hs1ToHg38 conversion chain
hs1ToHg38_chain_file <- file.path(IN_DIR, "hs1ToHg38.over.chain")
if (!file.exists(hs1ToHg38_chain_file)) {
    download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz", destfile = paste0(hs1ToHg38_chain_file, ".gz"))
    system(paste("gunzip -f", paste0(hs1ToHg38_chain_file, ".gz")))
}
hs1ToHg38_chain <- import.chain(hs1ToHg38_chain_file)

xG <- "CCND1"
Oncgr <- CGC %Q% (Gene_Symbol == xG)
window_oi <- Oncgr + 20e5
out_dir <- paste0(OUT_DIR, '/')



# 	CCND1 (H3K27ac) BRCAHR+ ####
# ZR751####
ZR751_H3K27ac_FE <- fmt_gt(gTrack(data = file.path(IN_DIR, "ZR751_FE.bw"), bars = T, name = "ZR751"))
ZR751_H3K27ac_peaks <- fread(file.path(IN_DIR, "ZR751_peaks.narrowPeak"))
colnames(ZR751_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
ZR751_H3K27ac_peaks_gt <- fmt_gt(gTrack(dt2gr(ZR751_H3K27ac_peaks), height = 0.5))
# CAMA1####
CAMA1_H3K27ac_FE_bdg <- readRDS(file.path(IN_DIR, "CAMA1_H3K27ac_FE_bdg.rds"))
CAMA1_H3K27ac_FE <- fmt_gt(gTrack(data = CAMA1_H3K27ac_FE_bdg, y.field = "score", bars = T, name = "CAMA1"))
rm(CAMA1_H3K27ac_FE_bdg)
CAMA1_H3K27ac_peaks <- fread(file.path(IN_DIR, "CAMA1_peaks.narrowPeak"))
colnames(CAMA1_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
CAMA1_H3K27ac_peaks_gt <- fmt_gt(gTrack(dt2gr(CAMA1_H3K27ac_peaks), height = 0.5))
# MDAMB134VI####
MDAMB134VI_H3K27ac_FE_bdg <- readRDS(file.path(IN_DIR, "MDAMB134VI_H3K27ac_FE_bdg.rds"))
MDAMB134VI_H3K27ac_FE <- fmt_gt(gTrack(data = MDAMB134VI_H3K27ac_FE_bdg, y.field = "score", bars = T, name = "MDAMB134VI"))
MDAMB134VI_H3K27ac_peaks <- fread(file.path(IN_DIR, "MDAMB134VI_peaks.narrowPeak"))
colnames(MDAMB134VI_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
MDAMB134VI_H3K27ac_peaks_gt <- fmt_gt(gTrack(dt2gr(MDAMB134VI_H3K27ac_peaks), height = 0.5))
# MDAMB361####
MDAMB361_H3K27ac_FE_bdg <- readRDS(file.path(IN_DIR, "MDAMB361_H3K27ac_FE_bdg.rds"))
MDAMB361_H3K27ac_FE <- fmt_gt(gTrack(data = MDAMB361_H3K27ac_FE_bdg, y.field = "score", bars = T, name = "MDAMB361"))
rm(MDAMB361_H3K27ac_FE_bdg)
MDAMB361_H3K27ac_peaks <- fread(file.path(IN_DIR, "MDAMB361_peaks.narrowPeak"))
colnames(MDAMB361_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
MDAMB361_H3K27ac_peaks_gt <- fmt_gt(gTrack(dt2gr(MDAMB361_H3K27ac_peaks), height = 0.5))



# 	CCND1 (H3K27ac) LUSC####
# HCC95####
HCC95_H3K27ac_bw_hg38 <- readRDS(file.path(IN_DIR, "GSM2784723_H3K27ac-HCC95-WT-A.bw_hg38.rds"))
HCC95_H3K27ac_gt <- fmt_gt(gTrack(data = gr.nochr(HCC95_H3K27ac_bw_hg38), y.field = "score", bars = T))
HCC95_H3K27ac_bed_hg19 <- fread(file.path(IN_DIR, "GSM4079747_HCC95_H3K27ac_peaks_hg19.bed"))
colnames(HCC95_H3K27ac_bed_hg19)[1:3] <- c("seqnames", "start", "end")
HCC95_H3K27ac_bed_hg38 <- unlist(liftOver(dt2gr(HCC95_H3K27ac_bed_hg19), chain = conversion_chain))
HCC95_H3K27ac_bed_gt <- fmt_gt(gTrack(data = gr.nochr(HCC95_H3K27ac_bed_hg38), height = 0.5))

# ATAC shared meta ####
T_normalizedCounts_Intgenes <- fread(file.path(IN_DIR, "20220104RNAseqwithATACbreastlungNormalizedCountsIntgenes.txt"))
atacMeta <- fread(file.path(IN_DIR, "The chromatin accessibility landscape of primary human cancers_SuppTbl1c.txt"))
atacMeta_B <- atacMeta[cohort %in% c("LUAD", "BRCA", "LUSC", "UCEC")]
atacMeta_B[, TCGA_ID := paste(unlist(strsplit(x = Tissue_Barcode, split = "[-]"))[1:4], collapse = "-"), by = Tissue_Barcode]
atacMeta_Bme <- as.data.table(merge(x = atacMeta_B, y = T_normalizedCounts_Intgenes, by.x = "TCGA_ID", by.y = "TCGA_Sample_ID", all.x = T)) # 4x diff -> one missing sample, check dropout
lof <- list.files(file.path(IN_DIR, "atac_bigwigs"))
lof_dt <- as.data.table(lof[lof != "Icon\r"])
lof_dt[, LibID := gsub(pattern = "_P[[:digit:]]*.insertions.bw", replacement = "", x = V1)]
lof_dt[, LibID := gsub(pattern = "_PMRG.insertions.bw", replacement = "", x = LibID)]
lof_dts <- lof_dt[LibID %in% atacMeta_Bme$Library_Name]



# LUSC ####
LUSC_H1 <- "TCGA-MF-A522-01A"
atacMeta_Bme_LUSC_H1 <- atacMeta_Bme[TCGA_ID == LUSC_H1]
LibLUSC_H1 <- atacMeta_Bme_LUSC_H1[which.max(Final_DeDup_reads), Library_Name]
goiLUSC_H1gt <- fmt_gt(gTrack(data = paste0(
    IN_DIR, "/atac_bigwigs/",
    lof_dts[LibID == LibLUSC_H1, V1]
), name = "LUSC_CCND1_h1", y1 = 150, bars = T))
goiLUSC_H1bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", LUSC_H1, "/", LUSC_H1, "_peaks.narrowPeak"))
colnames(goiLUSC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H1bed$seqnames)
goiLUSC_H1bed_gr <- dt2gr(goiLUSC_H1bed)
goiLUSC_H1bed_gt <- fmt_gt(gTrack(data = goiLUSC_H1bed_gr, height = 0.5))
LUSC_H2 <- "TCGA-98-A53J-01A"
atacMeta_Bme_LUSC_H2 <- atacMeta_Bme[TCGA_ID == LUSC_H2]
LibLUSC_H2 <- atacMeta_Bme_LUSC_H2[which.max(Final_DeDup_reads), Library_Name]
goiLUSC_H2gt <- fmt_gt(gTrack(data = paste0(
    IN_DIR, "/atac_bigwigs/",
    lof_dts[LibID == LibLUSC_H2, V1]
), name = "LUSC_CCND1_h2", y1 = 150, bars = T))
goiLUSC_H2bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", LUSC_H2, "/", LUSC_H2, "_peaks.narrowPeak"))
colnames(goiLUSC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H2bed$seqnames)
goiLUSC_H2bed_gr <- dt2gr(goiLUSC_H2bed)
goiLUSC_H2bed_gt <- fmt_gt(gTrack(data = goiLUSC_H2bed_gr, height = 0.5))
LUSC_H3 <- "TCGA-77-A5G6-01A"
atacMeta_Bme_LUSC_H3 <- atacMeta_Bme[TCGA_ID == LUSC_H3]
LibLUSC_H3 <- atacMeta_Bme_LUSC_H3[which.max(Final_DeDup_reads), Library_Name]
goiLUSC_H3gt <- fmt_gt(gTrack(data = paste0(
    IN_DIR, "/atac_bigwigs/",
    lof_dts[LibID == LibLUSC_H3, V1]
), name = "LUSC_CCND1_h3", y1 = 150, bars = T))
goiLUSC_H3bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", LUSC_H3, "/", LUSC_H3, "_peaks.narrowPeak"))
colnames(goiLUSC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H3bed$seqnames)
goiLUSC_H3bed_gr <- dt2gr(goiLUSC_H3bed)
goiLUSC_H3bed_gt <- fmt_gt(gTrack(data = goiLUSC_H3bed_gr, height = 0.5))
LUSC_H4 <- "TCGA-XC-AA0X-01A"
atacMeta_Bme_LUSC_H4 <- atacMeta_Bme[TCGA_ID == LUSC_H4]
LibLUSC_H4 <- atacMeta_Bme_LUSC_H4[which.max(Final_DeDup_reads), Library_Name]
goiLUSC_H4gt <- fmt_gt(gTrack(data = paste0(
    IN_DIR, "/atac_bigwigs/",
    lof_dts[LibID == LibLUSC_H4, V1]
), name = "LUSC_CCND1_H4", y1 = 150, bars = T))
goiLUSC_H4bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", LUSC_H4, "/", LUSC_H4, "_peaks.narrowPeak"))
colnames(goiLUSC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H4bed$seqnames)
goiLUSC_H4bed_gr <- dt2gr(goiLUSC_H4bed)
goiLUSC_H4bed_gt <- fmt_gt(gTrack(data = goiLUSC_H4bed_gr, height = 0.5))

# BRCA ####
lofBRCABRCA <- list.files(file.path(IN_DIR, "bigwigs"))
lofBRCA_dt <- as.data.table(lofBRCABRCA[lofBRCABRCA != "Icon\r"])
lofBRCA_dt[, LibID := gsub(pattern = "_P[[:digit:]]*.insertions.bw", replacement = "", x = V1)]
lofBRCA_dt[, LibID := gsub(pattern = "_PMRG.insertions.bw", replacement = "", x = LibID)]
lofBRCA_dts <- lofBRCA_dt[LibID %in% atacMeta_Bme$Library_Name]



# BRCA_HR+ ####
HRBC_H1 <- "TCGA-C8-A12M-01A" # CN(1.5)
atacMeta_Bme_HRBC_H1 <- atacMeta_Bme[TCGA_ID == HRBC_H1]
LibHRBC_H1 <- atacMeta_Bme_HRBC_H1[which.max(Final_DeDup_reads), Library_Name]
goiHRBC_H1gt <- fmt_gt(gTrack(data = paste0(
    IN_DIR, "/bigwigs/",
    lofBRCA_dts[LibID == LibHRBC_H1, V1]
), name = "HRBC_HR+_CCND1_h1", y1 = 300, bars = T))
goiHRBC_H1bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", HRBC_H1, "/", HRBC_H1, "_peaks.narrowPeak"))
colnames(goiHRBC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H1bed$seqnames)
goiHRBC_H1bed_gr <- dt2gr(goiHRBC_H1bed[V9 > -log10(0.05)])
goiHRBC_H1bed_gt <- fmt_gt(gTrack(data = reduce(goiHRBC_H1bed_gr), height = 0.5))
HRBC_H2 <- "TCGA-C8-A12Y-01A" # CN(1.76)
atacMeta_Bme_HRBC_H2 <- atacMeta_Bme[TCGA_ID == HRBC_H2]
LibHRBC_H2 <- atacMeta_Bme_HRBC_H2[which.max(Final_DeDup_reads), Library_Name]
goiHRBC_H2gt <- fmt_gt(gTrack(data = paste0(
    IN_DIR, "/bigwigs/",
    lofBRCA_dts[LibID == LibHRBC_H2, V1]
), name = "HRBC_HR+_CCND1_HRBC_H2", y1 = 300, bars = T))
goiHRBC_H2bed <- fread(paste0(IN_DIR, "/ATACbeds/out/", HRBC_H2, "/", HRBC_H2, "_peaks.narrowPeak"))
colnames(goiHRBC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H2bed$seqnames)
goiHRBC_H2bed_gr <- dt2gr(goiHRBC_H2bed[V9 > -log10(0.05)])
goiHRBC_H2bed_gt <- fmt_gt(gTrack(data = reduce(goiHRBC_H2bed_gr), height = 0.5))
HRBC_H3 <- "TCGA-C8-A12O-01A" # high exp & high CN (2.2)
atacMeta_Bme_HRBC_H3 <- atacMeta_Bme[TCGA_ID == HRBC_H3]
LibHRBC_H3 <- atacMeta_Bme_HRBC_H3[which.max(Final_DeDup_reads), Library_Name]
goiHRBC_H3gt <- fmt_gt(gTrack(data = paste0(
    IN_DIR, "/bigwigs/",
    lofBRCA_dts[LibID == LibHRBC_H3, V1]
), name = "HRBC_HR+_CCND1_HRBC_H3", y1 = 300, bars = T))
goiHRBC_H3bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", HRBC_H3, "/", HRBC_H3, "_peaks.narrowPeak"))
colnames(goiHRBC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H3bed$seqnames)
goiHRBC_H3bed_gr <- dt2gr(goiHRBC_H3bed[V9 > -log10(0.05)])
goiHRBC_H3bed_gt <- fmt_gt(gTrack(data = reduce(goiHRBC_H3bed_gr), height = 0.5))
HRBC_H4 <- "TCGA-A2-A0YF-01A" #
atacMeta_Bme_HRBC_H4 <- atacMeta_Bme[TCGA_ID == HRBC_H4]
LibHRBC_H4 <- atacMeta_Bme_HRBC_H4[which.max(Final_DeDup_reads), Library_Name]
goiHRBC_H4gt <- fmt_gt(gTrack(data = paste0(
    IN_DIR, "/bigwigs/",
    lofBRCA_dts[LibID == LibHRBC_H4, V1]
), name = "HRBC_HR+_CCND1_HRBC_H4", y1 = 300, bars = T))
goiHRBC_H4bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", HRBC_H4, "/", HRBC_H4, "_peaks.narrowPeak"))
colnames(goiHRBC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H4bed$seqnames)
goiHRBC_H4bed_gr <- dt2gr(goiHRBC_H4bed[V9 > -log10(0.05)])
goiHRBC_H4bed_gt <- fmt_gt(gTrack(data = reduce(goiHRBC_H4bed_gr), height = 0.5))

# functional val####
GERCARres <- as.data.table(read.csv2(file.path(IN_DIR, "GERCARres.csv")))
colnames(GERCARres)
GERCARres_gr <- dt2gr(GERCARres[chromosome == "chr11"])
GERCARres_prolifd21_gt <- fmt_gt(gTrack(GERCARres_gr, y.field = "logFC_prolif_21", bars = T))
GERCARres_gr_sig <- dt2gr(GERCARres[FDR_prolif_21 < 0.1])
GERCARres_prolifd21_q0.1_gt <- fmt_gt(gTrack(GERCARres_gr_sig))

# sigdiff ####
ATACfisher <- fread(file.path(IN_DIR, "20250614_CCND1_H3K27ac&ATAC_linDiffs_stats.csv"))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = "[:]"))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[3]), by = peakID]
fishres_A2gr <- dt2gr(ATACfisher)
fishres_A_sigs_BRCA_HR <- fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "BRCA_HR+")
fishres_sigs_BRCA_HR_gt <- fmt_gt(gTrack(data = reduce(fishres_A_sigs_BRCA_HR), ygap =0.5))
ATACfisher_LUSC_gt <- fmt_gt(gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "LUSC")), ygap =0.5))


xD <- 0.4
dmrs_d4 <- fread(paste0(IN_DIR, '/20260512_CCND1_TM_dmrs_d.csv'))
colnames(dmrs_d4)[6:7] = c("meanMethy_NSCLC", "meanMethy_BRCA")
# DMRs are called on hs1 -> lift to hg38 for this figure
dmrs_d4 <- as.data.table(unlist(liftOver(gr.chr(dt2gr(dmrs_d4)), hs1ToHg38_chain)))
# one greyscale-coded track per tissue instead of the three dmrs_d4/d5/d6 cutoff tracks:
# each DMR is shaded by the highest |diff.Methy| cutoff (0.4/0.5/0.6) it passes;
# categorical colormap as for the tissue-colored Fisher track (a continuous per-row grey() color column does not render in gTrack)
dmrs_d4[, methyl_level := ifelse(abs(diff.Methy) > 0.6, "d6", ifelse(abs(diff.Methy) > 0.5, "d5", "d4"))]
colormap_dmrs_d_grey <- list(methyl_level = c(d4 = grey(0.7), d5 = grey(0.4), d6 = grey(0.1)))
dmrs_d_NSCLC_grey_gt <- fmt_gt(gTrack(dt2gr(dmrs_d4[diff.Methy < -xD]), gr.colorfield = "color", colormaps = colormap_dmrs_d_grey))
dmrs_d_BRCA_grey_gt <- fmt_gt(gTrack(dt2gr(dmrs_d4[diff.Methy > xD]), gr.colorfield = "color", colormaps = colormap_dmrs_d_grey))


BRCA_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, "/20260907_Fig.2C_BRCA_CCND1_AmpFreq_track.rds"))
BRCA_xGtad_per_amp_gt_gr <- unlist(liftOver(BRCA_xGtad_per_amp_gt@data[[1]], hs1ToHg38_chain))
BRCA_xGtad_per_amp_gt <- fmt_gt(gTrack(data = BRCA_xGtad_per_amp_gt_gr, y.field = "AmpFreq", name = "%amp", col = "red"))
BRCA_xGtad_per_amp_gt$bars <- F

BRCA_CCND1_LvsR_sig_LR_gt <- readRDS(paste0(IN_DIR, "/20260907c_Fig.2C_BRCA_CCND1_LvsR_sig_LR_cons80.rds")) # nolint: line_length_linter.
BRCA_CCND1_LvsR_sig_LR_gt_gr <- unlist(liftOver(BRCA_CCND1_LvsR_sig_LR_gt@data[[1]], hs1ToHg38_chain))
BRCA_CCND1_LvsR_sig_LR_gt <- fmt_gt(gTrack(data = reduce(BRCA_CCND1_LvsR_sig_LR_gt_gr), name = "nlgp", col = "red", height = 0.5))

BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- readRDS(paste0(IN_DIR, "/20260907c_Fig.2C_BRCAvsNSCLC_CCND1_chr11_sig_track_cons80.rds")) # nolint: line_length_linter.
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr <- unlist(liftOver(BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt@data[[1]], hs1ToHg38_chain))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- fmt_gt(gTrack(data = BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr, name = "nlgp", col = "blue", height = 0.5))

NSCLC_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, "/20260907_Fig.2C_NSCLC_CCND1_AmpFreq_track.rds"))
NSCLC_xG_amps_grc <- NSCLC_xGtad_per_amp_gt@data[[1]]
NSCLC_xGtad_amps_hg38gr <- unlist(liftOver(NSCLC_xG_amps_grc, hs1ToHg38_chain))
NSCLC_xGtad_per_amp_gt <- fmt_gt(gTrack(data = NSCLC_xGtad_amps_hg38gr, y.field = "AmpFreq", name = "%amp", col = "blue", y0 = 0, y1 = 1))
NSCLC_xGtad_per_amp_gt$bars <- F

# plot ####
window_oi2 <- Oncgr + 4e5
window_oi3 <- Oncgr + 6e5

pdf(
    file = paste0(out_dir, "20260929b_Fig_3B_CCND1_BRCAvsNSCLC.pdf"),
    paper = "legal", height = 80
)
plot(c(
    gencode_hg38,
    dmrs_d_BRCA_grey_gt, fishres_sigs_BRCA_HR_gt,
    goiHRBC_H4bed_gt, goiHRBC_H4gt, goiHRBC_H3bed_gt, goiHRBC_H3gt, goiHRBC_H2bed_gt, goiHRBC_H2gt, goiHRBC_H1bed_gt, goiHRBC_H1gt,
    MDAMB361_H3K27ac_peaks_gt, MDAMB361_H3K27ac_FE, MDAMB134VI_H3K27ac_peaks_gt, MDAMB134VI_H3K27ac_FE, CAMA1_H3K27ac_peaks_gt, CAMA1_H3K27ac_FE, ZR751_H3K27ac_peaks_gt, ZR751_H3K27ac_FE,
    BRCA_xGtad_per_amp_gt, BRCA_CCND1_LvsR_sig_LR_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
    dmrs_d_NSCLC_grey_gt, ATACfisher_LUSC_gt,
    HCC95_H3K27ac_bed_gt, HCC95_H3K27ac_gt,
    goiLUSC_H4bed_gt, goiLUSC_H4gt, goiLUSC_H3bed_gt, goiLUSC_H3gt, goiLUSC_H2bed_gt, goiLUSC_H2gt, goiLUSC_H1bed_gt, goiLUSC_H1gt, NSCLC_xGtad_per_amp_gt
), window_oi2, cex.label = 0.1)
plot(c(
    GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt,
    dmrs_d_BRCA_grey_gt, fishres_sigs_BRCA_HR_gt,
    goiHRBC_H4bed_gt, goiHRBC_H4gt, goiHRBC_H3bed_gt, goiHRBC_H3gt, goiHRBC_H2bed_gt, goiHRBC_H2gt, goiHRBC_H1bed_gt, goiHRBC_H1gt,
    MDAMB361_H3K27ac_peaks_gt, MDAMB361_H3K27ac_FE, MDAMB134VI_H3K27ac_peaks_gt, MDAMB134VI_H3K27ac_FE, CAMA1_H3K27ac_peaks_gt, CAMA1_H3K27ac_FE, ZR751_H3K27ac_peaks_gt, ZR751_H3K27ac_FE,
    BRCA_xGtad_per_amp_gt, BRCA_CCND1_LvsR_sig_LR_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
    dmrs_d_NSCLC_grey_gt, ATACfisher_LUSC_gt,
    HCC95_H3K27ac_bed_gt, HCC95_H3K27ac_gt,
    goiLUSC_H4bed_gt, goiLUSC_H4gt, goiLUSC_H3bed_gt, goiLUSC_H3gt, goiLUSC_H2bed_gt, goiLUSC_H2gt, goiLUSC_H1bed_gt, goiLUSC_H1gt, NSCLC_xGtad_per_amp_gt
), window_oi2, cex.label = 0.1)
dev.off()
