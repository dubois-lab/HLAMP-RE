# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(gGnome)
library(gTrack)
library(rtracklayer)
# hg38
# Per-sample ATAC bigwig coverage and per-cell-line H3K27ac/ATAC FE (bigwig + bdgcmp .rds) tracks are replaced by
# per-tumor-type MEDIAN tracks:
#   - tumor ATAC medians:        .wig files in the mediantracks folder
#   - cell-line H3K27ac medians: {tumortype}_EGFR_CL.rds in the same folder
#     (HRpos has no cell-line median; HCC827 ATAC excluded from the LUAD median)
# Individual-sample ATAC peak beds and cell-line H3K27ac/ATAC peaks are kept.
CGC <- readRDS("FigS3A/source_data/CGC_hg38.rds")
gencode_hg38 <- readRDS("FigS3A/source_data/gencode_hg38.rds")
conversion_chain <- readRDS(file.path(IN_DIR, "hg19ToHg38.over.chain.rds"))

# Load hs1ToHg38 conversion chain (for T2T diffmethyl amp/sig tracks)
hs1ToHg38_chain_file <- file.path(IN_DIR, "hs1ToHg38.over.chain")
if (!file.exists(hs1ToHg38_chain_file)) {
    download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz", destfile = paste0(hs1ToHg38_chain_file, ".gz"))
    system(paste("gunzip -f", paste0(hs1ToHg38_chain_file, ".gz")))
}
hs1ToHg38_chain <- import.chain(hs1ToHg38_chain_file)

xG <- "EGFR"
Oncgr <- CGC %Q% (Gene_Symbol == xG)
window_oi <- Oncgr + 20e5

med_dir <- file.path(IN_DIR, "mediantracks")

#	EGFR (H3K27ac) BRCAHR- cell-line peaks #####
# MDAMB468
MDAMB468_H3K27ac_peaks <- fread(file.path(IN_DIR, "MDAMB468_REP5_peaks.narrowPeak"))
colnames(MDAMB468_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
MDAMB468_H3K27ac_peaks_gt <- gTrack(dt2gr(MDAMB468_H3K27ac_peaks), height = 0.5)
# HDQP1
HDQP1_H3K27ac_peaks <- fread(file.path(IN_DIR, "HDQP1_peaks.narrowPeak"))
colnames(HDQP1_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
HDQP1_H3K27ac_peaks_gt <- gTrack(dt2gr(HDQP1_H3K27ac_peaks), height = 0.5)
# SUM229PE
SUM229PE_H3K27ac_peaks <- fread(file.path(IN_DIR, "Sum229PE_H3K27ac_peaks.narrowPeak"))
colnames(SUM229PE_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
SUM229PE_H3K27ac_peaks_gt <- gTrack(dt2gr(SUM229PE_H3K27ac_peaks), height = 0.5)
# HCC1954
HCC1954_H3K27ac_peaks <- fread("FigS3A/source_data/HCC1954_peaks.narrowPeak")
colnames(HCC1954_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
HCC1954_H3K27ac_peaks_gt <- gTrack(dt2gr(HCC1954_H3K27ac_peaks), height = 0.5)
# BRCAHR- cell-line H3K27ac median
HRnegCL_gr <- readRDS(file.path(med_dir, "BRCA_TN_EGFR_CL.rds"))
HRnegCL_gt <- gTrack(HRnegCL_gr, y.field = "score", name = "HRneg_CL", bars = T, y1 = 15)

#	EGFR (H3K27ac) BRCAHR+ cell-line peaks ####
# T47D
T47D_H3K27ac_FE_bdg_gr = readRDS(file.path(IN_DIR, 'T47D_H3K27ac_FE_bdg.rds'))
T47D_H3K27ac_FE <- gTrack(data = T47D_H3K27ac_FE_bdg_gr, y.field = 'score', bars = T, name = 'T47D')
rm(T47D_H3K27ac_FE_bdg_gr)
T47D_H3K27ac_peaks <- fread(file.path(IN_DIR, "T47D_H3K27ac_peaks.narrowPeak"))
colnames(T47D_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
T47D_H3K27ac_peaks_gt <- gTrack(dt2gr(T47D_H3K27ac_peaks), height = 0.5)
# (HRpos has no cell-line H3K27ac median track)

#	EGFR (H3K27ac/ATAC) LUAD cell-line peaks ####
# HCC827 (ATAC)
HCC827_ATAC_peaks <- fread(file.path(IN_DIR, "GSM5777958_X1c_3562AZ_HCC827_CTRL-NF_peaks.narrowPeak"))
colnames(HCC827_ATAC_peaks)[1:3] <- c("chr", "start", "end")
HCC827_ATAC_peaks_gt <- gTrack(dt2gr(HCC827_ATAC_peaks[V9 > 5 & V7 > 10]), height = 0.5)
# PC9
PC9_H3K27ac_peaks <- fread(file.path(IN_DIR, "PC9_H3K27ac_ENCFF214UQW.bed"))
colnames(PC9_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
PC9_H3K27ac_peaks_gt <- gTrack(dt2gr(PC9_H3K27ac_peaks[V9 > 2]), height = 0.5)
# PC14
PC14_H3K27ac_peaks <- fread(file.path(IN_DIR, "PC14_H3K27ac_peaks.narrowPeak"))
colnames(PC14_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
PC14_H3K27ac_peaks_gt <- gTrack(dt2gr(PC14_H3K27ac_peaks), height = 0.5)
# NCIH1975
NCIH1975_H3K27ac_peaks <- fread(file.path(IN_DIR, "NCIH1975_peaks.narrowPeak"))
colnames(NCIH1975_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
NCIH1975_H3K27ac_peaks_gt <- gTrack(dt2gr(NCIH1975_H3K27ac_peaks), height = 0.5)
# LUAD cell-line H3K27ac median (HCC827 ATAC excluded)
LUADCL_gr <- readRDS(file.path(med_dir, "LUAD_EGFR_CL.rds"))
LUADCL_gt <- gTrack(LUADCL_gr, y.field = "score", name = "LUAD_CL", bars = T, y1 = 15)

#	EGFR (H3K27ac) LUSC cell-line peaks ####
# HCC2814_PE
HCC2814_PE_H3K27ac_peaks <- fread(file.path(IN_DIR, "HCC2814_PE_H3K27ac_peaks.narrowPeak"))
colnames(HCC2814_PE_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
HCC2814_PE_H3K27ac_peaks_gt <- gTrack(dt2gr(HCC2814_PE_H3K27ac_peaks), height = 0.5)
# HCC95
HCC95_H3K27ac_bed_hg38 <- fread(file.path(IN_DIR, "GSM4079747_HCC95_H3K27ac_peaks_hg38.bed"))
colnames(HCC95_H3K27ac_bed_hg38)[1:3] <- c("seqnames", "start", "end")
HCC95_H3K27ac_bed_gt <- gTrack(data = dt2gr(HCC95_H3K27ac_bed_hg38), height = 0.5)
# LUSC cell-line H3K27ac median
LUSCCL_gr <- readRDS(file.path(med_dir, "LUSC_EGFR_CL.rds"))
LUSCCL_gt <- gTrack(LUSCCL_gr, y.field = "score", name = "LUSC_CL", bars = T, y1 = 15)

# Median ATAC tumor tracks (per tumor type) ####
LUAD_gt       <- gTrack(import(file.path(med_dir, "LUAD_median.wig")),    y.field = "score", name = "LUAD",       bars = T, y1 = 150)
LUSC_gt       <- gTrack(import(file.path(med_dir, "LUSC_median.wig")),    y.field = "score", name = "LUSC",       bars = T, y1 = 150)
UCEC_gt       <- gTrack(import(file.path(med_dir, "UCEC_median.wig")),    y.field = "score", name = "UCEC",       bars = T, y1 = 150)
HNSC_gt       <- gTrack(import(file.path(med_dir, "HNSC_median.wig")),    y.field = "score", name = "HNSC",       bars = T, y1 = 150)
BRCA_HRneg_gt <- gTrack(import(file.path(med_dir, "BRCA_TN_median.wig")), y.field = "score", name = "BRCA_HRneg", bars = T, y1 = 150)
BRCA_HRpos_gt <- gTrack(import(file.path(med_dir, "HRpos_median.wig")),   y.field = "score", name = "BRCA_HRpos", bars = T, y1 = 150)

# BRCA ATAC peak-to-gene links (EGFR) ####
BRCApeakToGeneLinks <- fread(file.path(IN_DIR, "BRCAtcgaATACpeakToGeneLinks.txt"))
BRCApeakToGeneLinks_goi <- BRCApeakToGeneLinks[Linked_Gene == "EGFR"]
BRCApeakToGeneLinks_goigr <- dt2gr(BRCApeakToGeneLinks_goi)
BRCApeakToGeneLinks_goigt <- gTrack(data = BRCApeakToGeneLinks_goigr, name = "BRCApeakToGeneLinks")

# Individual-sample ATAC peak beds ####
atac_beds <- paste0(IN_DIR, "/ATACbeds/out2/")
ucec_hnsc_beds <- paste0(IN_DIR, "/atac_peaks/")

# LUAD ####
H1 <- "TCGA-MP-A4TF-01A"
goiH1bed <- fread(paste0(atac_beds, H1, "/", H1, "_peaks.narrowPeak"))
colnames(goiH1bed)[1:3] <- c("seqnames", "start", "end")
goiH1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH1bed$seqnames)
goiH1bed_gr <- dt2gr(goiH1bed)
goiH1bed_gt <- gTrack(data = goiH1bed_gr, height = 0.5)
H2 <- "TCGA-49-AARO-01A"
goiH2bed <- fread(paste0(IN_DIR, "/20250707new_macs3_atac/", H2, "_atac__peaks.narrowPeak"))
colnames(goiH2bed)[1:3] <- c("seqnames", "start", "end")
goiH2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH2bed$seqnames)
goiH2bed_gr <- dt2gr(goiH2bed)
goiH2bed_gt <- gTrack(data = goiH2bed_gr, height = 0.5)
H3 <- "TCGA-44-3918-01A"
goiH3bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", H3, "/", H3, "_peaks.narrowPeak"))
colnames(goiH3bed)[1:3] <- c("seqnames", "start", "end")
goiH3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH3bed$seqnames)
goiH3bed_gr <- dt2gr(goiH3bed)
goiH3bed_gt <- gTrack(data = goiH3bed_gr, height = 0.5)
H4 <- "TCGA-86-A4P7-01A"
goiH4bed <- fread(paste0(IN_DIR, "/ATACbeds/out/", H4, "_atac__peaks.narrowPeak"))
colnames(goiH4bed)[1:3] <- c("seqnames", "start", "end")
goiH4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH4bed$seqnames)
goiH4bed_gr <- dt2gr(goiH4bed)
goiH4bed_gt <- gTrack(data = goiH4bed_gr, height = 0.5)

# LUSC ####
LUSC_H1 <- "TCGA-O2-A52N-01A"
goiLUSC_H1bed <- fread(paste0(atac_beds, LUSC_H1, "/", LUSC_H1, "_peaks.narrowPeak"))
colnames(goiLUSC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H1bed$seqnames)
goiLUSC_H1bed_gr <- dt2gr(goiLUSC_H1bed)
goiLUSC_H1bed_gt <- gTrack(data = goiLUSC_H1bed_gr, height = 0.5)
LUSC_H2 <- "TCGA-85-A50M-01A"
goiLUSC_H2bed <- fread(paste0(atac_beds, LUSC_H2, "/", LUSC_H2, "_peaks.narrowPeak"))
colnames(goiLUSC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H2bed$seqnames)
goiLUSC_H2bed_gr <- dt2gr(goiLUSC_H2bed)
goiLUSC_H2bed_gt <- gTrack(data = goiLUSC_H2bed_gr, height = 0.5)
LUSC_H3 <- "TCGA-77-A5G6-01A"
goiLUSC_H3bed <- fread(paste0(atac_beds, LUSC_H3, "/", LUSC_H3, "_peaks.narrowPeak"))
colnames(goiLUSC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H3bed$seqnames)
goiLUSC_H3bed_gr <- dt2gr(goiLUSC_H3bed)
goiLUSC_H3bed_gt <- gTrack(data = goiLUSC_H3bed_gr, height = 0.5)
LUSC_H4 <- "TCGA-98-A53J-01A"
goiLUSC_H4bed <- fread(paste0(atac_beds, LUSC_H4, "/", LUSC_H4, "_peaks.narrowPeak"))
colnames(goiLUSC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H4bed$seqnames)
goiLUSC_H4bed_gr <- dt2gr(goiLUSC_H4bed)
goiLUSC_H4bed_gt <- gTrack(data = goiLUSC_H4bed_gr, height = 0.5)

# BRCA_TN ####
BRCA_TN_H1 <- "TCGA-AR-A0U4-01A"
goiBRCA_TN_H1bed <- fread(paste0(atac_beds, BRCA_TN_H1, "/", BRCA_TN_H1, "_peaks.narrowPeak"))
colnames(goiBRCA_TN_H1bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H1bed$seqnames)
goiBRCA_TN_H1bed_gr <- dt2gr(goiBRCA_TN_H1bed[V9 > -log10(0.05)])
goiBRCA_TN_H1bed_gt <- gTrack(data = reduce(goiBRCA_TN_H1bed_gr), height = 0.5)
BRCA_TN_H2 <- "TCGA-BH-A0E0-01A"
goiBRCA_TN_H2bed <- fread(paste0(atac_beds, BRCA_TN_H2, "/", BRCA_TN_H2, "_peaks.narrowPeak"))
colnames(goiBRCA_TN_H2bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H2bed$seqnames)
goiBRCA_TN_H2bed_gr <- dt2gr(goiBRCA_TN_H2bed[V9 > -log10(0.05)])
goiBRCA_TN_H2bed_gt <- gTrack(data = reduce(goiBRCA_TN_H2bed_gr), height = 0.5)
BRCA_TN_H3 <- "TCGA-AR-A0U0-01A"
goiBRCA_TN_H3bed <- fread(paste0(atac_beds, BRCA_TN_H3, "/", BRCA_TN_H3, "_summits.bed"))
colnames(goiBRCA_TN_H3bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H3bed$seqnames)
goiBRCA_TN_H3bed_gr <- dt2gr(goiBRCA_TN_H3bed)
goiBRCA_TN_H3bed_gt <- gTrack(data = reduce(goiBRCA_TN_H3bed_gr), height = 0.5)
BRCA_TN_H4 <- "TCGA-C8-A12V-01A"
goiBRCA_TN_H4bed <- fread(paste0(IN_DIR, "/20250707new_macs3_atac/", BRCA_TN_H4, "_atac__peaks.narrowPeak"))
colnames(goiBRCA_TN_H4bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H4bed$seqnames)
goiBRCA_TN_H4bed_gr <- dt2gr(goiBRCA_TN_H4bed[V9 > -log10(0.05)])
goiBRCA_TN_H4bed_gt <- gTrack(data = reduce(goiBRCA_TN_H4bed_gr), height = 0.5)

# BRCA_HR+ ####
HRBC_H1 <- "TCGA-BH-A0DL-01A"
goiHRBC_H1bed <- fread(paste0(IN_DIR, "/20250707new_macs3_atac/", HRBC_H1, "_atac__peaks.narrowPeak"))
colnames(goiHRBC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H1bed$seqnames)
goiHRBC_H1bed_gr <- dt2gr(goiHRBC_H1bed[V9 > -log10(0.05)])
goiHRBC_H1bed_gt <- gTrack(data = reduce(goiHRBC_H1bed_gr), height = 0.5)
HRBC_H2 <- "TCGA-A2-A0YJ-01A"
goiHRBC_H2bed <- fread(paste0(atac_beds, HRBC_H2, "/", HRBC_H2, "_peaks.narrowPeak"))
colnames(goiHRBC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H2bed$seqnames)
goiHRBC_H2bed_gr <- dt2gr(goiHRBC_H2bed[V9 > -log10(0.05)])
goiHRBC_H2bed_gt <- gTrack(data = reduce(goiHRBC_H2bed_gr), height = 0.5)
HRBC_H3 <- "TCGA-A7-A13E-01A"
goiHRBC_H3bed <- fread(paste0(atac_beds, HRBC_H3, "/", HRBC_H3, "_peaks.narrowPeak"))
colnames(goiHRBC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H3bed$seqnames)
goiHRBC_H3bed_gr <- dt2gr(goiHRBC_H3bed[V9 > -log10(0.05)])
goiHRBC_H3bed_gt <- gTrack(data = reduce(goiHRBC_H3bed_gr), height = 0.5)
HRBC_H4 <- "TCGA-S3-AA0Z-01A"
goiHRBC_H4bed <- fread(paste0(atac_beds, HRBC_H4, "/", HRBC_H4, "_peaks.narrowPeak"))
colnames(goiHRBC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H4bed$seqnames)
goiHRBC_H4bed_gr <- dt2gr(goiHRBC_H4bed[V9 > -log10(0.05)])
goiHRBC_H4bed_gt <- gTrack(data = reduce(goiHRBC_H4bed_gr), height = 0.5)

# UCEC ####
UCEC_H1 <- "TCGA-B5-A0JN-01A"
goiUCEC_H1bed <- fread(paste0(ucec_hnsc_beds, UCEC_H1, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H1bed$seqnames)
goiUCEC_H1bed_gr <- dt2gr(goiUCEC_H1bed)
goiUCEC_H1bed_gt <- gTrack(data = goiUCEC_H1bed_gr, height = 0.5)
UCEC_H2 <- "TCGA-B5-A11Y-01A"
goiUCEC_H2bed <- fread(paste0(ucec_hnsc_beds, UCEC_H2, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H2bed$seqnames)
goiUCEC_H2bed_gr <- dt2gr(goiUCEC_H2bed)
goiUCEC_H2bed_gt <- gTrack(data = goiUCEC_H2bed_gr, height = 0.5)
UCEC_H3 <- "TCGA-BG-A0RY-01A"
goiUCEC_H3bed <- fread(paste0(ucec_hnsc_beds, UCEC_H3, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H3bed$seqnames)
goiUCEC_H3bed_gr <- dt2gr(goiUCEC_H3bed)
goiUCEC_H3bed_gt <- gTrack(data = goiUCEC_H3bed_gr, height = 0.5)
UCEC_H4 <- "TCGA-BG-A0MS-01A"
goiUCEC_H4bed <- fread(paste0(ucec_hnsc_beds, UCEC_H4, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H4bed$seqnames)
goiUCEC_H4bed_gr <- dt2gr(goiUCEC_H4bed)
goiUCEC_H4bed_gt <- gTrack(data = goiUCEC_H4bed_gr, height = 0.5)

# HNSC ####
HNSC_H1 <- "TCGA-BA-A8YP-01A"
goiHNSC_H1bed <- fread(paste0(ucec_hnsc_beds, HNSC_H1, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H1bed$seqnames)
goiHNSC_H1bed_gr <- dt2gr(goiHNSC_H1bed)
goiHNSC_H1bed_gt <- gTrack(data = goiHNSC_H1bed_gr, height = 0.5)
HNSC_H2 <- "TCGA-BA-A4IH-01A"
goiHNSC_H2bed <- fread(paste0(ucec_hnsc_beds, HNSC_H2, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H2bed$seqnames)
goiHNSC_H2bed_gr <- dt2gr(goiHNSC_H2bed)
goiHNSC_H2bed_gt <- gTrack(data = goiHNSC_H2bed_gr, height = 0.5)
HNSC_H3 <- "TCGA-D6-A6EP-01A"
goiHNSC_H3bed <- fread(paste0(ucec_hnsc_beds, HNSC_H3, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H3bed$seqnames)
goiHNSC_H3bed_gr <- dt2gr(goiHNSC_H3bed)
goiHNSC_H3bed_gt <- gTrack(data = goiHNSC_H3bed_gr, height = 0.5)
HNSC_H4 <- "TCGA-CN-A49B-01A"
goiHNSC_H4bed <- fread(paste0(ucec_hnsc_beds, HNSC_H4, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H4bed$seqnames)
goiHNSC_H4bed_gr <- dt2gr(goiHNSC_H4bed)
goiHNSC_H4bed_gt <- gTrack(data = goiHNSC_H4bed_gr, height = 0.5)

#functional val####
SID_MCF7_E2 <- as.data.table(read.csv2(file.path(IN_DIR, "cd-23-1157_supplementary_table_s3_suppst3_MCF7_E2.csv")))
SID_MCF7_E2_gr <- dt2gr(SID_MCF7_E2)
SID_MCF7_E2_prolifd21_gt <- gTrack(SID_MCF7_E2_gr, y.field = "log2FC_day21_vs_day07", bars = T)
SID_MCF7_E2_sig <- SID_MCF7_E2[FDR_day21_vs_day07 < 0.1]
SID_MCF7_E2_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_MCF7_E2_sig))

SID_MCF7_LTED <- as.data.table(read.csv2(file.path(IN_DIR, "cd-23-1157_supplementary_table_s9_MCF7_LTED.csv")))
SID_MCF7_LTED[, seqnamesL := unlist(strsplit(identifier, "[:]"))[1], by = identifier]
SID_MCF7_LTED[, seqnames := unlist(strsplit(seqnamesL, "[_]"))[2], by = seqnamesL]
SID_MCF7_LTED[, pos12 := unlist(strsplit(identifier, "[:]"))[2], by = identifier]
SID_MCF7_LTED[, start := as.numeric(unlist(strsplit(pos12, "[-]"))[1]), by = pos12]
SID_MCF7_LTED[, end := as.numeric(unlist(strsplit(pos12, "[-]"))[2]), by = pos12]
SID_MCF7_LTED_gr <- dt2gr(SID_MCF7_LTED)
SID_MCF7_LTED_prolifd21_gt <- gTrack(SID_MCF7_LTED_gr, y.field = "log2FC_day21_vs_day07", bars = T)
SID_MCF7_LTED_sig <- SID_MCF7_LTED[FDR_day21_vs_day07 < 0.1]
SID_MCF7_LTED_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_MCF7_LTED_sig))

SID_T47D_E2 <- as.data.table(read.csv2(file.path(IN_DIR, "cd-23-1157_supplementary_table_s4_suppst4_T47D_E2.csv")))
SID_T47D_E2_gr <- dt2gr(SID_T47D_E2)
SID_T47D_E2_prolifd21_gt <- gTrack(SID_T47D_E2_gr, y.field = "log2FC_day21_vs_lib", bars = T, name = "SID_T47D")
SID_T47D_E2_sig <- SID_T47D_E2[FDR_day21_vs_lib < 0.1]
SID_T47D_E2_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_T47D_E2_sig))

GERCARres <- as.data.table(read.csv2(file.path(IN_DIR, "GERCARres.csv")))
GERCARres_gr <- dt2gr(GERCARres[chromosome == as.character(seqnames(gr.chr(Oncgr)))])
GERCARres_prolifd21_gt <- gTrack(gr.nochr(GERCARres_gr), y.field = "logFC_prolif_21", bars = T, name = "GERCAR_MCF7")
GERCARres_gr_sig <- dt2gr(GERCARres[FDR_prolif_21 < 0.1 & chromosome == as.character(seqnames(gr.chr(Oncgr)))])
GERCARres_prolifd21_q0.1_gt <- gTrack(gr.nochr(GERCARres_gr_sig))

# #sigdiff (cell-line H3K27ac/ATAC Fisher lineage differences) ####
fishres_A2 <- fread(file.path(IN_DIR, "20260623q_EGFR_H3K27ac_ATAC_linDiffs_stats.csv"))
fishres_A2[, seqnames := unlist(strsplit(x = peakID, split = "[:]"))[1], by = peakID]
fishres_A2[, start := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[2]), by = peakID]
fishres_A2[, end := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[3]), by = peakID]
fishres_A2_gr <- dt2gr(fishres_A2)
ATACfisher_LUAD_gt  <- gTrack(data = reduce(fishres_A2_gr %Q% (q_val < 0.15 & tstdttype == "LUAD")))
ATACfisher_HRBC_gt  <- gTrack(data = reduce(fishres_A2_gr %Q% (q_val < 0.15 & tstdttype == "BRCA_HR+")))
ATACfisher_HRNBC_gt <- gTrack(data = reduce(fishres_A2_gr %Q% (q_val < 0.15 & tstdttype == "BRCA_HR-")))
ATACfisher_LUSC_gt  <- gTrack(data = reduce(fishres_A2_gr %Q% (q_val < 0.15 & tstdttype == "LUSC")))
ATACfisher_UCEC_gt  <- gTrack(data = reduce(fishres_A2_gr %Q% (q_val < 0.15 & tstdttype == "UCEC")))
ATACfisher_HNSC_gt  <- gTrack(data = reduce(fishres_A2_gr %Q% (q_val < 0.15 & tstdttype == "HNSC")))

# diffmethyl ####
xD <- 0.4
dmrs_d4 <- fread(paste0(IN_DIR, "/20251003_EGFR_TM_dmrs_d_delt_0.4.csv"))
colnames(dmrs_d4)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d4_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs_d4)), hs1ToHg38_chain))
dmrs_d4 <- as.data.table(dmrs_d4_gr)
dmrs_d4_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -xD])))
dmrs_d4_BRCAgt <- gTrack(dt2gr(dmrs_d4[diff.Methy > xD]))

# amp/sig tracks (T2T/hs1 -> hg38) ####
# Built from the CNMatm/Fisher tables; each resulting GRanges is lifted
# hs1->hg38 before being wrapped in gTrack.
# BRCA/NSCLC EGFR amp-frequency tracks
BRCA_EGFR_CNMatm <- readRDS(file.path(IN_DIR, "20260910_BRCA_EGFR_CNMatm_combined_hs1.rds"))
NSCLC_EGFR_CNMatm <- readRDS(file.path(IN_DIR, "20260910_NSCLC_LUSC_EGFR_CNMatm_combined_hs1.rds"))
BRCA_EGFR_CNMatm_hg38gr <- unlist(liftOver(gr.chr(dt2gr(BRCA_EGFR_CNMatm[, 1:10])), hs1ToHg38_chain))
NSCLC_EGFR_CNMatm_hg38gr <- unlist(liftOver(gr.chr(dt2gr(NSCLC_EGFR_CNMatm[, 1:10])), hs1ToHg38_chain))
# region with >80% amp frequency in BRCA or NSCLC (hg38)
window_oi <- reduce(dt2gr(rbind(gr2dt(BRCA_EGFR_CNMatm_hg38gr)[AmpPercent > 0.8], gr2dt(NSCLC_EGFR_CNMatm_hg38gr)[AmpPercent > 0.8])))

BRCA_xGtad_per_amp_gt <- gTrack(BRCA_EGFR_CNMatm_hg38gr, y.field = "AmpPercent", col = "brown3", bars = F, y0 = 0, y1 = 1)
NSCLC_xGtad_per_amp_gt <- gTrack(NSCLC_EGFR_CNMatm_hg38gr, y.field = "AmpPercent", col = "royalblue3", bars = F, y0 = 0, y1 = 1)

# BRCA vs NSCLC EGFR differential amp-frequency significance track (chr7)
BRCAvsNSCLC_EGFR_fisher_gr <- readRDS(file.path(IN_DIR, "20260910_BRCAvsNSCLC_EGFR_chr7_FisherStats_CN10_cons80_gr.rds"))
BRCAvsNSCLC_EGFR_sig <- BRCAvsNSCLC_EGFR_fisher_gr[BRCAvsNSCLC_EGFR_fisher_gr$Fisher_pval < 0.05]
mcols(BRCAvsNSCLC_EGFR_sig)$tissue <- ifelse(mcols(BRCAvsNSCLC_EGFR_sig)$BRCA_AmpFreq > mcols(BRCAvsNSCLC_EGFR_sig)$NSCLC_AmpFreq, "BRCA", "NSCLC")
colormap_BRCAvsNSCLC_EGFR <- list(tissue = c(NSCLC = "royalblue3", BRCA = "brown3"))
BRCAvsNSCLC_EGFR_sig_hg38gr <- unlist(liftOver(gr.chr(BRCAvsNSCLC_EGFR_sig), hs1ToHg38_chain))
BRCAvsNSCLC_EGFR_chr7_sig_track_cons80_gt <- gTrack(BRCAvsNSCLC_EGFR_sig_hg38gr, gr.colorfield = "tissue", colormaps = colormap_BRCAvsNSCLC_EGFR, height = 0.5)

# plot ####
window_oi2 <- Oncgr + 3e5
window_oi3 <- Oncgr + 1e5

pdf(file = paste0(OUT_DIR, "/20260929_Fig_S3G_EGFR_H3K27acATAC_BRCAvsNSCLC_with_peaks_and_medianwigs.pdf")
    , paper = "legal", height = 80)
plot(c(gencode_hg38,
    ATACfisher_UCEC_gt,
    goiUCEC_H4bed_gt, goiUCEC_H3bed_gt, goiUCEC_H2bed_gt, goiUCEC_H1bed_gt, UCEC_gt,
    ATACfisher_HRBC_gt,
    goiHRBC_H4bed_gt, goiHRBC_H3bed_gt, goiHRBC_H2bed_gt, goiHRBC_H1bed_gt, BRCA_HRpos_gt,
    T47D_H3K27ac_peaks_gt,T47D_H3K27ac_FE,
    ATACfisher_HRNBC_gt, 
    MDAMB468_H3K27ac_peaks_gt, HDQP1_H3K27ac_peaks_gt, SUM229PE_H3K27ac_peaks_gt, HCC1954_H3K27ac_peaks_gt, HRnegCL_gt,
    goiBRCA_TN_H4bed_gt, goiBRCA_TN_H3bed_gt, goiBRCA_TN_H2bed_gt, goiBRCA_TN_H1bed_gt, BRCA_HRneg_gt,
    dmrs_d4_BRCAgt,
    BRCA_xGtad_per_amp_gt, BRCAvsNSCLC_EGFR_chr7_sig_track_cons80_gt,NSCLC_xGtad_per_amp_gt,
    dmrs_d4_NSCLCgt, ATACfisher_LUSC_gt, HCC2814_PE_H3K27ac_peaks_gt, HCC95_H3K27ac_bed_gt, LUSCCL_gt,
    goiLUSC_H4bed_gt, goiLUSC_H3bed_gt, goiLUSC_H2bed_gt, goiLUSC_H1bed_gt, LUSC_gt,
    ATACfisher_LUAD_gt,
    PC14_H3K27ac_peaks_gt, PC9_H3K27ac_peaks_gt, HCC827_ATAC_peaks_gt, NCIH1975_H3K27ac_peaks_gt, LUADCL_gt,
    goiH4bed_gt, goiH3bed_gt, goiH2bed_gt, goiH1bed_gt, LUAD_gt, 
    ATACfisher_HNSC_gt,
    goiHNSC_H4bed_gt, goiHNSC_H3bed_gt, goiHNSC_H2bed_gt, goiHNSC_H1bed_gt, HNSC_gt
), window_oi2, cex.label = 0.1)
plot(c(gencode_hg38,
    ATACfisher_UCEC_gt,
    goiUCEC_H4bed_gt, goiUCEC_H3bed_gt, goiUCEC_H2bed_gt, goiUCEC_H1bed_gt, UCEC_gt,
    SID_MCF7_LTED_sig_prolifd21_q0.1_gt, SID_MCF7_LTED_prolifd21_gt,
    dmrs_d4_BRCAgt, ATACfisher_HRNBC_gt, 
    MDAMB468_H3K27ac_peaks_gt, HDQP1_H3K27ac_peaks_gt, SUM229PE_H3K27ac_peaks_gt, HCC1954_H3K27ac_peaks_gt, HRnegCL_gt,
    goiBRCA_TN_H4bed_gt, goiBRCA_TN_H3bed_gt, goiBRCA_TN_H2bed_gt, goiBRCA_TN_H1bed_gt, BRCA_HRneg_gt,
    GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt, SID_MCF7_E2_sig_prolifd21_q0.1_gt, SID_MCF7_E2_prolifd21_gt, SID_T47D_E2_sig_prolifd21_q0.1_gt, SID_T47D_E2_prolifd21_gt, ATACfisher_HRBC_gt,
    goiHRBC_H4bed_gt, goiHRBC_H3bed_gt, goiHRBC_H2bed_gt, goiHRBC_H1bed_gt, BRCA_HRpos_gt,
    T47D_H3K27ac_peaks_gt, T47D_H3K27ac_FE,
    BRCApeakToGeneLinks_goigt, BRCA_xGtad_per_amp_gt, BRCAvsNSCLC_EGFR_chr7_sig_track_cons80_gt,
    dmrs_d4_NSCLCgt, ATACfisher_LUSC_gt, HCC2814_PE_H3K27ac_peaks_gt, HCC95_H3K27ac_bed_gt, LUSCCL_gt,
    goiLUSC_H4bed_gt, goiLUSC_H3bed_gt, goiLUSC_H2bed_gt, goiLUSC_H1bed_gt, LUSC_gt,
    ATACfisher_LUAD_gt,
    PC14_H3K27ac_peaks_gt, PC9_H3K27ac_peaks_gt, HCC827_ATAC_peaks_gt, NCIH1975_H3K27ac_peaks_gt, LUADCL_gt,
    goiH4bed_gt, goiH3bed_gt, goiH2bed_gt, goiH1bed_gt, LUAD_gt, NSCLC_xGtad_per_amp_gt,
    ATACfisher_HNSC_gt,
    goiHNSC_H4bed_gt, goiHNSC_H3bed_gt, goiHNSC_H2bed_gt, goiHNSC_H1bed_gt, HNSC_gt
), window_oi2, cex.label = 0.1)
plot(c(gencode_hg38,
    ATACfisher_UCEC_gt,
    goiUCEC_H4bed_gt, goiUCEC_H3bed_gt, goiUCEC_H2bed_gt, goiUCEC_H1bed_gt, UCEC_gt,
    SID_MCF7_LTED_sig_prolifd21_q0.1_gt, SID_MCF7_LTED_prolifd21_gt,
    dmrs_d4_BRCAgt, ATACfisher_HRNBC_gt, 
    MDAMB468_H3K27ac_peaks_gt, HDQP1_H3K27ac_peaks_gt, SUM229PE_H3K27ac_peaks_gt, HCC1954_H3K27ac_peaks_gt, HRnegCL_gt,
    goiBRCA_TN_H4bed_gt, goiBRCA_TN_H3bed_gt, goiBRCA_TN_H2bed_gt, goiBRCA_TN_H1bed_gt, BRCA_HRneg_gt,
    GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt, SID_MCF7_E2_sig_prolifd21_q0.1_gt, SID_MCF7_E2_prolifd21_gt, SID_T47D_E2_sig_prolifd21_q0.1_gt, SID_T47D_E2_prolifd21_gt, ATACfisher_HRBC_gt,
    goiHRBC_H4bed_gt, goiHRBC_H3bed_gt, goiHRBC_H2bed_gt, goiHRBC_H1bed_gt, BRCA_HRpos_gt,
    T47D_H3K27ac_peaks_gt, T47D_H3K27ac_FE,
    BRCApeakToGeneLinks_goigt, BRCA_xGtad_per_amp_gt, BRCAvsNSCLC_EGFR_chr7_sig_track_cons80_gt,
    dmrs_d4_NSCLCgt, ATACfisher_LUSC_gt, HCC2814_PE_H3K27ac_peaks_gt, HCC95_H3K27ac_bed_gt, LUSCCL_gt,
    goiLUSC_H4bed_gt, goiLUSC_H3bed_gt, goiLUSC_H2bed_gt, goiLUSC_H1bed_gt, LUSC_gt,
    ATACfisher_LUAD_gt,
    PC14_H3K27ac_peaks_gt, PC9_H3K27ac_peaks_gt, HCC827_ATAC_peaks_gt, NCIH1975_H3K27ac_peaks_gt, LUADCL_gt,
    goiH4bed_gt, goiH3bed_gt, goiH2bed_gt, goiH1bed_gt, LUAD_gt, NSCLC_xGtad_per_amp_gt,
    ATACfisher_HNSC_gt,
    goiHNSC_H4bed_gt, goiHNSC_H3bed_gt, goiHNSC_H2bed_gt, goiHNSC_H1bed_gt, HNSC_gt
), window_oi3, cex.label = 0.1)
dev.off()
