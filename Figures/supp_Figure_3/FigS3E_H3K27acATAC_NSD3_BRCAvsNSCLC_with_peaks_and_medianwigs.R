# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(gGnome)
library(gTrack)
library(rtracklayer)
# hg38
# Per-sample coverage is replaced by per-tumor-type MEDIAN tracks (no individual-sample coverage bigwig tracks):
#   - tumor ATAC medians: .wig files in the mediantracks folder
#   - cell-line H3K27ac/ATAC medians: {tumortype}_NSD3_CL.rds in the same folder
# Individual-sample ATAC peak beds and cell-line H3K27ac peaks are kept.
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

xG <- "NSD3"
Oncgr <- CGC %Q% (Gene_Symbol == xG)

med_dir <- file.path(IN_DIR, "mediantracks")

#	NSD3 (H3K27ac) BRCAHR- cell-line peaks #####
# MFM223
MFM223_H3K27ac_peaks <- fread(file.path(IN_DIR, "MFM223_peaks.narrowPeak"))
colnames(MFM223_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
MFM223_H3K27ac_peaks_gt <- gTrack(dt2gr(MFM223_H3K27ac_peaks), height = 0.5)
# CAL120
CAL120_H3K27ac_peaks <- fread(file.path(IN_DIR, "CAL120_peaks.narrowPeak"))
colnames(CAL120_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
CAL120_H3K27ac_peaks_gt <- gTrack(dt2gr(CAL120_H3K27ac_peaks), height = 0.5)
# MDAMB157
MDAMB157_H3K27ac_peaks <- fread(file.path(IN_DIR, "MDAMB157_peaks.narrowPeak"))
colnames(MDAMB157_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
MDAMB157_H3K27ac_peaks_gt <- gTrack(dt2gr(MDAMB157_H3K27ac_peaks), height = 0.5)
# BRCAHR- cell-line H3K27ac/ATAC median
HRnegCL_gr <- readRDS(file.path(med_dir, "BRCA_TN_NSD3_CL.rds"))
HRnegCL_gt <- gTrack(HRnegCL_gr, y.field = "score", name = "HRneg_CL", bars = T, y1 = 15)

#	NSD3 (H3K27ac) BRCAHR+ cell-line peaks ####
# ZR751
ZR751_H3K27ac_peaks <- fread(file.path(IN_DIR, "ZR751_peaks.narrowPeak"))
colnames(ZR751_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
ZR751_H3K27ac_peaks_gt <- gTrack(dt2gr(ZR751_H3K27ac_peaks), height = 0.5)
# CAMA1
CAMA1_H3K27ac_peaks <- fread(file.path(IN_DIR, "CAMA1_peaks.narrowPeak"))
colnames(CAMA1_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
CAMA1_H3K27ac_peaks_gt <- gTrack(dt2gr(CAMA1_H3K27ac_peaks), height = 0.5)
# MDAMB134VI
MDAMB134VI_H3K27ac_peaks <- fread(file.path(IN_DIR, "MDAMB134VI_peaks.narrowPeak"))
colnames(MDAMB134VI_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
MDAMB134VI_H3K27ac_peaks_gt <- gTrack(dt2gr(MDAMB134VI_H3K27ac_peaks), height = 0.5)
# BRCAHR+ cell-line H3K27ac/ATAC median
HRposCL_gr <- readRDS(file.path(med_dir, "HRpos_NSD3_CL.rds"))
HRposCL_gt <- gTrack(HRposCL_gr, y.field = "score", name = "HRpos_CL", bars = T, y1 = 10)

#	NSD3 (H3K27ac) LUAD cell-line peaks ####
# NCIH520
NCIH520_H3K27ac_peaks <- fread(file.path(IN_DIR, "NCIH520_peaks.narrowPeak"))
colnames(NCIH520_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
NCIH520_H3K27ac_peaks_gt <- gTrack(dt2gr(NCIH520_H3K27ac_peaks), height = 0.5)
# NCIH2087
NCIH2087_H3K27ac_peaks <- fread(file.path(IN_DIR, "NCIH2087_peaks.narrowPeak"))
colnames(NCIH2087_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
NCIH2087_H3K27ac_peaks_gt <- gTrack(dt2gr(NCIH2087_H3K27ac_peaks), height = 0.5)
# NCIH1703
NCIH1703_H3K27ac_peaks <- fread(file.path(IN_DIR, "NCIH1703_peaks.narrowPeak"))
colnames(NCIH1703_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
NCIH1703_H3K27ac_peaks_gt <- gTrack(dt2gr(NCIH1703_H3K27ac_peaks), height = 0.5)
# LUAD cell-line H3K27ac/ATAC median
LUADCL_gr <- readRDS(file.path(med_dir, "LUAD_NSD3_CL.rds"))
LUADCL_gt <- gTrack(LUADCL_gr, y.field = "score", name = "LUAD_CL", bars = T, y1 = 15)

# HCC95 (LUSC cell line) peaks ####
HCC95_H3K27ac_bed_hg38 <- fread(file.path(IN_DIR, "GSM4079747_HCC95_H3K27ac_peaks_hg38.bed"))
colnames(HCC95_H3K27ac_bed_hg38)[1:3] <- c("seqnames", "start", "end")
HCC95_H3K27ac_bed_gt <- gTrack(data = dt2gr(HCC95_H3K27ac_bed_hg38), height = 0.5)
# LUSC cell-line H3K27ac/ATAC median (HCC95 only -> median of one = itself)
LUSCCL_gr <- readRDS(file.path(med_dir, "LUSC_NSD3_CL.rds"))
LUSCCL_gt <- gTrack(LUSCCL_gr, y.field = "score", name = "LUSC_CL", bars = T, y1 = 15)

# Median ATAC tumor tracks (per tumor type) ####
LUAD_gt       <- gTrack(import(file.path(med_dir, "LUAD_median.wig")),    y.field = "score", name = "LUAD",       bars = T, y1 = 150)
LUSC_gt       <- gTrack(import(file.path(med_dir, "LUSC_median.wig")),    y.field = "score", name = "LUSC",       bars = T, y1 = 150)
UCEC_gt       <- gTrack(import(file.path(med_dir, "UCEC_median.wig")),    y.field = "score", name = "UCEC",       bars = T, y1 = 150)
HNSC_gt       <- gTrack(import(file.path(med_dir, "HNSC_median.wig")),    y.field = "score", name = "HNSC",       bars = T, y1 = 150)
BRCA_HRneg_gt <- gTrack(import(file.path(med_dir, "BRCA_TN_median.wig")), y.field = "score", name = "BRCA_HRneg", bars = T, y1 = 150)
BRCA_HRpos_gt <- gTrack(import(file.path(med_dir, "HRpos_median.wig")),   y.field = "score", name = "BRCA_HRpos", bars = T, y1 = 150)

# Individual-sample ATAC peak beds ####
atac_beds <- paste0(IN_DIR, "/ATACbeds/out2/")

# LUAD ####
H1 <- "TCGA-86-A4D0-01A"
goiH1bed <- fread(paste0(atac_beds, H1, "/", H1, "_peaks.narrowPeak"))
colnames(goiH1bed)[1:3] <- c("seqnames", "start", "end")
goiH1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH1bed$seqnames)
goiH1bed_gr <- dt2gr(goiH1bed)
goiH1bed_gt <- gTrack(data = goiH1bed_gr, height = 0.5)
H2 <- "TCGA-44-6146-01A"
goiH2bed <- fread(paste0(atac_beds, H2, "/", H2, "_peaks.narrowPeak"))
colnames(goiH2bed)[1:3] <- c("seqnames", "start", "end")
goiH2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH2bed$seqnames)
goiH2bed_gr <- dt2gr(goiH2bed)
goiH2bed_gt <- gTrack(data = goiH2bed_gr, height = 0.5)
H3 <- "TCGA-44-3918-01A"
goiH3bed <- fread(paste0(atac_beds, H3, "/", H3, "_peaks.narrowPeak"))
colnames(goiH3bed)[1:3] <- c("seqnames", "start", "end")
goiH3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH3bed$seqnames)
goiH3bed_gr <- dt2gr(goiH3bed)
goiH3bed_gt <- gTrack(data = goiH3bed_gr, height = 0.5)
H4 <- "TCGA-MP-A4TE-01A"
goiH4bed <- fread(paste0(atac_beds, H4, "/", H4, "_peaks.narrowPeak"))
colnames(goiH4bed)[1:3] <- c("seqnames", "start", "end")
goiH4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH4bed$seqnames)
goiH4bed_gr <- dt2gr(goiH4bed)
goiH4bed_gt <- gTrack(data = goiH4bed_gr, height = 0.5)

# LUSC ####
LUSC_H1 <- "TCGA-85-A4QR-01A"
goiLUSC_H1bed <- fread(paste0(atac_beds, LUSC_H1, "_atac__peaks.narrowPeak"))
colnames(goiLUSC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H1bed$seqnames)
goiLUSC_H1bed_gr <- dt2gr(goiLUSC_H1bed)
goiLUSC_H1bed_gt <- gTrack(data = goiLUSC_H1bed_gr, height = 0.5)
LUSC_H2 <- "TCGA-33-A4WN-01A"
goiLUSC_H2bed <- fread(paste0(atac_beds, LUSC_H2, "/", LUSC_H2, "_peaks.narrowPeak"))
colnames(goiLUSC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H2bed$seqnames)
goiLUSC_H2bed_gr <- dt2gr(goiLUSC_H2bed)
goiLUSC_H2bed_gt <- gTrack(data = goiLUSC_H2bed_gr, height = 0.5)
LUSC_H3 <- "TCGA-85-A4JB-01A"
goiLUSC_H3bed <- fread(paste0(atac_beds, LUSC_H3, "/", LUSC_H3, "_peaks.narrowPeak"))
colnames(goiLUSC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H3bed$seqnames)
goiLUSC_H3bed_gr <- dt2gr(goiLUSC_H3bed)
goiLUSC_H3bed_gt <- gTrack(data = goiLUSC_H3bed_gr, height = 0.5)
LUSC_H4 <- "TCGA-J1-A4AH-01A"
goiLUSC_H4bed <- fread(paste0(atac_beds, LUSC_H4, "/", LUSC_H4, "_peaks.narrowPeak"))
colnames(goiLUSC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H4bed$seqnames)
goiLUSC_H4bed_gr <- dt2gr(goiLUSC_H4bed)
goiLUSC_H4bed_gt <- gTrack(data = goiLUSC_H4bed_gr, height = 0.5)

# BRCA_TN ####
BRCA_TN_H1 <- "TCGA-C8-A12V-01A"
goiBRCA_TN_H1bed <- fread(paste0(IN_DIR, "/atac_peaks/", BRCA_TN_H1, "_atac__peaks.narrowPeak"))
colnames(goiBRCA_TN_H1bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H1bed$seqnames)
goiBRCA_TN_H1bed_gr <- dt2gr(goiBRCA_TN_H1bed[V9 > -log10(0.05)])
goiBRCA_TN_H1bed_gt <- gTrack(data = reduce(goiBRCA_TN_H1bed_gr), height = 0.5)
BRCA_TN_H2 <- "TCGA-D8-A13Z-01A"
goiBRCA_TN_H2bed <- fread(paste0(atac_beds, BRCA_TN_H2, "/", BRCA_TN_H2, "_peaks.narrowPeak"))
colnames(goiBRCA_TN_H2bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H2bed$seqnames)
goiBRCA_TN_H2bed_gr <- dt2gr(goiBRCA_TN_H2bed[V9 > -log10(0.05)])
goiBRCA_TN_H2bed_gt <- gTrack(data = reduce(goiBRCA_TN_H2bed_gr), height = 0.5)
BRCA_TN_H3 <- "TCGA-AR-A0U4-01A"
goiBRCA_TN_H3bed <- fread(paste0(atac_beds, BRCA_TN_H3, "/", BRCA_TN_H3, "_peaks.narrowPeak"))
colnames(goiBRCA_TN_H3bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H3bed$seqnames)
goiBRCA_TN_H3bed_gr <- dt2gr(goiBRCA_TN_H3bed[V9 > -log10(0.05)])
goiBRCA_TN_H3bed_gt <- gTrack(data = reduce(goiBRCA_TN_H3bed_gr), height = 0.5)
BRCA_TN_H4 <- "TCGA-C8-A12Q-01A"
goiBRCA_TN_H4bed <- fread(paste0(atac_beds, BRCA_TN_H4, "/", BRCA_TN_H4, "_peaks.narrowPeak"))
colnames(goiBRCA_TN_H4bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H4bed$seqnames)
goiBRCA_TN_H4bed_gr <- dt2gr(goiBRCA_TN_H4bed[V9 > -log10(0.05)])
goiBRCA_TN_H4bed_gt <- gTrack(data = reduce(goiBRCA_TN_H4bed_gr), height = 0.5)

# BRCA_HR+ ####
HRBC_H1 <- "TCGA-A7-A0D9-01A"
goiHRBC_H1bed <- fread(paste0(atac_beds, HRBC_H1, "_atac__peaks.narrowPeak"))
colnames(goiHRBC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H1bed$seqnames)
goiHRBC_H1bed_gr <- dt2gr(goiHRBC_H1bed[V9 > -log10(0.05)])
goiHRBC_H1bed_gt <- gTrack(data = reduce(goiHRBC_H1bed_gr), height = 0.5)
HRBC_H2 <- "TCGA-A2-A0YH-01A"
goiHRBC_H2bed <- fread(paste0(atac_beds, HRBC_H2, "_atac__peaks.narrowPeak"))
colnames(goiHRBC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H2bed$seqnames)
goiHRBC_H2bed_gr <- dt2gr(goiHRBC_H2bed[V9 > -log10(0.05)])
goiHRBC_H2bed_gt <- gTrack(data = reduce(goiHRBC_H2bed_gr), height = 0.5)
HRBC_H3 <- "TCGA-A2-A0EX-01A"
goiHRBC_H3bed <- fread(paste0(atac_beds, HRBC_H3, "/", HRBC_H3, "_peaks.narrowPeak"))
colnames(goiHRBC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H3bed$seqnames)
goiHRBC_H3bed_gr <- dt2gr(goiHRBC_H3bed[V9 > -log10(0.05)])
goiHRBC_H3bed_gt <- gTrack(data = reduce(goiHRBC_H3bed_gr), height = 0.5)
HRBC_H4 <- "TCGA-BH-A0B1-01A"
goiHRBC_H4bed <- fread(paste0(atac_beds, HRBC_H4, "/", HRBC_H4, "_peaks.narrowPeak"))
colnames(goiHRBC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H4bed$seqnames)
goiHRBC_H4bed_gr <- dt2gr(goiHRBC_H4bed[V9 > -log10(0.05)])
goiHRBC_H4bed_gt <- gTrack(data = reduce(goiHRBC_H4bed_gr), height = 0.5)

# UCEC ####
UCEC_H1 <- "TCGA-AP-A051-01A"
goiUCEC_H1bed <- fread(paste0(IN_DIR, "/atac_peaks/", UCEC_H1, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H1bed$seqnames)
goiUCEC_H1bed_gr <- dt2gr(goiUCEC_H1bed)
goiUCEC_H1bed_gt <- gTrack(data = goiUCEC_H1bed_gr, height = 0.5)
UCEC_H2 <- "TCGA-B5-A0JT-01A"
goiUCEC_H2bed <- fread(paste0(IN_DIR, "/atac_peaks/", UCEC_H2, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H2bed$seqnames)
goiUCEC_H2bed_gr <- dt2gr(goiUCEC_H2bed)
goiUCEC_H2bed_gt <- gTrack(data = goiUCEC_H2bed_gr, height = 0.5)
UCEC_H3 <- "TCGA-BK-A13B-01A"
goiUCEC_H3bed <- fread(paste0(IN_DIR, "/atac_peaks/", UCEC_H3, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H3bed$seqnames)
goiUCEC_H3bed_gr <- dt2gr(goiUCEC_H3bed)
goiUCEC_H3bed_gt <- gTrack(data = goiUCEC_H3bed_gr, height = 0.5)
UCEC_H4 <- "TCGA-B5-A11Y-01A"
goiUCEC_H4bed <- fread(paste0(IN_DIR, "/atac_peaks/", UCEC_H4, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H4bed$seqnames)
goiUCEC_H4bed_gr <- dt2gr(goiUCEC_H4bed)
goiUCEC_H4bed_gt <- gTrack(data = goiUCEC_H4bed_gr, height = 0.5)

# HNSC ####
HNSC_H1 <- "TCGA-KU-A66S-01A"
goiHNSC_H1bed <- fread(paste0(IN_DIR, "/atac_peaks/", HNSC_H1, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H1bed$seqnames)
goiHNSC_H1bed_gr <- dt2gr(goiHNSC_H1bed)
goiHNSC_H1bed_gt <- gTrack(data = goiHNSC_H1bed_gr, height = 0.5)
HNSC_H2 <- "TCGA-BA-A4IH-01A"
goiHNSC_H2bed <- fread(paste0(IN_DIR, "/atac_peaks/", HNSC_H2, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H2bed$seqnames)
goiHNSC_H2bed_gr <- dt2gr(goiHNSC_H2bed)
goiHNSC_H2bed_gt <- gTrack(data = goiHNSC_H2bed_gr, height = 0.5)
HNSC_H3 <- "TCGA-BA-A8YP-01A"
goiHNSC_H3bed <- fread(paste0(IN_DIR, "/atac_peaks/", HNSC_H3, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H3bed$seqnames)
goiHNSC_H3bed_gr <- dt2gr(goiHNSC_H3bed)
goiHNSC_H3bed_gt <- gTrack(data = goiHNSC_H3bed_gr, height = 0.5)
HNSC_H4 <- "TCGA-H7-A76A-01A"
goiHNSC_H4bed <- fread(paste0(IN_DIR, "/atac_peaks/", HNSC_H4, "_atac__peaks.narrowPeak"))
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

# #sigdiff ####
fishres_A2 <- fread(file.path(IN_DIR, "20250619_NSD3_H3K27ac_ATAC_linDiffs_stats.csv"))
fishres_A2[, seqnames := unlist(strsplit(x = peakID, split = "[:]"))[1], by = peakID]
fishres_A2[, start := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[2]), by = peakID]
fishres_A2[, end := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[3]), by = peakID]
fishres_A2gr <- dt2gr(fishres_A2)
fishres_A_sigs_LUAD <- fishres_A2gr %Q% c(q_val < 0.1 & tstdttype == "LUAD")
fishres_sigs_LUADgt <- gTrack(data = reduce(fishres_A_sigs_LUAD))
fishres_A_sigs_BRCA_HR <- fishres_A2gr %Q% (q_val < 0.1 & tstdttype == "BRCA_HR+")
fishres_sigs_BRCA_HR_gt <- gTrack(data = reduce(fishres_A_sigs_BRCA_HR))
fishres_A_sigs_BRCA_HRN <- fishres_A2gr %Q% (q_val < 0.1 & tstdttype == "BRCA_HR-")
fishres_sigs_BRCA_HRNgt <- gTrack(data = reduce(fishres_A_sigs_BRCA_HRN))
fishres_sigs_LUSCgt <- gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.1 & tstdttype == "LUSC")))
ATACfisher_UCEC_gt <- gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.1 & tstdttype == "UCEC")))
ATACfisher_HNSC_gt <- gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.1 & tstdttype == "HNSC")))

# diffmethyl
xD <- 0.4
dmrs_d4 <- fread(paste0(IN_DIR, "/20260615_NSD3_TM_dmrs_d.csv"))
colnames(dmrs_d4)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d4_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs_d4)), hs1ToHg38_chain))
dmrs_d4 <- as.data.table(dmrs_d4_gr)
dmrs_d4_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -xD])))
dmrs_d4_BRCAgt <- gTrack(dt2gr(dmrs_d4[diff.Methy > xD]))
dmrs_d5_NSCLCgt <- gTrack(dt2gr(dmrs_d4[diff.Methy < -0.5]))
dmrs_d5_BRCAgt <- gTrack(dt2gr(dmrs_d4[diff.Methy > 0.5]))

dmrs <- fread(paste0(IN_DIR, "/20260615_NSD3_TM_dmrs.csv"))
colnames(dmrs)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs)), hs1ToHg38_chain))
dmrs <- gr2dt(dmrs_gr)
dmrs3_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > 0.3]))
dmrs3_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -0.3]))
dmrs4_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > 0.4]))
dmrs4_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -0.4]))
dmrs5_BRCAgt <- gTrack(dt2gr(dmrs[diff.Methy > 0.5]))
dmrs5_NSCLCgt <- gTrack(dt2gr(dmrs[diff.Methy < -0.5]))

# amp/sig tracks (T2T/hs1 -> hg38) ####
# Built from the CNMatm/Fisher tables; each resulting GRanges is lifted
# hs1->hg38 before being wrapped in gTrack, same as the dmrs/dmrs_d4 tracks above.
# BRCA/NSCLC NSD3 amp-frequency tracks
BRCA_NSD3_CNMatm <- readRDS(file.path(IN_DIR, "20260910_BRCA_NSD3_CNMatm_combined_hs1.rds"))
NSCLC_NSD3_CNMatm <- readRDS(file.path(IN_DIR, "20260910_NSCLC_LUSC_NSD3_CNMatm_combined_hs1.rds"))
BRCA_NSD3_CNMatm_hg38gr <- unlist(liftOver(gr.chr(dt2gr(BRCA_NSD3_CNMatm[, 1:10])), hs1ToHg38_chain))
NSCLC_NSD3_CNMatm_hg38gr <- unlist(liftOver(gr.chr(dt2gr(NSCLC_NSD3_CNMatm[, 1:10])), hs1ToHg38_chain))
window_oi <- reduce(dt2gr(rbind(gr2dt(BRCA_NSD3_CNMatm_hg38gr)[AmpPercent > 0.8], gr2dt(NSCLC_NSD3_CNMatm_hg38gr)[AmpPercent > 0.8])))

BRCA_xGtad_per_amp_gt <- gTrack(BRCA_NSD3_CNMatm_hg38gr, y.field = "AmpPercent", col = "brown3", bars = F, y1=1, y0=0)
NSCLC_xGtad_per_amp_gt <- gTrack(NSCLC_NSD3_CNMatm_hg38gr, y.field = "AmpPercent", col = "royalblue3", bars = F, y1=1, y0=0)

# BRCA vs NSCLC NSD3 differential amp-frequency significance track (chr8)
BRCAvsNSCLC_NSD3_fisher_gr <- readRDS(file.path(IN_DIR, "20260910_BRCAvsNSCLC_NSD3_chr8_FisherStats_CN10_cons80_gr.rds"))
BRCAvsNSCLC_NSD3_sig <- BRCAvsNSCLC_NSD3_fisher_gr[BRCAvsNSCLC_NSD3_fisher_gr$Fisher_pval < 0.05]
mcols(BRCAvsNSCLC_NSD3_sig)$tissue <- ifelse(mcols(BRCAvsNSCLC_NSD3_sig)$BRCA_AmpFreq > mcols(BRCAvsNSCLC_NSD3_sig)$NSCLC_AmpFreq, "BRCA", "NSCLC")
colormap_BRCAvsNSCLC_NSD3 <- list(tissue = c(NSCLC = "royalblue3", BRCA = "brown3"))
BRCAvsNSCLC_NSD3_sig_hg38gr <- unlist(liftOver(gr.chr(BRCAvsNSCLC_NSD3_sig), hs1ToHg38_chain))
BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt <- gTrack(BRCAvsNSCLC_NSD3_sig_hg38gr, gr.colorfield = "color", colormaps = colormap_BRCAvsNSCLC_NSD3, height = 0.5)

# NSCLC NSD3 left-vs-right (amp boundary) significance track
NSCLC_NSD3_LvsR_fisher <- fread(file.path(IN_DIR, "20260910_NSCLC_NSD3_LeftVsRight_FisherStats_CN10_cons80.txt"))
NSCLC_NSD3_LvsR_fisher_gr <- dt2gr(NSCLC_NSD3_LvsR_fisher[, ':='(seqnames = fifelse(right_AmpFreq > left_AmpFreq, right_chr, left_chr),
                                             start = fifelse(right_AmpFreq > left_AmpFreq, right_start, left_start),
                                             end = fifelse(right_AmpFreq > left_AmpFreq, right_end, left_end)
)])
NSCLC_NSD3_LvsR_sig <- NSCLC_NSD3_LvsR_fisher_gr[NSCLC_NSD3_LvsR_fisher_gr$Fisher_pval < 0.05]
NSCLC_NSD3_LvsR_sig_hg38gr <- unlist(liftOver(gr.chr(NSCLC_NSD3_LvsR_sig), hs1ToHg38_chain))
NSCLC_NSD3_LvsR_sig_LR_gt <- gTrack(NSCLC_NSD3_LvsR_sig_hg38gr, col = "royalblue3", height = 0.5)

# plot ####
window_oi2 <- window_oi + 1e5
window_oi <- reduce(window_oi2) - 1e5
window_oi3 <- window_oi - 1e5


pdf(file = paste0(OUT_DIR, "/20260919_Fig_S3E_NSD3_H3K27acATAC_BRCAvsNSCLC_with_peaks_and_medianwigs.pdf")
    , paper = "legal", height = 80)
    plot(c(gencode_hg38,
        goiUCEC_H4bed_gt, goiUCEC_H3bed_gt, goiUCEC_H2bed_gt, goiUCEC_H1bed_gt, UCEC_gt,
        MDAMB157_H3K27ac_peaks_gt, CAL120_H3K27ac_peaks_gt, MFM223_H3K27ac_peaks_gt, HRnegCL_gt,
        goiBRCA_TN_H4bed_gt, goiBRCA_TN_H3bed_gt, goiBRCA_TN_H2bed_gt, goiBRCA_TN_H1bed_gt, BRCA_HRneg_gt,
        fishres_sigs_BRCA_HR_gt,
        goiHRBC_H4bed_gt, goiHRBC_H3bed_gt, goiHRBC_H2bed_gt, goiHRBC_H1bed_gt, BRCA_HRpos_gt,
        MDAMB134VI_H3K27ac_peaks_gt, CAMA1_H3K27ac_peaks_gt, ZR751_H3K27ac_peaks_gt, HRposCL_gt,
        dmrs3_BRCAgt, dmrs4_BRCAgt, 
        BRCA_xGtad_per_amp_gt, BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt,NSCLC_NSD3_LvsR_sig_LR_gt, NSCLC_xGtad_per_amp_gt,
        dmrs3_NSCLCgt, dmrs4_NSCLCgt, 
        fishres_sigs_LUSCgt,
        HCC95_H3K27ac_bed_gt, LUSCCL_gt,
        goiLUSC_H4bed_gt, goiLUSC_H3bed_gt, goiLUSC_H2bed_gt, goiLUSC_H1bed_gt, LUSC_gt,
        NCIH1703_H3K27ac_peaks_gt, NCIH2087_H3K27ac_peaks_gt, NCIH520_H3K27ac_peaks_gt, LUADCL_gt,
        goiH4bed_gt, goiH3bed_gt, goiH2bed_gt, goiH1bed_gt, LUAD_gt,
        goiHNSC_H4bed_gt, goiHNSC_H3bed_gt, goiHNSC_H2bed_gt, goiHNSC_H1bed_gt, HNSC_gt
        ), window_oi, cex.label = 0.1)
    plot(c(gencode_hg38,
        goiUCEC_H4bed_gt, goiUCEC_H3bed_gt, goiUCEC_H2bed_gt, goiUCEC_H1bed_gt, UCEC_gt,
        MDAMB157_H3K27ac_peaks_gt, CAL120_H3K27ac_peaks_gt, MFM223_H3K27ac_peaks_gt, HRnegCL_gt,
        goiBRCA_TN_H4bed_gt, goiBRCA_TN_H3bed_gt, goiBRCA_TN_H2bed_gt, goiBRCA_TN_H1bed_gt, BRCA_HRneg_gt,
        fishres_sigs_BRCA_HR_gt,
        goiHRBC_H4bed_gt, goiHRBC_H3bed_gt, goiHRBC_H2bed_gt, goiHRBC_H1bed_gt, BRCA_HRpos_gt,
        MDAMB134VI_H3K27ac_peaks_gt, CAMA1_H3K27ac_peaks_gt, ZR751_H3K27ac_peaks_gt, HRposCL_gt,
        dmrs3_BRCAgt, dmrs4_BRCAgt, dmrs_d4_BRCAgt, 
        BRCA_xGtad_per_amp_gt, BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt,NSCLC_NSD3_LvsR_sig_LR_gt, NSCLC_xGtad_per_amp_gt,
        dmrs3_NSCLCgt, dmrs4_NSCLCgt, dmrs_d4_NSCLCgt, 
        fishres_sigs_LUSCgt,
        HCC95_H3K27ac_bed_gt, LUSCCL_gt,
        goiLUSC_H4bed_gt, goiLUSC_H3bed_gt, goiLUSC_H2bed_gt, goiLUSC_H1bed_gt, LUSC_gt,
        NCIH1703_H3K27ac_peaks_gt, NCIH2087_H3K27ac_peaks_gt, NCIH520_H3K27ac_peaks_gt, LUADCL_gt,
        goiH4bed_gt, goiH3bed_gt, goiH2bed_gt, goiH1bed_gt, LUAD_gt,
        goiHNSC_H4bed_gt, goiHNSC_H3bed_gt, goiHNSC_H2bed_gt, goiHNSC_H1bed_gt, HNSC_gt
    ), window_oi2, cex.label = 0.1)
    plot(c(gencode_hg38,
        goiUCEC_H4bed_gt, goiUCEC_H3bed_gt, goiUCEC_H2bed_gt, goiUCEC_H1bed_gt, UCEC_gt,
        MDAMB157_H3K27ac_peaks_gt, CAL120_H3K27ac_peaks_gt, MFM223_H3K27ac_peaks_gt, HRnegCL_gt,
        goiBRCA_TN_H4bed_gt, goiBRCA_TN_H3bed_gt, goiBRCA_TN_H2bed_gt, goiBRCA_TN_H1bed_gt, BRCA_HRneg_gt,
        fishres_sigs_BRCA_HR_gt,
        goiHRBC_H4bed_gt, goiHRBC_H3bed_gt, goiHRBC_H2bed_gt, goiHRBC_H1bed_gt, BRCA_HRpos_gt,
        MDAMB134VI_H3K27ac_peaks_gt, CAMA1_H3K27ac_peaks_gt, ZR751_H3K27ac_peaks_gt, HRposCL_gt,
        dmrs3_BRCAgt,# dmrs4_BRCAgt, dmrs_d4_BRCAgt, dmrs_d5_BRCAgt,
        BRCA_xGtad_per_amp_gt, BRCAvsNSCLC_NSD3_chr8_sig_track_cons80_gt,NSCLC_NSD3_LvsR_sig_LR_gt, NSCLC_xGtad_per_amp_gt,
        dmrs3_NSCLCgt, #dmrs4_NSCLCgt, dmrs_d4_NSCLCgt, dmrs_d5_NSCLCgt,
        fishres_sigs_LUSCgt,
        HCC95_H3K27ac_bed_gt, LUSCCL_gt,
        goiLUSC_H4bed_gt, goiLUSC_H3bed_gt, goiLUSC_H2bed_gt, goiLUSC_H1bed_gt, LUSC_gt,
        NCIH1703_H3K27ac_peaks_gt, NCIH2087_H3K27ac_peaks_gt, NCIH520_H3K27ac_peaks_gt, LUADCL_gt,
        goiH4bed_gt, goiH3bed_gt, goiH2bed_gt, goiH1bed_gt, LUAD_gt,
        goiHNSC_H4bed_gt, goiHNSC_H3bed_gt, goiHNSC_H2bed_gt, goiHNSC_H1bed_gt, HNSC_gt
    ), window_oi3, cex.label = 0.1)

dev.off()
