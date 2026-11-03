# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(gGnome)
library(rtracklayer)
# library(BSgenome.Hsapiens.UCSC.hg38)
# hg38
CGC <- readRDS("FigS3A/source_data/CGC_hg38.rds")
gencode_hg38 <- readRDS("FigS3A/source_data/gencode_hg38.rds")
conversion_chain <- import.chain(file.path(IN_DIR, "hg38ToHs1.over.chain")) # hg38ToHs1 chain
conversion_chain_hg19ToHg38 <- readRDS(file.path(IN_DIR, 'hg19ToHg38.over.chain.rds'))

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

# 	CCND1 (H3K27ac) BRCAHR- #####
# HCC1143
HCC1143_H3K27ac_peaks <- fread(file.path(IN_DIR, "HCC1143_peaks.narrowPeak"))
colnames(HCC1143_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
HCC1143_H3K27ac_peaks_gt <- gTrack(dt2gr(HCC1143_H3K27ac_peaks), height = 0.5)
# HCC1954
HCC1954_H3K27ac_peaks <- fread("FigS3A/source_data/HCC1954_peaks.narrowPeak")
colnames(HCC1954_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
HCC1954_H3K27ac_peaks_gt <- gTrack(dt2gr(HCC1954_H3K27ac_peaks), height = 0.5)
# MFM223
MFM223_H3K27ac_peaks <- fread(file.path(IN_DIR, "MFM223_peaks.narrowPeak"))
colnames(MFM223_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
MFM223_H3K27ac_peaks_gt <- gTrack(dt2gr(MFM223_H3K27ac_peaks), height = 0.5)
# HDQP1
HDQP1_H3K27ac_peaks <- fread(file.path(IN_DIR, "HDQP1_peaks.narrowPeak"))
colnames(HDQP1_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
HDQP1_H3K27ac_peaks_gt <- gTrack(dt2gr(HDQP1_H3K27ac_peaks), height = 0.5)
# BRCAHR- ChIP median ####
HRnegH3K27ac_avg <- readRDS(file.path(IN_DIR, 'HRneg_ChIPavg_CCND1.rds'))
HRnegH3K27ac_gt <- gTrack(HRnegH3K27ac_avg, y.field = 'score', name = "HRneg", bars = T, y1 = 15)

# 	CCND1 (H3K27ac) BRCAHR+ ####
# ZR751####
ZR751_H3K27ac_peaks <- fread(file.path(IN_DIR, "ZR751_peaks.narrowPeak"))
colnames(ZR751_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
ZR751_H3K27ac_peaks_gt <- gTrack(dt2gr(ZR751_H3K27ac_peaks), height = 0.5)
# CAMA1####
CAMA1_H3K27ac_peaks <- fread(file.path(IN_DIR, "CAMA1_peaks.narrowPeak"))
colnames(CAMA1_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
CAMA1_H3K27ac_peaks_gt <- gTrack(dt2gr(CAMA1_H3K27ac_peaks), height = 0.5)
# MDAMB134VI####
MDAMB134VI_H3K27ac_peaks <- fread(file.path(IN_DIR, "MDAMB134VI_peaks.narrowPeak"))
colnames(MDAMB134VI_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
MDAMB134VI_H3K27ac_peaks_gt <- gTrack(dt2gr(MDAMB134VI_H3K27ac_peaks), height = 0.5)
# MDAMB361####
MDAMB361_H3K27ac_peaks <- fread(file.path(IN_DIR, "MDAMB361_peaks.narrowPeak"))
colnames(MDAMB361_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
MDAMB361_H3K27ac_peaks_gt <- gTrack(dt2gr(MDAMB361_H3K27ac_peaks), height = 0.5)
# BRCAHR+ ChIP median ####
HRposH3K27ac_avg <- readRDS(file.path(IN_DIR, 'HRpos_ChIPavg_CCND1.rds'))
HRposH3K27ac_gt <- gTrack(HRposH3K27ac_avg, y.field = 'score', name = "HRpos", bars = T, y1 = 10)

# 	CCND1 (H3K27ac) LUAD####
# NCIH23####
NCIH23_H3K27ac_peaks <- fread(file.path(IN_DIR, "NCIH23_peaks.narrowPeak"))
colnames(NCIH23_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
NCIH23_H3K27ac_peaks_gt <- gTrack(dt2gr(NCIH23_H3K27ac_peaks), height = 0.5)
# EKVX####
EKVX_H3K27ac_peaks <- fread(file.path(IN_DIR, "EKVX_peaks.narrowPeak"))
colnames(EKVX_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
EKVX_H3K27ac_peaks_gt <- gTrack(dt2gr(EKVX_H3K27ac_peaks), height = 0.5)
# NCIH3122
NCIH3122_H3K27ac_peaks <- fread(file.path(IN_DIR, "NCIH3122_H3K27ac_peaks.narrowPeak"))
colnames(NCIH3122_H3K27ac_peaks)[1:3] <- c("chr", "start", "end")
NCIH3122_H3K27ac_peaks_gt <- gTrack(dt2gr(NCIH3122_H3K27ac_peaks), height = 0.5)
# LUAD ChIP median ####
LUAD_H3K27ac_avg <- readRDS(file.path(IN_DIR, 'LUAD_ChIPavg_CCND1.rds'))
LUAD_H3K27ac_gt <- gTrack(LUAD_H3K27ac_avg, y.field = 'score', name = "LUAD", bars = T, y1 = 15)

# 	CCND1 (H3K27ac) LUSC####
# HCC95####
HCC95_H3K27ac_bw_hg38 <- readRDS(file.path(IN_DIR, "GSM2784723_H3K27ac-HCC95-WT-A.bw_hg38.rds"))
HCC95_H3K27ac_gt <- gTrack(data = gr.nochr(HCC95_H3K27ac_bw_hg38), y.field = "score", bars = T)
HCC95_H3K27ac_bed_hg19 <- fread(file.path(IN_DIR, "GSM4079747_HCC95_H3K27ac_peaks_hg19.bed"))
colnames(HCC95_H3K27ac_bed_hg19)[1:3] <- c("seqnames", "start", "end")
HCC95_H3K27ac_bed_hg38 <- unlist(liftOver(dt2gr(HCC95_H3K27ac_bed_hg19), chain = conversion_chain_hg19ToHg38))
write.csv(gr2dt(HCC95_H3K27ac_bed_hg38), file.path(OUT_DIR, "GSM4079747_HCC95_H3K27ac_peaks_hg38.bed"), row.names = F)
HCC95_H3K27ac_bed_gt <- gTrack(data = gr.nochr(HCC95_H3K27ac_bed_hg38), height = 0.5)

# Median ATAC tracks ####
# LUAD ####
LUAD_wig <- import(file.path(IN_DIR, 'LUAD_median.wig'))
LUAD_gt <- gTrack(LUAD_wig, y.field = 'score', name = "LUAD", bars = T, y1 = 150)
# LUSC ####
LUSC_wig <- import(file.path(IN_DIR, 'LUSC_median.wig'))
LUSC_gt <- gTrack(LUSC_wig, y.field = 'score', name = "LUSC", bars = T, y1 = 150)
# UCEC ####
UCEC_wig <- import(file.path(IN_DIR, 'UCEC_median.wig'))
UCEC_gt <- gTrack(UCEC_wig, y.field = 'score', name = "UCEC", bars = T, y1 = 150)
# HNSC ####
HNSC_wig <- import(file.path(IN_DIR, 'HNSC_median.wig'))
HNSC_gt <- gTrack(HNSC_wig, y.field = 'score', name = "HNSC", bars = T, y1 = 150)
# BRCA_HRneg ####
BRCA_HRneg_wig <- import(file.path(IN_DIR, 'BRCA_TN_median.wig'))
BRCA_HRneg_gt <- gTrack(BRCA_HRneg_wig, y.field = 'score', name = "BRCA_HRneg", bars = T, y1 = 150)
# BRCA_HRpos ####
BRCA_HRpos_wig <- import(file.path(IN_DIR, 'HRpos_median.wig'))
BRCA_HRpos_gt <- gTrack(BRCA_HRpos_wig, y.field = 'score', name = "BRCA_HRpos", bars = T, y1 = 150)

# NSCLC ATAC peak beds ####
# LUAD ####
H1 <- "TCGA-MP-A4TF-01A"
goiH1bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", H1, "/", H1, "_peaks.narrowPeak"))
colnames(goiH1bed)[1:3] <- c("seqnames", "start", "end")
goiH1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH1bed$seqnames)
goiH1bed_gr <- dt2gr(goiH1bed)
goiH1bed_gt <- gTrack(data = goiH1bed_gr, height = 0.5)
H2 <- "TCGA-MP-A4TE-01A"
goiH2bed <- fread(paste0(IN_DIR, "/ATACbeds/out/", H2, "/", H2, "_peaks.narrowPeak")) # TCGA-MP-A4TE-01A bed missing
colnames(goiH2bed)[1:3] <- c("seqnames", "start", "end")
goiH2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH2bed$seqnames)
goiH2bed_gr <- dt2gr(goiH2bed)
goiH2bed_gt <- gTrack(data = goiH2bed_gr, height = 0.5)
H3 <- "TCGA-L9-A8F4-01A"
goiH3bed <- fread(paste0(IN_DIR, "/ATACbeds/out/", H3, "/", H3, "_peaks.narrowPeak")) # TCGA-L9-A8F4-01A bed missing
colnames(goiH3bed)[1:3] <- c("seqnames", "start", "end")
goiH3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH3bed$seqnames)
goiH3bed_gr <- dt2gr(goiH3bed)
goiH3bed_gt <- gTrack(data = goiH3bed_gr, height = 0.5)
H4 <- "TCGA-44-6146-01A"
goiH4bed <- fread(paste0(IN_DIR, "/ATACbeds/out/", H4, "/", H4, "_peaks.narrowPeak")) # TCGA-L9-A8F4-01A bed missing
colnames(goiH4bed)[1:3] <- c("seqnames", "start", "end")
goiH4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiH4bed$seqnames)
goiH4bed_gr <- dt2gr(goiH4bed)
goiH4bed_gt <- gTrack(data = goiH4bed_gr, height = 0.5)

# UCEC ####
UCEC_H1 <- "TCGA-AP-A051-01A"
goiUCEC_H1bed <- fread(paste0(IN_DIR, "/atac_peaks/", UCEC_H1, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H1bed$seqnames)
goiUCEC_H1bed_gr <- dt2gr(goiUCEC_H1bed)
goiUCEC_H1bed_gt <- gTrack(data = goiUCEC_H1bed_gr, height = 0.5)

UCEC_H2 <- "TCGA-BG-A0M8-01A"
goiUCEC_H2bed <- fread(paste0(IN_DIR, "/atac_peaks/", UCEC_H2, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H2bed$seqnames)
goiUCEC_H2bed_gr <- dt2gr(goiUCEC_H2bed)
goiUCEC_H2bed_gt <- gTrack(data = goiUCEC_H2bed_gr, height = 0.5)

UCEC_H3 <- "TCGA-BG-A0MS-01A"
goiUCEC_H3bed <- fread(paste0(IN_DIR, "/atac_peaks/", UCEC_H3, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H3bed$seqnames)
goiUCEC_H3bed_gr <- dt2gr(goiUCEC_H3bed)
goiUCEC_H3bed_gt <- gTrack(data = goiUCEC_H3bed_gr, height = 0.5)

UCEC_H4 <- "TCGA-B5-A0K9-01A"
goiUCEC_H4bed <- fread(paste0(IN_DIR, "/atac_peaks/", UCEC_H4, "_atac__peaks.narrowPeak"))
colnames(goiUCEC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiUCEC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiUCEC_H4bed$seqnames)
goiUCEC_H4bed_gr <- dt2gr(goiUCEC_H4bed)
goiUCEC_H4bed_gt <- gTrack(data = goiUCEC_H4bed_gr, height = 0.5)

# HNSC ####
HNSC_H1 <- "TCGA-CN-A49B-01A"
goiHNSC_H1bed <- fread(paste0(IN_DIR, "/atac_peaks/", HNSC_H1, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H1bed$seqnames)
goiHNSC_H1bed_gr <- dt2gr(goiHNSC_H1bed)
goiHNSC_H1bed_gt <- gTrack(data = goiHNSC_H1bed_gr, height = 0.5)

HNSC_H2 <- "TCGA-D6-A6EP-01A"
goiHNSC_H2bed <- fread(paste0(IN_DIR, "/atac_peaks/", HNSC_H2, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H2bed$seqnames)
goiHNSC_H2bed_gr <- dt2gr(goiHNSC_H2bed)
goiHNSC_H2bed_gt <- gTrack(data = goiHNSC_H2bed_gr, height = 0.5)

HNSC_H3 <- "TCGA-F7-A624-01A"
goiHNSC_H3bed <- fread(paste0(IN_DIR, "/atac_peaks/", HNSC_H3, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H3bed$seqnames)
goiHNSC_H3bed_gr <- dt2gr(goiHNSC_H3bed)
goiHNSC_H3bed_gt <- gTrack(data = goiHNSC_H3bed_gr, height = 0.5)

HNSC_H4 <- "TCGA-KU-A66S-01A"
goiHNSC_H4bed <- fread(paste0(IN_DIR, "/atac_peaks/", HNSC_H4, "_atac__peaks.narrowPeak"))
colnames(goiHNSC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiHNSC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHNSC_H4bed$seqnames)
goiHNSC_H4bed_gr <- dt2gr(goiHNSC_H4bed)
goiHNSC_H4bed_gt <- gTrack(data = goiHNSC_H4bed_gr, height = 0.5)

# LUSC ####
LUSC_H1 <- "TCGA-MF-A522-01A"
goiLUSC_H1bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", LUSC_H1, "/", LUSC_H1, "_peaks.narrowPeak"))
colnames(goiLUSC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H1bed$seqnames)
goiLUSC_H1bed_gr <- dt2gr(goiLUSC_H1bed)
goiLUSC_H1bed_gt <- gTrack(data = goiLUSC_H1bed_gr, height = 0.5)
LUSC_H2 <- "TCGA-98-A53J-01A"
goiLUSC_H2bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", LUSC_H2, "/", LUSC_H2, "_peaks.narrowPeak"))
colnames(goiLUSC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H2bed$seqnames)
goiLUSC_H2bed_gr <- dt2gr(goiLUSC_H2bed)
goiLUSC_H2bed_gt <- gTrack(data = goiLUSC_H2bed_gr, height = 0.5)
LUSC_H3 <- "TCGA-77-A5G6-01A"
goiLUSC_H3bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", LUSC_H3, "/", LUSC_H3, "_peaks.narrowPeak"))
colnames(goiLUSC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H3bed$seqnames)
goiLUSC_H3bed_gr <- dt2gr(goiLUSC_H3bed)
goiLUSC_H3bed_gt <- gTrack(data = goiLUSC_H3bed_gr, height = 0.5)
LUSC_H4 <- "TCGA-XC-AA0X-01A"
goiLUSC_H4bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", LUSC_H4, "/", LUSC_H4, "_peaks.narrowPeak"))
colnames(goiLUSC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiLUSC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiLUSC_H4bed$seqnames)
goiLUSC_H4bed_gr <- dt2gr(goiLUSC_H4bed)
goiLUSC_H4bed_gt <- gTrack(data = goiLUSC_H4bed_gr, height = 0.5)

# BRCA ####
# BRCA_TN ####
BRCA_TN_H1 <- "TCGA-AR-A0U4-01A" # TPM: 16928.66  TCGA-AR-A0U0-01A ->TPM: 19 425.7388
goiBRCA_TN_H1bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", BRCA_TN_H1, "/", BRCA_TN_H1, "_peaks.narrowPeak"))
colnames(goiBRCA_TN_H1bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H1bed$seqnames)
goiBRCA_TN_H1bed_gr <- dt2gr(goiBRCA_TN_H1bed[V9 > -log10(0.05)])
goiBRCA_TN_H1bed_gt <- gTrack(data = reduce(goiBRCA_TN_H1bed_gr), height = 0.5)
BRCA_TN_H2 <- "TCGA-A8-A094-01A" # TPM 11 416.3
goiBRCA_TN_H2bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", BRCA_TN_H2, "/", BRCA_TN_H2, "_peaks.narrowPeak")) # rerun MACS2
colnames(goiBRCA_TN_H2bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H2bed$seqnames)
goiBRCA_TN_H2bed_gr <- dt2gr(goiBRCA_TN_H2bed[V9 > -log10(0.05)])
goiBRCA_TN_H2bed_gt <- gTrack(data = reduce(goiBRCA_TN_H2bed_gr), height = 0.5)
BRCA_TN_H3 <- "TCGA-C8-A12Q-01A" # 22 851.3
goiBRCA_TN_H3bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", BRCA_TN_H3, "/", BRCA_TN_H3, "_peaks.narrowPeak"))
colnames(goiBRCA_TN_H3bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H3bed$seqnames)
goiBRCA_TN_H3bed_gr <- dt2gr(goiBRCA_TN_H3bed[V9 > -log10(0.05)])
goiBRCA_TN_H3bed_gt <- gTrack(data = reduce(goiBRCA_TN_H3bed_gr), height = 0.5)
BRCA_TN_H4 <- "TCGA-AO-A12F-01A" # TPM: 19 610.7061
goiBRCA_TN_H4bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", BRCA_TN_H4, "/", BRCA_TN_H4, "_peaks.narrowPeak"))
colnames(goiBRCA_TN_H4bed)[1:3] <- c("seqnames", "start", "end")
goiBRCA_TN_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiBRCA_TN_H4bed$seqnames)
goiBRCA_TN_H4bed_gr <- dt2gr(goiBRCA_TN_H4bed[V9 > -log10(0.05)])
goiBRCA_TN_H4bed_gt <- gTrack(data = reduce(goiBRCA_TN_H4bed_gr), height = 0.5)

# BRCA_HR+ ####
HRBC_H1 <- "TCGA-C8-A12M-01A" # CN(1.5)
goiHRBC_H1bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", HRBC_H1, "/", HRBC_H1, "_peaks.narrowPeak"))
colnames(goiHRBC_H1bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H1bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H1bed$seqnames)
goiHRBC_H1bed_gr <- dt2gr(goiHRBC_H1bed[V9 > -log10(0.05)])
goiHRBC_H1bed_gt <- gTrack(data = reduce(goiHRBC_H1bed_gr), height = 0.5)
HRBC_H2 <- "TCGA-C8-A12Y-01A" # CN(1.76)
goiHRBC_H2bed <- fread(paste0(IN_DIR, "/ATACbeds/out/", HRBC_H2, "/", HRBC_H2, "_peaks.narrowPeak"))
colnames(goiHRBC_H2bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H2bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H2bed$seqnames)
goiHRBC_H2bed_gr <- dt2gr(goiHRBC_H2bed[V9 > -log10(0.05)])
goiHRBC_H2bed_gt <- gTrack(data = reduce(goiHRBC_H2bed_gr), height = 0.5)
HRBC_H3 <- "TCGA-C8-A12O-01A" # high exp & high CN (2.2)
goiHRBC_H3bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", HRBC_H3, "/", HRBC_H3, "_peaks.narrowPeak"))
colnames(goiHRBC_H3bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H3bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H3bed$seqnames)
goiHRBC_H3bed_gr <- dt2gr(goiHRBC_H3bed[V9 > -log10(0.05)])
goiHRBC_H3bed_gt <- gTrack(data = reduce(goiHRBC_H3bed_gr), height = 0.5)
HRBC_H4 <- "TCGA-A2-A0YF-01A" #
goiHRBC_H4bed <- fread(paste0(IN_DIR, "/ATACbeds/out2/", HRBC_H4, "/", HRBC_H4, "_peaks.narrowPeak"))
colnames(goiHRBC_H4bed)[1:3] <- c("seqnames", "start", "end")
goiHRBC_H4bed$seqnames <- gsub(pattern = "chr", replacement = "", goiHRBC_H4bed$seqnames)
goiHRBC_H4bed_gr <- dt2gr(goiHRBC_H4bed[V9 > -log10(0.05)])
goiHRBC_H4bed_gt <- gTrack(data = reduce(goiHRBC_H4bed_gr), height = 0.5)

# functional val####
GERCARres <- as.data.table(read.csv2(file.path(IN_DIR, "GERCARres.csv")))
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
fishres_A_sigs_LUAD <- fishres_A2gr %Q% c(q_val < 0.05 & tstdttype == "LUAD")
fishres_sigs_LUADgt <- gTrack(data = reduce(fishres_A_sigs_LUAD))
fishres_A_sigs_BRCA_HR <- fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "BRCA_HR+")
fishres_sigs_BRCA_HR_gt <- gTrack(data = reduce(fishres_A_sigs_BRCA_HR))
fishres_A_sigs_BRCA_HRN <- fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "BRCA_HR-")
fishres_sigs_BRCA_HRNgt <- gTrack(data = reduce(fishres_A_sigs_BRCA_HRN))
ATACfisher_LUSC_gt <- gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "LUSC")))
ATACfisher_UCEC_gt <- gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "UCEC")))
ATACfisher_HNSC_gt <- gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "HNSC")))

xD <- 0.4
dmrs_d4 <- fread(paste0(IN_DIR, "/20260512_CCND1_TM_dmrs_d.csv"))
colnames(dmrs_d4)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d4_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs_d4)), hs1ToHg38_chain))
dmrs_d4 <- as.data.table(dmrs_d4_gr)
dmrs_d4_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -xD])), ygap =0.5)
dmrs_d4_BRCAgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy > xD])), ygap =0.5)
dmrs_d5_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -0.5])), ygap =0.5)
dmrs_d5_BRCAgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy > 0.5])), ygap =0.5)
dmrs_d6_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -0.6])), ygap =0.5)
dmrs_d6_BRCAgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy > 0.6])), ygap =0.5)


BRCA_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, "/20260907_Fig.2C_BRCA_CCND1_AmpFreq_track.rds"))
BRCA_xGtad_per_amp_gt_gr <- unlist(liftOver(BRCA_xGtad_per_amp_gt@data[[1]], hs1ToHg38_chain))
BRCA_xGtad_per_amp_gt <- gTrack(data = BRCA_xGtad_per_amp_gt_gr, y.field = "AmpFreq", name = "%amp", col = "red")
BRCA_xGtad_per_amp_gt$bars <- F

BRCA_CCND1_LvsR_sig_LR_gt <- readRDS(paste0(IN_DIR, "/20260907c_Fig.2C_BRCA_CCND1_LvsR_sig_LR_cons80.rds")) 
BRCA_CCND1_LvsR_sig_LR_gt_gr <- unlist(liftOver(BRCA_CCND1_LvsR_sig_LR_gt@data[[1]], hs1ToHg38_chain))
BRCA_CCND1_LvsR_sig_LR_gt <- gTrack(data = BRCA_CCND1_LvsR_sig_LR_gt_gr, name = "nlgp", col = "red", y0 = 0, y1 = 1, , height = 0.5)

BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- readRDS(paste0(IN_DIR, "/20260907c_Fig.2C_BRCAvsNSCLC_CCND1_chr11_sig_track_cons80.rds")) 
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr <- unlist(liftOver(BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt@data[[1]], hs1ToHg38_chain))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- gTrack(data = BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr, name = "nlgp", col = "blue", height = 0.5)

NSCLC_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, "/20260907_Fig.2C_NSCLC_CCND1_AmpFreq_track.rds"))
NSCLC_xG_amps_grc <- NSCLC_xGtad_per_amp_gt@data[[1]]
NSCLC_xGtad_amps_hg38gr <- unlist(liftOver(NSCLC_xG_amps_grc, hs1ToHg38_chain))
NSCLC_xGtad_per_amp_gt <- gTrack(data = NSCLC_xGtad_amps_hg38gr, y.field = "AmpFreq", name = "%amp", col = "blue", y0 = 0, y1 = 1)
NSCLC_xGtad_per_amp_gt$bars <- F

# plot ####
window_oi2 <- Oncgr + 5e5
window_oi3 <- Oncgr + 6e5

pdf(
    file = paste0(OUT_DIR, "/20260911_Fig_S3C_CCND1_with_peaks_and_medianwigs.pdf"),
    paper = "legal", height = 80
)
plot(c(gencode_hg38,
    ATACfisher_UCEC_gt, goiUCEC_H4bed_gt, goiUCEC_H3bed_gt, goiUCEC_H2bed_gt, goiUCEC_H1bed_gt, UCEC_gt,
    fishres_sigs_BRCA_HRNgt,
    MFM223_H3K27ac_peaks_gt, HCC1954_H3K27ac_peaks_gt, HCC1143_H3K27ac_peaks_gt, HDQP1_H3K27ac_peaks_gt, HRnegH3K27ac_gt,
    goiBRCA_TN_H4bed_gt, goiBRCA_TN_H3bed_gt, goiBRCA_TN_H2bed_gt, goiBRCA_TN_H1bed_gt, BRCA_HRneg_gt,
    fishres_sigs_BRCA_HR_gt,
    goiHRBC_H4bed_gt, goiHRBC_H3bed_gt, goiHRBC_H2bed_gt, goiHRBC_H1bed_gt, BRCA_HRpos_gt,
    MDAMB361_H3K27ac_peaks_gt, MDAMB134VI_H3K27ac_peaks_gt, CAMA1_H3K27ac_peaks_gt, ZR751_H3K27ac_peaks_gt, HRposH3K27ac_gt,
    dmrs_d5_BRCAgt,
    BRCA_xGtad_per_amp_gt, BRCA_CCND1_LvsR_sig_LR_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
    NSCLC_xGtad_per_amp_gt,
    dmrs_d5_NSCLCgt, ATACfisher_LUSC_gt, HCC95_H3K27ac_bed_gt, HCC95_H3K27ac_gt,
    goiLUSC_H4bed_gt, goiLUSC_H3bed_gt, goiLUSC_H2bed_gt, goiLUSC_H1bed_gt, LUSC_gt,
    fishres_sigs_LUADgt,
    NCIH3122_H3K27ac_peaks_gt, EKVX_H3K27ac_peaks_gt, NCIH23_H3K27ac_peaks_gt, LUAD_H3K27ac_gt,
    goiH4bed_gt, goiH3bed_gt, goiH2bed_gt, goiH1bed_gt, LUAD_gt,
    ATACfisher_HNSC_gt,
    goiHNSC_H4bed_gt, goiHNSC_H3bed_gt, goiHNSC_H2bed_gt, goiHNSC_H1bed_gt, HNSC_gt
), window_oi2, cex.label = 0.1)
plot(c(gencode_hg38,
    ATACfisher_UCEC_gt, goiUCEC_H4bed_gt, goiUCEC_H3bed_gt, goiUCEC_H2bed_gt, goiUCEC_H1bed_gt, UCEC_gt,
    fishres_sigs_BRCA_HRNgt,
    MFM223_H3K27ac_peaks_gt, HCC1954_H3K27ac_peaks_gt, HCC1143_H3K27ac_peaks_gt, HDQP1_H3K27ac_peaks_gt, HRnegH3K27ac_gt,
    goiBRCA_TN_H4bed_gt, goiBRCA_TN_H3bed_gt, goiBRCA_TN_H2bed_gt, goiBRCA_TN_H1bed_gt, BRCA_HRneg_gt,
    fishres_sigs_BRCA_HR_gt,
    goiHRBC_H4bed_gt, goiHRBC_H3bed_gt, goiHRBC_H2bed_gt, goiHRBC_H1bed_gt, BRCA_HRpos_gt,
    MDAMB361_H3K27ac_peaks_gt, MDAMB134VI_H3K27ac_peaks_gt, CAMA1_H3K27ac_peaks_gt, ZR751_H3K27ac_peaks_gt, HRposH3K27ac_gt,
    dmrs_d4_BRCAgt,dmrs_d5_BRCAgt,dmrs_d6_BRCAgt,
    BRCA_xGtad_per_amp_gt, BRCA_CCND1_LvsR_sig_LR_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
    NSCLC_xGtad_per_amp_gt,
    dmrs_d4_NSCLCgt, dmrs_d5_NSCLCgt,dmrs_d6_NSCLCgt,
    ATACfisher_LUSC_gt, HCC95_H3K27ac_bed_gt, HCC95_H3K27ac_gt,
    goiLUSC_H4bed_gt, goiLUSC_H3bed_gt, goiLUSC_H2bed_gt, goiLUSC_H1bed_gt, LUSC_gt,
    fishres_sigs_LUADgt,
    NCIH3122_H3K27ac_peaks_gt, EKVX_H3K27ac_peaks_gt, NCIH23_H3K27ac_peaks_gt, LUAD_H3K27ac_gt,
    goiH4bed_gt, goiH3bed_gt, goiH2bed_gt, goiH1bed_gt, LUAD_gt,
    ATACfisher_HNSC_gt,
    goiHNSC_H4bed_gt, goiHNSC_H3bed_gt, goiHNSC_H2bed_gt, goiHNSC_H1bed_gt, HNSC_gt
), window_oi2, cex.label = 0.1)
plot(c(gencode_hg38,
    ATACfisher_UCEC_gt, goiUCEC_H4bed_gt, goiUCEC_H3bed_gt, goiUCEC_H2bed_gt, goiUCEC_H1bed_gt, UCEC_gt,
    fishres_sigs_BRCA_HRNgt,
    MFM223_H3K27ac_peaks_gt, HCC1954_H3K27ac_peaks_gt, HCC1143_H3K27ac_peaks_gt, HDQP1_H3K27ac_peaks_gt, HRnegH3K27ac_gt,
    goiBRCA_TN_H4bed_gt, goiBRCA_TN_H3bed_gt, goiBRCA_TN_H2bed_gt, goiBRCA_TN_H1bed_gt, BRCA_HRneg_gt,
    GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt, fishres_sigs_BRCA_HR_gt,
    goiHRBC_H4bed_gt, goiHRBC_H3bed_gt, goiHRBC_H2bed_gt, goiHRBC_H1bed_gt, BRCA_HRpos_gt,
    MDAMB361_H3K27ac_peaks_gt, MDAMB134VI_H3K27ac_peaks_gt, CAMA1_H3K27ac_peaks_gt, ZR751_H3K27ac_peaks_gt, HRposH3K27ac_gt,
    dmrs_d4_BRCAgt,
    BRCA_xGtad_per_amp_gt, BRCA_CCND1_LvsR_sig_LR_gt, BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt,
    NSCLC_xGtad_per_amp_gt,
    dmrs_d4_NSCLCgt, ATACfisher_LUSC_gt, HCC95_H3K27ac_bed_gt, HCC95_H3K27ac_gt,
    goiLUSC_H4bed_gt, goiLUSC_H3bed_gt, goiLUSC_H2bed_gt, goiLUSC_H1bed_gt, LUSC_gt,
    fishres_sigs_LUADgt,
    NCIH3122_H3K27ac_peaks_gt, EKVX_H3K27ac_peaks_gt, NCIH23_H3K27ac_peaks_gt, LUAD_H3K27ac_gt,
    goiH4bed_gt, goiH3bed_gt, goiH2bed_gt, goiH1bed_gt, LUAD_gt,
    ATACfisher_HNSC_gt,
    goiHNSC_H4bed_gt, goiHNSC_H3bed_gt, goiHNSC_H2bed_gt, goiHNSC_H1bed_gt, HNSC_gt
), window_oi2, cex.label = 0.1)
dev.off()
