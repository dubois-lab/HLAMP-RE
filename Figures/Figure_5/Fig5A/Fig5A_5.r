# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(GenomicRanges)
library(GxG)
library(parallel)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)


## Load core annotation
CGC <- readRDS('../../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
gencode_hg38 <- readRDS('../../supp_Figure_3/FigS3A/source_data/gencode_hg38.rds')
hs1ToHg38_chain <- import.chain(file.path(IN_DIR, 'hs1ToHg38.over.chain'))

# HMF&PCAWG ampdist
BRCA_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, "/20260907_Fig.2C_BRCA_CCND1_AmpFreq_track.rds"))
BRCA_xGtad_per_amp_gt_gr <- unlist(liftOver(BRCA_xGtad_per_amp_gt@data[[1]], hs1ToHg38_chain))
BRCA_xGtad_per_amp_gt <- gTrack(data = BRCA_xGtad_per_amp_gt_gr, y.field = "AmpFreq", name = "%amp", col = "red")
BRCA_xGtad_per_amp_gt$bars <- F

BRCA_CCND1_LvsR_sig_LR_gt <- readRDS(paste0(IN_DIR, "/20260907c_Fig.2C_BRCA_CCND1_LvsR_sig_LR_cons80.rds")) # nolint: line_length_linter.
BRCA_CCND1_LvsR_sig_LR_gt_gr <- unlist(liftOver(BRCA_CCND1_LvsR_sig_LR_gt@data[[1]], hs1ToHg38_chain))
BRCA_CCND1_LvsR_sig_LR_gt <- gTrack(data = reduce(BRCA_CCND1_LvsR_sig_LR_gt_gr), name = "nlgp", col = "red", height = 0.5)

BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- readRDS(paste0(IN_DIR, "/20260907c_Fig.2C_BRCAvsNSCLC_CCND1_chr11_sig_track_cons80.rds")) # nolint: line_length_linter.
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr <- unlist(liftOver(BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt@data[[1]], hs1ToHg38_chain))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- gTrack(data = BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr, name = "nlgp", col = "black")

NSCLC_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, "/20260907_Fig.2C_NSCLC_CCND1_AmpFreq_track.rds"))
NSCLC_xG_amps_grc <- NSCLC_xGtad_per_amp_gt@data[[1]]
NSCLC_xGtad_amps_hg38gr <- unlist(liftOver(NSCLC_xG_amps_grc, hs1ToHg38_chain))
NSCLC_xGtad_per_amp_gt <- gTrack(data = NSCLC_xGtad_amps_hg38gr, y.field = "AmpFreq", name = "%amp", col = "blue", y0 = 0, y1 = 1)
NSCLC_xGtad_per_amp_gt$bars <- F


xG <- "CCND1"
Oncgr <- CGC %Q% (Gene_Symbol == xG)
window_oi <- Oncgr + 12e5


## TF binding data
MCF7_ESR1_hg38 <- gTrack(file.path(IN_DIR, 'MCF7_ESR1_hg38_FE_ENCFF237WTX.bigWig')
                         , bars = T, y1 = 100, name = 'ESR1')
MCF7_FOXA1_hg38 <- gTrack(file.path(IN_DIR, 'MCF7_FOXA1_hg38_FE_ENCFF795BHZ.bigWig')
                          , bars = T, y1 = 40, name = 'FOXA1')
MCF7_GATA3_hg38 <- gTrack(file.path(IN_DIR, 'MCF7_GATA3_hg38_FE_ENCFF384CPN.bigWig'),
                          bars = T, y1 = 40, name = 'GATA3')
MCF7_TRPS1_hg38 <- gTrack(file.path(IN_DIR, 'MCF7_TRPS1_hg38_FE.bw'),
                          bars = T, y1 = 40, name = 'TRPS1')
MCF7_SPDEF_hg38 <- gTrack(file.path(IN_DIR, 'MCF7-SPDEF-ENCFF114YMT-hg38.bigWig'),
                          bars = T, y1 = 40, name = 'SPDEF')
MCF7_TFAP2C_hg38 <- gTrack(file.path(IN_DIR, 'MCF7_TFAP2C_hg38_FE.bw'),
                          bars = T, y1 = 40, name = 'TFAP2C')
MCF7_ESR1_hg38_bed <-fread(file.path(IN_DIR, 'MCF7_ESR1_ENCFF138XTJ.bed'))
colnames(MCF7_ESR1_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_ESR1_hg38_bed_gt <- gTrack(dt2gr(MCF7_ESR1_hg38_bed), height = 0.5)
MCF7_FOXA1_hg38_bed <- fread(file.path(IN_DIR, 'MCF7_FOXA1_ENCFF112JVK.bed'))
colnames(MCF7_FOXA1_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_FOXA1_hg38_bed_gt <- gTrack(dt2gr(MCF7_FOXA1_hg38_bed), height = 0.5)
MCF7_GATA3_hg38_bed <- fread(file.path(IN_DIR, 'MCF7-GATA3-ENCFF217KJR-hg38.bed'))
colnames(MCF7_GATA3_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_GATA3_hg38_bed_gt <- gTrack(dt2gr(MCF7_GATA3_hg38_bed), height = 0.5)
MCF7_TRPS1_hg38_bed <- fread(file.path(IN_DIR, 'MCF7_TRPS1_hg38_peaks.narrowPeak'))
colnames(MCF7_TRPS1_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_TRPS1_hg38_bed_gt <- gTrack(dt2gr(MCF7_TRPS1_hg38_bed), height = 0.5)
MCF7_TFAP2C_hg38_bed <- fread(file.path(IN_DIR, 'MCF7_TFAP2C_hg38_peaks.narrowPeak'))
colnames(MCF7_TFAP2C_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_TFAP2C_hg38_bed_gt <- gTrack(dt2gr(MCF7_TFAP2C_hg38_bed), height = 0.5)
MCF7_SPDEF_hg38_bed <- fread(file.path(IN_DIR, 'MCF7-SPDEF-ENCFF464GFC-hg38.bed'))
colnames(MCF7_SPDEF_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_SPDEF_hg38_bed_gt <- gTrack(dt2gr(MCF7_SPDEF_hg38_bed), height = 0.5)

# sigdiff ####
ATACfisher <- fread(file.path(IN_DIR, "20260614_CCND1_H3K27ac&ATAC_linDiffs_stats.csv"))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = "[:]"))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[3]), by = peakID]
fishres_A2gr <- dt2gr(ATACfisher)
fishres_A_sigs_BRCA_HR <- fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "BRCA_HR+")
fishres_sigs_BRCA_HR_gt <- gTrack(data = reduce(fishres_A_sigs_BRCA_HR), ygap = 0.5)
ATACfisher_LUSC_gt <- gTrack(data = reduce(fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "LUSC")), ygap = 0.5)

xD <- 0.4
dmrs_d4 <- fread(paste0(IN_DIR, "/20260512_CCND1_TM_dmrs_d.csv"))
colnames(dmrs_d4)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d4_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs_d4)), hs1ToHg38_chain))
dmrs_d4 <- as.data.table(dmrs_d4_gr)
dmrs_d4_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -xD])), ygap = 0.5)
dmrs_d4_BRCAgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy > xD])), ygap = 0.5)
dmrs_d5_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -0.5])), ygap = 0.5)
dmrs_d5_BRCAgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy > 0.5])), ygap = 0.5)
dmrs_d6_NSCLCgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy < -0.6])), ygap = 0.5)
dmrs_d6_BRCAgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy > 0.6])), ygap = 0.5)

# CCND1 BRCA virtual-4C combined (Fisher) p-values (hs1 -> hg38)
comb_Pvals_v4C_CCND1_BRCA_gr <- readRDS(file.path(IN_DIR, '20260702_comb_Pvals_v4C_CCND1_BRCA_gr.rds'))
comb_Pvals_v4C_CCND1_BRCA_gr_hg38 <- grl.unlist(liftOver(comb_Pvals_v4C_CCND1_BRCA_gr, hs1ToHg38_chain))
comb_Pvals_v4C_CCND1_BRCA_gt <- gTrack(comb_Pvals_v4C_CCND1_BRCA_gr_hg38, y.field = 'nlg_comb_q', bars = T, col = 'red')

## Functional validation datasets
GERCARres <- as.data.table(read.csv2(file.path(IN_DIR, 'GERCARres.csv')))
GERCARres_gr <- dt2gr(GERCARres[chromosome == as.character(seqnames(gr.chr(Oncgr)))])
GERCARres_prolifd21_gt <- gTrack(gr.nochr(GERCARres_gr), y.field = 'logFC_prolif_21', bars = TRUE)

GERCARres_gr_sig <- dt2gr(GERCARres[FDR_prolif_21 < 0.1 & chromosome == as.character(seqnames(gr.chr(Oncgr)))])
GERCARres_prolifd21_q0.1_gt <- gTrack(gr.nochr(GERCARres_gr_sig))

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


CCND1_TSS <- fread(file.path(IN_DIR, 'CCND1_TSS.csv'))
colnames(CCND1_TSS)[4] = 'seqnames'
colnames(CCND1_TSS)[8] = 'start'
CCND1_TSS[, end := start +1]
CCND1_TSSgr <- dt2gr(CCND1_TSS)

window_oiz <- gr.nochr(GRanges("chr11:69000000-70000000") - 2e5)
window_oiz2 <- gr.nochr(GRanges("chr11:69000000-70000000") - 3e5)

pdf(file = paste0(OUT_DIR, '/20260914_Fig5a_CCND1_TF_BRCA_zo_5e5.pdf')
    ,paper = 'a4', height = 70)
  plot(c(gencode_hg38, GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt,
         dmrs_d6_BRCAgt, comb_Pvals_v4C_CCND1_BRCA_gt,
         MCF7_TFAP2C_hg38_bed_gt, MCF7_TFAP2C_hg38, MCF7_TRPS1_hg38_bed_gt, MCF7_TRPS1_hg38, MCF7_SPDEF_hg38_bed_gt, MCF7_SPDEF_hg38, MCF7_GATA3_hg38_bed_gt, MCF7_GATA3_hg38, MCF7_FOXA1_hg38_bed_gt, MCF7_FOXA1_hg38, MCF7_ESR1_hg38_bed_gt, MCF7_ESR1_hg38,
         BRCA_CCND1_LvsR_sig_LR_gt,
         BRCA_xGtad_per_amp_gt)
       , window_oiz, cex.label = 0.1)
  plot(c(gencode_hg38, GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt,
         dmrs_d6_BRCAgt, comb_Pvals_v4C_CCND1_BRCA_gt,
         MCF7_TFAP2C_hg38_bed_gt, MCF7_TFAP2C_hg38, MCF7_TRPS1_hg38_bed_gt, MCF7_TRPS1_hg38, MCF7_SPDEF_hg38_bed_gt, MCF7_SPDEF_hg38, MCF7_GATA3_hg38_bed_gt, MCF7_GATA3_hg38, MCF7_FOXA1_hg38_bed_gt, MCF7_FOXA1_hg38, MCF7_ESR1_hg38_bed_gt, MCF7_ESR1_hg38,
         BRCA_CCND1_LvsR_sig_LR_gt,
         BRCA_xGtad_per_amp_gt)
       , window_oiz2, cex.label = 0.1)
dev.off()
# #
