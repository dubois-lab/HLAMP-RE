# CCND1 TCGA-C8-A12M HiChIP (BRCA) locus visualisation
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Load packages tolerantly. A transitive dependency (biomaRt, pulled in via gGnome/khtools)
## runs an Ensembl connectivity check at load time that can fail. biomaRt is not used, so suppress the
## load-time messages/warnings and ignore any such load-time error.
suppressWarnings(suppressMessages({
  try(library(parallel), silent = TRUE)
  try(library(GxG), silent = TRUE)
  try(library(rtracklayer), silent = TRUE)
}))
CGC <- readRDS('../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
gencode_hg38 = readRDS('../supp_Figure_3/FigS3A/source_data/gencode_hg38.rds')

# hs1 -> hg38 liftOver chain (CCND1 amp / lindiff / methylation / v4C tracks are in hs1/T2T coordinates)
hs1ToHg38_chain <- import.chain(file.path(IN_DIR, 'hs1ToHg38.over.chain'))

xG = 'CCND1'
xT = 'TCGA_C8_A12M'
Oncgr <- CGC %Q%  (Gene_Symbol == xG)
window_oi = Oncgr + 12e5

# CCND1 BRCA amplification frequency + left-vs-right and BRCA-vs-NSCLC significance tracks (hs1 -> hg38)
BRCA_xGtad_per_amp_gt <- readRDS(paste0(IN_DIR, "/20260907_Fig.2C_BRCA_CCND1_AmpFreq_track.rds"))
BRCA_xGtad_per_amp_gt_gr <- unlist(liftOver(BRCA_xGtad_per_amp_gt@data[[1]], hs1ToHg38_chain))
BRCA_xGtad_per_amp_gt <- gTrack(data = BRCA_xGtad_per_amp_gt_gr, y.field = "AmpFreq", name = "%amp", col = "red")
BRCA_xGtad_per_amp_gt$bars <- F

BRCA_CCND1_LvsR_sig_LR_gt <- readRDS(paste0(IN_DIR, "/20260907c_Fig.2C_BRCA_CCND1_LvsR_sig_LR_cons80.rds"))
BRCA_CCND1_LvsR_sig_LR_gt_gr <- unlist(liftOver(BRCA_CCND1_LvsR_sig_LR_gt@data[[1]], hs1ToHg38_chain))
BRCA_CCND1_LvsR_sig_LR_gt <- gTrack(data = BRCA_CCND1_LvsR_sig_LR_gt_gr, name = "nlgp", col = "red", height = 0.5)


TCGA_C8_A12M_HiChiP = straw(
  file.path(IN_DIR, 'BRCA-C9C8D426-A3FD-4455-89A9-768BC01D66A9-X003-S02-B1-T1_H3K27ac.allValidPairs.hic'), 
  gr = window_oi, res = 5e3)
TCGA_C8_A12M_seg_dt <- fread(file.path(IN_DIR, '6b20f4a5-91c3-482f-a923-6be878ec942f.wgs.ASCAT.copy_number_variation.seg.txt'))
TCGA_C8_A12M_seg_gt <- gTrack(data = dt2gr(TCGA_C8_A12M_seg_dt), y.field = 'Copy_Number')
TCGA_HiChIP_metadata <- fread(file.path(IN_DIR, 'TCGA_HiChIP_metadata.txt'))

TCGA_C8_A12M_peaks_dt <- fread(paste0(IN_DIR, '/TCGA_HiChIP_1Dsignal_peaks/', 
                                      TCGA_HiChIP_metadata[case_submitter_id == 'TCGA-C8-A12M' & data_format == 'narrowPeak', file_name]))
colnames(TCGA_C8_A12M_peaks_dt)[1:3] <- c('seqnames',  'start',    'end')
TCGA_C8_A12M_peaks_gt <- gTrack(data = dt2gr(TCGA_C8_A12M_peaks_dt))
TCGA_C8_A12M_loops_dt <- fread(paste0(IN_DIR, '/TCGA_HiChIP_FitHiChIP_loop_calls/', 
                                      TCGA_HiChIP_metadata[case_submitter_id == 'TCGA-C8-A12M' & data_format == 'bed', file_name]))
colnames(TCGA_C8_A12M_loops_dt)[1:3] = c('seqnames', 'start', 'end')
TCGA_C8_A12M_loops_dt_b1gr <- dt2gr(TCGA_C8_A12M_loops_dt)
TCGA_C8_A12M_loops_dt_b1grOnc <- gr.nochr(TCGA_C8_A12M_loops_dt_b1gr) %&% Oncgr
# bin2
colnames(TCGA_C8_A12M_loops_dt)[1:3] = c('seqnames_b1', 'start_b1', 'end_b1')
colnames(TCGA_C8_A12M_loops_dt)[4:6] = c('seqnames', 'start', 'end')
TCGA_C8_A12M_loops_dt_b2gr <- dt2gr(TCGA_C8_A12M_loops_dt)
TCGA_C8_A12M_loops_dt_b2grOnc <- gr.nochr(TCGA_C8_A12M_loops_dt_b2gr) %&% Oncgr
TCGA_C8_A12M_loops_dt_Onc_dt <- rbind(gr2dt(TCGA_C8_A12M_loops_dt_b2grOnc), gr2dt(TCGA_C8_A12M_loops_dt_b1grOnc), use.names=F)
colnames(TCGA_C8_A12M_loops_dt_Onc_dt)[1:3] = c('seqnames_b3', 'start_b3', 'end_b3')
colnames(TCGA_C8_A12M_loops_dt_Onc_dt)[6:8] = c('seqnames', 'start', 'end')

TCGA_C8_A12M_loops_dt_OncLinkedPeaks_gr <- dt2gr(TCGA_C8_A12M_loops_dt_Onc_dt[Q.Value_Bias < 0.001])
TCGA_C8_A12M_loops_dt_OncLinkedPeaks_gt <- gTrack(reduce(TCGA_C8_A12M_loops_dt_OncLinkedPeaks_gr), name = 'TCGA_C8_A12M_loops')

H3K27ac_1Dsignal.norm <- gTrack(data = file.path(IN_DIR, 'BRCA-C9C8D426-A3FD-4455-89A9-768BC01D66A9-X003-S02-B1-T1_H3K27ac_1D-signal.norm.bw')
                                ,bars = T )

# make 4C track from the HiChIP matrix at the CCND1 TSS
# (hg38; 'chromosome_name' column is non-'chr' -> matches the straw HiChIP gr directly)
CCND1_TSS <- fread(file.path(IN_DIR, 'CCND1_TSS_hg38.csv'))
colnames(CCND1_TSS)[4] = 'seqnames'
colnames(CCND1_TSS)[8] = 'start'
CCND1_TSS[, end := start + 1]
CCND1_TSS[, strand := '*']
CCND1_TSSgr_hg38 <- dt2gr(CCND1_TSS)
TCGA_xT_HiChiP_4C <- TCGA_C8_A12M_HiChiP[TCGA_C8_A12M_HiChiP$gr %^% CCND1_TSSgr_hg38, drop = T]
TCGA_xT_HiChiP_4C_gt <- gTrack(data = TCGA_xT_HiChiP_4C, y.field = 'value', bars = T, col = 'red', y1 = 200)

# our tumor methylation data: DMRs (hs1 -> hg38)
xD = 0.4
dmrs_d4 <- readRDS(file.path(IN_DIR, '20260512_CCND1_TM_dmrs_d.rds'))
colnames(dmrs_d4)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d4_gr <- unlist(liftOver(gr.chr(dt2gr(dmrs_d4)), hs1ToHg38_chain))
dmrs_d4 <- as.data.table(dmrs_d4_gr)
dmrs_d4_BRCAgt <- gTrack(dt2gr(dmrs_d4[diff.Methy > xD]))
dmrs_d5_BRCAgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy > 0.5])), ygap =0.5)
dmrs_d6_BRCAgt <- gTrack(reduce(dt2gr(dmrs_d4[diff.Methy > 0.6])), ygap =0.5)


#functional val####
GERCARres <- as.data.table(read.csv2(file.path(IN_DIR, 'GERCARres.csv')))
GERCARres_gr <- dt2gr(GERCARres[chromosome == as.character(seqnames(gr.chr(Oncgr)))])
GERCARres_prolifd21_gt <- gTrack(gr.nochr(GERCARres_gr), y.field = 'logFC_prolif_21', bars = T)
GERCARres_gr_sig <- dt2gr(GERCARres[FDR_prolif_21 < 0.1 & chromosome == as.character(seqnames(gr.chr(Oncgr)))])
GERCARres_prolifd21_q0.1_gt <- gTrack(gr.nochr(GERCARres_gr_sig))

# significant lineage-differential H3K27ac & ATAC peaks (BRCA_HR+)
ATACfisher <- fread(file.path(IN_DIR, '20260614_CCND1_H3K27ac&ATAC_linDiffs_stats.csv'))
ATACfisher[, chr := unlist(strsplit(x = peakID, split = "[:]"))[1]]
ATACfisher[, start := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[2]), by = peakID]
ATACfisher[, end := as.numeric(unlist(strsplit(x = peakID, split = "[:]"))[3]), by = peakID]
fishres_A2gr <- dt2gr(ATACfisher)
fishres_A_sigs_BRCA_HR <- fishres_A2gr %Q% (q_val < 0.05 & tstdttype == "BRCA_HR+")
fishres_sigs_BRCA_HR_gt <- gTrack(data = reduce(fishres_A_sigs_BRCA_HR))

# CCND1 BRCA virtual-4C combined (Fisher) p-values (hs1 -> hg38)
comb_Pvals_v4C_CCND1_BRCA_gr <- readRDS(file.path(IN_DIR, '20260702_comb_Pvals_v4C_CCND1_BRCA_gr.rds'))
comb_Pvals_v4C_CCND1_BRCA_gr_hg38 <- grl.unlist(liftOver(comb_Pvals_v4C_CCND1_BRCA_gr, hs1ToHg38_chain))
comb_Pvals_v4C_CCND1_BRCA_gt <- gTrack(comb_Pvals_v4C_CCND1_BRCA_gr_hg38, y.field = 'nlg_comb_q', bars = T, y1 =30, col = 'red', sep.draw =F)

# plot window (window_oiz)
# stripped to non-'chr' seqnames to match the straw HiChIP matrix
window_oi <- gr.nochr(GRanges("chr11:69000000-70000000"))
window_oiz <- gr.nochr(GRanges("chr11:69000000-70000000") - 2e5)

pdf(file = paste0(OUT_DIR, '/20260908_', xT, '_HiChIP_BRCA_CCND1.pdf')
    ,paper = 'a4', height = 70)
  plot(c(gencode_hg38,
       GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt,
       #  fishres_sigs_BRCA_HR_gt,
       dmrs_d4_BRCAgt, dmrs_d5_BRCAgt, dmrs_d6_BRCAgt,
       BRCA_CCND1_LvsR_sig_LR_gt, BRCA_xGtad_per_amp_gt, 
       comb_Pvals_v4C_CCND1_BRCA_gt,
       TCGA_C8_A12M_seg_gt,H3K27ac_1Dsignal.norm,
       TCGA_xT_HiChiP_4C_gt,
       TCGA_C8_A12M_HiChiP$gtrack(name = paste0(xT, "_HiChiP"), height = 90, clim = c(0, 50))),
     window_oi, cex.label = 0.1)
  plot(c(gencode_hg38,
         GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt,
         #  fishres_sigs_BRCA_HR_gt,
         dmrs_d4_BRCAgt, dmrs_d5_BRCAgt, dmrs_d6_BRCAgt,
         BRCA_CCND1_LvsR_sig_LR_gt, BRCA_xGtad_per_amp_gt, 
         comb_Pvals_v4C_CCND1_BRCA_gt,
         TCGA_C8_A12M_seg_gt,H3K27ac_1Dsignal.norm, TCGA_C8_A12M_loops_dt_OncLinkedPeaks_gt,
         TCGA_xT_HiChiP_4C_gt,
         TCGA_C8_A12M_HiChiP$gtrack(name = paste0(xT, "_HiChiP"), height = 90, clim = c(0, 50))),
       window_oi, cex.label = 0.1)
  plot(c(gencode_hg38,
         GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt,
        #  fishres_sigs_BRCA_HR_gt,
         dmrs_d4_BRCAgt, dmrs_d5_BRCAgt, dmrs_d6_BRCAgt,
         BRCA_CCND1_LvsR_sig_LR_gt, BRCA_xGtad_per_amp_gt, 
         comb_Pvals_v4C_CCND1_BRCA_gt,
        TCGA_C8_A12M_seg_gt,H3K27ac_1Dsignal.norm,
         TCGA_xT_HiChiP_4C_gt,
         TCGA_C8_A12M_HiChiP$gtrack(name = paste0(xT, "_HiChiP"), height = 90, clim = c(0, 50))),
       window_oiz, cex.label = 0.1)
dev.off()
