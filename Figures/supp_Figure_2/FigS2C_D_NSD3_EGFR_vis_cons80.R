# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(gGnome)
library(gTrack)

xG <- "EGFR" # NSD3, EGFR
chrom <- "chr7" # chr8,7

# gene track
CGC <- readRDS('../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')
gene <- CGC %Q% (Gene_Symbol == xG)
plot_window <- gene + 2e6
gt_ge <- readRDS(paste0('../../common/data/refGenome/',xG,'_plotting_genetrack_T2T.rds'))

ccnd1_l <- readRDS(paste0(IN_DIR, '/20260825_NSCLC_LUSC_',xG,'_CNMatm_combined_hs1.rds'))
ccnd1_b <- readRDS(paste0(IN_DIR, '/20260825_BRCA_',xG,'_CNMatm_combined_hs1.rds'))

ccnd1_b_gr <- dt2gr(ccnd1_b)
ccnd1_l_gr <- dt2gr(ccnd1_l)
ccnd1_b_gt <- gTrack(ccnd1_b_gr, y.field="AmpPercent", col="brown3",bars=T)
ccnd1_l_gt <- gTrack(ccnd1_l_gr, y.field="AmpPercent", col="royalblue3",bars=T)

fisher_bVl <- readRDS(paste0(IN_DIR, '/20260825_BRCAvsNSCLC_',xG,'_',chrom,'_FisherStats_CN10_cons80_gr.rds'))
sign_BvL <- fisher_bVl[fisher_bVl$Fisher_pval<0.05]
#sign_BvL <- fisher_bVl[fisher_bVl$qval<0.05]
mcols(sign_BvL)$tissue <- ifelse(mcols(sign_BvL)$BRCA_AmpFreq > mcols(sign_BvL)$NSCLC_AmpFreq, "BRCA", "NSCLC")
colormap_BvL <- list(tissue=c(NSCLC="royalblue3", BRCA="brown3"))
sign_BvL_gt <- gTrack(sign_BvL,gr.colorfield = "color",colormaps = colormap_BvL, height=5)

fisher_LvR_l <- fread(paste0(IN_DIR, '/20260825_NSCLC_',xG,'_LeftVsRight_FisherStats_CN10_cons80.txt'))
fisher_LvR_b <- fread(paste0(IN_DIR, '/20260825_BRCA_',xG,'_LeftVsRight_FisherStats_CN10_cons80.txt'))

fisher_LvR_l_gr <- dt2gr(fisher_LvR_l[, ':='(seqnames = fifelse(right_AmpFreq > left_AmpFreq, right_chr, left_chr),
                                             start = fifelse(right_AmpFreq > left_AmpFreq, right_start, left_start),
                                             end = fifelse(right_AmpFreq > left_AmpFreq, right_end, left_end)
)])
fisher_LvR_b_gr <- dt2gr(fisher_LvR_b[, ':='(seqnames = fifelse(right_AmpFreq > left_AmpFreq, right_chr, left_chr),
                                             start = fifelse(right_AmpFreq > left_AmpFreq, right_start, left_start),
                                             end = fifelse(right_AmpFreq > left_AmpFreq, right_end, left_end)
)])

sign_LvR_l <- fisher_LvR_l_gr[fisher_LvR_l_gr$Fisher_pval<0.05]
sign_LvR_b <- fisher_LvR_b_gr[fisher_LvR_b_gr$Fisher_pval<0.05]
#sign_LvR_l <- fisher_LvR_l_gr[fisher_LvR_l_gr$qval<0.05]
#sign_LvR_b <- fisher_LvR_b_gr[fisher_LvR_b_gr$qval<0.05]
sign_LvR_l_gt <- gTrack(sign_LvR_l, col="royalblue3", height=5)
sign_LvR_b_gt <- gTrack(sign_LvR_b, col="brown3", height=5)


pdf(file = paste0(OUT_DIR, '/FigS2C_D_',xG,'_vis_cons80.pdf'), paper = "a4", height = 24)
plot(c(gt_ge,
       ccnd1_l_gt,
       sign_LvR_l_gt,
       sign_BvL_gt,
       sign_LvR_b_gt,
       ccnd1_b_gt),plot_window)
dev.off()
