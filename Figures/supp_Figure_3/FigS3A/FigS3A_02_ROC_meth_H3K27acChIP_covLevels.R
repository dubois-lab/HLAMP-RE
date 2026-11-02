library(gGnome)
library(gTrack)
library(rtracklayer)

gencode_hg38 = readRDS('source_data/gencode_hg38.rds')
CGC <- readRDS('source_data/CGC_hg38.rds')

gene <- CGC %Q%  (Gene_Symbol == 'CCND1')
window = gene + 10e5

## methylation data
HCC1954_bedmethyl_ungrouped <- readRDS('source_data/HCC1954.ungrouped.bedmethyl.rds')
# filter for meth. (CpG)
all_cpg <- dt2gr(HCC1954_bedmethyl_ungrouped[modified_base_code == 'm' & Nvalid_cov > 5])
all_cpg <- gr.chr(all_cpg)
rm(HCC1954_bedmethyl_ungrouped)

## coverage
cnvkit_cov <- read.table("source_data/HCC1954_Tumor_ONT.GRCh38.sorted.targetcoverage.cnn", sep="\t",header=T,stringsAsFactors=F)
cnvkit_cov_gr <- dt2gr(cnvkit_cov)
cnvkit_cov_gr <- cnvkit_cov_gr[seqnames(cnvkit_cov_gr)=="chr11"]

## H3K27ac ChIP
bw_FE <- import('source_data/HCC1954_H3K27ac_ChIP_FE.bw')
chip_FE_gt <- gTrack(data = bw_FE, y.field = 'score', name = 'FE_score',bars = T,yaxis.cex=0.7,y1=8)
rm(bw_FE)
H3K27ac_peaks <- fread('source_data/HCC1954_peaks.narrowPeak')
colnames(H3K27ac_peaks)[1:3] = c('chr', 'start', 'end')
H3K27ac_peaks_gr <- dt2gr(H3K27ac_peaks)

## define coverage ranges
amp_range_high <- cnvkit_cov_gr[cnvkit_cov_gr$depth >= 350 & cnvkit_cov_gr$depth <= 550]
amp_range_high <- reduce(amp_range_high)
seqlengths(amp_range_high) <- NA
amp_range_medium <- cnvkit_cov_gr[cnvkit_cov_gr$depth >= 125 & cnvkit_cov_gr$depth <= 250]
amp_range_medium <- reduce(amp_range_medium)
seqlengths(amp_range_medium) <- NA

cpg_range_high <- subsetByOverlaps(all_cpg,amp_range_high)
cpg_range_medium <- subsetByOverlaps(all_cpg,amp_range_medium)

## baseline truth
seqlengths(H3K27ac_peaks_gr) <- NA
H3K27ac_peaks_gr_sub <- H3K27ac_peaks_gr[H3K27ac_peaks_gr$V7 >= 5]
chip_cpg_range_high <- subsetByOverlaps(cpg_range_high,H3K27ac_peaks_gr_sub)
chip_cpg_range_high0 <- chip_cpg_range_high[, FALSE]
chip_cpg_range_medium <- subsetByOverlaps(cpg_range_medium,H3K27ac_peaks_gr_sub)
chip_cpg_range_medium0 <- chip_cpg_range_medium[, FALSE]

cov_chip_range_high <- coverage(chip_cpg_range_high0)
cov_chip_range_medium <- coverage(chip_cpg_range_medium0)

## ROC analysis
for (typ in c("medium","high")) {

cov_chip <- get(paste0("cov_chip_range_",typ))
cpg_amp <- get(paste0("cpg_range_",typ))
sens_v <- c()
spec_v <- c()

for (k in 0:100) {
  meth_amp <- cpg_amp[mcols(cpg_amp)$percent_modified <= k]
  meth_amp0 <- meth_amp[, FALSE]
  cov_meth <- coverage(meth_amp0)
  
  TP <- FP <- FN <- TN <- 0
  for (i in seq_along(cpg_amp)) {
    r <- cpg_amp[i]
    chr <- as.character(seqnames(r))
    start_r <- start(r)
    end_r <- end(r)
    meth_vec <- as.vector(cov_meth[[chr]][start_r:end_r] > 0)
    chip_vec <- as.vector(cov_chip[[chr]][start_r:end_r] > 0)
    TP <- TP + sum(chip_vec & meth_vec)
    FP <- FP + sum(chip_vec & !meth_vec)
    FN <- FN + sum(!chip_vec & meth_vec)
    TN <- TN + sum(!chip_vec & !meth_vec)
  }
  conf_mat <- matrix(
    c(TP, FP, FN, TN),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      "ChIP" = c("TRUE", "FALSE"),
      "de-meth" = c("TRUE", "FALSE"))
    )
    TP <- conf_mat["TRUE","TRUE"]
    FP <- conf_mat["FALSE","TRUE"]
    FN <- conf_mat["TRUE","FALSE"]
    TN <- conf_mat["FALSE","FALSE"]
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    precision   <- TP / (TP + FP)

    sens_v <- c(sens_v,sensitivity)
    spec_v <- c(spec_v,specificity)
  }
saveRDS(sens_v,paste0('source_data/HCC1954_ChIPvsMeth_AUC_sensitivity_covRange_',typ,'.rds'))
saveRDS(spec_v,paste0('source_data/HCC1954_ChIPvsMeth_AUC_specificity_covRange_',typ,'.rds'))
}
