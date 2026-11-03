# Tests which H3K27ac and ATAC peaks in the consensus amplified regions of one gene are enriched in one tumor type,
# compared with all other samples (Fisher's exact test, q-values across all tumor types).
# Run from this folder: Rscript 01_peak_enrichment_fisher.R <gene>, with <gene> a name in gene_settings.R

# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
N_CORES <- as.integer(Sys.getenv('N_CORES', unset = '10'))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(gUtils)
library(parallel)
library(qvalue)
library(rtracklayer)

source('gene_settings.R')

xG <- commandArgs(trailingOnly = TRUE)[1]
if (!xG %in% names(gene_settings)) stop('Usage: Rscript 01_peak_enrichment_fisher.R <', paste(names(gene_settings), collapse = '|'), '>')
settings <- gene_settings[[xG]]

# Load hs1ToHg38 conversion chain
hs1ToHg38_chain_file <- file.path(IN_DIR, 'hs1ToHg38.over.chain')
if (!file.exists(hs1ToHg38_chain_file)) {
    download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz', destfile = paste0(hs1ToHg38_chain_file, '.gz'))
    system(paste('gunzip -f', paste0(hs1ToHg38_chain_file, '.gz')))
}
hs1ToHg38_chain <- import.chain(hs1ToHg38_chain_file)

# Consensus amplified regions (hs1) in hg38, the coordinates of the peaks: amplified in more than 80% of the walks
cons_amp_gr <- lapply(settings$fisher_amp_tracks, function(track_file) {
  amp_gt <- readRDS(file.path(IN_DIR, track_file))
  unlist(liftOver(amp_gt@data[[1]], hs1ToHg38_chain)) %Q% (AmpFreq > 0.8)
})

# Peaks of each tumor type in the consensus amplified regions. Column 4 is the number of samples with the peak
peaks_co_amped <- rbindlist(lapply(names(tumor_types), function(xTT) {
  peaks <- fread(file.path(IN_DIR, settings$peak_files[[xTT]]))
  peaks <- peaks[V4 >= settings$min_samples[[xTT]]]
  colnames(peaks)[1:3] <- c('seqnames', 'start', 'end')
  peaks_co_amped_tt <- gr2dt(dt2gr(peaks) %&% cons_amp_gr[[tumor_types[[xTT]]]])
  peaks_co_amped_tt[, tt := xTT]
}), fill = TRUE)
peaks_co_amped[, PeakID := paste0(seqnames, ':', start, ':', end)]

# The same intervals across all samples
all_peaks <- fread(file.path(IN_DIR, settings$all_peaks_file))
colnames(all_peaks)[1:5] <- c('seqnames', 'start', 'end', 'samples_with_peaks_sum', 'samples_with_peaks_IDn')
all_peaks_gr <- dt2gr(all_peaks)

fisher_stats <- rbindlist(lapply(unique(peaks_co_amped$tt), function(xTT) {
  print(xTT)
  peaks_tt <- peaks_co_amped[tt == xTT]
  rbindlist(mclapply(peaks_tt$PeakID, function(xPeak) {
    peak <- peaks_tt[PeakID == xPeak]
    peak_gr <- dt2gr(peak)
    # the widest interval across all samples that shares more than 5 bp with the peak
    all_peaks_ovlp_gr <- all_peaks_gr %&% peak_gr
    all_peaks_ovlp <- gr2dt(all_peaks_ovlp_gr)[(all_peaks_ovlp_gr %o% peak_gr) > 5]
    all_peaks_ovlp <- all_peaks_ovlp[which.max(width)]
    # samples with and without the peak, in this tumor type and in all others. The group sizes are the highest sample
    # counts in the tables, which equal the number of samples if some interval is shared by all of them.
    # A peak with no such interval has an empty count_others, so matrix()
    # repeats the first row: the test then compares the tumor type with itself (p = 1) and compttype_count is NA
    n_tt <- max(peaks_tt$V4)
    n_others <- max(all_peaks$samples_with_peaks_sum) - n_tt
    count_tt <- peak$V4
    count_others <- all_peaks_ovlp$samples_with_peaks_sum - count_tt
    contMat <- matrix(c(count_tt, n_tt - count_tt,
                        count_others, n_others - count_others),
                      nrow = 2, ncol = 2, byrow = TRUE,
                      dimnames = list(c(xTT, 'others'), c('with_peaks', 'no_peaks')))
    fisher_res <- fisher.test(x = contMat)
    data.table(p_val = fisher_res$p.value, OR = fisher_res$estimate, peakID = xPeak, tstdttype = xTT,
               tstdttype_count = count_tt, compttype_count = count_others)
  }, mc.cores = N_CORES))
}))
# fisher.test can return a p-value just above 1 from rounding, which qvalue rejects
fisher_stats[, p_val := pmin(p_val, 1)]
fisher_stats[, q_val := qvalue(p = p_val)$qvalue]
write.csv(fisher_stats, file.path(OUT_DIR, settings$stats_file), row.names = FALSE)
