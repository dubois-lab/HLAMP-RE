# Virtual 4C tracks of one sample: the row of the HiCDCPlus contact matrix at the bin of each gene's TSS.
# Run from this folder: Rscript 06_virtual_4C_tracks.R <sample> <gene> [<gene> ...]
# Writes one track per gene and value: <sample>_covCovariate_v4C_<gene>_<value>.rds, with value nlg_p, nlg_q, counts or normcounts.

# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'results')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'results')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(gGnome)
library(GxG)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('Usage: Rscript 06_virtual_4C_tracks.R <sample> <gene> [<gene> ...]')
sid <- args[1]
genes <- args[-1]

binsize <- 5e3

# Bins of the contact matrix, 0-based
bins <- fread(file.path(IN_DIR, 'normalized_abs.bed'), select = 1:3, col.names = c('chr', 'start', 'end'))
bins[, start := start + 1]

interactions_dt <- fread(file.path(IN_DIR, paste0(sid, '_covCovariate_interactions_dt.csv')))
interactions_dt[, nlg_p := -log10(pvalue)]
interactions_dt[, nlg_q := -log10(qvalue)]
interactions_dt[, i := end1 / binsize]
interactions_dt[, j := end2 / binsize]

for (xG in genes) {
  TSS_gr <- readRDS(file.path(IN_DIR, paste0(xG, '_TSSgr_T2T.rds')))
  xChr <- as.character(seqnames(TSS_gr))[1]
  for (datafield in c('nlg_p', 'nlg_q', 'counts', 'normcounts')) {
    interactions_chr <- interactions_dt[seqnames1 == xChr & seqnames2 == xChr, .(i, j, value = get(datafield))]
    gmat <- gM(gr = dt2gr(bins[chr == xChr]), dat = interactions_chr, fill = 0)
    v4C_gr <- gmat[gmat$gr %^% TSS_gr, drop = TRUE]
    saveRDS(v4C_gr, file.path(OUT_DIR, paste0(sid, '_covCovariate_v4C_', xG, '_', datafield, '.rds')))
  }
}
