## =============================================================================
## Walks co-amplification profiles per cohort
## -----------------------------------------------------------------------------
## Four pages per cohort, three tracks each (%amp, AmpCount, MedianCN; rolling
## mean over SMOOTH_K bins within a chromosome):
##   1  gene +/- 5 Mb   2  gene +/- 2.5 Mb   3  whole gene chromosome   4  all chromosomes
##
## IN   <OUT_DIR>/<GENE>/{BRCA,NSCLC_LUSC}_<GENE>_CNMatm_combined_hs1.txt.gz
## OUT  <OUT_DIR>/<GENE>/{BRCA,NSCLC_LUSC}_<GENE>_CNMatm_combined_hs1_coamps.pdf
##
## Usage (env): GENE=CCND1|NSD3|EGFR [OUT_DIR=output] [REF_DIR=...] [SMOOTH_K=5]
##              Rscript 03_plot_coamps.R
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(gGnome)
  library(gUtils)
  library(gTrack)
})
options(scipen = 999)

xG       <- toupper(Sys.getenv('GENE', unset = 'CCND1'))
stopifnot(xG %in% c('CCND1', 'NSD3', 'EGFR'))
SMOOTH_K <- as.integer(Sys.getenv('SMOOTH_K', unset = '5'))
REF_DIR  <- Sys.getenv('REF_DIR', unset = '../../../data/refGenome')
OUT_DIR  <- Sys.getenv('OUT_DIR', unset = 'output')
in_dir   <- file.path(OUT_DIR, xG)

gt_ge  <- readRDS(file.path(REF_DIR, 'T2T_gencode_gt_UCSC.rds'))
CGC    <- readRDS(file.path(REF_DIR, 'cgc_T2T_gr_consAmpWalks.rds'))
window <- CGC[CGC$Gene_Symbol == xG][1]
intchr <- as.character(seqnames(window))
hs1_sizes <- fread(file.path(REF_DIR, 'hs1.chrom.sizes'),
                   header = FALSE, col.names = c('seqnames', 'seqlength'))
hs1_sizes <- hs1_sizes[seqnames %in% paste0('chr', c(1:22, 'X', 'Y', 'M'))]
hs1_seqlengths <- setNames(hs1_sizes$seqlength, hs1_sizes$seqnames)

roll_mean <- function(x, k = SMOOTH_K) {
  if (length(x) < k) return(x)
  y <- as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
  ix <- is.na(y); y[ix] <- x[ix]; y
}
keep_cols <- c('binID', 'seqnames', 'start', 'end', 'n_samples_total', 'AmpPercent', 'AmpSum', 'MedianCN')
for (tt in c('BRCA', 'NSCLC_LUSC')) {
  in_file <- file.path(in_dir, sprintf('%s_%s_CNMatm_combined_hs1.txt.gz', tt, xG))
  out_pdf <- file.path(in_dir, sprintf('%s_%s_CNMatm_combined_hs1_coamps.pdf', tt, xG))
  d <- fread(in_file, select = keep_cols)
  d[is.na(MedianCN), MedianCN := 0]
  n_tot <- max(d$n_samples_total)
  setorder(d, seqnames, start)
  d[, AmpFraction_smoothed  := roll_mean(AmpPercent),           by = seqnames]
  d[, mean_support_smoothed := roll_mean(as.numeric(AmpSum)),   by = seqnames]
  d[, median_CN_smoothed    := roll_mean(as.numeric(MedianCN)), by = seqnames]
  gr <- dt2gr(d[, .(seqnames, start, end, AmpFraction_smoothed, mean_support_smoothed, median_CN_smoothed)])
  seqlevels(gr) <- intersect(seqlevels(gr), names(hs1_seqlengths))
  seqlengths(gr) <- hs1_seqlengths[seqlevels(gr)]
  gt <- c(gTrack(gr, y.field = 'AmpFraction_smoothed',  name = '%amp',     col = 'blue',    bars = TRUE, y0 = 0, y1 = 1),
          gTrack(gr, y.field = 'mean_support_smoothed', name = 'AmpCount', col = 'blue',    bars = TRUE, y0 = 0, y1 = max(1, n_tot)),
          gTrack(gr, y.field = 'median_CN_smoothed',    name = 'MedianCN', col = 'darkred', bars = TRUE, y0 = 0))
  pdf(out_pdf, paper = 'a4', height = 24)
  try(plot(c(gt_ge, gt), window + 50e5, cex.label = 0.1, sep.draw = FALSE, y0 = 0))
  try(plot(c(gt_ge, gt), window + 25e5, cex.label = 0.1, sep.draw = FALSE, y0 = 0))
  try(plot(gt, tileGenome(seqlengths = hs1_seqlengths[intchr], tilewidth = 1e6, cut.last.tile.in.chrom = TRUE),
           cex.label = 0.3, sep.draw = FALSE, y0 = 0))
  try(plot(gt, tileGenome(seqlengths = hs1_seqlengths[paste0('chr', c(1:22, 'X', 'Y'))], tilewidth = 1e6,
                          cut.last.tile.in.chrom = TRUE), cex.label = 0.3, sep.draw = FALSE, y0 = 0))
  invisible(dev.off())
  cat(sprintf('[%s] %d bins drawn -> %s\n', tt, nrow(d), basename(out_pdf)))
}
cat('Done.\n')
