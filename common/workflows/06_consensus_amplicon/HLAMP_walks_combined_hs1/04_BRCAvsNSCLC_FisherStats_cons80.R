## =============================================================================
## BRCA vs NSCLC per 10 kb bin on the gene chromosome -- cons80 Fisher tests
## -----------------------------------------------------------------------------
## 2x2 Fisher test of (amplified, not amplified) x (BRCA, NSCLC) per bin, only
## where at least one cohort reaches AmpFreq >= 80% (cons80).  A sample counts
## as amplified in a bin at CN >= 10 (AmpSum of the combined matrix); the
## denominator is the cohort size (n_samples_total).  q-values by qvalue() over
## the tested bins, BH fallback.
##
## IN   <OUT_DIR>/<GENE>/{BRCA,NSCLC_LUSC}_<GENE>_CNMatm_combined_hs1.txt.gz
## OUT  <OUT_DIR>/<GENE>/stats/cons80/<GENE>_BRCAvsNSCLC_*
##
## Usage (env): GENE=CCND1|NSD3|EGFR [OUT_DIR=output] [REF_DIR=...]
##              [AMP_FREQ_MIN=0.80] [CN_CUTOFF=10, used in file names only]
##              Rscript 04_BRCAvsNSCLC_FisherStats_cons80.R
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(qvalue)
  library(gTrack)
  library(gUtils)
})
options(scipen = 999)

xG           <- toupper(Sys.getenv('GENE', unset = 'CCND1'))
stopifnot(xG %in% c('CCND1', 'NSD3', 'EGFR'))
CNcutoff     <- as.numeric(Sys.getenv('CN_CUTOFF', unset = '10'))
AMP_FREQ_MIN <- as.numeric(Sys.getenv('AMP_FREQ_MIN', unset = '0.80'))
REF_DIR      <- Sys.getenv('REF_DIR', unset = '../../../data/refGenome')
OUT_DIR      <- Sys.getenv('OUT_DIR', unset = 'output')
in_dir     <- file.path(OUT_DIR, xG)
brca_path  <- file.path(in_dir, sprintf('BRCA_%s_CNMatm_combined_hs1.txt.gz', xG))
nsclc_path <- file.path(in_dir, sprintf('NSCLC_LUSC_%s_CNMatm_combined_hs1.txt.gz', xG))
out_dir    <- file.path(in_dir, 'stats', 'cons80')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(brca_path), file.exists(nsclc_path))

gt_ge   <- readRDS(file.path(REF_DIR, 'T2T_gencode_gt_UCSC.rds'))
CGC     <- readRDS(file.path(REF_DIR, 'cgc_T2T_gr_consAmpWalks.rds'))
gene_gr <- CGC[CGC$Gene_Symbol == xG][1]
stopifnot(length(gene_gr) == 1)
hs1_sizes <- fread(file.path(REF_DIR, 'hs1.chrom.sizes'),
                   header = FALSE, col.names = c('seqnames', 'seqlength'))
hs1_sizes <- hs1_sizes[seqnames %in% paste0('chr', c(1:22, 'X', 'Y', 'M'))]
hs1_seqlengths <- setNames(hs1_sizes$seqlength, hs1_sizes$seqnames)
intchr <- as.character(seqnames(gene_gr))
cat(sprintf('%s on %s:%d-%d\n', xG, intchr, start(gene_gr), end(gene_gr)))

keep_cols <- c('binID', 'seqnames', 'start', 'end', 'n_samples_total', 'AmpSum')
brca  <- fread(brca_path,  select = keep_cols)
nsclc <- fread(nsclc_path, select = keep_cols)
brca_chr  <- brca[seqnames == intchr]
nsclc_chr <- nsclc[seqnames == intchr]
cat(sprintf('%s bins: BRCA=%d, NSCLC=%d\n', intchr, nrow(brca_chr), nrow(nsclc_chr)))
setnames(brca_chr,  c('n_samples_total', 'AmpSum'), c('BRCA_n_total',  'BRCA_AmpSum'))
setnames(nsclc_chr, c('n_samples_total', 'AmpSum'), c('NSCLC_n_total', 'NSCLC_AmpSum'))

merged <- merge(brca_chr[,  .(binID, seqnames, start, end, BRCA_n_total,  BRCA_AmpSum)],
                nsclc_chr[, .(binID,                       NSCLC_n_total, NSCLC_AmpSum)],
                by = 'binID', all = FALSE)
setorder(merged, seqnames, start)
cat(sprintf('%s bins shared between BRCA and NSCLC: %d\n', intchr, nrow(merged)))

merged[, BRCA_Samples_amped      := as.integer(BRCA_AmpSum)]
merged[, NSCLC_Samples_amped     := as.integer(NSCLC_AmpSum)]
merged[, BRCA_Samples_non_amped  := pmax(0L, as.integer(BRCA_n_total)  - BRCA_Samples_amped)]
merged[, NSCLC_Samples_non_amped := pmax(0L, as.integer(NSCLC_n_total) - NSCLC_Samples_amped)]
merged[, frac_amp_BRCA  := ifelse(BRCA_n_total  > 0, BRCA_Samples_amped  / BRCA_n_total,  NA_real_)]
merged[, frac_amp_NSCLC := ifelse(NSCLC_n_total > 0, NSCLC_Samples_amped / NSCLC_n_total, NA_real_)]

testable <- !is.na(merged$frac_amp_BRCA) & !is.na(merged$frac_amp_NSCLC) &
            (merged$frac_amp_BRCA >= AMP_FREQ_MIN | merged$frac_amp_NSCLC >= AMP_FREQ_MIN)
cat(sprintf('cons80 (>=%.0f%% in at least one cohort): %d / %d bins\n',
            AMP_FREQ_MIN * 100, sum(testable), length(testable)))
merged_test <- merged[testable]

cat("running Fisher's exact tests...\n")
merged_test[, Fisher_pval := vapply(seq_len(.N), function(i) {
  m <- matrix(c(BRCA_Samples_amped[i],  BRCA_Samples_non_amped[i],
                NSCLC_Samples_amped[i], NSCLC_Samples_non_amped[i]), nrow = 2, byrow = TRUE)
  if (any(rowSums(m) == 0)) NA_real_ else fisher.test(m)$p.value
}, numeric(1))]
ok <- !is.na(merged_test$Fisher_pval)
qv <- rep(NA_real_, nrow(merged_test))
qres <- try(qvalue(merged_test$Fisher_pval[ok])$qvalues, silent = TRUE)
if (inherits(qres, 'try-error')) {
  cat('qvalue() failed (likely flat p-value distribution); falling back to BH.\n')
  qv[ok] <- p.adjust(merged_test$Fisher_pval[ok], method = 'BH')
} else qv[ok] <- qres
merged_test[, qval := qv]

out_chr <- merged_test[, .(chr = seqnames, start, end,
                           BRCA_Samples_amped, NSCLC_Samples_amped,
                           BRCA_Samples_non_amped, NSCLC_samples_non_amped = NSCLC_Samples_non_amped,
                           frac_amp_BRCA, frac_amp_NSCLC, Fisher_pval, qval)]
chr_path <- file.path(out_dir, sprintf('%s_BRCAvsNSCLC_FisherStats_CN%s_cons80.txt', xG, CNcutoff))
fwrite(out_chr, chr_path, sep = '\t', quote = FALSE)
cat(sprintf('Wrote %s (%d rows)\n', basename(chr_path), nrow(out_chr)))

gene_window <- gene_gr + 25e5
in_win <- overlapsAny(GRanges(merged_test$seqnames, IRanges(merged_test$start, merged_test$end)), gene_window)
fwrite(out_chr[in_win], file.path(out_dir, sprintf('%s_BRCAvsNSCLC_25e5window_FisherStats_CN%s_cons80.txt', xG, CNcutoff)), sep = '\t', quote = FALSE)

gr_full_dt <- merged_test[, .(seqnames, start, end,
                              BRCA_n_total, BRCA_Samples_amped, BRCA_AmpFreq = frac_amp_BRCA,
                              NSCLC_n_total, NSCLC_Samples_amped, NSCLC_AmpFreq = frac_amp_NSCLC,
                              BRCA_Samples_non_amped, NSCLC_Samples_non_amped, Fisher_pval, qval)]
for (cl in c('BRCA_AmpFreq', 'NSCLC_AmpFreq', 'BRCA_Samples_amped', 'NSCLC_Samples_amped'))
  gr_full_dt[is.na(get(cl)), (cl) := 0]
gr_full_gr <- dt2gr(gr_full_dt)
seqlevels(gr_full_gr)  <- intersect(seqlevels(gr_full_gr), names(hs1_seqlengths))
seqlengths(gr_full_gr) <- hs1_seqlengths[seqlevels(gr_full_gr)]
gr_full_gr$neglog10_Fisher_pval <- ifelse(is.na(gr_full_gr$Fisher_pval), 0, -log10(gr_full_gr$Fisher_pval))
saveRDS(gr_full_gr, file.path(out_dir, sprintf('%s_BRCAvsNSCLC_FisherStats_CN%s_cons80_gr.rds', xG, CNcutoff)))

sig_gr <- gr_full_gr[!is.na(gr_full_gr$Fisher_pval) & gr_full_gr$Fisher_pval < 0.05]
cat(sprintf('Significant cons80 bins (p<0.05): %d\n', length(sig_gr)))

chr_dt <- copy(merged)[, `:=`(BRCA_AmpFreq = frac_amp_BRCA, NSCLC_AmpFreq = frac_amp_NSCLC)]
for (cl in c('BRCA_AmpFreq', 'NSCLC_AmpFreq', 'BRCA_Samples_amped', 'NSCLC_Samples_amped'))
  chr_dt[is.na(get(cl)), (cl) := 0]
chr_gr <- dt2gr(chr_dt[, .(seqnames, start, end, BRCA_AmpFreq, NSCLC_AmpFreq,
                           BRCA_Samples_amped, NSCLC_Samples_amped)])
seqlevels(chr_gr)  <- intersect(seqlevels(chr_gr), names(hs1_seqlengths))
seqlengths(chr_gr) <- hs1_seqlengths[seqlevels(chr_gr)]

BRCA_freq_gt  <- gTrack(chr_gr, y.field = 'BRCA_AmpFreq',        bars = TRUE, col = 'red',  name = 'BRCA AmpFreq',  y0 = 0, y1 = 1)
NSCLC_freq_gt <- gTrack(chr_gr, y.field = 'NSCLC_AmpFreq',       bars = TRUE, col = 'blue', name = 'NSCLC AmpFreq', y0 = 0, y1 = 1)
BRCA_cnt_gt   <- gTrack(chr_gr, y.field = 'BRCA_Samples_amped',  bars = TRUE, col = 'red',  name = 'BRCA AmpCount', y0 = 0)
NSCLC_cnt_gt  <- gTrack(chr_gr, y.field = 'NSCLC_Samples_amped', bars = TRUE, col = 'blue', name = 'NSCLC AmpCount', y0 = 0)
sig_track <- if (length(sig_gr)) gTrack(sig_gr, y.field = 'neglog10_Fisher_pval', bars = TRUE, col = 'black',
                                        name = '-log10(p) p<0.05 cons80', y0 = 1, height = 5) else NULL
gt_concat <- function(...) Reduce(c, Filter(Negate(is.null), list(...)))

pdf(file.path(out_dir, sprintf('%s_BRCAvsNSCLC_FisherStats_CN%s_cons80_p0.05.pdf', xG, CNcutoff)), paper = 'a4', height = 12)
if (length(sig_gr)) {
  try(plot(gt_concat(gt_ge, sig_track, NSCLC_freq_gt, BRCA_freq_gt), sig_gr + 5e5, cex.label = 0.1, sep.draw = FALSE, y0 = 0))
  try(plot(gt_concat(gt_ge, sig_track, NSCLC_cnt_gt,  BRCA_cnt_gt),  sig_gr + 5e5, cex.label = 0.1, sep.draw = FALSE, y0 = 0))
}
chr_tiles <- tileGenome(seqlengths = hs1_seqlengths[intchr], tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
try(plot(gt_concat(sig_track, NSCLC_freq_gt, BRCA_freq_gt), chr_tiles, cex.label = 0.3, sep.draw = FALSE, y0 = 0, y1 = 2))
try(plot(gt_concat(sig_track, BRCA_cnt_gt, NSCLC_cnt_gt),   chr_tiles, cex.label = 0.3, sep.draw = FALSE, y0 = 0))
try(plot(gt_concat(gt_ge, sig_track, NSCLC_freq_gt, BRCA_freq_gt, NSCLC_cnt_gt, BRCA_cnt_gt),
         gene_window, cex.label = 0.1, sep.draw = FALSE))
invisible(dev.off())
cat('Done.\n')
