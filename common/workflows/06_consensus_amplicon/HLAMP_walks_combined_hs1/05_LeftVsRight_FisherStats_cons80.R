## =============================================================================
## Left-vs-Right mirror-pair Fisher tests around the gene -- cons80
## -----------------------------------------------------------------------------
## In gene +/- 2.5 Mb, the i-th bin right of the gene is paired with the i-th
## bin left of it.  Within one cohort, a 2x2 Fisher test compares the amplified
## counts (CN >= 10) of the two bins; only pairs where one member reaches
## AmpFreq >= 80% (cons80) are tested.  q-values by qvalue(), BH fallback.
## Pairing is by rank among the bins present in the combined matrix.
##
## IN   <OUT_DIR>/<GENE>/{BRCA,NSCLC_LUSC}_<GENE>_CNMatm_combined_hs1.txt.gz
##      <OUT_DIR>/<GENE>/stats/cons80/<GENE>_BRCAvsNSCLC_FisherStats_CN10_cons80.txt
##        (from 04_; optional, only drawn as an extra track)
## OUT  <OUT_DIR>/<GENE>/stats/cons80/<GENE>_{BRCA,NSCLC}_LeftVsRight_*
##
## Usage (env): GENE=CCND1|NSD3|EGFR [OUT_DIR=output] [REF_DIR=...]
##              [AMP_FREQ_MIN=0.80] [CN_CUTOFF=10, used in file names only]
##              Rscript 05_LeftVsRight_FisherStats_cons80.R
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
intchr  <- as.character(seqnames(gene_gr))
gene_st <- start(gene_gr); gene_en <- end(gene_gr)
gene_window <- gene_gr + 25e5
win_st <- start(gene_window); win_en <- end(gene_window)
cat(sprintf('%s = %s:%d-%d ; window = %s:%d-%d\n', xG, intchr, gene_st, gene_en, intchr, win_st, win_en))

keep_cols <- c('binID', 'seqnames', 'start', 'end', 'n_samples_total', 'AmpSum')
brca  <- fread(brca_path,  select = keep_cols)
nsclc <- fread(nsclc_path, select = keep_cols)

mirror_pairs <- function(d, tag) {
  bw <- d[seqnames == intchr & start >= win_st & end <= win_en]
  right_bins <- bw[start > gene_en]; setorder(right_bins, start)
  left_bins  <- bw[end   < gene_st]; setorder(left_bins, -end)
  npair <- min(nrow(right_bins), nrow(left_bins))
  cat(sprintf('[%s] window bins: right=%d, left=%d, candidate pairs=%d\n',
              tag, nrow(right_bins), nrow(left_bins), npair))
  if (npair == 0) return(data.table())
  rows <- lapply(seq_len(npair), function(i) {
    rb <- right_bins[i]; lb <- left_bins[i]
    r_a <- as.integer(rb$AmpSum); l_a <- as.integer(lb$AmpSum)
    r_tot <- as.integer(rb$n_samples_total); l_tot <- as.integer(lb$n_samples_total)
    r_n <- pmax(0L, r_tot - r_a); l_n <- pmax(0L, l_tot - l_a)
    r_freq <- if (r_tot > 0) r_a / r_tot else NA_real_
    l_freq <- if (l_tot > 0) l_a / l_tot else NA_real_
    pass <- !is.na(r_freq) & !is.na(l_freq) & ((r_freq >= AMP_FREQ_MIN) | (l_freq >= AMP_FREQ_MIN))
    pv <- NA_real_
    if (pass) {
      cm <- matrix(c(r_a, r_n, l_a, l_n), nrow = 2, byrow = TRUE)
      if (!any(rowSums(cm) == 0)) pv <- fisher.test(cm)$p.value
    }
    data.table(pair_index = i,
               right_chr = as.character(rb$seqnames), right_start = as.integer(rb$start), right_end = as.integer(rb$end),
               left_chr  = as.character(lb$seqnames), left_start  = as.integer(lb$start), left_end  = as.integer(lb$end),
               right_samples_amped = r_a, right_samples_non_amped = r_n,
               left_samples_amped  = l_a, left_samples_non_amped  = l_n,
               right_AmpFreq = r_freq, left_AmpFreq = l_freq,
               cons80_pass = pass, Fisher_pval = pv)
  })
  out <- rbindlist(rows)[cons80_pass == TRUE]
  cat(sprintf('[%s] cons80 (>=%.0f%% in at least one paired bin): %d / %d\n',
              tag, AMP_FREQ_MIN * 100, nrow(out), npair))
  out[, cons80_pass := NULL]
  ok <- !is.na(out$Fisher_pval); qv <- rep(NA_real_, nrow(out))
  if (sum(ok)) {
    qres <- try(qvalue(out$Fisher_pval[ok])$qvalues, silent = TRUE)
    if (inherits(qres, 'try-error')) {
      cat(sprintf('[%s] qvalue() failed; falling back to BH.\n', tag))
      qv[ok] <- p.adjust(out$Fisher_pval[ok], method = 'BH')
    } else qv[ok] <- qres
  }
  out[, qval := qv]
  out
}
brca_pairs  <- mirror_pairs(brca,  'BRCA')
nsclc_pairs <- mirror_pairs(nsclc, 'NSCLC')
fwrite(brca_pairs,  file.path(out_dir, sprintf('%s_BRCA_LeftVsRight_FisherStats_CN%s_cons80.txt', xG, CNcutoff)), sep = '\t', quote = FALSE)
fwrite(nsclc_pairs, file.path(out_dir, sprintf('%s_NSCLC_LeftVsRight_FisherStats_CN%s_cons80.txt', xG, CNcutoff)), sep = '\t', quote = FALSE)
cat(sprintf('[BRCA] p<0.05: %d | [NSCLC] p<0.05: %d\n',
            sum(brca_pairs$Fisher_pval < 0.05, na.rm = TRUE), sum(nsclc_pairs$Fisher_pval < 0.05, na.rm = TRUE)))

make_cohort_tracks <- function(d, tag, col) {
  d <- d[seqnames == intchr]
  d[, AmpFreq := ifelse(n_samples_total > 0, AmpSum / n_samples_total, 0)]
  d[, AmpCount := as.integer(AmpSum)]
  gr <- dt2gr(d[, .(seqnames, start, end, AmpFreq, AmpCount)])
  list(freq = gTrack(gr, y.field = 'AmpFreq',  bars = TRUE, col = col, name = paste0(tag, ' AmpFreq'),  y0 = 0, y1 = 1),
       cnt  = gTrack(gr, y.field = 'AmpCount', bars = TRUE, col = col, name = paste0(tag, ' AmpCount'), y0 = 0))
}
brca_tracks  <- make_cohort_tracks(brca,  'BRCA',  'red')
nsclc_tracks <- make_cohort_tracks(nsclc, 'NSCLC', 'blue')

LR_sig_track <- function(pairs, tag, col) {
  if (!nrow(pairs)) return(NULL)
  sig <- pairs[!is.na(Fisher_pval) & Fisher_pval < 0.05]
  if (!nrow(sig)) return(NULL)
  right <- sig$right_samples_amped >= sig$left_samples_amped
  gr <- GRanges(ifelse(right, sig$right_chr, sig$left_chr),
                IRanges(ifelse(right, sig$right_start, sig$left_start), ifelse(right, sig$right_end, sig$left_end)))
  gr$neglog10_Fisher_pval <- -log10(sig$Fisher_pval)
  gTrack(gr, y.field = 'neglog10_Fisher_pval', bars = TRUE, col = col,
         name = paste0(tag, ' L-vs-R -log10(p)'), y0 = 1, height = 2)
}
chr_path <- file.path(out_dir, sprintf('%s_BRCAvsNSCLC_FisherStats_CN%s_cons80.txt', xG, CNcutoff))
bn_track <- NULL
if (file.exists(chr_path)) {
  cs <- fread(chr_path)[!is.na(Fisher_pval) & Fisher_pval < 0.05]
  if (nrow(cs)) {
    cgr <- GRanges(cs$chr, IRanges(cs$start, cs$end)); cgr$neglog10_Fisher_pval <- -log10(cs$Fisher_pval)
    bn_track <- gTrack(cgr, y.field = 'neglog10_Fisher_pval', bars = TRUE, col = 'black',
                       name = 'BRCAvsNSCLC -log10(p) cons80', y0 = 1, height = 2)
  }
}
gt_concat <- function(...) Reduce(c, Filter(Negate(is.null), list(...)))
pdf(file.path(out_dir, sprintf('%s_BRCA_NSCLC_LeftVsRight_FisherStats_CN%s_cons80_p0.05.pdf', xG, CNcutoff)),
    paper = 'a4', height = 12)
plot(gt_concat(gt_ge, bn_track, nsclc_tracks$freq, brca_tracks$freq, nsclc_tracks$cnt, brca_tracks$cnt,
               LR_sig_track(brca_pairs, 'BRCA', 'red'), LR_sig_track(nsclc_pairs, 'NSCLC', 'blue')),
     gene_window, cex.label = 0.1, sep.draw = FALSE)
invisible(dev.off())
cat('Done.\n')
