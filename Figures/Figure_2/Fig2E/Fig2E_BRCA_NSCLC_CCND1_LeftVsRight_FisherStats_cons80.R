## =====================================================================
## Fig.2E -- BRCA vs NSCLC CCND1 Left-vs-Right mirror-pair Fisher (cons80)
##
## Window = gene +/- 20e5, cons80 mirror-pairing.
## =====================================================================
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(qvalue)
  library(gTrack)
  library(gUtils)
})
options(scipen = 999)

xG       <- "CCND1"
CNcutoff <- 10
AMP_FREQ_MIN <- 0.80   # cons80: test only if at least one paired bin has AmpFreq >= 80%
date_tag <- "20260907c"

# combined matrices
data_dir   <- IN_DIR
brca_path  <- file.path(data_dir, "20260824_BRCA_CCND1_CNMatm_combined_hs1.txt")
nsclc_path <- file.path(data_dir, "20260824_NSCLC_LUSC_CCND1_CNMatm_combined_hs1.txt")

# output folders for the saved tracks and the pdf
src_dir  <- file.path(OUT_DIR, "source_data")
pdf_dir  <- OUT_DIR
dir.create(src_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)

# ---- T2T reference tracks --------------------------------------------------
gt_ge   <- readRDS(file.path(IN_DIR, "gencodeTrack_t2t.rds"))
CGC     <- readRDS(file.path(IN_DIR, "CGC_gr_T2T.rds"))
gene_gr <- CGC[CGC$Gene_Symbol == xG]
stopifnot(length(gene_gr) >= 1)
gene_gr <- gene_gr[1]

intchr   <- as.character(seqnames(gene_gr))
gene_st  <- start(gene_gr)
gene_en  <- end(gene_gr)
ccnd1_window <- GRanges(seqnames = intchr,
                        ranges   = IRanges(start = 67657962, end = 71669679),
                        seqinfo  = seqinfo(gene_gr))
win_st   <- start(ccnd1_window)
win_en   <- end(ccnd1_window)
cat(sprintf("CCND1 = %s:%d-%d ; window = %s:%d-%d\n",
            intchr, gene_st, gene_en,
            as.character(seqnames(ccnd1_window)), win_st, win_en))

# BRCAvsNSCLC chr-level cons80 file
brcaNSCLC_chr_path <- file.path(IN_DIR, '20260824_BRCAvsNSCLC_CCND1_chr11_FisherStats_CN10_cons80.txt')

keep_cols <- c("binID", "seqnames", "start", "end", "width", "strand",
               "n_samples_total", "n_samples_obs", "n_samples_CNge10",
               "AmpPercent", "fraction_CNge10_obs", "AmpSum", "MedianCN", "AmpMean")

brca  <- fread(brca_path,  select = keep_cols)
nsclc <- fread(nsclc_path, select = keep_cols)

window_bins <- function(d) {
  d <- d[seqnames == intchr]
  d[start >= win_st & end <= win_en]
}


brca_out_path  <- file.path(src_dir, sprintf("%s_BRCA_%s_LeftVsRight_FisherStats_CN%s_cons80.txt",  date_tag, xG, CNcutoff))
nsclc_out_path <- file.path(src_dir, sprintf("%s_NSCLC_%s_LeftVsRight_FisherStats_CN%s_cons80.txt", date_tag, xG, CNcutoff))

# ---- Cohort frequency / count tracks (full window) --------------
make_cohort_tracks <- function(d, tag, col) {
  d <- d[seqnames == intchr]
  d[, AmpFreq  := ifelse(n_samples_total > 0, AmpSum / n_samples_total, 0)]
  d[is.na(AmpFreq), AmpFreq := 0]
  d[, AmpCount := as.integer(AmpSum)]
  gr <- dt2gr(d[, .(seqnames, start, end, AmpFreq, AmpCount)])
  list(
    freq = gTrack(gr, y.field = "AmpFreq",  bars = TRUE, col = col,
                  name = paste0(tag, " AmpFreq"),  y0 = 0, y1 = 1, ygap =0.5),
    cnt  = gTrack(gr, y.field = "AmpCount", bars = TRUE, col = col,
                  name = paste0(tag, " AmpCount"), y0 = 0)
  )
}

brca_tracks  <- make_cohort_tracks(brca,  "BRCA",  "red")
nsclc_tracks <- make_cohort_tracks(nsclc, "NSCLC", "blue")

# ---- L-vs-R sig_tracks (cons80) --------------------------------------------

brca_pairs  <- fread(file.path(IN_DIR, '20260824_BRCA_CCND1_LeftVsRight_FisherStats_CN10_cons80.txt'))
nsclc_pairs  <- fread(file.path(IN_DIR, '20260824_NSCLC_CCND1_LeftVsRight_FisherStats_CN10_cons80.txt'))

# ---- L-vs-R sig_tracks (cons80) --------------------------------------------
build_LR_sig_track <- function(pairs, tag, col, p_thresh = 0.1) {
  if (!nrow(pairs)) return(NULL)
  sig <- pairs[!is.na(qval) & qval < p_thresh]
  if (!nrow(sig)) return(NULL)
  sig[, higher_side := ifelse(right_samples_amped >= left_samples_amped, "right", "left")]
  sig[, sig_chr     := ifelse(higher_side == "right", right_chr,   left_chr)]
  sig[, sig_start   := ifelse(higher_side == "right", right_start, left_start)]
  sig[, sig_end     := ifelse(higher_side == "right", right_end,   left_end)]
  sig[, neglog10_Fisher_pval := -log10(Fisher_pval)]
  gr <- GRanges(sig$sig_chr, IRanges(sig$sig_start, sig$sig_end))
  mcols(gr)$neglog10_Fisher_pval <- sig$neglog10_Fisher_pval
  cat(sprintf("[%s L-vs-R cons80] sig bins (p<%.2g): %d\n", tag, p_thresh, length(gr)))
  gTrack(gr, bars = TRUE, col = col,
         name = paste0(tag, " L-vs-R -log10(p)"), height = 0.5)
}

brca_sig_LR  <- build_LR_sig_track(brca_pairs,  "BRCA",  "red")
nsclc_sig_LR <- build_LR_sig_track(nsclc_pairs, "NSCLC", "blue")


# ---- BRCAvsNSCLC chr11 sig_track from cons80 file (carried into page) -----
brcaNSCLC_sig_track <- NULL
if (file.exists(brcaNSCLC_chr_path)) {
  BRCAvsNSCLC_dt <- fread(brcaNSCLC_chr_path)
  BRCAvsNSCLC_dt[, colr := ifelse(frac_amp_BRCA > frac_amp_NSCLC, 'red', 'blue')]
  csig <- BRCAvsNSCLC_dt[!is.na(qval) & qval < 0.05]
  if (nrow(csig) > 0) {
    csig[, neglog10_Fisher_pval := -log10(Fisher_pval)]
    cgr <- GRanges(csig$chr, IRanges(csig$start, csig$end))
    mcols(cgr)$neglog10_Fisher_pval <- csig$neglog10_Fisher_pval
    brcaNSCLC_sig_track <- gTrack(cgr, 
                                  bars = TRUE, col = csig$colr[1],
                                  name = "BRCAvsNSCLC -log10(p) cons80",
                                  height = 0.5)
    cat(sprintf("Loaded BRCAvsNSCLC cons80 chr11 sig bins (p<0.05): %d\n", length(cgr)))
  }
} else {
  cat(sprintf("Note: cons80 chr11 stats file not found at %s\n", brcaNSCLC_chr_path))
}

gt_concat <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  Reduce(c, parts)
}

# ---- Save gTracks as source data ------------------------------------------
save_track <- function(obj, fname) {
  if (is.null(obj)) {
    cat(sprintf("Skipped (NULL): %s\n", fname))
    return(invisible(NULL))
  }
  fpath <- file.path(src_dir, fname)
  saveRDS(obj, fpath)
  cat(sprintf("Wrote: %s\n", fpath))
}

save_track(nsclc_sig_LR,        sprintf("%s_Fig.2C_NSCLC_%s_LvsR_sig_LR_cons80.rds",       date_tag, xG))
save_track(nsclc_tracks$freq,   sprintf("%s_Fig.2C_NSCLC_%s_AmpFreq_track.rds",            date_tag, xG))
save_track(brca_sig_LR,         sprintf("%s_Fig.2C_BRCA_%s_LvsR_sig_LR_cons80.rds",        date_tag, xG))
save_track(brca_tracks$freq,    sprintf("%s_Fig.2C_BRCA_%s_AmpFreq_track.rds",             date_tag, xG))
save_track(brcaNSCLC_sig_track, sprintf("%s_Fig.2C_BRCAvsNSCLC_%s_chr11_sig_track_cons80.rds", date_tag, xG))

# ---- Plot (CCND1 +/- 20e5 window, single page) -----------------------------
pdf_path <- file.path(pdf_dir, sprintf("%s_BRCA_NSCLC_%s_LeftVsRight_FisherStats_CN%s_cons80_p0.05.pdf",
                                       date_tag, xG, CNcutoff))
pdf(file = pdf_path, paper = "a4", height = 12)
# gTrack: last in c() is plotted at the top
plot(gt_concat(gt_ge,
               nsclc_tracks$freq, nsclc_sig_LR,brcaNSCLC_sig_track, brca_tracks$freq, brca_sig_LR),
     ccnd1_window, cex.label = 0.1, sep.draw = FALSE)
dev.off()
cat(sprintf("Wrote PDF: %s\n", pdf_path))
