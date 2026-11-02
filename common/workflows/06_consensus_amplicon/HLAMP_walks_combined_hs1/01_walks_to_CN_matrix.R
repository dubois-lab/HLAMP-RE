## =============================================================================
## JaBbA walks -> per-sample co-amplification CN matrix (one genome build)
## -----------------------------------------------------------------------------
## For one gene and one build (hg19 or hs1), from every walk file listed in the
## samplesheet:
##
##   1. A walk QUALIFIES if its segments at CN >= MIN_CN cover >= MIN_COV of the
##      gene target interval (below).
##   2. All qualifying walks of a sample, over all of its walk files, are
##      collapsed into one segment set (union of segments).
##   3. A sample enters the matrix if it has at least one qualifying walk.  If a
##      `patient` is given, only one sample per patient is kept: the one with the
##      highest gene coverage, then the highest CN at the gene of its best walk,
##      then the sample ID.
##   4. Each 10 kb bin of the build's genome gets, per sample, the MAX CN over the
##      collapsed segments hitting it; bins hit by no segment of a sample are CN 2.
##      Rows are the bins hit by at least one sample.
##
##
## Gene target intervals (1-based, inclusive):
##   hg19  CCND1  11:69,455,873-69,469,242   NM_053056 coding sequence
##         NSD3    8:38,127,215-38,239,790   CGC gene element (WHSC1L1)
##         EGFR    7:55,086,714-55,324,313   CGC gene element
##   hs1   from cgc_T2T_gr_consAmpWalks.rds (native T2T gene element)
##
## SAMPLESHEET  csv, one row per walk file (gGnome gWalk .rds):
##   sample     sample ID (becomes the matrix column)
##   cohort     BRCA / NSCLC (hg19), BRCA / LUSC (hs1); used in output names
##   patient    patient ID for one-sample-per-patient dedup, or empty
##   walk_file  path to the walk .rds
##   priority   integer tie-break between walk sources, lower wins (optional)
##
## OUT  <OUT_DIR>/<GENE>/<BUILD>/<cohort>_<GENE>_CNMatm_JabbaWlks[_hs1].txt
##      <OUT_DIR>/<GENE>/<BUILD>/<cohort>_<GENE>_walks_collapsed_segments.rds
##      <OUT_DIR>/<GENE>/<BUILD>/sample_status.csv
##
## Usage (env): GENE=CCND1|NSD3|EGFR BUILD=hg19|hs1 SAMPLESHEET=... ARM_DIR=...
##              [REF_DIR=...] [MIN_GENE_COV=0.30] [MIN_CN=10] [MC_CORES=8]
##              Rscript 01_walks_to_CN_matrix.R
## =============================================================================

suppressPackageStartupMessages({
  library(gGnome)
  library(parallel)
  library(data.table)
  library(GenomicRanges)
})
options(scipen = 999)
## gGnome attaches BiocGenerics, whose S4 `merge` bypasses data.table's method
mrg <- function(...) data.table::merge.data.table(...)

GENE  <- toupper(Sys.getenv('GENE', unset = 'CCND1'))
stopifnot(GENE %in% c('CCND1', 'NSD3', 'EGFR'))
BUILD <- tolower(Sys.getenv('BUILD', unset = 'hg19'))
stopifnot(BUILD %in% c('hg19', 'hs1'))
SAMPLESHEET <- Sys.getenv('SAMPLESHEET'); stopifnot(nzchar(SAMPLESHEET))
ARM_DIR  <- Sys.getenv('ARM_DIR'); stopifnot(nzchar(ARM_DIR))
REF_DIR  <- Sys.getenv('REF_DIR', unset = '../../../data/refGenome')
MIN_COV  <- as.numeric(Sys.getenv('MIN_GENE_COV', unset = '0.30'))
MIN_CN   <- as.numeric(Sys.getenv('MIN_CN', unset = '10'))
MC_CORES <- as.integer(Sys.getenv('MC_CORES', unset = '8'))
CHUNK    <- 400L
out_dir  <- file.path(ARM_DIR, GENE, BUILD)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- genome bins and gene target --------------------------------------------
if (BUILD == 'hg19') {
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  seql <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
  TGT <- list(CCND1 = list(chr = '11', start = 69455873L, end = 69469242L),
              NSD3  = list(chr = '8',  start = 38127215L, end = 38239790L),
              EGFR  = list(chr = '7',  start = 55086714L, end = 55324313L))[[GENE]]
} else {
  sz <- fread(file.path(REF_DIR, 'hs1.chrom.sizes'), header = FALSE, col.names = c('seqnames', 'seqlength'))
  sz <- sz[seqnames %in% paste0('chr', c(1:22, 'X', 'Y', 'M'))]
  seql <- setNames(sz$seqlength, sz$seqnames)
  g <- readRDS(file.path(REF_DIR, 'cgc_T2T_gr_consAmpWalks.rds'))
  g <- g[g$Gene_Symbol == GENE][1]
  TGT <- list(chr = sub('^chr', '', as.character(seqnames(g))),
              start = as.integer(start(g)), end = as.integer(end(g)))
}
bins <- gr.sub(tileGenome(seqlengths = seql, tilewidth = 10000, cut.last.tile.in.chrom = TRUE), 'chr', '')
bins_dt <- gr2dt(bins)[, .(seqnames = as.character(seqnames), start, end)]
TGT_LEN    <- TGT$end - TGT$start + 1L
MIN_COV_BP <- as.integer(ceiling(MIN_COV * TGT_LEN))
cat(sprintf('=== %s %s | target %s:%d-%d (%d bp) | qualifying walk: >=%d bp at CN>=%g ===\n',
            GENE, BUILD, TGT$chr, TGT$start, TGT$end, TGT_LEN, MIN_COV_BP, MIN_CN))

## ---- samplesheet ------------------------------------------------------------
ss <- fread(SAMPLESHEET, colClasses = 'character')
stopifnot(all(c('sample', 'cohort', 'walk_file') %in% names(ss)))
if (!'patient'  %in% names(ss)) ss[, patient := NA_character_]
if (!'priority' %in% names(ss)) ss[, priority := '1']
ss[!nzchar(patient), patient := NA_character_]
ss[, priority := as.integer(priority)]
if (nrow(unique(ss[, .(sample, cohort)])) != uniqueN(ss$sample)) stop('a sample has more than one cohort')
ss[, base := basename(walk_file)]
cm <- regexpr('cnLrg[0-9]+(\\.[0-9]+)?', ss$base)
fcn <- rep(0, nrow(ss)); fcn[cm > 0] <- as.numeric(sub('cnLrg', '', regmatches(ss$base, cm)))
ss[, file_cn := fcn]
cat(sprintf('samplesheet: %d walk files, %d samples\n', nrow(ss), uniqueN(ss$sample)))

## ---- walk segments ----------------------------------------------------------
## Node table + per-walk node lists when the gWalk carries its graph (fast), the
## walk GRangesList otherwise.  Returns one row per (walk, segment).
walk_segments <- function(w) {
  nodes <- tryCatch(as.data.table(w$graph$nodes$gr), error = function(e) NULL)
  sn    <- tryCatch(w$dt$snode.id, error = function(e) NULL)
  if (!is.null(nodes) && !is.null(sn) && is.list(sn) &&
      all(c('node.id', 'start', 'end') %in% names(nodes))) {
    nd <- nodes[, .(node.id  = as.integer(node.id),
                    seqnames = sub('^chr', '', as.character(seqnames)),
                    start    = as.integer(start), end = as.integer(end),
                    cn       = if ('cn' %in% names(nodes)) as.numeric(cn) else NA_real_)]
    flat <- abs(as.integer(unlist(sn, use.names = FALSE)))
    wvec <- rep.int(seq_along(sn), vapply(sn, length, integer(1)))
    out  <- nd[match(flat, nd$node.id), .(seqnames, start, end, cn)]
    out[, wix := wvec]
    return(out)
  }
  d <- tryCatch(gr2dt(w$grl), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  if (!'cn' %in% names(d)) d[, cn := NA_real_]
  d[, .(seqnames = sub('^chr', '', as.character(seqnames)), start = as.integer(start),
        end = as.integer(end), cn = as.numeric(cn), wix = as.integer(group))]
}

## bp of the target covered by a set of segments
tgt_bp <- function(s, e) {
  s <- pmax(s, TGT$start); e <- pmin(e, TGT$end)
  if (!length(s)) return(0L)
  as.integer(sum(IRanges::width(IRanges::reduce(IRanges::IRanges(s, e)))))
}

## per file: the union of the segments of ALL its qualifying walks (-> matrix),
## and its best walk (-> only to rank samples of the same patient)
scan_file <- function(i) {
  w <- tryCatch(readRDS(ss$walk_file[i]), error = function(e) NULL)
  if (is.null(w)) return(list(info = data.table(i = i, ok = FALSE)))
  if (!isTRUE(length(w) > 0)) return(list(info = data.table(i = i, ok = TRUE, cov_bp = 0L)))
  segs <- walk_segments(w)
  if (is.null(segs)) return(list(info = data.table(i = i, ok = FALSE)))
  h <- segs[seqnames == TGT$chr & start <= TGT$end & end >= TGT$start]
  if (!nrow(h)) return(list(info = data.table(i = i, ok = TRUE, cov_bp = 0L)))
  cov <- h[, {
    amp <- !is.na(cn) & cn >= MIN_CN
    .(cov_bp = tgt_bp(start[amp], end[amp]), cov_any = tgt_bp(start, end),
      gene_cn = suppressWarnings(max(cn, na.rm = TRUE)))
  }, by = wix]
  spans <- segs[wix %in% cov$wix, .(span = sum(as.numeric(end - start + 1))), by = wix]
  cov[spans, on = 'wix', span := i.span]
  setorder(cov, -cov_bp, -cov_any, -span)
  q <- cov[cov_bp >= MIN_COV_BP]$wix
  list(info = data.table(i = i, ok = TRUE, cov_bp = cov$cov_bp[1], span = cov$span[1],
                         gene_cn = cov$gene_cn[1]),
       segs = if (length(q)) unique(segs[wix %in% q, .(seqnames, start, end, cn)])[, i := i])
}

cat(sprintf('Scanning %d walk files on %d cores ...\n', nrow(ss), MC_CORES))
t0 <- Sys.time(); res <- list()
for (ii in split(seq_len(nrow(ss)), ceiling(seq_len(nrow(ss)) / CHUNK))) {
  r <- mclapply(ii, function(i) tryCatch(scan_file(i), error = function(e)
         list(info = data.table(i = i, ok = FALSE))), mc.cores = MC_CORES)
  bad <- !vapply(r, is.list, logical(1))
  r[bad] <- lapply(ii[bad], function(i) list(info = data.table(i = i, ok = FALSE)))
  res <- c(res, r)
}
fi <- rbindlist(lapply(res, `[[`, 'info'), fill = TRUE)
fi[is.na(cov_bp), cov_bp := 0L]
seg <- rbindlist(lapply(res, `[[`, 'segs'))
rm(res); invisible(gc(verbose = FALSE))
cat(sprintf('Scan done in %.1f min | unreadable files: %d\n',
            as.numeric(Sys.time() - t0, units = 'mins'), sum(!fi$ok)))
if (any(!fi$ok)) print(ss[fi[ok == FALSE, i], .(sample, walk_file)])

## ---- samples: qualification and one sample per patient ----------------------
## gene_cov of the best walk = max over all walks, so `qualifies` is the same as
## "has at least one qualifying walk".  best_walk_cn only breaks ties between
## samples of the same patient; it does not affect what goes into the matrix.
fi <- cbind(fi, ss[fi$i, .(sample, cohort, patient, priority, base, file_cn)])
hits <- fi[cov_bp > 0]
setorder(hits, sample, -cov_bp, -span, -file_cn, priority, base)
status <- hits[!duplicated(sample), .(sample, gene_cov = cov_bp / TGT_LEN,
                                      best_walk_cn = gene_cn, best_walk_file = base)]
status <- mrg(unique(ss[, .(sample, cohort, patient)]), status, by = 'sample', all.x = TRUE)
status[is.na(gene_cov), gene_cov := 0]
status[, qualifies := gene_cov >= MIN_COV]
setorder(status, -gene_cov, -best_walk_cn, sample)
status[, duplicate_patient := FALSE]
status[qualifies & !is.na(patient), duplicate_patient := duplicated(patient)]
status[, in_matrix := qualifies & !duplicate_patient]
setorder(status, cohort, -in_matrix, sample)
fwrite(status, file.path(out_dir, 'sample_status.csv'))
print(status[, .(samples = .N, qualifying = sum(qualifies),
                 duplicate_patient = sum(duplicate_patient), in_matrix = sum(in_matrix)), by = cohort])

seg <- cbind(seg, ss[seg$i, .(sample, cohort)])[, i := NULL]
seg <- unique(seg, by = c('sample', 'seqnames', 'start', 'end', 'cn'))

## ---- CN matrix per cohort ---------------------------------------------------
sfx <- if (BUILD == 'hs1') '_hs1' else ''
for (tt in sort(unique(status$cohort))) {
  keep <- status[cohort == tt & in_matrix == TRUE, sample]
  if (!length(keep)) { cat(sprintf('[%s] no samples\n', tt)); next }
  s_tt <- seg[sample %in% keep & !is.na(cn)]
  saveRDS(s_tt, file.path(out_dir, sprintf('%s_%s_walks_collapsed_segments.rds', tt, GENE)))
  bin_cn <- rbindlist(mclapply(keep, function(fx) {
    sx <- s_tt[sample == fx]
    ov <- findOverlaps(bins, GRanges(sx$seqnames, IRanges(sx$start, sx$end)))
    if (!length(ov)) return(NULL)
    d <- data.table(bin_i = queryHits(ov), cn = sx$cn[subjectHits(ov)])[, .(cn = max(cn)), by = bin_i]
    cbind(bins_dt[d$bin_i], d[, .(cn, sample = fx)])
  }, mc.cores = MC_CORES))
  bin_cn[, binID := paste0(seqnames, ':', start, ':', end)]
  CNMat <- dcast(bin_cn, binID ~ sample, value.var = 'cn', fun.aggregate = max, fill = 2)
  CNMat <- mrg(unique(bin_cn[, .(binID, seqnames, start, end)]), CNMat, by = 'binID')
  write.table(CNMat, file.path(out_dir, sprintf('%s_%s_CNMatm_JabbaWlks%s.txt', tt, GENE, sfx)),
              quote = FALSE, row.names = FALSE, sep = '\t')
  cat(sprintf('[%s] CN matrix: %d bins x %d samples\n', tt, nrow(CNMat), length(keep)))
}
cat('Done.\n')
