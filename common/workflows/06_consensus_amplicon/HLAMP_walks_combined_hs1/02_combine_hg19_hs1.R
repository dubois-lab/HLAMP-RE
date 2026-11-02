## =============================================================================
## Combine the hg19 and hs1 (T2T-CHM13 v2.0) walk matrices on one hs1 10 kb grid
## -----------------------------------------------------------------------------
## Input per gene and cohort: two per-sample copy-number matrices built from the
## JaBbA walks (10 kb bins, max CN of any walk covering the bin):
##   hg19 arm  samples reconstructed on hg19 (HMF / PCAWG)
##   hs1  arm  samples reconstructed on hs1 (ONT)
##
## The hg19 arm is projected onto hs1: every hs1 10 kb tile inherits its CN from
## the hg19 10 kb bin that covers most of it after liftOver.  The winning bin is
## chosen against the COMPLETE hg19 10 kb grid (the same tileGenome() the walk
## scripts build their bins from), and only then looked up in the walks matrix.
## If the winning bin carries no walk, the cell is NA.
##
## Rows: the hs1 arm tiles + the hs1 tiles spanned by the lifted hg19 walk bins.
##
## Per tile: n_samples_total, n_samples_obs (non-NA), n_samples_CNge10 / AmpSum
## (CN >= 10), AmpPercent (= AmpSum / n_samples_total), fraction_CNge10_obs,
## MedianCN, AmpMean, followed by one column per sample.
##
## IN   <ARM_DIR>/<GENE>/hg19/{BRCA,NSCLC}_<GENE>_CNMatm_JabbaWlks.txt
##      <ARM_DIR>/<GENE>/hs1/{BRCA,LUSC}_<GENE>_CNMatm_JabbaWlks_hs1.txt
## OUT  <OUT_DIR>/<GENE>/{BRCA,NSCLC_LUSC}_<GENE>_CNMatm_combined_hs1.{txt.gz,rds}
##
## Usage (env): GENE=CCND1|NSD3|EGFR ARM_DIR=... [OUT_DIR=output] [REF_DIR=...]
##              Rscript 02_combine_hg19_hs1.R
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(IRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
})
options(scipen = 999)

GENE <- toupper(Sys.getenv('GENE', unset = 'CCND1'))
stopifnot(GENE %in% c('CCND1', 'NSD3', 'EGFR'))
ARM_DIR <- Sys.getenv('ARM_DIR'); stopifnot(nzchar(ARM_DIR))
REF_DIR <- Sys.getenv('REF_DIR', unset = '../../../data/refGenome')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
out_dir <- file.path(OUT_DIR, GENE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
chain   <- import.chain(file.path(REF_DIR, 'hg19ToHs1.over.chain'))

arm <- function(a) file.path(ARM_DIR, GENE, a)

meta_cols <- c('binID','seqnames','start','end','strand','width',
               'geneName','Gene_Symbol','Name','Entrez_GeneId','Genome_Location','Tier',
               'Hallmark','Chr_Band','Somatic','Germline',
               'Tumour_Types.Somatic.','Tumour_Types.Germline.','Cancer_Syndrome',
               'Tissue_Type','Molecular_Genetics','Role_in_Cancer','Mutation_Types',
               'Translocation_Partner','Other_Germline_Mut','Other_Syndrome',
               'Synonyms','locC','grl.ix','grl.iix','source',
               'AmpPercent','AmpMean','AmpSum','amped',
               'n_samples_total','n_samples_obs','n_samples_CNge10',
               'fraction_CNge10_obs','MedianCN')

cat(sprintf('=== %s | out=%s ===\n', GENE, OUT_DIR))

## -----------------------------------------------------------------------------
## The complete hg19 10 kb grid, lifted once and reused for both cohorts
## -----------------------------------------------------------------------------
## Same construction as the walk scripts' bins, so every bin a walks matrix can
## contain is present here and is matched by binID.  Bins that do not lift, or
## lift onto more than one chromosome, are dropped.
cat('Lifting the complete hg19 10 kb grid ...\n')
ref <- tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),
                  tilewidth = 10000, cut.last.tile.in.chrom = TRUE)
ref_chr <- sub('^chr', '', as.character(seqnames(ref)))
ref_binID <- paste(ref_chr, start(ref), end(ref), sep = ':')
lifted_ref <- suppressMessages(suppressWarnings(liftOver(ref, chain)))
nh  <- elementNROWS(lifted_ref)
nuc <- lengths(unique(as(seqnames(lifted_ref), 'CharacterList')))
keep_ref <- nh >= 1 & nuc == 1
lr  <- lifted_ref[keep_ref]
rseg      <- unlist(lr, use.names = FALSE)
rseg_ref  <- rep(which(keep_ref), elementNROWS(lr))
rseg_grr  <- GRanges(sub('^chr', '', as.character(seqnames(rseg))),
                     IRanges(start(rseg), end(rseg)))
cat(sprintf('reference grid: %d bins, %d lift (%.1f%%), %d lifted segments\n',
            length(ref), sum(keep_ref), 100 * sum(keep_ref) / length(ref), length(rseg_grr)))
rm(lifted_ref, lr, rseg); invisible(gc(verbose = FALSE))

combine_one <- function(tag, hg19_path, hs1_path, out_prefix) {
  cat(sprintf('\n=== %s %s ===\n', GENE, tag))

  hg19 <- fread(hg19_path, sep = '\t', header = TRUE)
  hs1  <- fread(hs1_path,  sep = '\t', header = TRUE)
  hg19_samples <- setdiff(colnames(hg19), meta_cols)
  hs1_samples  <- setdiff(colnames(hs1),  meta_cols)
  cat(sprintf('hg19 samples: %d, hs1 samples: %d\n', length(hg19_samples), length(hs1_samples)))
  stopifnot(length(intersect(hg19_samples, hs1_samples)) == 0)

  ## ---- tile universe --------------------------------------------------------
  ## hs1 arm tiles, unioned with the hs1 tiles spanned by the lifted walk bins
  hp <- tstrsplit(hs1$binID, ':', fixed = TRUE)
  hs1[, `:=`(hs1_chr = hp[[1]], hs1_start = as.integer(hp[[2]]), hs1_end = as.integer(hp[[3]]))]
  hs1_tiles <- unique(hs1[, .(seqnames = hs1_chr, start = hs1_start, end = hs1_end, binID)])

  wi <- match(hg19$binID, ref_binID)
  if (anyNA(wi)) stop(sum(is.na(wi)), ' walks hg19 binIDs are not on the reference grid')
  inw <- rep(FALSE, length(ref)); inw[wi] <- TRUE
  wseg <- rseg_grr[inw[rseg_ref]]
  ts <- as.integer(floor((start(wseg) - 1L) / 10000L) * 10000L + 1L)
  te <- as.integer(floor((end(wseg)   - 1L) / 10000L) * 10000L + 1L)
  nt <- as.integer((te - ts) / 10000L) + 1L
  ix <- rep.int(seq_along(ts), nt)
  st <- ts[ix] + (sequence(nt) - 1L) * 10000L
  hg19_tiles <- data.table(seqnames = as.character(seqnames(wseg))[ix],
                           start = st, end = st + 9999L)
  hg19_tiles[, binID := paste(seqnames, start, end, sep = ':')]
  all_tiles <- unique(rbindlist(list(hs1_tiles, hg19_tiles)), by = 'binID')
  cat(sprintf('Tiles: %d (hs1 arm=%d, hg19-derived=%d)\n',
              nrow(all_tiles), nrow(hs1_tiles), uniqueN(hg19_tiles$binID)))
  rm(hg19_tiles, wseg); invisible(gc(verbose = FALSE))

  setorder(all_tiles, seqnames, start)
  if (any(grepl('e+', all_tiles$binID, fixed = TRUE))) stop('scientific-notation binID')
  if (anyDuplicated(all_tiles[, .(seqnames, start, end)])) stop('duplicate tile coordinates')
  tile_gr <- GRanges(all_tiles$seqnames, IRanges(all_tiles$start, all_tiles$end))

  ## ---- tile -> hg19 bin: max overlap against the complete hg19 grid ---------
  hits <- findOverlaps(tile_gr, rseg_grr)
  ov   <- pintersect(tile_gr[queryHits(hits)], rseg_grr[subjectHits(hits)])
  hd   <- data.table(q = queryHits(hits), refrow = rseg_ref[subjectHits(hits)],
                     w = width(ov))
  setorder(hd, q, -w)
  hd <- hd[!duplicated(q)]
  rm(hits, ov); invisible(gc(verbose = FALSE))
  src_binID <- rep(NA_character_, nrow(all_tiles))
  src_binID[hd$q] <- ref_binID[hd$refrow]
  src_row <- match(src_binID, hg19$binID)          # NA = source bin has no walk
  cat(sprintf('tiles with a covering hg19 bin: %d / %d\n', sum(!is.na(src_binID)), nrow(all_tiles)))
  cat(sprintf('  of those, the covering bin carries a walk: %d\n', sum(!is.na(src_row))))

  ## ---- fill ------------------------------------------------------------------
  all_samples <- c(hg19_samples, hs1_samples)
  M <- matrix(NA_integer_, nrow = nrow(all_tiles), ncol = length(all_samples),
              dimnames = list(NULL, all_samples))
  for (s in hg19_samples) M[, s] <- as.integer(hg19[[s]])[src_row]
  m_hs1 <- match(all_tiles$binID, hs1$binID); hit <- !is.na(m_hs1)
  for (s in hs1_samples) M[hit, s] <- as.integer(hs1[[s]])[m_hs1[hit]]

  ## sanity check: each hg19 sample's maximum CN survives the projection
  mism <- character(0)
  for (s in hg19_samples) {
    a <- suppressWarnings(max(as.integer(hg19[[s]]), na.rm = TRUE))
    b <- suppressWarnings(max(M[, s], na.rm = TRUE))
    if (!isTRUE(all.equal(a, b))) mism <- c(mism, sprintf('%s in=%g out=%g', s, a, b))
  }
  cat(sprintf('Input==projected max: %d / %d hg19 samples OK\n',
              length(hg19_samples) - length(mism), length(hg19_samples)))
  for (x in mism) cat(sprintf('  MISMATCH %s\n', x))

  n_total <- ncol(M)
  n_amp   <- rowSums(M >= 10L, na.rm = TRUE)
  n_obs   <- rowSums(!is.na(M))
  out_dt <- data.table(binID = all_tiles$binID,
                       seqnames = paste0('chr', all_tiles$seqnames),
                       start = all_tiles$start, end = all_tiles$end,
                       width = 10000L, strand = '*')
  out_dt[, `:=`(n_samples_total     = n_total,
                n_samples_obs       = n_obs,
                n_samples_CNge10    = n_amp,
                AmpPercent          = n_amp / n_total,
                fraction_CNge10_obs = ifelse(n_obs > 0, n_amp / n_obs, NA_real_),
                AmpSum              = n_amp,
                MedianCN            = if (requireNamespace('matrixStats', quietly = TRUE))
                                        matrixStats::rowMedians(M, na.rm = TRUE)
                                      else apply(M, 1, median, na.rm = TRUE),
                AmpMean             = rowMeans(M, na.rm = TRUE))]
  out_dt <- cbind(out_dt, as.data.table(M))

  fwrite(out_dt,  file.path(out_dir, paste0(out_prefix, '_combined_hs1.txt.gz')), sep = '\t')
  saveRDS(out_dt, file.path(out_dir, paste0(out_prefix, '_combined_hs1.rds')))
  cat(sprintf('Wrote %s (%d rows, %d samples)\n', paste0(out_prefix, '_combined_hs1.txt.gz'),
              nrow(out_dt), n_total))
  rm(M, out_dt, hg19, hs1); invisible(gc(verbose = FALSE))
  invisible(NULL)
}

combine_one('BRCA',
  hg19_path  = file.path(arm('hg19'), sprintf('BRCA_%s_CNMatm_JabbaWlks.txt', GENE)),
  hs1_path   = file.path(arm('hs1'),  sprintf('BRCA_%s_CNMatm_JabbaWlks_hs1.txt', GENE)),
  out_prefix = sprintf('BRCA_%s_CNMatm', GENE))

combine_one('NSCLC',
  hg19_path  = file.path(arm('hg19'), sprintf('NSCLC_%s_CNMatm_JabbaWlks.txt', GENE)),
  hs1_path   = file.path(arm('hs1'),  sprintf('LUSC_%s_CNMatm_JabbaWlks_hs1.txt', GENE)),
  out_prefix = sprintf('NSCLC_LUSC_%s_CNMatm', GENE))

cat('\nDone.\n')
