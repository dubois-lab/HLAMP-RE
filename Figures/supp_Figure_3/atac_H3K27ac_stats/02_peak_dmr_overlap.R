# Overlaps the peaks tested by 01_peak_enrichment_fisher.R with the differentially methylated regions (DMRs) of the
# matching amplification group inside its consensus amplified regions.
# Run from this folder: Rscript 02_peak_dmr_overlap.R <gene>, with <gene> a name in gene_settings.R
# Writes <gene>_peak_dmr_overlaps.csv with one row per peak with OR > 1 and DMR, <gene>_peak_dmr_overlap_summary.csv
# with the number of DMRs at each peak cutoff, and <gene>_NSCLC_dmrs_in_preferential_amp.csv if the gene has a BRCA vs NSCLC track.

# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(gUtils)
library(rtracklayer)

source('gene_settings.R')

xG <- commandArgs(trailingOnly = TRUE)[1]
if (!xG %in% names(gene_settings)) stop('Usage: Rscript 02_peak_dmr_overlap.R <', paste(names(gene_settings), collapse = '|'), '>')
settings <- gene_settings[[xG]]
if (is.null(settings$overlap_amp_tracks)) stop('gene_settings.R has no overlap settings for ', xG)

hs1ToHg38_chain <- import.chain(file.path(IN_DIR, 'hs1ToHg38.over.chain'))

# Consensus amplified regions: amplified in more than 80% of the walks. Size in bp is reported in hs1
amp_tracks <- lapply(settings$overlap_amp_tracks, function(track_file) readRDS(file.path(IN_DIR, track_file))@data[[1]])
cons_amp_bp <- sapply(amp_tracks, function(amp_gr) sum(width(reduce((amp_gr %Q% (AmpFreq > 0.8)) + 1))))
cons_amp_gr <- lapply(amp_tracks, function(amp_gr) unlist(liftOver(amp_gr, hs1ToHg38_chain)) %Q% (AmpFreq > 0.8))

# DMRs in hg38. diff.Methy is the first DSS group (LUSC) minus the second (BRCA), so a positive value means lower methylation in BRCA
dmrs <- as.data.table(unlist(liftOver(gr.chr(dt2gr(fread(file.path(IN_DIR, settings$dmr_file)))), hs1ToHg38_chain)))
dmrs_by_group <- list(NSCLC = dt2gr(dmrs[diff.Methy < -settings$min_diff_methyl]),
                      BRCA = dt2gr(dmrs[diff.Methy > settings$min_diff_methyl]))
dmrs_cons <- lapply(names(dmrs_by_group), function(xGroup) dmrs_by_group[[xGroup]] %&% cons_amp_gr[[xGroup]])
names(dmrs_cons) <- names(dmrs_by_group)

# Peak coordinates are stored in peakID as chr:start:end
stats <- fread(file.path(IN_DIR, settings$stats_file))
stats[, c('seqnames', 'start', 'end') := tstrsplit(peakID, ':', type.convert = TRUE)]

overlaps <- rbindlist(lapply(names(tumor_types), function(xTT) {
  dmrs_tt <- dmrs_cons[[tumor_types[[xTT]]]]
  peaks_tt <- stats[tstdttype == xTT & OR > 1]
  hits <- gr.findoverlaps(dmrs_tt, dt2gr(peaks_tt), ignore.strand = TRUE)
  if (length(hits) == 0) return(NULL)
  dmr_hits <- gr2dt(dmrs_tt[hits$query.id])[, .(dmr_chr = seqnames, dmr_start = start, dmr_end = end, diff.Methy, areaStat)]
  dmr_hits[, dmr_id := paste0(dmr_chr, ':', dmr_start, ':', dmr_end)]
  cbind(tumor_type = xTT, amp_group = tumor_types[[xTT]], dmr_hits,
        peaks_tt[hits$subject.id, .(peakID, OR, p_val, q_val, tstdttype_count, compttype_count)])
}))
write.csv(overlaps, file.path(OUT_DIR, paste0(xG, '_peak_dmr_overlaps.csv')), row.names = FALSE)

# Number of DMRs that overlap a peak passing each cutoff, always with OR > 1
overlap_summary <- rbindlist(lapply(names(tumor_types), function(xTT) {
  xGroup <- tumor_types[[xTT]]
  overlaps_tt <- overlaps[tumor_type == xTT]
  data.table(tumor_type = xTT, amp_group = xGroup, cons_amp_bp = cons_amp_bp[[xGroup]],
             dmrs_in_cons_amp = length(dmrs_cons[[xGroup]]),
             q_below_0.05 = uniqueN(overlaps_tt[q_val < 0.05, dmr_id]),
             q_below_gene_cutoff = uniqueN(overlaps_tt[q_val < settings$q_cutoff, dmr_id]),
             p_below_0.05 = uniqueN(overlaps_tt[p_val < 0.05, dmr_id]),
             at_least_4_samples = uniqueN(overlaps_tt[tstdttype_count > 3, dmr_id]),
             any_peak = uniqueN(overlaps_tt$dmr_id))
}))
overlap_summary[, gene_q_cutoff := settings$q_cutoff]
write.csv(overlap_summary, file.path(OUT_DIR, paste0(xG, '_peak_dmr_overlap_summary.csv')), row.names = FALSE)

# NSCLC DMRs in regions amplified more often in NSCLC than in BRCA
if (!is.null(settings$preferential_amp_file)) {
  pref_amp_gr <- unlist(liftOver(readRDS(file.path(IN_DIR, settings$preferential_amp_file)), hs1ToHg38_chain))
  pref_amp_gr <- pref_amp_gr %Q% (NSCLC_AmpFreq > BRCA_AmpFreq & qval < 0.1)
  write.csv(gr2dt(dmrs_by_group$NSCLC %&% pref_amp_gr), file.path(OUT_DIR, paste0(xG, '_NSCLC_dmrs_in_preferential_amp.csv')), row.names = FALSE)
}
