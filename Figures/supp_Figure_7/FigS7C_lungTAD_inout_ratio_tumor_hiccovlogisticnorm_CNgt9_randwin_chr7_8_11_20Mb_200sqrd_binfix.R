## TAD ratio analysis for the tumors with 0-based contact maps: the bin index of the dumped txt map is detected per map (1-based in the nf-core maps, 0-based in the maps dumped with cooler: 1IDEVD, K0HCJ3, NODWTC);
## CN models refitted with this bin index
## TAD ratio analysis with the 200 x 200 cell background models (joint chr7/8/11 logistic HiC coverage CN model, pairs up to 20Mb, df 8;
## CN model file as 7th argument), partners and random window buffer up to 20Mb
## Within- vs outside-TAD contacts of the normal lung TAD overlapping a gene (CCND1 or NSD3) in one HLAMP tumor (hs1), normalized by the CN
## effect of the logistic HiC coverage CN model fitted jointly on chr7, chr8 and chr11 up to 20Mb, restricted to TAD tiles with JaBbA CN > 9,
## plus two sets of 200 random windows of the tested TAD's size across chr7, chr8 and chr11 as background; partners up to 20Mb apart
## - CN model: fitted by the Fig. S4 script (df 8; raw 5kb contact map, intrachromosomal
##   chr7/chr8/chr11 pairs up to 20Mb): observed_ij = E(d_ij) * (CN_i*CN_j)^gamma(d, CN_i*CN_j), CN_i = HiC coverage of tile i (row sum of the
##   raw map over its chromosome) / median QC tile of the three chromosomes; one model for TAD and random windows
## - CN_i for the normalization: the HiC coverage and its median stored with the CN model, matched to the tiles by tile midpoint
## - normal lung TAD overlapping the gene (hg38, Lung_norm_16f_TAD_1_ENCFF581SYY) lifted to hs1 with the UCSC hg38ToHs1 chain
## - tile counts from the raw 5kb contact map of the nf-core/hic run (<hic_dir>/contact_maps/txt/<map_id>.5000_raw.txt, intrachromosomal pairs
##   of chr7, chr8, chr11), normalized to the model CN effect: norm_count_ij = count_ij / (CN_i*CN_j)^gamma(d_ij, CN_i*CN_j); the distance decay
##   stays in; self contacts take the CN effect at the shortest fitted distance (one tile step, 5kb)
## - HiC coverage of a tile (row sum of the raw map over its chromosome, self-contact once, as in the CN model) for the tile QC and the
##   track: tiles with HiC coverage <= min_hiccov_frac * median of the chromosome (or without a CN of the model) drop out of all sums (as
##   anchors and as partners), for the TAD and the random windows
## - TAD tiles with JaBbA CN <= 9 (or without a JaBbA CN) are dropped from the TAD: no ratio for them and they count neither as within- nor
##   as outside-TAD partners
## - per kept 5kb TAD tile: ratio = sum(norm. contacts with kept TAD tiles) / sum(norm. contacts with tiles outside the TAD); partners on the
##   same chromosome up to Dmax = 20Mb (the range of the CN model fit)
## - random windows: as many tiles as kept TAD tiles; same per-tile ratios with the window in place of the TAD (no CN filter on the tiles);
##   candidate window starts on chr7, chr8 and chr11 where the window +-Dmax lies on the chromosome, the window itself has a JaBbA CN for
##   every tile and >= min_good_frac of its tiles pass the tile QC, and
##   - set 'randomwin': no tile of the gene's TAD and no tile with JaBbA CN > 9 within window +- Dmax
##   - set 'randomwin_noampfilter': amplified tiles allowed; the window itself holds no tile of the gene's TAD
##   candidates of the three chromosomes pooled and drawn in random order (seed 20261018 for each set), a window is kept if it does not
##   overlap a kept window of its chromosome, up to n_win = 200
## - diagnostic plot with the HiC coverage track, the lifted TAD and kept TAD tiles: pages 1/2 window_oi, pages 3/4 window_wide, pages 5/6
##   window_wide2 (CCND1: chr11:69-70Mb, 68.5-70.5Mb and 68-71Mb; other genes: gene midpoint +-0.5Mb, +-1Mb and +-1.5Mb), each first the
##   CN-model-normalized matrix, its HiC coverage (row sum of counts / CN effect over the partners up to Dmax) and v4C of the gene TSS, then the
##   raw count matrix (all tiles) and raw v4C; HiC colour scale and v4C y axis 0 to 4x the median v4C over the lifted TAD
## hic_dir: results folder of an nf-core/hic run or a .hic file (read with straw, raw counts, hs1 5kb tiles of the nf-core run as reference);
## methylation calls from either bedmethyl folder; no NA window row for a chromosome without candidates
## run: FigS7C_job-script_binfix_TADratio_logistic.sh; args: sampleid ampChr tissue hic_dir gene outdir cn_model_file [map_id] (map_id: stem of the raw map file if it differs from sampleid)
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')

library(GxG)
library(gGnome)
library(rtracklayer)

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop('usage: Rscript <script> sampleid ampChr tissue hic_dir gene outdir cn_model_dir [map_id]', call. = FALSE)
sampleid <- args[1]
ampChr <- args[2]
tissue <- args[3]
hic_dir <- args[4]
xG <- args[5]
outdir <- sub('/*$', '/', args[6])
cn_model_dir <- sub('/*$', '/', args[7])
map_id <- if (length(args) >= 8) args[8] else sampleid
model <- 'logistic'
print(paste0('Starting sample: ', sampleid, ' (', hic_dir, ', map ', map_id, '); gene ', xG, '; CN model: ', model, ' (chr7/chr8/chr11 joint fit, 20Mb)'))
# jabba results sit under the sample sheet's amplicon (C17 is listed as chr8/FGFR1)
amp <- c(chr8 = 'FGFR1', chr11 = 'CCND1')[[ampChr]]

norm_tag <- paste0('hiccov_', model, 'norm')
cn_model_file <- if (grepl('\\.rds$', args[7])) args[7] else paste0(cn_model_dir, '20261023_', sampleid, '_chr7_8_11_20Mb_HiCcov_logistic_df8.rds')
chrs <- c('chr7', 'chr8', 'chr11')
binsize = 5e3
Dmax = 2e7
cn_cutoff = 9
n_win = 200
min_hiccov_frac = 0.1
min_good_frac = 0.9
chunk_size = 2e6
clim_factor = 4
win_sets <- c('randomwin', 'randomwin_noampfilter')

hg38ToHs1_chain <- import.chain(file.path(IN_DIR, 'hg38ToHs1.over.chain'))
gencode_T2T_gt <- readRDS(file.path(IN_DIR, 'final_T2T_gencode_gt_USCS.rds'))
CGC_gr_T2T <- readRDS('../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')
CGC_gr_hg38 <- readRDS('../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
Oncgr <- gr.chr(CGC_gr_T2T[mcols(CGC_gr_T2T)$Gene_Symbol == xG])
xchr <- as.character(seqnames(Oncgr))[1]
stopifnot(xchr %in% chrs)
# plot windows and the TSS of the gene for the v4C (hs1)
if (xG == 'CCND1') {
    window_oi = GRanges("chr11:69000000-70000000")
    window_wide = GRanges("chr11:68500000-70500000")
    window_wide2 = GRanges("chr11:68000000-71000000")
    xG_TSSgr <- readRDS(file.path(IN_DIR, 'CCND1_TSSgr_T2T.rds'))
} else {
    xG_mid <- round((min(start(Oncgr)) + max(end(Oncgr))) / 2)
    window_oi <- GRanges(xchr, IRanges(xG_mid - 5e5, xG_mid + 5e5))
    window_wide <- GRanges(xchr, IRanges(xG_mid - 1e6, xG_mid + 1e6))
    window_wide2 <- GRanges(xchr, IRanges(xG_mid - 1.5e6, xG_mid + 1.5e6))
    xG_TSS <- fread(paste0(IN_DIR, '/', xG, '_TSS_hg38.csv'))
    colnames(xG_TSS)[4] = 'seqnames'
    colnames(xG_TSS)[8] = 'start'
    xG_TSS[, end := start + 1]
    xG_TSSgr <- unlist(liftOver(gr.chr(dt2gr(xG_TSS)), hg38ToHs1_chain))
}
xG_TSSgr <- xG_TSSgr[1]
print(paste('TSS of', xG, '(hs1):', as.character(seqnames(xG_TSSgr)), start(xG_TSSgr), '; plot windows', paste(gr.string(c(window_oi, window_wide, window_wide2)), collapse = ', ')))

# normal lung TADs (hg38) lifted to hs1; a lifted TAD spans the range of its lifted pieces on xchr
lift_hs1 <- function(gr) {
    gr_lift <- unlist(range(liftOver(gr, hg38ToHs1_chain)))
    gr_lift[seqnames(gr_lift) == xchr]
}
Lung_norm_16f_TAD_1 <- fread(file.path(IN_DIR, 'Lung_norm_16f_TAD_1_ENCFF581SYY.bedpe'))
colnames(Lung_norm_16f_TAD_1)[1:3] = c('chr', 'start', 'end')
Lung_TAD_gr <- gr.chr(dt2gr(Lung_norm_16f_TAD_1[, 1:3]))
TAD_hg38 <- reduce(Lung_TAD_gr[Lung_TAD_gr %^% gr.chr(CGC_gr_hg38 %Q% (Gene_Symbol == xG))], ignore.strand = TRUE)
TAD_oi <- lift_hs1(TAD_hg38)
if (length(TAD_oi) != 1) stop(paste('lifted TAD does not map to a single range on', xchr))
TAD_str <- paste0(xchr, ':', start(TAD_oi), '-', end(TAD_oi))
print(paste('TAD used (hs1):', TAD_str, '; width hg38', width(TAD_hg38), ', hs1', width(TAD_oi)))
Lung_norm_16f_TAD_1_hs1_gt <- gTrack(data = lift_hs1(Lung_TAD_gr[Lung_TAD_gr %^% gr.chr(CGC_gr_hg38)]), name = 'Lung_norm_16f_TAD_1_hs1', col = 'grey', height = 0.5)

# contact map: raw 5kb map of the nf-core/hic run (hic_dir = results dir) or a .hic file (hic_dir = path to <id>_1.hic, linkprep runs; read
# with straw, norm NONE = raw counts); 5kb tiles per chromosome from the bins of the raw contact map (genome-wide cooler bin index, 1- or
# 0-based in the dumped txt, map_id_off; for a .hic file the hs1 5kb bins of the nf-core run as reference); tile k of a chromosome ends at
# k*binsize (the last one at the chromosome end), so tile index = ceiling(end/binsize)
# genome-wide bin index of the dumped 5kb txt map: 1-based in the nf-core maps (last id = number of bins), 0-based in the maps dumped with
# cooler (id 0 present: 1IDEVD, K0HCJ3, NODWTC); returns id_off = cooler bin row - map id
map_id_off <- function(map_file, n_bins) {
    ids <- function(cmd) as.numeric(strsplit(trimws(system(paste(cmd, map_file), intern = TRUE)), '\\s+')[[1]][1:2])
    first <- min(ids('head -n 1'))
    last <- max(ids('tail -n 1'))
    off <- if (first == 0 && last < n_bins) 1L else if (first >= 1 && last == n_bins) 0L else
        stop(paste0('bin index of ', map_file, ' undetermined (first id ', first, ', last id ', last, ', ', n_bins, ' bins)'))
    print(paste0(basename(map_file), ': ', 1L - off, '-based bin index (first id ', first, ', last id ', last, ', ', n_bins, ' cooler bins)'))
    off
}
is_hic <- grepl('\\.hic$', hic_dir)
bins_file <- if (is_hic) file.path(IN_DIR, 'hic/nfcore_results/contact_maps/bins/cooler_bins_5000.bed') else paste0(hic_dir, '/contact_maps/bins/cooler_bins_5000.bed')
cbins <- fread(bins_file, header = FALSE, select = 1:3, col.names = c('chr', 'start', 'end'))
n_bins <- nrow(cbins)
chr_dt <- rbindlist(lapply(chrs, function(ch) {
    idx <- which(cbins$chr == ch)
    stopifnot(length(idx) > 0, all(diff(idx) == 1))
    data.table(chr = ch, lo = idx[1], hi = idx[length(idx)], len = cbins$end[idx[length(idx)]])
}))
tiles <- lapply(setNames(seq_along(chrs), chrs), function(q) {
    gr <- GRanges(chrs[q], IRanges(cbins$start[chr_dt$lo[q]:chr_dt$hi[q]] + 1, cbins$end[chr_dt$lo[q]:chr_dt$hi[q]]))
    stopifnot(all(ceiling(end(gr) / binsize) == seq_along(gr)))
    gr
})
rm(cbins)

# joint HiC coverage CN model: CN per tile (HiC coverage of the fit / its median QC tile, matched by tile midpoint) and the CN effect
cn_model <- readRDS(cn_model_file)
th <- cn_model$fits[[model]]$par[cn_model$k + seq_len(cn_model$forms[[model]]$n)]
cn_effect <- function(ld, lx) exp(cn_model$forms[[model]]$eta(th, cn_model$design(ld, lx))$eta)
print(paste0('CN model ', cn_model_file, ' (chromosomes ', paste(cn_model$chrs, collapse = ', '), ', distance spline df ', cn_model$ns_df, ', pairs up to ',
    cn_model$Dmax, '): ', cn_model$forms[[model]]$label))
print(data.table(parameter = cn_model$forms[[model]]$par_names, value = th))
print(cn_model$comparison[, .(form, model, n_par, loglik, delta_AIC, delta_BIC, nlminb)])
stopifnot(Dmax <= cn_model$Dmax)

# JaBbA CN per tile, matched by tile midpoint (TAD tile filter and random windows only)
jabba_res_gg = gG(jabba=paste0(IN_DIR, '/jabba/',amp,'/',tissue,'/',sampleid,'/jabba.gg.rds'))
jabba_nodes_gr <- gr.chr(gr.stripstrand(jabba_res_gg$nodes$gr))

# raw contact map, intrachromosomal pairs of the three chromosomes in one pass (upper triangle: i in the chromosome and j <= its last bin);
# local tile index per chromosome
print('reading in the raw contact map')
if (is_hic) {
    # .hic: whole chromosome at 5kb, raw counts; straw bins carry 0-based starts, tile index = start %/% binsize + 1
    map_dt <- rbindlist(lapply(seq_len(nrow(chr_dt)), function(q) {
        m <- straw(hic_dir, gr = GRanges(chr_dt$chr[q], IRanges(1, chr_dt$len[q])), res = binsize, norm = 'NONE')
        t_map <- as.integer(start(m$gr) %/% binsize + 1)
        stopifnot(max(t_map) <= chr_dt$hi[q] - chr_dt$lo[q] + 1)
        d <- m$dat
        data.table(i = pmin(t_map[d$i], t_map[d$j]), j = pmax(t_map[d$i], t_map[d$j]), counts = as.integer(round(d$value)), q = q)
    }))
} else {
    # map id = cooler bin row - id_off
    map_file <- paste0(hic_dir, '/contact_maps/txt/', map_id, '.5000_raw.txt')
    id_off <- map_id_off(map_file, n_bins)
    awk_prog <- paste(sprintf('($1 >= %d && $1 <= %d && $2 <= %d) { print $1, $2, $3, %d; next }', chr_dt$lo - id_off, chr_dt$hi - id_off, chr_dt$hi - id_off,
        seq_len(nrow(chr_dt))), collapse = ' ')
    map_dt <- fread(cmd = paste0("awk '", awk_prog, "' ", map_file), col.names = c('i', 'j', 'counts', 'q'))
    stopifnot(map_dt[, all(i <= j)])
    map_dt[, c('i', 'j') := .(i + id_off - chr_dt$lo[q] + 1L, j + id_off - chr_dt$lo[q] + 1L)]
}
print(map_dt[, .(nonzero_tile_pairs = .N, pairs_up_to_Dmax = sum(j - i <= Dmax / binsize), counts = sum(as.numeric(counts))), keyby = .(chr = chrs[q])])
wide_t <- which(tiles[[xchr]] %^% window_wide2)

# per chromosome: HiC coverage, CN of the model, tile QC, JaBbA CN and the tile pairs up to Dmax normalized to the model CN effect (raw
# counts of the widest plot window kept for xchr)
prep_chr <- function(ch) {
    xT_gr <- tiles[[ch]]
    n_t <- length(xT_gr)
    it <- map_dt[q == match(ch, chrs), .(i, j, counts)]
    # HiC coverage: sum of all counts of a tile on the chromosome (row sum of the symmetric matrix, self-contact once)
    total <- numeric(n_t)
    s1 <- it[, .(t = sum(counts)), by = i]
    total[s1$i] <- s1$t
    s2 <- it[i != j, .(t = sum(counts)), by = j]
    total[s2$j] <- total[s2$j] + s2$t
    it <- it[j - i <= Dmax / binsize]
    cn_t <- (cn_model$hic_cov / cn_model$hiccov_median)[gr.match(gr.mid(xT_gr), cn_model$tile_gr)]
    both_pos <- is.finite(cn_t) & cn_t > 0 & total > 0
    valid <- total > 0 & is.finite(cn_t) & cn_t > 0
    raw_wide <- if (ch == xchr) it[i %in% wide_t & j %in% wide_t, .(i, j, value = counts)] else NULL
    it <- it[valid[i] & valid[j]]
    it[, c('ld', 'lx') := .(log10(pmax(abs(j - i), 1) * binsize), log10(cn_t[i] * cn_t[j]))]
    eff <- numeric(nrow(it))
    for (s in seq(1, nrow(it), by = chunk_size)) {
        ch_rows <- s:min(nrow(it), s + chunk_size - 1)
        eff[ch_rows] <- cn_effect(it$ld[ch_rows], it$lx[ch_rows])
    }
    it[, value := counts / eff]
    min_hiccov <- min_hiccov_frac * median(total[total > 0])
    good <- valid & total > min_hiccov
    tile_cn <- jabba_nodes_gr$cn[gr.match(gr.mid(xT_gr), jabba_nodes_gr)]
    print(paste0(ch, ': ', n_t, ' tiles, ', sum(total == 0), ' without HiC counts, ', sum(total > 0 & !valid), ' more without a CN of the model; ',
        sum(good), ' pass tile QC (HiC coverage > ', min_hiccov, '); CN of the model vs raw map HiC coverage, Spearman ',
        round(cor(cn_t[both_pos], total[both_pos], method = 'spearman'), 3), '; CN effect over the pairs (median / max) ', signif(median(eff), 3), ' / ',
        signif(max(eff), 3), '; tiles with JaBbA CN > ', cn_cutoff, ': ', sum(tile_cn > cn_cutoff, na.rm = TRUE), ', without JaBbA CN: ', sum(is.na(tile_cn))))
    list(xT_gr = xT_gr, n_t = n_t, total = total, cn_t = cn_t, good = good, tile_cn = tile_cn, pairs = it[, .(i, j, value)], raw_wide = raw_wide)
}
cd <- lapply(setNames(chrs, chrs), prep_chr)
rm(map_dt)
gc()

# lung TAD of the gene on xchr: tiles kept with JaBbA CN > cn_cutoff; excluded TAD tiles drop out as partners too
cx <- cd[[xchr]]
TAD_tiles <- which(gr.mid(cx$xT_gr) %^% TAD_oi)
keep_tile <- cx$good & !is.na(cx$tile_cn) & cx$tile_cn > cn_cutoff
TAD_kept <- TAD_tiles[keep_tile[TAD_tiles]]
in_TAD <- seq_len(cx$n_t) %in% TAD_kept
use_tile <- cx$good & !(seq_len(cx$n_t) %in% setdiff(TAD_tiles, TAD_kept))
print(paste(length(TAD_tiles), 'tiles in TAD,', length(TAD_kept), 'kept with JaBbA CN >', cn_cutoff, 'and passing tile QC'))
print('JaBbA CN of TAD tiles:')
print(table(cx$tile_cn[TAD_tiles], useNA = 'ifany'))
fwrite(data.table(sample = sampleid, seqnames = xchr, start = start(cx$xT_gr)[TAD_tiles], end = end(cx$xT_gr)[TAD_tiles],
    hic_cov = cx$total[TAD_tiles], cn_model = cx$cn_t[TAD_tiles], jabba_cn = cx$tile_cn[TAD_tiles], kept = keep_tile[TAD_tiles]),
    paste0(outdir, '20261024_', sampleid, '_', xG, '_lungTAD_hs1_tiles_jabbaCN.tsv'), sep = '\t')

# random windows of the tested (CN-filtered) TAD's size across the chromosomes (same seed for every tumor and set): candidate starts where
# the window +- Dmax lies on the chromosome, the window has a JaBbA CN for every tile and >= min_good_frac of its tiles pass tile QC, and
# (randomwin) no TAD tile and no tile with JaBbA CN > cn_cutoff within window +- Dmax, or (randomwin_noampfilter) no TAD tile in the window;
# candidates pooled over the chromosomes and drawn in random order; a window is kept if it does not overlap a kept window of its chromosome
win_n <- length(TAD_kept)
D_tiles <- Dmax / binsize
win_cand <- function(set) rbindlist(lapply(chrs, function(ch) {
    cc <- cd[[ch]]
    if (win_n == 0 || cc$n_t < win_n + 2 * D_tiles) return(NULL)
    TAD_ch <- seq_len(cc$n_t) %in% (if (ch == xchr) TAD_tiles else integer(0))
    na_cs <- c(0, cumsum(is.na(cc$tile_cn)))
    good_cs <- c(0, cumsum(cc$good))
    s <- (D_tiles + 1):(cc$n_t - win_n + 1 - D_tiles)
    ok <- na_cs[s + win_n] - na_cs[s] == 0 & good_cs[s + win_n] - good_cs[s] >= min_good_frac * win_n
    if (set == 'randomwin') {
        amp_cs <- c(0, cumsum((!is.na(cc$tile_cn) & cc$tile_cn > cn_cutoff) | TAD_ch))
        ok <- ok & amp_cs[s + win_n + D_tiles] - amp_cs[s - D_tiles] == 0
    } else {
        TAD_cs <- c(0, cumsum(TAD_ch))
        ok <- ok & TAD_cs[s + win_n] - TAD_cs[s] == 0
    }
    # no candidate on this chromosome: no row (data.table() would recycle chr and fill s with NA)
    if (!any(ok)) return(NULL)
    data.table(chr = ch, s = s[ok])
}))
draw_windows <- function(cand) {
    set.seed(20261018)
    win_dt <- data.table(chr = character(0), s = integer(0))
    if (nrow(cand) > 0) {
        cand <- cand[sample.int(nrow(cand))]
        chosen <- lapply(setNames(chrs, chrs), function(ch) integer(0))
        keep_r <- integer(0)
        for (r in seq_len(nrow(cand))) {
            ch <- cand$chr[r]
            if (all(abs(cand$s[r] - chosen[[ch]]) >= win_n)) {
                chosen[[ch]] <- c(chosen[[ch]], cand$s[r])
                keep_r <- c(keep_r, r)
            }
            if (length(keep_r) == n_win) break
        }
        win_dt <- cand[keep_r]
    }
    win_dt[, window := seq_len(.N)]
    win_dt
}
wins <- list()
for (set in win_sets) {
    cand <- win_cand(set)
    wins[[set]] <- draw_windows(cand)
    w <- wins[[set]]
    # amplified tiles per window (JaBbA CN > cn_cutoff)
    w[, n_amp_tiles := sapply(seq_len(.N), function(r) sum(cd[[chr[r]]]$tile_cn[s[r] + 0:(win_n - 1)] > cn_cutoff, na.rm = TRUE))]
    wins[[set]] <- w
    print(paste0(set, ': ', nrow(cand), ' candidate window starts (', paste(chrs, sapply(chrs, function(ch) sum(cand$chr == ch)), collapse = ', '), '), ', nrow(w),
        ' random windows of ', win_n, ' tiles (', paste(chrs, sapply(chrs, function(ch) sum(w$chr == ch)), collapse = ', '), '); ', sum(w$n_amp_tiles > 0),
        ' windows with a tile of JaBbA CN > ', cn_cutoff))
    if (nrow(w) < n_win) warning(paste(sampleid, set, ': only', nrow(w), 'random windows'))
}

add_ratios <- function(dt) {
    dt[, ratio := within / outside]
    dt[, ratio_incl_self := (within + self) / outside]
    # length normalized interaction count sums (self = 1 tile)
    dt[, within_lennorm := within / n_within]
    dt[, outside_lennorm := outside / n_outside]
    dt[, ratio_lennorm := within_lennorm / outside_lennorm]
    dt[, ratio_lennorm_incl_self := ((within + self) / (n_within + 1)) / outside_lennorm]
    dt
}
ratio_cols <- c('ratio', 'ratio_incl_self', 'ratio_lennorm', 'ratio_lennorm_incl_self')

# v4C of the gene TSS and HiC colour scale (0 to clim_factor x the median v4C over the lifted lung TAD) of a matrix
v4C_of <- function(gmat) gmat[gmat$gr %^% xG_TSSgr, drop = T]
clim_of <- function(v4C_gr, label) {
    clim_max <- clim_factor * median(v4C_gr$value[gr.mid(v4C_gr) %^% TAD_oi], na.rm = T)
    print(paste0('HiC colour scale max, ', label, ' (', clim_factor, 'x median ', xG, ' v4C over the TAD): ', clim_max))
    stopifnot(clim_max > 0)
    clim_max
}

# per chromosome: normalized matrix, TAD (xchr) and random window sums, diagnostic plots (xchr); one matrix in memory at a time
ratio_dt <- NULL
rand_list <- lapply(setNames(win_sets, win_sets), function(set) list())
for (ch in chrs) {
    cc <- cd[[ch]]
    w_ch <- lapply(wins, function(w) w[chr == ch])
    if (ch != xchr && all(sapply(w_ch, nrow) == 0)) next
    tile_idx <- seq_len(cc$n_t)
    gmat <- gM(gr = cc$xT_gr, dat = cc$pairs, fill = 0)
    # within/outside sums of tile k for the region in_region (logical per tile); partners: tiles in use, same chromosome, <= Dmax away
    tile_sums <- function(k, in_region, use) {
        # virtual 4C of tile k against every tile on the chromosome
        xT_map_4C <- gmat[tile_idx == k, drop = T]
        p <- ceiling(end(xT_map_4C) / binsize)
        p_ok <- use[p] & abs(p - k) * binsize <= Dmax
        partner <- use & abs(tile_idx - k) * binsize <= Dmax & tile_idx != k
        data.table(sample = sampleid, tissue = tissue, seqnames = ch, start = start(cc$xT_gr)[k], end = end(cc$xT_gr)[k],
            hic_cov = cc$total[k], cn_model = cc$cn_t[k], jabba_cn = cc$tile_cn[k],
            self = sum(xT_map_4C$value[p_ok & p == k], na.rm = T),
            within = sum(xT_map_4C$value[p_ok & p != k & in_region[p]], na.rm = T),
            outside = sum(xT_map_4C$value[p_ok & !in_region[p]], na.rm = T),
            n_within = sum(partner & in_region), n_outside = sum(partner & !in_region))
    }
    # the 4C row of a tile (kept TAD tile or window tile passing tile QC) must hold its contacts from both triangles of the tile-pair table
    check_t <- c(if (ch == xchr) TAD_kept, unlist(lapply(w_ch, function(w) unlist(lapply(w$s, function(s) s + 0:(win_n - 1))))))
    check_t <- check_t[cc$good[check_t]]
    if (length(check_t) > 0) {
        k <- check_t[1]
        xT_map_4C <- gmat[tile_idx == k, drop = T]
        sum_4C <- sum(xT_map_4C$value[ceiling(end(xT_map_4C) / binsize) != k], na.rm = T)
        sum_pairs <- cc$pairs[xor(i == k, j == k), sum(value)]
        print(paste(ch, ': check tile', k, ': off-diagonal sum from 4C', sum_4C, ', from tile pairs', sum_pairs))
        if (!isTRUE(all.equal(sum_4C, sum_pairs))) stop('virtual 4C does not match the tile-pair table')
    }

    if (ch == xchr) {
        if (length(TAD_kept) == 0) {
            warning(paste(sampleid, ': no TAD tiles with JaBbA CN >', cn_cutoff, '- no ratio table written'))
        } else {
            ratio_dt <- add_ratios(rbindlist(lapply(TAD_kept, tile_sums, in_region = in_TAD, use = use_tile)))
            fwrite(ratio_dt, paste0(outdir, '20261024_', sampleid, '_', xG, '_lungTAD_hs1_inout_ratio_', norm_tag, '_CNgt', cn_cutoff, '.tsv'), sep = '\t')
            print('TAD median:')
            print(ratio_dt[, c(.(n_tiles = .N), lapply(.SD, function(x) median(x[is.finite(x)]))), .SDcols = ratio_cols])
        }
        print('diagnostic plots')
        # methylation calls: HiCDCPlus inputs, else the bedmethyl rds folder used for the linkprep samples
        bedmeth_files <- paste0(IN_DIR, c('/bedmeth/', '/bedmethylRDS/'), sampleid, '.wf_mods.ungrouped.bedmethyl.rds')
        mods.ungrouped.bedmethyl <- readRDS(bedmeth_files[file.exists(bedmeth_files)][1])
        mods.ungrouped.bedmethyl_logmCe3_gt <- gTrack(dt2gr(mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3)
        rm(mods.ungrouped.bedmethyl)
        hic_cov_gr <- cc$xT_gr
        hic_cov_gr$hic_cov <- cc$total
        hic_cov_gt <- gTrack(data = hic_cov_gr, y.field = 'hic_cov', name = 'HiC_cov', bars = T, y0 = 0)
        # HiC coverage of the normalized matrix: row sum of counts / CN effect over the partners up to Dmax (self-contact once; tiles
        # without a CN of the model 0)
        cov_norm <- numeric(cc$n_t)
        s1 <- cc$pairs[, .(t = sum(value)), by = i]
        cov_norm[s1$i] <- s1$t
        s2 <- cc$pairs[i != j, .(t = sum(value)), by = j]
        cov_norm[s2$j] <- cov_norm[s2$j] + s2$t
        hic_cov_gr$hic_cov_norm <- cov_norm
        hic_cov_norm_gt <- gTrack(data = hic_cov_gr, y.field = 'hic_cov_norm', name = paste0('HiC_cov_CNnorm_', Dmax / 1e6, 'Mb'), bars = T, y0 = 0)
        TAD_gts <- Lung_norm_16f_TAD_1_hs1_gt
        if (length(TAD_kept) > 0) TAD_gts <- c(TAD_gts, gTrack(data = reduce(cc$xT_gr[TAD_kept]), name = paste0('TAD_CN>', cn_cutoff), col = 'black', height = 0.5))
        # v4C y axis capped at the HiC colour scale max
        v4C_xG_gr <- v4C_of(gmat)
        clim_norm <- clim_of(v4C_xG_gr, 'normalized')
        v4C_xG_gt <- gTrack(data = v4C_xG_gr, y.field = 'value', bars = T, name = 'counts/CNeffect', y0 = 0, y1 = clim_norm, col = 'red')
        # raw counts within the widest plot window
        gmat_raw <- gM(gr = cc$xT_gr, dat = cc$raw_wide, fill = 0)
        v4C_raw_gr <- v4C_of(gmat_raw)
        clim_raw <- clim_of(v4C_raw_gr, 'raw counts')
        v4C_raw_gt <- gTrack(data = v4C_raw_gr, y.field = 'value', bars = T, name = 'raw counts', y0 = 0, y1 = clim_raw, col = 'red')
        base_gts <- c(gencode_T2T_gt, hic_cov_gt, jabba_res_gg$gt, mods.ungrouped.bedmethyl_logmCe3_gt)
        pdf(file = paste0(outdir, '20261024_', sampleid, '_', xG, '_', norm_tag, '_matrix_v4C_lungTAD_hs1_CNgt', cn_cutoff, '.pdf'), paper = 'a4', height = 24)
        for (win in list(window_oi, window_wide, window_wide2)) {
            plot(c(base_gts, hic_cov_norm_gt, v4C_xG_gt, TAD_gts,
                gmat$gtrack(name = paste0('counts/(CN_i*CN_j)^gamma, ', model, ' HiCcov CN model (chr7/8/11)'), height = 90, clim = c(0, clim_norm))), win, cex.label = 0.1)
            plot(c(base_gts, v4C_raw_gt, TAD_gts,
                gmat_raw$gtrack(name = 'raw counts (raw 5kb contact map)', height = 90, clim = c(0, clim_raw))), win, cex.label = 0.1)
        }
        dev.off()
        rm(gmat_raw)
    }

    for (set in win_sets) {
        w <- w_ch[[set]]
        if (nrow(w) == 0) next
        rand_list[[set]][[ch]] <- rbindlist(lapply(seq_len(nrow(w)), function(r) {
            win_tiles <- w$s[r] + 0:(win_n - 1)
            # window tiles failing tile QC drop out like the excluded TAD tiles
            win_good <- win_tiles[cc$good[win_tiles]]
            win_str <- paste0(ch, ':', start(cc$xT_gr)[win_tiles[1]], '-', end(cc$xT_gr)[win_tiles[win_n]])
            data.table(window = w$window[r], window_range = win_str, window_n_amp_tiles = w$n_amp_tiles[r],
                rbindlist(lapply(win_good, tile_sums, in_region = tile_idx %in% win_good, use = cc$good)))
        }))
        print(paste(ch, ':', set, ':', nrow(w), 'random windows done'))
    }
    rm(gmat)
    gc()
}

for (set in win_sets) {
    if (length(rand_list[[set]]) == 0) next
    rand_dt <- add_ratios(rbindlist(rand_list[[set]]))
    setorder(rand_dt, window)
    fwrite(rand_dt, paste0(outdir, '20261024_', sampleid, '_', xG, '_', set, '_hs1_inout_ratio_', norm_tag, '.tsv'), sep = '\t')
    print(paste(set, ': random window medians:'))
    print(rand_dt[, c(.(n_tiles = .N, n_amp_tiles = window_n_amp_tiles[1]), lapply(.SD, function(x) median(x[is.finite(x)]))), by = .(window, window_range),
        .SDcols = ratio_cols])
}
print(paste0('Done: ', sampleid, ' (', xG, ', ', model, ')'))
