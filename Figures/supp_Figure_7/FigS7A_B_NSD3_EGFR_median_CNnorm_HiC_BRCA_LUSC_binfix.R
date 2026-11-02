## Fig S7A/B: median CN-normalized HiC contact maps of the NSD3 and the EGFR HLAMP tumors, BRCA and LUSC on separate pages, with the gencode
## track; the gene is the argument:
## - NSD3, chr8:37,900,000-39,400,000 (hs1): nf-core HiC tumors (4 BRCA: C12, C132, C17, Cb5; 4 LUSC: 1IDEVD, 1E2AVA, NODTWC (map
##   NODWTC), K0HCJ3; the linkprep tumor C168 not included)
## - EGFR, chr7:54,500,000-56,000,000 (hs1): the linkprep tumors (.hic read with straw, raw counts, hs1 5kb tiles of the
##   nf-core run as reference; 4 BRCA: C109, C15, C168, C208; 4 LUSC: 1HGXBV, 1WNRWL, AJZ5JO, XO7JTV)
## - per tumor: raw 5kb counts of the window, normalized to the CN effect of the tumor's 20Mb background model (200 x 200 cells,
##   C17 fitted in the CCND1 set): count_ij / (CN_i*CN_j)^gamma(d_ij, CN_i*CN_j), CN_i = HiC coverage /
##   median QC tile of the model; self contacts take the CN effect at one tile step; tile pairs without counts 0; tiles without a CN of the model NA
## - sequencing depth: each tumor's normalized counts divided by its E(d) of the model (geometric mean over the distances of the window, 5kb
##   steps) relative to the median of that over the eight tumors of the gene, i.e. counts at the median depth
## - median per tile pair over the tumors with a CN for both tiles, separately for the BRCA (page 1) and the LUSC tumors (page 2); one HiC
##   colour scale for both pages: 0 to 4x the median v4C of the gene TSS (hg38 TSS lifted to hs1) over the lifted normal lung TAD (part in the
##   window) in the median map of all eight tumors
## - outputs per gene: _median_matrix.tsv (tile pairs: median and number of tumors per group: BRCA, LUSC, all), _depth_scale.tsv, .pdf (median
##   contact map and gencode track, BRCA and LUSC page)
## - the bin index of the dumped txt map is detected per map (1-based in the nf-core maps, 0-based in the maps dumped with cooler: 1IDEVD, K0HCJ3, NODWTC); CN models of 1IDEVD, NODTWC and K0HCJ3 refitted with this bin index
## run: FigS7A_B_job-script_median_CNnorm_HiC.sh (array: 0 NSD3, 1 EGFR); args: gene
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(GxG)
library(gGnome)
library(rtracklayer)
library(parallel)

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop('usage: Rscript <script> gene (NSD3 or EGFR)', call. = FALSE)
xG <- args[1]
outdir <- paste0(OUT_DIR, '/')
binsize = 5e3
clim_factor = 4
model <- 'logistic'
# tumors per gene; model_set: folder of the tumor's background model; map_id: stem of the raw map file
egfr_tumors <- data.table(sample = c('C109', 'C15', 'C168', 'C208', '1HGXBV', '1WNRWL', 'AJZ5JO', 'XO7JTV'), tissue = rep(c('breast', 'lung'), each = 4),
    model_set = 'EGFR')
egfr_tumors[, c('map_id', 'hic_dir') := .(sample, paste0(IN_DIR, '/hic/linkprep/', sample, '_1.hic'))]
cfg <- list(
    NSD3 = list(window = 'chr8:37900000-39400000', tumors = data.table(sample = c('C12', 'C132', 'C17', 'Cb5', '1IDEVD', '1E2AVA', 'NODTWC', 'K0HCJ3'),
        map_id = c('C12', 'C132', 'C17', 'Cb5', '1IDEVD', '1E2AVA', 'NODWTC', 'K0HCJ3'), tissue = rep(c('breast', 'lung'), each = 4),
        hic_dir = file.path(IN_DIR, 'hic/nfcore_results'), model_set = c('NSD3', 'NSD3', 'CCND1', 'NSD3', 'NSD3', 'NSD3', 'NSD3', 'NSD3'))),
    EGFR = list(window = 'chr7:54500000-56000000', tumors = egfr_tumors))
if (!xG %in% names(cfg)) stop(paste('gene must be one of', paste(names(cfg), collapse = ', ')))
window_oi <- GRanges(cfg[[xG]]$window)
xchr <- as.character(seqnames(window_oi))
tumors <- cfg[[xG]]$tumors
# background models refitted with the bin index fix for the tumors with 0-based maps
binfix_tumors <- c('1IDEVD', 'NODTWC', 'K0HCJ3')
model_file <- function(s, set) paste0(IN_DIR, '/', if (s %in% binfix_tumors) 'CNmodels_binfix/' else 'CNmodels/',
    set, '/20261023_', s, '_chr7_8_11_20Mb_HiCcov_logistic_df8_x200_d200.rds')
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
out_prefix <- paste0(outdir, 'FigS7A_B_', xG, '_median_CNnorm_HiC_BRCA_vs_LUSC_', xchr, '_', start(window_oi) / 1e6, '_', end(window_oi) / 1e6, 'Mb')
# tumor groups of the pages (tissue column of tumors)
groups <- c(BRCA = 'breast', LUSC = 'lung')
print(paste0(xG, ': window ', gr.string(window_oi), '; tumors ', paste(tumors$sample, collapse = ', ')))

# one tumor: CN-normalized counts of all tile pairs of the window (upper triangle incl. the diagonal) and the depth of its model
norm_window <- function(r) {
    s <- tumors$sample[r]
    hic_dir <- tumors$hic_dir[r]
    is_hic <- grepl('\\.hic$', hic_dir)
    # 5kb tiles of xchr from the bins of the raw contact map (genome-wide bin index, 1- or 0-based in the dumped txt, map_id_off; .hic: hs1
    # 5kb bins of the nf-core run); window tiles by midpoint
    bins_file <- if (is_hic) file.path(IN_DIR, 'hic/nfcore_results/contact_maps/bins/cooler_bins_5000.bed') else paste0(hic_dir, '/contact_maps/bins/cooler_bins_5000.bed')
    cbins <- fread(bins_file, header = FALSE, select = 1:3, col.names = c('chr', 'start', 'end'))
    idx <- which(cbins$chr == xchr)
    tiles <- GRanges(xchr, IRanges(cbins$start[idx] + 1, cbins$end[idx]))
    stopifnot(all(ceiling(end(tiles) / binsize) == seq_along(tiles)))
    win_k <- which(gr.mid(tiles) %^% window_oi)
    stopifnot(length(win_k) > 0, all(diff(win_k) == 1))
    wt <- tiles[win_k]
    if (is_hic) {
        # .hic: raw counts over the window tiles with straw (norm NONE); straw bins carry 0-based starts, tile index = start %/% binsize + 1
        mm <- straw(hic_dir, gr = GRanges(xchr, IRanges(start(wt)[1], end(wt)[length(wt)])), res = binsize, norm = 'NONE')
        t_map <- as.integer(start(mm$gr) %/% binsize + 1)
        d <- mm$dat
        pr <- data.table(i = pmin(t_map[d$i], t_map[d$j]), j = pmax(t_map[d$i], t_map[d$j]), counts = as.integer(round(d$value)))
        pr <- pr[i >= win_k[1] & j <= win_k[length(win_k)]]
        pr[, c('i', 'j') := .(i - win_k[1] + 1L, j - win_k[1] + 1L)]
    } else {
        # map ids of the window (map id = cooler bin row - id_off)
        map_file <- paste0(hic_dir, '/contact_maps/txt/', tumors$map_id[r], '.5000_raw.txt')
        id_off <- map_id_off(map_file, nrow(cbins))
        glo <- idx[win_k[1]] - id_off
        ghi <- idx[win_k[length(win_k)]] - id_off
        pr <- fread(cmd = sprintf("awk '$1 >= %d && $1 <= %d && $2 >= %d && $2 <= %d { print $1, $2, $3 }' %s",
            glo, ghi, glo, ghi, map_file), col.names = c('i', 'j', 'counts'))
        pr[, c('i', 'j') := .(i - glo + 1L, j - glo + 1L)]
    }
    stopifnot(pr[, all(i <= j)], !anyDuplicated(pr[, .(i, j)]))
    # CN effect of the tumor's background model, CN per tile matched by tile midpoint
    m <- readRDS(model_file(s, tumors$model_set[r]))
    th <- m$fits[[model]]$par[m$k + seq_len(m$forms[[model]]$n)]
    cn_effect <- function(ld, lx) exp(m$forms[[model]]$eta(th, m$design(ld, lx))$eta)
    cn_t <- (m$hic_cov / m$hiccov_median)[gr.match(gr.mid(wt), m$tile_gr)]
    valid <- is.finite(cn_t) & cn_t > 0
    n <- length(wt)
    grid <- CJ(i = seq_len(n), j = seq_len(n))[i <= j]
    grid[, counts := 0]
    grid[pr, counts := as.numeric(i.counts), on = .(i, j)]
    ok <- valid[grid$i] & valid[grid$j]
    ld <- log10(pmax(grid$j[ok] - grid$i[ok], 1) * binsize)
    lx <- log10(cn_t[grid$i[ok]] * cn_t[grid$j[ok]])
    grid[, value := NA_real_]
    grid[ok, value := counts / cn_effect(ld, lx)]
    # depth: E(d) of the model (CN 1 tile pairs), geometric mean over the distance steps of the window
    depth <- exp(mean(m$logE_by_dist[seq_len(n - 1)]))
    print(paste0(s, ': ', n, ' window tiles (', sum(valid), ' with a CN of the model), ', nrow(pr), ' nonzero tile pairs, ', sum(pr$counts), ' counts; CN of the model ',
        paste(round(range(cn_t[valid]), 2), collapse = '-'), '; depth (E(d), geometric mean over 5kb-', (n - 1) * binsize / 1e6, 'Mb) ', signif(depth, 4)))
    list(tiles = wt, grid = grid, depth = depth, n_valid = sum(valid), counts = sum(pr$counts))
}
res <- mclapply(seq_len(nrow(tumors)), function(r) try(norm_window(r)), mc.cores = nrow(tumors))
failed <- tumors$sample[sapply(res, function(x) inherits(x, 'try-error'))]
if (length(failed) > 0) stop(paste('failed:', paste(failed, collapse = ', '), '-', paste(sapply(res[sapply(res, inherits, 'try-error')], as.character), collapse = '; ')))
wt <- res[[1]]$tiles
stopifnot(all(sapply(res, function(x) identical(start(x$tiles), start(wt)) && identical(end(x$tiles), end(wt)))))

depth_dt <- data.table(tumors[, .(sample, tissue, model_set)], window_counts = sapply(res, `[[`, 'counts'), window_tiles_with_CN = sapply(res, `[[`, 'n_valid'),
    depth_Ed = sapply(res, `[[`, 'depth'))
depth_dt[, depth_scale := depth_Ed / median(depth_Ed)]
fwrite(depth_dt, paste0(out_prefix, '_depth_scale.tsv'), sep = '\t')
print(depth_dt)

# median over the tumors of a group per tile pair (depth-scaled normalized counts; tumors without a CN for a tile left out): BRCA, LUSC and
# all eight (colour scale)
all_dt <- rbindlist(lapply(seq_along(res), function(r) res[[r]]$grid[, .(sample = tumors$sample[r], tissue = tumors$tissue[r], i, j,
    value = value / depth_dt$depth_scale[r])]))
group_median <- function(dt) dt[, .(median_value = if (any(!is.na(value))) median(value, na.rm = TRUE) else NA_real_, n_tumors = sum(!is.na(value))), keyby = .(i, j)]
med_list <- c(lapply(groups, function(tis) group_median(all_dt[tissue == tis])), list(all = group_median(all_dt)))
med_dt <- rbindlist(med_list, idcol = 'group')
fwrite(data.table(med_dt[, .(group)], seqnames = xchr, start_i = start(wt)[med_dt$i], end_i = end(wt)[med_dt$i], start_j = start(wt)[med_dt$j],
    end_j = end(wt)[med_dt$j], med_dt[, !'group']), paste0(out_prefix, '_median_matrix.tsv'), sep = '\t')
for (g in names(med_list)) print(paste0(g, ': tile pairs ', nrow(med_list[[g]]), '; median from ', max(med_list[[g]]$n_tumors), ' tumors: ',
    med_list[[g]][, sum(n_tumors == max(n_tumors))], '; without any tumor: ', med_list[[g]][, sum(n_tumors == 0)]))
gmats <- lapply(med_list, function(md) gM(gr = wt, dat = md[!is.na(median_value), .(i, j, value = median_value)], fill = 0))

# colour scale: 4x the median v4C of the gene TSS over the normal lung TAD lifted to hs1
hg38ToHs1_chain <- import.chain(file.path(IN_DIR, 'hg38ToHs1.over.chain'))
CGC_gr_hg38 <- readRDS('../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
Lung_norm_16f_TAD_1 <- fread(file.path(IN_DIR, 'Lung_norm_16f_TAD_1_ENCFF581SYY.bedpe'))
colnames(Lung_norm_16f_TAD_1)[1:3] = c('chr', 'start', 'end')
Lung_TAD_gr <- gr.chr(dt2gr(Lung_norm_16f_TAD_1[, 1:3]))
TAD_hg38 <- reduce(Lung_TAD_gr[Lung_TAD_gr %^% gr.chr(CGC_gr_hg38 %Q% (Gene_Symbol == xG))], ignore.strand = TRUE)
TAD_oi <- unlist(range(liftOver(TAD_hg38, hg38ToHs1_chain)))
TAD_oi <- TAD_oi[seqnames(TAD_oi) == xchr]
stopifnot(length(TAD_oi) == 1)
xG_TSS <- fread(paste0(IN_DIR, '/', xG, '_TSS_hg38.csv'))
colnames(xG_TSS)[4] = 'seqnames'
colnames(xG_TSS)[8] = 'start'
xG_TSS[, end := start + 1]
xG_TSSgr <- unlist(liftOver(gr.chr(dt2gr(xG_TSS)), hg38ToHs1_chain))[1]
stopifnot(xG_TSSgr %^% window_oi)
print(paste0('TSS of ', xG, ' (hs1): ', gr.string(xG_TSSgr), '; lung TAD (hs1): ', gr.string(TAD_oi), ', ', round(100 * sum(width(GenomicRanges::intersect(TAD_oi, window_oi))) / width(TAD_oi)),
    '% in the window'))
clim_of <- function(gmat) clim_factor * median({v4C_gr <- gmat[gmat$gr %^% xG_TSSgr, drop = T]; v4C_gr$value[gr.mid(v4C_gr) %^% TAD_oi]}, na.rm = T)
clims <- sapply(gmats, clim_of)
print(paste0(clim_factor, 'x median ', xG, ' v4C over the TAD per median map: ', paste(names(clims), signif(clims, 4), collapse = ', ')))
# one colour scale for the BRCA and LUSC pages: that of the median map of all eight tumors
clim_max <- clims[['all']]
stopifnot(clim_max > 0)

gencode_T2T_gt <- readRDS(file.path(IN_DIR, 'final_T2T_gencode_gt_USCS.rds'))
pdf(file = paste0(out_prefix, '.pdf'), width = 8, height = 10)
for (g in names(groups)) {
    smp <- tumors[tissue == groups[[g]], sample]
    plot(c(gencode_T2T_gt, gmats[[g]]$gtrack(name = paste0('median counts/(CN_i*CN_j)^gamma, ', length(smp), ' ', g, ' ', xG, ' HLAMP tumors (',
        paste(smp, collapse = ', '), ')'), height = 90, clim = c(0, clim_max))), window_oi, cex.label = 0.1)
}
dev.off()
print('Done')
