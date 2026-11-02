## Fig S5A: median CN-normalized HiC contact map of the CCND1 HLAMP tumors with nf-core HiC (3 BRCA: C169, C17, C3; 4 LUSC: H3D44S, KDAMIF,
## Y57U8R, F57J9F; the linkprep tumors C15 and 1WNRWL not included) in chr11:69,000,000-70,500,000 (hs1), with the gencode track
## - per tumor: raw 5kb counts of the window (raw contact map of the nf-core run), normalized to the CN effect of the tumor's 20Mb background
##   model (200 x 200 cells): count_ij / (CN_i*CN_j)^gamma(d_ij, CN_i*CN_j), CN_i = HiC coverage /
##   median QC tile of the model; self contacts take the CN effect at one tile step; tile pairs without counts 0; tiles without a CN of the model NA
## - sequencing depth: each tumor's normalized counts divided by its E(d) of the model (geometric mean over the distances of the window, 5kb
##   steps) relative to the median of that over the seven tumors, i.e. counts at the median depth
## - median per tile pair over the tumors with a CN for both tiles, separately for the 3 BRCA (page 1) and the 4 LUSC tumors (page 2); one HiC
##   colour scale for both pages: 0 to 4x the median v4C of the CCND1 TSS over the lifted normal lung TAD in the median map of all seven
##   tumors
## - outputs: _median_matrix.tsv (tile pairs: median and number of tumors per group: BRCA, LUSC, all), _depth_scale.tsv, .pdf (median contact
##   map and gencode track, BRCA and LUSC page)
## run: FigS5A_job-script_CCND1_median_CNnorm_HiC.sh
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(GxG)
library(gGnome)
library(rtracklayer)
library(parallel)

outdir <- paste0(OUT_DIR, '/')
window_oi <- GRanges('chr11:69000000-70500000')
xchr <- 'chr11'
xG <- 'CCND1'
binsize = 5e3
clim_factor = 4
model <- 'logistic'
tumors <- data.table(sample = c('C169', 'C17', 'C3', 'H3D44S', 'KDAMIF', 'Y57U8R', 'F57J9F'),
    tissue = c('breast', 'breast', 'breast', 'lung', 'lung', 'lung', 'lung'),
    hic_dir = c(file.path(IN_DIR, 'hic/nfcore_results'), file.path(IN_DIR, 'hic/nfcore_results'), file.path(IN_DIR, 'hic/arima_results'),
        file.path(IN_DIR, 'hic/H3D44S_results'), file.path(IN_DIR, 'hic/arima_results'), file.path(IN_DIR, 'hic/arima_results'),
        file.path(IN_DIR, 'hic/arima_results')))
model_file <- function(s) paste0(IN_DIR, '/CNmodels/CCND1/20261023_', s, '_chr7_8_11_20Mb_HiCcov_logistic_df8_x200_d200.rds')
out_prefix <- paste0(outdir, 'FigS5A_', xG, '_median_CNnorm_HiC_BRCA_vs_LUSC_', xchr, '_', start(window_oi) / 1e6, '_', end(window_oi) / 1e6, 'Mb')
# tumor groups of the pages (tissue column of tumors)
groups <- c(BRCA = 'breast', LUSC = 'lung')

# one tumor: CN-normalized counts of all tile pairs of the window (upper triangle incl. the diagonal) and the depth of its model
norm_window <- function(r) {
    s <- tumors$sample[r]
    hic_dir <- tumors$hic_dir[r]
    # 5kb tiles of xchr from the bins of the raw contact map (genome-wide bin index, 1-based in the dumped txt); window tiles by midpoint
    cbins <- fread(paste0(hic_dir, '/contact_maps/bins/cooler_bins_5000.bed'), header = FALSE, select = 1:3, col.names = c('chr', 'start', 'end'))
    idx <- which(cbins$chr == xchr)
    tiles <- GRanges(xchr, IRanges(cbins$start[idx] + 1, cbins$end[idx]))
    win_k <- which(gr.mid(tiles) %^% window_oi)
    stopifnot(length(win_k) > 0, all(diff(win_k) == 1))
    wt <- tiles[win_k]
    glo <- idx[win_k[1]]
    ghi <- idx[win_k[length(win_k)]]
    pr <- fread(cmd = sprintf("awk '$1 >= %d && $1 <= %d && $2 >= %d && $2 <= %d { print $1, $2, $3 }' %s/contact_maps/txt/%s.5000_raw.txt",
        glo, ghi, glo, ghi, hic_dir, s), col.names = c('i', 'j', 'counts'))
    stopifnot(pr[, all(i <= j)])
    pr[, c('i', 'j') := .(i - glo + 1L, j - glo + 1L)]
    # CN effect of the tumor's background model, CN per tile matched by tile midpoint
    m <- readRDS(model_file(s))
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

depth_dt <- data.table(tumors[, .(sample, tissue)], window_counts = sapply(res, `[[`, 'counts'), window_tiles_with_CN = sapply(res, `[[`, 'n_valid'),
    depth_Ed = sapply(res, `[[`, 'depth'))
depth_dt[, depth_scale := depth_Ed / median(depth_Ed)]
fwrite(depth_dt, paste0(out_prefix, '_depth_scale.tsv'), sep = '\t')
print(depth_dt)

# median over the tumors of a group per tile pair (depth-scaled normalized counts; tumors without a CN for a tile left out): BRCA, LUSC and
# all seven (colour scale)
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

# colour scale: 4x the median v4C of the CCND1 TSS over the normal lung TAD lifted to hs1
hg38ToHs1_chain <- import.chain(file.path(IN_DIR, 'hg38ToHs1.over.chain'))
CGC_gr_hg38 <- readRDS('../../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
Lung_norm_16f_TAD_1 <- fread(file.path(IN_DIR, 'Lung_norm_16f_TAD_1_ENCFF581SYY.bedpe'))
colnames(Lung_norm_16f_TAD_1)[1:3] = c('chr', 'start', 'end')
Lung_TAD_gr <- gr.chr(dt2gr(Lung_norm_16f_TAD_1[, 1:3]))
TAD_hg38 <- reduce(Lung_TAD_gr[Lung_TAD_gr %^% gr.chr(CGC_gr_hg38 %Q% (Gene_Symbol == xG))], ignore.strand = TRUE)
TAD_oi <- unlist(range(liftOver(TAD_hg38, hg38ToHs1_chain)))
TAD_oi <- TAD_oi[seqnames(TAD_oi) == xchr]
stopifnot(length(TAD_oi) == 1)
xG_TSSgr <- readRDS(file.path(IN_DIR, 'CCND1_TSSgr_T2T.rds'))[1]
clim_of <- function(gmat) clim_factor * median({v4C_gr <- gmat[gmat$gr %^% xG_TSSgr, drop = T]; v4C_gr$value[gr.mid(v4C_gr) %^% TAD_oi]}, na.rm = T)
clims <- sapply(gmats, clim_of)
print(paste0('lung TAD (hs1): ', gr.string(TAD_oi), '; ', clim_factor, 'x median ', xG, ' v4C over the TAD per median map: ',
    paste(names(clims), signif(clims, 4), collapse = ', ')))
# one colour scale for the BRCA and LUSC pages: that of the median map of all seven tumors
clim_max <- clims[['all']]
stopifnot(clim_max > 0)

gencode_T2T_gt <- readRDS(file.path(IN_DIR, 'final_T2T_gencode_gt_USCS.rds'))
pdf(file = paste0(out_prefix, '.pdf'), width = 8, height = 10)
for (g in names(groups)) {
    smp <- tumors[tissue == groups[[g]], sample]
    plot(c(gencode_T2T_gt, gmats[[g]]$gtrack(name = paste0('median counts/(CN_i*CN_j)^gamma, ', length(smp), ' ', g, ' CCND1 HLAMP tumors (',
        paste(smp, collapse = ', '), ')'), height = 90, clim = c(0, clim_max))), window_oi, cex.label = 0.1)
}
dev.off()
print('Done')
