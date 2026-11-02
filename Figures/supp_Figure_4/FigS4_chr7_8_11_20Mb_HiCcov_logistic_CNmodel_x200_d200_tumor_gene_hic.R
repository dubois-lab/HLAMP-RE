## Continuous copy number (HiC coverage) model with a logistic gamma for the HiC contact counts of one HLAMP tumor, fitted jointly on
## chr7, chr8 and chr11 (intrachromosomal tile pairs, 5kb, up to 20Mb apart, from the raw contact map):
##   observed_ij = E(d_ij) * (CN_i*CN_j)^gamma(d_ij, CN_i*CN_j)
## - CN_i: HiC coverage of tile i (row sum over its own chromosome, self-contact once) relative to the median QC tile of the three
##   chromosomes (tile QC: HiC coverage > 10% of the median of the three chromosomes); one CN scale for all three, so chromosome-wide
##   gains and losses stay in the CN
## - E(d): one natural spline in log10(distance) for the three chromosomes, df = ns_df (8; monotonically decreasing for all df 4-10 up to
##   10Mb)
## - gamma(d, x), x = CN_i*CN_j: g_lo + (g_hi - g_lo) * plogis(a + b z_d + c z_x), z_d / z_x: standardized log10(distance) / log10(x)
## - fit: negative binomial likelihood (NB2, dispersion fitted jointly; nlminb, analytic gradient, two starts) on pairs sampled per cell of
##   200 log-spaced x bins x 200 log-spaced distance bins (up to 400 pairs per cell over the three chromosomes, weighted back to all pairs);
##   compared with gamma = 0.5, gamma = 1 and a constant gamma (MASS::glm.nb) by log-likelihood, AIC and BIC; the fit set does not depend
##   on ns_df, so the fits of different ns_df compare on the same pairs
## - CN-corrected count: count_ij / (CN_i*CN_j)^gamma(d_ij, CN_i*CN_j)
## - pdf (<prefix>.pdf): 3D per cell (x: log10(CN_i*CN_j), y: log10(distance)): mean observed count; model with the observed means as
##   points; gamma with the lung TAD tile pairs of the gene (purple); residuals; mean CN-corrected count. 2D: distance decay E(d); gamma vs
##   log10(CN_i*CN_j) at fixed distances (range of the TAD tile pairs as purple bar)
## - 2D per distance band (<prefix>_distbands{,_linY,_linXY}.pdf; 0-50kb, 50-250kb, 250kb-1Mb, 1-2Mb, 2-5Mb, 5-10Mb, 10-20Mb): density of
##   count / E(d) (grey), mean count / E(d) per log-spaced bin of CN_i*CN_j (100 bins, >= 100 pairs, zeros included; three chromosomes in
##   red, each chromosome as thin line), slope 1 and slope 0.5 lines through the binned means, model CN effect (CN_i*CN_j)^gamma at the
##   lower edge, geometric middle and upper edge of the band (blue), range of the TAD tile pairs (purple bar); log10 x / log10 y (count > 0),
##   log10 x / linear y and linear x / linear y (y cut at the 99.9th percentile)
## - rds: the model for correcting other analyses (beta, theta, gamma_fun, logE_by_dist, z_d, z_x, hic_cov, hiccov_median, good_tile, tile_gr;
##   tiles of all three chromosomes; logistic form: fits$logistic$par, k, forms$logistic, design, comparison), read by the TAD ratio scripts
## hic_dir: results folder of an nf-core/hic run or a .hic file (read with GxG::straw, raw counts; hs1 5kb tiles of the nf-core run as reference)
## run: FigS4_job-script_HiCcov_logistic_chr7_8_11_20Mb_x200_d200.sh; args: sampleid hic_dir ns_df outdir [gene] [map_id]
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')

library(data.table)
library(gUtils)
library(splines)
library(rtracklayer)
library(MASS)

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop('usage: Rscript <script> sampleid hic_dir ns_df outdir [gene, default CCND1] [map_id, default sampleid]', call. = FALSE)
sampleid <- args[1]
hic_dir <- args[2]
ns_df <- as.integer(args[3])
outdir <- sub('/*$', '/', args[4])
# gene whose normal lung TAD is marked in the plots; stem of the raw map file name if it differs from the sample id (NODTWC: NODWTC)
xG <- if (length(args) >= 5) args[5] else 'CCND1'
map_id <- if (length(args) >= 6) args[6] else sampleid
print(paste0('Starting sample: ', sampleid, ' (', hic_dir, '); distance spline df ', ns_df, '; output in ', outdir))

out_prefix <- paste0(outdir, '20261023_', sampleid, '_chr7_8_11_20Mb_HiCcov_logistic_df', ns_df, '_x200_d200')
chrs <- c('chr7', 'chr8', 'chr11')
binsize = 5e3
Dmax = 2e7
n_dstep <- Dmax / binsize
min_hiccov_frac = 0.1
n_xbins = 200
n_dbins = 200
cell_cap = 400
min_cell_pairs = 20
dist_bands <- c(0, 5e4, 2.5e5, 1e6, 2e6, 5e6, 1e7, 2e7)
n_band_bins = 100
min_bin_pairs = 100
chr_cols <- c(chr7 = 'forestgreen', chr8 = 'saddlebrown', chr11 = 'magenta3')
curve_cols <- c('lightskyblue', 'deepskyblue3', 'navy')

# contact map: raw 5kb map of the nf-core/hic run (hic_dir = results dir) or a .hic file (hic_dir = path to <id>_1.hic, linkprep runs; read
# with GxG::straw, norm NONE = raw counts); tiles of the three chromosomes in the genome-wide cooler bin index (1-based in the dumped txt; for
# a .hic file the hs1 5kb bins of the nf-core run as reference); global tile index over the three
is_hic <- grepl('\\.hic$', hic_dir)
bins_file <- if (is_hic) file.path(IN_DIR, 'hic/nfcore_results/contact_maps/bins/cooler_bins_5000.bed') else paste0(hic_dir, '/contact_maps/bins/cooler_bins_5000.bed')
cbins <- fread(bins_file, header = FALSE, select = 1:3, col.names = c('chr', 'start', 'end'))
chr_dt <- rbindlist(lapply(chrs, function(ch) {
    idx <- which(cbins$chr == ch)
    stopifnot(length(idx) > 0, all(diff(idx) == 1))
    data.table(chr = ch, lo = idx[1], hi = idx[length(idx)], len = cbins$end[idx[length(idx)]])
}))
chr_dt[, n := hi - lo + 1L]
chr_dt[, goff := cumsum(c(0L, n[-.N]))]
tile_dt <- rbindlist(lapply(seq_len(nrow(chr_dt)), function(q) data.table(q = q, start = cbins$start[chr_dt$lo[q]:chr_dt$hi[q]] + 1L,
    end = cbins$end[chr_dt$lo[q]:chr_dt$hi[q]])))
rm(cbins)
tile_gr <- GRanges(chr_dt$chr[tile_dt$q], IRanges(tile_dt$start, tile_dt$end))
n_tiles <- length(tile_gr)
tile_q <- tile_dt$q
chr_last <- (chr_dt$goff + chr_dt$n)[tile_q]
rm(tile_dt)

if (is_hic) {
    # .hic: whole chromosome at 5kb, raw counts; straw bins carry 0-based starts, tile index = start %/% binsize + 1
    map_dt <- rbindlist(lapply(seq_len(nrow(chr_dt)), function(q) {
        m <- GxG::straw(hic_dir, gr = GRanges(chr_dt$chr[q], IRanges(1, chr_dt$len[q])), res = binsize, norm = 'NONE')
        t_map <- as.integer(start(m$gr) %/% binsize + 1)
        stopifnot(max(t_map) <= chr_dt$n[q])
        d <- m$dat
        print(paste(chr_dt$chr[q], ':', nrow(d), 'nonzero tile pairs from', hic_dir))
        data.table(i = pmin(t_map[d$i], t_map[d$j]) + chr_dt$goff[q], j = pmax(t_map[d$i], t_map[d$j]) + chr_dt$goff[q], counts = as.integer(round(d$value)))
    }))
} else {
    # raw contact map, intrachromosomal pairs of the three chromosomes in one pass (upper triangle: i in the chromosome and j <= its last bin)
    awk_prog <- paste(sprintf('($1 >= %d && $1 <= %d && $2 <= %d) { print $1, $2, $3, %d; next }', chr_dt$lo, chr_dt$hi, chr_dt$hi, seq_len(nrow(chr_dt))),
        collapse = ' ')
    map_dt <- fread(cmd = paste0("awk '", awk_prog, "' ", hic_dir, '/contact_maps/txt/', map_id, '.5000_raw.txt'), col.names = c('i', 'j', 'counts', 'q'))
    stopifnot(map_dt[, all(i <= j)])
    map_dt[, c('i', 'j') := .(i - chr_dt$lo[q] + 1L + chr_dt$goff[q], j - chr_dt$lo[q] + 1L + chr_dt$goff[q])]
    map_dt[, q := NULL]
}

# HiC coverage: sum of all counts of a tile on its chromosome (row sum, self-contact once); tile QC; CN relative to the median QC tile
hic_cov <- numeric(n_tiles)
s1 <- map_dt[, .(t = sum(counts)), by = i]
hic_cov[s1$i] <- s1$t
s2 <- map_dt[i != j, .(t = sum(counts)), by = j]
hic_cov[s2$j] <- hic_cov[s2$j] + s2$t
rm(s1, s2)
map_dt <- map_dt[j - i >= 1L & j - i <= n_dstep]
good_tile <- hic_cov > min_hiccov_frac * median(hic_cov[hic_cov > 0])
hiccov_median <- median(hic_cov[good_tile])
cn_tile <- hic_cov / hiccov_median
print(data.table(chr = chrs, tiles = tabulate(tile_q, length(chrs)), QC_tiles = tabulate(tile_q[good_tile], length(chrs)),
    median_CN = sapply(seq_along(chrs), function(q) round(median(cn_tile[good_tile & tile_q == q]), 3)),
    max_CN = sapply(seq_along(chrs), function(q) round(max(cn_tile[good_tile & tile_q == q]), 1)),
    tiles_CN_gt_4 = sapply(seq_along(chrs), function(q) sum(good_tile & tile_q == q & cn_tile > 4))))

# all pairs of QC tiles up to Dmax apart on the same chromosome; pairs without contacts count as 0
good_idx <- which(good_tile)
pair_dt <- rbindlist(lapply(seq_len(n_dstep), function(d) {
    i <- good_idx[good_idx + d <= chr_last[good_idx]]
    i <- i[good_tile[i + d]]
    data.table(i = i, j = i + d)
}))
setkey(pair_dt, i, j)
setkey(map_dt, i, j)
pair_dt[map_dt, counts := i.counts]
pair_dt[is.na(counts), counts := 0L]
rm(map_dt)
gc()
pair_dt[, c('ds', 'q') := .(j - i, tile_q[i])]
pair_dt[, c('ld', 'lx') := .(log10(ds * binsize), log10(cn_tile[i] * cn_tile[j]))]

# normal lung TAD of the gene xG (hg38) lifted to hs1; TAD tile pairs: both tiles in the lifted TAD
hg38ToHs1_chain <- import.chain(file.path(IN_DIR, 'hg38ToHs1.over.chain'))
CGC_gr_hg38 <- readRDS('../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
Lung_norm_16f_TAD_1 <- fread(file.path(IN_DIR, 'Lung_norm_16f_TAD_1_ENCFF581SYY.bedpe'))
colnames(Lung_norm_16f_TAD_1)[1:3] = c('chr', 'start', 'end')
Lung_TAD_gr <- gr.chr(dt2gr(Lung_norm_16f_TAD_1[, 1:3]))
xG_gr_hg38 <- gr.chr(CGC_gr_hg38 %Q% (Gene_Symbol == xG))
TAD_hg38 <- reduce(Lung_TAD_gr[Lung_TAD_gr %^% xG_gr_hg38], ignore.strand = TRUE)
TAD_lift <- unlist(range(liftOver(TAD_hg38, hg38ToHs1_chain)))
TAD_oi <- TAD_lift[seqnames(TAD_lift) == as.character(seqnames(xG_gr_hg38))[1]]
in_TAD_t <- seq_len(n_tiles) %in% which(gr.mid(tile_gr) %^% TAD_oi)
pair_dt[, in_TAD := in_TAD_t[i] & in_TAD_t[j]]
TAD_lx <- pair_dt[in_TAD == TRUE, range(lx)]
print(paste(sum(good_tile), 'of', n_tiles, 'tiles pass QC (median HiC coverage', hiccov_median, ');', nrow(pair_dt), 'tile pairs up to', Dmax / 1e6, 'Mb apart (',
    paste(paste(chrs, tabulate(pair_dt$q, length(chrs))), collapse = ', '), '),', round(100 * mean(pair_dt$counts == 0), 1), '% with count 0;',
    pair_dt[, sum(in_TAD)], xG, 'lung TAD tile pairs, log10(CN_i*CN_j)', round(TAD_lx[1], 2), '-', round(TAD_lx[2], 2)))

# cells: 100 log-spaced bins of CN_i*CN_j x 50 log-spaced distance bins
x_breaks <- seq(min(pair_dt$lx), max(pair_dt$lx), length.out = n_xbins + 1)
d_breaks <- seq(log10(binsize), log10(Dmax), length.out = n_dbins + 1)
x_mid <- (x_breaks[-1] + x_breaks[-length(x_breaks)]) / 2
d_mid <- (d_breaks[-1] + d_breaks[-length(d_breaks)]) / 2
pair_dt[, c('xbin', 'dbin') := .(cut(lx, x_breaks, include.lowest = TRUE, labels = FALSE), cut(ld, d_breaks, include.lowest = TRUE, labels = FALSE))]

# fit set: up to cell_cap pairs per cell, weighted back to all pairs of the cell
set.seed(20261022)
pair_dt[, idx := .I]
samp_dt <- pair_dt[, .(idx = if (.N <= cell_cap) idx else sample(idx, cell_cap), w = max(1, .N / cell_cap)), by = .(xbin, dbin)]
fit_dt <- pair_dt[samp_dt$idx, .(counts, ld, lx, q)]
fit_dt[, w := samp_dt$w / mean(samp_dt$w)]
rm(samp_dt)
pair_dt[, idx := NULL]
print(paste(nrow(fit_dt), 'pairs in the fit set from', uniqueN(pair_dt, by = c('xbin', 'dbin')), 'cells (', paste(paste(chrs, tabulate(fit_dt$q, length(chrs))),
    collapse = ', '), ')'))

# gamma and the CN effect in a function, so the closures carry what they need when saved in the rds and loaded elsewhere
make_model <- function(z_d, z_x) {
    for (a in list(z_d, z_x)) force(a)
    gamma_fun <- function(th, zd, zx) th[1] + (th[2] - th[1]) * plogis(th[3] + th[4] * zd + th[5] * zx)
    design <- function(ld, lx) list(zd = (ld - z_d[1]) / z_d[2], zx = (lx - z_x[1]) / z_x[2], L = lx * log(10))
    forms <- list(logistic = list(label = 'logistic: eta = (g_lo + (g_hi - g_lo) * plogis(a + b z_d + c z_x)) * log(x)',
        par_names = c('g_lo', 'g_hi', 'a', 'b', 'c'), n = 5L,
        eta = function(th, D, deriv = FALSE) list(eta = gamma_fun(th, D$zd, D$zx) * D$L)))
    list(gamma_fun = gamma_fun, design = design, forms = forms)
}

# model: log E[count] = B(d) beta + gamma(d, x) * log(x); B: intercept + natural spline in log10(distance)
ns_d <- ns(fit_dt$ld, df = ns_df)
B <- cbind(1, ns_d)
k <- ncol(B)
y <- fit_dt$counts
w <- fit_dt$w
L <- fit_dt$lx * log(10)
zd_c <- c(mean(fit_dt$ld), sd(fit_dt$ld))
zx_c <- c(mean(fit_dt$lx), sd(fit_dt$lx))
zd <- (fit_dt$ld - zd_c[1]) / zd_c[2]
zx <- (fit_dt$lx - zx_c[1]) / zx_c[2]
mf <- make_model(zd_c, zx_c)
gamma_fun <- mf$gamma_fun
# matrix products with base functions: IRanges masks drop() and gUtils masks %*%, together they return a placeholder symbol
lin_pred <- function(X, b) base::drop(base::`%*%`(X, b))
# weighted negative binomial (NB2, size r) log-likelihood for means mu
nb_ll <- function(mu, r) sum(w * (lgamma(y + r) - lgamma(r) - lgamma(y + 1) + r * (log(r) - log(r + mu)) + y * (log(mu) - log(r + mu))))
# parameters p: beta (k), gamma (g_lo, g_hi, a, b, c), log(r)
nb_nll <- function(p) {
    eta <- lin_pred(B, p[1:k]) + gamma_fun(p[k + 1:5], zd, zx) * L
    -nb_ll(exp(eta), exp(p[k + 6]))
}
nb_nll_grad <- function(p) {
    th <- p[k + 1:5]
    r <- exp(p[k + 6])
    s <- plogis(th[3] + th[4] * zd + th[5] * zx)
    eta <- lin_pred(B, p[1:k]) + (th[1] + (th[2] - th[1]) * s) * L
    mu <- exp(eta)
    g_eta <- w * r * (y - mu) / (r + mu)
    ds <- (th[2] - th[1]) * s * (1 - s) * L
    g_logr <- sum(w * r * (digamma(y + r) - digamma(r) + log(r) - log(r + mu) + 1 - (y + r) / (r + mu)))
    -c(base::drop(base::crossprod(B, g_eta)), sum(g_eta * L * (1 - s)), sum(g_eta * L * s), sum(g_eta * ds), sum(g_eta * ds * zd), sum(g_eta * ds * zx), g_logr)
}

# reference models: gamma = 0.5, gamma = 1, constant gamma
nb_sqrt <- suppressWarnings(glm.nb(y ~ ns_d + offset(0.5 * L), weights = w))
nb_prod <- suppressWarnings(glm.nb(y ~ ns_d + offset(L), weights = w))
nb_const <- suppressWarnings(glm.nb(y ~ ns_d + L, weights = w))
gamma_const <- coef(nb_const)[['L']]

# logistic gamma from two starts (from the gamma = 0.5 and from the constant gamma fit), the better one kept
starts <- list(c(unname(coef(nb_sqrt)), 0.5, 1, 0, 1, 1, log(nb_sqrt$theta)),
    c(unname(coef(nb_const)[1:k]), max(0.05, gamma_const - 0.1), gamma_const + 0.3, 0, 1, 1, log(nb_const$theta)))
opts <- lapply(starts, function(p0) nlminb(p0, nb_nll, nb_nll_grad, lower = c(rep(-Inf, k), 0, 0, -30, -30, -30, -10),
    upper = c(rep(Inf, k), 3, 3, 30, 30, 30, 20), control = list(iter.max = 3000, eval.max = 6000, rel.tol = 1e-10)))
for (o in opts) print(paste('nlminb:', o$message, '; iterations', o$iterations, '; log-likelihood', round(-o$objective, 1)))
opt <- opts[[which.min(sapply(opts, `[[`, 'objective'))]]
beta <- opt$par[1:k]
theta <- setNames(opt$par[k + 1:5], c('g_lo', 'g_hi', 'a', 'b', 'c'))
nb_size <- exp(opt$par[k + 6])
mu_fit <- exp(lin_pred(B, beta) + gamma_fun(theta, zd, zx) * L)
n_fit <- nrow(fit_dt)
dev_dt <- data.table(sample = sampleid, ns_df = ns_df, form = c('gamma_0.5', 'gamma_1', 'constant', 'logistic'),
    model = c('gamma = 0.5', 'gamma = 1', paste0('gamma = ', round(gamma_const, 3), ' (constant)'), 'gamma(d, CN_i*CN_j)'),
    nb_size = c(nb_sqrt$theta, nb_prod$theta, nb_const$theta, nb_size), n_par = c(k + 1, k + 1, k + 2, k + 6),
    loglik = c(nb_ll(fitted(nb_sqrt), nb_sqrt$theta), nb_ll(fitted(nb_prod), nb_prod$theta), nb_ll(fitted(nb_const), nb_const$theta), nb_ll(mu_fit, nb_size)),
    nlminb = c(NA, NA, NA, opt$message))
dev_dt[, c('AIC', 'BIC') := .(2 * n_par - 2 * loglik, log(n_fit) * n_par - 2 * loglik)]
dev_dt[, c('delta_AIC', 'delta_BIC', 'delta_cv_loglik') := .(AIC - min(AIC), BIC - min(BIC), NA_real_)]
print(theta)
print(paste('negative binomial size (dispersion):', round(nb_size, 3)))
print(dev_dt[, .(model, nb_size, n_par, loglik, AIC, delta_AIC, BIC, delta_BIC)])

# decay and gamma for all pairs; model mean and CN-corrected count
logE_by_dist <- lin_pred(cbind(1, predict(ns_d, log10(seq_len(n_dstep) * binsize))), beta)
# E(d) must decrease monotonically with the distance
dle <- diff(logE_by_dist) / log(10)
dev_dt[, c('Ed_monotone_decreasing', 'Ed_n_increasing_steps', 'Ed_max_local_slope') := .(all(dle < 0), sum(dle >= 0),
    max(dle / diff(log10(seq_len(n_dstep) * binsize))))]
print(paste0('E(d) monotonically decreasing over 5kb-', Dmax / 1e6, 'Mb: ', all(dle < 0), ' (', sum(dle >= 0), ' non-decreasing steps of ', length(dle),
    ', from ', if (any(dle >= 0)) paste0(min(which(dle >= 0)) * binsize / 1e3, 'kb') else '-', '; flattest local slope d log10 E / d log10 d ',
    round(max(dle / diff(log10(seq_len(n_dstep) * binsize))), 3), ')'))
if (!all(dle < 0)) warning(paste(sampleid, ': E(d) not monotonically decreasing'))
pair_dt[, gam := gamma_fun(theta, (ld - zd_c[1]) / zd_c[2], (lx - zx_c[1]) / zx_c[2])]
pair_dt[, mu := exp(logE_by_dist[ds] + gam * lx * log(10))]
pair_dt[, count_cncorr := counts / 10^(gam * lx)]
print(paste('gamma over all pairs:', paste(round(quantile(pair_dt$gam, c(0, 0.01, 0.5, 0.99, 1)), 3), collapse = ' / '), '(min / 1% / median / 99% / max);',
    'TAD tile pairs:', paste(round(pair_dt[in_TAD == TRUE, quantile(gam, c(0, 0.5, 1))], 3), collapse = ' / '), '(min / median / max)'))

# cell means (all three chromosomes)
cell_dt <- pair_dt[, .(n_pairs = .N, mean_count = mean(counts), mean_model = mean(mu), mean_cncorr = mean(count_cncorr), mean_gamma = mean(gam),
    n_TAD = sum(in_TAD)), keyby = .(xbin, dbin)]
cell_dt[, c('log10_x', 'log10_dist') := .(x_mid[xbin], d_mid[dbin])]
fwrite(data.table(sample = sampleid, ns_df = ns_df, cell_dt), paste0(out_prefix, '_cells.tsv'), sep = '\t')
cell_ok <- cell_dt$n_pairs >= min_cell_pairs
cell_mat <- function(v) {
    m <- matrix(NA_real_, n_xbins, n_dbins)
    keep <- cell_ok & is.finite(v)
    m[cbind(cell_dt$xbin[keep], cell_dt$dbin[keep])] <- v[keep]
    m
}
log10_pos <- function(v) ifelse(v > 0, log10(v), NA_real_)

# 3D surface over the cells, facets coloured by height
persp_page <- function(zmat, zlab, main, pal, zlim = range(zmat, finite = TRUE), add_fun = NULL, theta_view = -40, phi_view = 25) {
    zf <- (zmat[-1, -1] + zmat[-1, -n_dbins] + zmat[-n_xbins, -1] + zmat[-n_xbins, -n_dbins]) / 4
    col_idx <- pmax(1, pmin(length(pal), ceiling((zf - zlim[1]) / diff(zlim) * length(pal))))
    facet_col <- pal[col_idx]
    facet_col[is.na(facet_col)] <- 'white'
    pm <- persp(x_mid, d_mid, zmat, theta = theta_view, phi = phi_view, col = facet_col, border = NA, ticktype = 'detailed', zlim = zlim,
        xlab = 'log10(CN_i*CN_j)', ylab = 'log10(distance)', zlab = zlab, main = main, expand = 0.75, cex.axis = 0.6, cex.lab = 0.8, cex.main = 0.75)
    if (!is.null(add_fun)) add_fun(pm)
    invisible(pm)
}
obs_mat <- cell_mat(log10_pos(cell_dt$mean_count))
model_mat <- cell_mat(log10_pos(cell_dt$mean_model))
count_zlim <- range(c(obs_mat, model_mat), finite = TRUE)
pal_count <- hcl.colors(100, 'viridis')
model_txt <- paste0('gamma = ', round(theta[['g_lo']], 2), ' + ', round(theta[['g_hi']] - theta[['g_lo']], 2), ' * plogis(', round(theta[['a']], 2), ' + ',
    round(theta[['b']], 2), ' z_d + ', round(theta[['c']], 2), ' z_x); E(d) spline df ', ns_df)
chr_txt <- paste(chrs, collapse = ', ')

print('plotting')
pdf(file = paste0(out_prefix, '.pdf'), width = 8, height = 8)
persp_page(obs_mat, 'log10(mean observed count)', paste0(sampleid, ': ', chr_txt, ' tile pairs (5kb, <= ', Dmax / 1e6, 'Mb), mean observed count per cell\n(',
    n_xbins, ' log-spaced CN_i*CN_j bins x ', n_dbins, ' log-spaced distance bins, cells >= ', min_cell_pairs, ' pairs, zeros incl.; CN = HiC coverage / median tile)'),
    pal_count, zlim = count_zlim)
persp_page(model_mat, 'log10(mean count)', paste0(sampleid, ': model E(d) * (CN_i*CN_j)^gamma(d, CN_i*CN_j) (surface) and observed cell means (points:\n',
    'red above, blue below the model); ', model_txt), pal_count, zlim = count_zlim,
    add_fun = function(pm) {
        keep <- cell_ok & cell_dt$mean_count > 0
        # trans3d fails on an empty point set
        if (any(keep)) {
            p <- trans3d(cell_dt$log10_x[keep], cell_dt$log10_dist[keep], log10(cell_dt$mean_count[keep]), pm)
            points(p, pch = 16, cex = 0.25, col = ifelse(cell_dt$mean_count[keep] > cell_dt$mean_model[keep], 'red', 'blue'))
        }
    })
gamma_mat <- cell_mat(cell_dt$mean_gamma)
persp_page(gamma_mat, 'gamma', paste0(sampleid, ': gamma(d, CN_i*CN_j) per cell; purple: cells with ', xG, ' lung TAD tile pairs\n', model_txt),
    hcl.colors(100, 'Plasma'), zlim = range(c(0.5, 1, gamma_mat), finite = TRUE), theta_view = 40,
    add_fun = function(pm) {
        # with 200 x 200 cells the TAD tile pairs can fall into cells below min_cell_pairs only (1E2AVA, NSD3): no TAD points then
        keep <- cell_ok & cell_dt$n_TAD > 0
        if (any(keep)) points(trans3d(cell_dt$log10_x[keep], cell_dt$log10_dist[keep], cell_dt$mean_gamma[keep], pm), pch = 16, cex = 0.5, col = 'purple')
    })
resid_mat <- cell_mat(log10_pos(cell_dt$mean_count) - log10_pos(cell_dt$mean_model))
resid_max <- max(abs(resid_mat), na.rm = T)
persp_page(resid_mat, 'log10(observed / model)', paste0(sampleid, ': residuals per cell, log10(mean observed count / mean model count)'),
    hcl.colors(100, 'Blue-Red 3'), zlim = c(-resid_max, resid_max))
persp_page(cell_mat(log10_pos(cell_dt$mean_cncorr)), 'log10(mean CN-corrected count)',
    paste0(sampleid, ': mean CN-corrected count per cell, count / (CN_i*CN_j)^gamma(d, CN_i*CN_j)\n(should depend on the distance only)'), pal_count)
# 2D: distance decay and gamma at fixed distances
d_grid <- seq(log10(binsize), log10(Dmax), length.out = 300)
plot(d_grid, lin_pred(cbind(1, predict(ns_d, d_grid)), beta) / log(10), type = 'l', lwd = 2, las = 1, xlab = 'log10(distance)',
    ylab = 'log10 E(d) (count of two median-coverage tiles)', main = paste0(sampleid, ': distance decay E(d) at CN_i*CN_j = 1 (natural spline, df ', ns_df,
        ')\n', chr_txt), cex.main = 0.85)
x_grid <- seq(min(pair_dt$lx), max(pair_dt$lx), length.out = 300)
d_show <- c(5e3, 5e4, 2.5e5, 1e6, 5e6, 1e7, 2e7)
d_cols <- hcl.colors(length(d_show), 'Dark 3')
plot(NA, xlim = range(x_grid), ylim = range(c(0.4, 1.1, theta[c('g_lo', 'g_hi')])), las = 1, xlab = 'log10(CN_i*CN_j)', ylab = 'gamma',
    main = paste0(sampleid, ': gamma(d, CN_i*CN_j) at fixed distances, ', chr_txt, '\n', model_txt), cex.main = 0.8)
abline(h = c(0.5, 1), lty = c(1, 2), col = c('darkorange', 'black'))
for (q in seq_along(d_show)) lines(x_grid, gamma_fun(theta, (log10(d_show[q]) - zd_c[1]) / zd_c[2], (x_grid - zx_c[1]) / zx_c[2]), col = d_cols[q], lwd = 2)
y_bar <- par('usr')[3] + 0.015 * diff(par('usr')[3:4])
segments(TAD_lx[1], y_bar, TAD_lx[2], y_bar, col = 'purple', lwd = 6, lend = 1)
legend('topleft', legend = c(paste0(d_show / 1e3, 'kb apart'), 'gamma 0.5 (sqrt)', 'gamma 1 (product)', paste0('range of the ', xG, ' lung TAD tile pairs')),
    col = c(d_cols, 'darkorange', 'black', 'purple'), lty = c(rep(1, length(d_show)), 1, 2, 1), lwd = c(rep(2, length(d_show)), 1, 1, 6), cex = 0.7,
    bg = adjustcolor('white', 0.8), box.col = NA)
dev.off()

# 2D per distance band: observed count / E(d) per log-spaced bin of CN_i*CN_j (all
# chromosomes and per chromosome) and the model CN effect (CN_i*CN_j)^gamma at three fixed distances
cn_effect <- function(ld, lx) 10^(gamma_fun(theta, (ld - zd_c[1]) / zd_c[2], (lx - zx_c[1]) / zx_c[2]) * lx)
band_labels <- paste0(dist_bands[-length(dist_bands)] / 1e3, '-', dist_bands[-1] / 1e3, 'kb apart')
band_bins <- function(b) {
    pdt <- pair_dt[ds > dist_bands[b] / binsize & ds <= dist_bands[b + 1] / binsize, .(counts, lx, gam, q, in_TAD, ds)]
    pdt[, c('yn', 'eff') := .(counts / exp(logE_by_dist[ds]), 10^(gam * lx))]
    d_lo <- max(binsize, dist_bands[b])
    d_edges <- c(d_lo, sqrt(d_lo * dist_bands[b + 1]), dist_bands[b + 1])
    xb <- seq(min(pdt$lx), max(pdt$lx), length.out = n_band_bins + 1)
    pdt[, xbin := cut(lx, xb, include.lowest = TRUE, labels = FALSE)]
    all_dt <- pdt[, .(chr = 'all', n_pairs = .N, mean_x = mean(10^lx), mean_count = mean(counts), mean_count_over_Ed = mean(yn), mean_model_over_Ed = mean(eff)),
        keyby = xbin]
    per_chr <- pdt[, .(n_pairs = .N, mean_x = mean(10^lx), mean_count = mean(counts), mean_count_over_Ed = mean(yn), mean_model_over_Ed = mean(eff)),
        keyby = .(q, xbin)][, chr := chrs[q]][, q := NULL]
    bin_dt <- rbind(all_dt, per_chr, use.names = TRUE)[n_pairs >= min_bin_pairs & mean_count_over_Ed > 0]
    lbx <- log10(bin_dt$mean_x)
    for (q in seq_along(d_edges)) bin_dt[, paste0('model_at_', signif(d_edges[q] / 1e3, 3), 'kb') := cn_effect(rep(log10(d_edges[q]), length(lbx)), lbx)]
    lgrid <- seq(min(pdt$lx), max(pdt$lx), length.out = 300)
    curves <- lapply(d_edges, function(d) log10(cn_effect(rep(log10(d), length(lgrid)), lgrid)))
    min_slope <- sapply(curves, function(cv) min(diff(cv) / diff(lgrid)))
    TAD_b <- if (pdt[, any(in_TAD)]) pdt[in_TAD == TRUE, range(lx)] else c(NA_real_, NA_real_)
    resid_chr <- bin_dt[, .(med = round(median(abs(log10(mean_count_over_Ed / mean_model_over_Ed))), 3)), keyby = chr]
    print(paste0(band_labels[b], ': ', nrow(pdt), ' pairs; median |log10(observed / model)| over the bins (both / E(d), model averaged over the pairs ',
        'of a bin): ', paste(resid_chr$chr, resid_chr$med, collapse = ', '), '; lowest local exponent of the model at ', paste(signif(d_edges / 1e3, 3),
        collapse = ' / '), 'kb: ', paste(round(min_slope, 3), collapse = ' / ')))
    list(label = band_labels[b], pdt = pdt, bin_dt = bin_dt, lgrid = lgrid, curves = curves, d_edges = d_edges, min_slope = min_slope, TAD_lx = TAD_b,
        bins = data.table(sample = sampleid, ns_df = ns_df, band = band_labels[b], bin_dt))
}
# density of the pairs in grey; falls back to a random subset of points without KernSmooth
dens_plot <- function(x, y, ...) {
    if (requireNamespace('KernSmooth', quietly = TRUE)) {
        smoothScatter(x, y, nrpoints = 0, colramp = colorRampPalette(c('white', 'grey25')), ...)
    } else {
        s <- sample.int(length(x), min(2e5, length(x)))
        plot(x[s], y[s], pch = 16, cex = 0.2, col = adjustcolor('grey30', 0.1), ...)
    }
}
band_page <- function(fit, x_log, y_log) {
    to_x <- function(l) if (x_log) l else 10^l
    to_y <- function(l) if (y_log) l else 10^l
    pdt <- fit$pdt
    ba <- fit$bin_dt[chr == 'all']
    lbx <- log10(ba$mean_x)
    lby <- log10(ba$mean_count_over_Ed)
    if (y_log) {
        sel <- pdt$counts > 0
        y <- log10(pdt$yn[sel])
        sub_title <- paste0('count > 0: ', sum(sel), ' of ', nrow(pdt), ' pairs')
    } else {
        y_max <- max(quantile(pdt$yn, 0.999), 10^lby)
        sel <- pdt$yn <= y_max
        y <- pdt$yn[sel]
        sub_title <- paste0('y axis cut at the 99.9th percentile (', signif(y_max, 3), '; ', sum(!sel), ' of ', nrow(pdt), ' pairs above)')
    }
    par(mar = c(6, 4.5, 4, 1))
    dens_plot(to_x(pdt$lx[sel]), y, xlim = range(to_x(pdt$lx)), ylim = range(c(y, to_y(lby), to_y(unlist(fit$curves)))),
        xlab = if (x_log) 'log10(CN_i*CN_j)' else 'CN_i*CN_j', ylab = if (y_log) 'log10(observed count / E(d))' else 'observed count / E(d)',
        main = paste0(sampleid, ': ', chr_txt, ' tile pairs (5kb), ', fit$label, '; CN = HiC coverage / median tile\nE(d): distance decay of the model (df ',
            ns_df, '); ', sub_title), cex.main = 0.75)
    for (ch in chrs) {
        bc <- fit$bin_dt[chr == ch]
        if (nrow(bc) > 1) lines(to_x(log10(bc$mean_x)), to_y(log10(bc$mean_count_over_Ed)), col = chr_cols[[ch]], lwd = 1)
    }
    points(to_x(lbx), to_y(lby), pch = 16, col = 'red', cex = 0.6)
    lines(to_x(lbx), to_y(lby), col = 'red')
    lgrid <- seq(min(lbx), max(lbx), length.out = 400)
    lines(to_x(lgrid), to_y(mean(lby - lbx) + lgrid), col = 'black', lty = 2, lwd = 2)
    lines(to_x(lgrid), to_y(mean(lby - 0.5 * lbx) + 0.5 * lgrid), col = 'darkorange', lwd = 2)
    for (q in seq_along(fit$curves)) lines(to_x(fit$lgrid), to_y(fit$curves[[q]]), col = curve_cols[q], lwd = 2.5)
    leg <- c(paste0('mean observed count / E(d) per log-spaced bin, ', chr_txt, ' (', n_band_bins, ' bins, zeros incl., >= ', min_bin_pairs, ' pairs)'),
        paste0('the same per chromosome: ', chrs), 'slope 1 (count ~ CN_i*CN_j)', 'slope 0.5 (count ~ sqrt(CN_i*CN_j))',
        paste0('model CN effect at ', signif(fit$d_edges / 1e3, 3), 'kb (lowest local exponent ', round(fit$min_slope, 2), ')'))
    leg_col <- c('red', chr_cols[chrs], 'black', 'darkorange', curve_cols)
    leg_lty <- c(1, rep(1, length(chrs)), 2, 1, 1, 1, 1)
    leg_lwd <- c(2, rep(1, length(chrs)), 2, 2, 2.5, 2.5, 2.5)
    if (all(is.finite(fit$TAD_lx))) {
        y_bar <- par('usr')[3] + 0.015 * diff(par('usr')[3:4])
        segments(to_x(fit$TAD_lx[1]), y_bar, to_x(fit$TAD_lx[2]), y_bar, col = 'purple', lwd = 6, lend = 1)
        leg <- c(leg, paste0('range of the ', xG, ' lung TAD tile pairs'))
    } else {
        leg <- c(leg, paste0('no ', xG, ' lung TAD tile pairs at these distances'))
    }
    legend('topleft', legend = leg, col = c(leg_col, 'purple'), lty = c(leg_lty, 1), lwd = c(leg_lwd, 6), cex = 0.55, bg = adjustcolor('white', 0.8), box.col = NA)
    mtext(paste0('logistic CN model: CN effect (CN_i*CN_j)^gamma, ', model_txt), side = 1, line = 4.5, cex = 0.5)
}
# the three pdfs stay open, every band page goes into each
pdf(file = paste0(out_prefix, '_distbands.pdf'), width = 7, height = 7.5)
dev_log <- dev.cur()
pdf(file = paste0(out_prefix, '_distbands_linY.pdf'), width = 7, height = 7.5)
dev_linY <- dev.cur()
pdf(file = paste0(out_prefix, '_distbands_linXY.pdf'), width = 7, height = 7.5)
dev_linXY <- dev.cur()
bins_list <- list()
for (b in seq_along(band_labels)) {
    fit <- band_bins(b)
    dev.set(dev_log)
    band_page(fit, x_log = TRUE, y_log = TRUE)
    dev.set(dev_linY)
    band_page(fit, x_log = TRUE, y_log = FALSE)
    dev.set(dev_linXY)
    band_page(fit, x_log = FALSE, y_log = FALSE)
    bins_list[[b]] <- fit$bins
    rm(fit)
    gc()
}
dev.off(dev_linXY)
dev.off(dev_linY)
dev.off(dev_log)
fwrite(rbindlist(bins_list, fill = TRUE), paste0(out_prefix, '_distbands_bins.tsv'), sep = '\t')

# parameters, deviances and the model for correcting other analyses
fwrite(data.table(sample = sampleid, ns_df = ns_df, parameter = c(names(theta), 'nb_size', 'z_d_mean', 'z_d_sd', 'z_x_mean', 'z_x_sd', 'hiccov_median', 'gamma_constant'),
    value = c(theta, nb_size, zd_c, zx_c, hiccov_median, gamma_const)), paste0(out_prefix, '_parameters.tsv'), sep = '\t')
fwrite(dev_dt, paste0(out_prefix, '_model_comparison.tsv'), sep = '\t')
saveRDS(list(sample = sampleid, chrs = chrs, binsize = binsize, Dmax = Dmax, ns_df = ns_df, beta = beta, theta = theta, nb_size = nb_size,
    ns_knots = attr(ns_d, 'knots'), ns_boundary_knots = attr(ns_d, 'Boundary.knots'), z_d = zd_c, z_x = zx_c, logE_by_dist = logE_by_dist,
    hiccov_median = hiccov_median, hic_cov = hic_cov, good_tile = good_tile, tile_gr = tile_gr, gamma_fun = mf$gamma_fun,
    k = k, fits = list(logistic = list(par = opt$par, loglik = -opt$objective, message = opt$message)), forms = mf$forms, design = mf$design,
    comparison = dev_dt), paste0(out_prefix, '.rds'))
print(paste0('Done: ', sampleid, ', df ', ns_df))
