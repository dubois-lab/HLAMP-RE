## TAD vs random window plots and permutation tests for the tumors of the 200 x 200 cell background models, partners and random
## window buffer up to 20Mb
## Within- vs outside-TAD contact ratios per 5kb bin of the normal lung TAD overlapping a gene (CCND1 or NSD3), next to two sets of 200 random
## windows of the tested TAD's size placed across chr7, chr8 and chr11 (one value per window: median of its per-bin ratios):
## three normal tissue Hi-C (hg38, raw counts) vs the HLAMP tumors (hs1, lifted TAD, counts normalized by the CN effect of the logistic HiC
## coverage CN model fitted jointly on chr7, chr8 and chr11 up to 20Mb, only TAD tiles with JaBbA CN > 9; per-tumor TSVs from the tumor script)
## Partners up to Dmax = 20Mb apart for normals and tumors (tumors from the raw 5kb contact maps; 20Mb = range of the CN model fit).
## Random windows: candidate starts on chr7, chr8 and chr11 where the window +-Dmax lies on the chromosome and the window itself has a HiC bin
## (normals) / a JaBbA CN (tumors) for every tile, and
## - set 'randomwin': no tile of the gene's TAD within window +- Dmax (tumors: and no tile with JaBbA CN > 9)
## - set 'randomwin_noampfilter': the window itself holds no tile of the gene's TAD (tumors: amplified tiles allowed; windows with a tile of
##   JaBbA CN > 9 marked orange in the plots)
## candidates of the three chromosomes pooled and drawn in random order (seed 20261018 per set), a window is kept if it does not overlap a
## kept window of its chromosome, up to n_win = 200
## Permutation test per sample, window set and in/out ratio (ratio, ratio_incl_self, ratio_lennorm, ratio_lennorm_incl_self): one-sided
## empirical p-value of the lung TAD median ratio against the random window medians of the same sample, p = (1 + #windows >= TAD) /
## (1 + #windows); written with the ratios to _permtest_*.csv
## - q-values: Benjamini-Hochberg over the samples (normals and tumors) of the gene per window set and ratio type (q_value in _permtest_*.csv);
##   in the plots a black star above a sample for q < 0.1, a grey star for p < 0.05 (and q >= 0.1)
## - normals vs tumors: per sample, lung TAD median ratio_lennorm_incl_self / median of the random window medians; three normal tissues vs
##   all tumors of the gene, two-sided Wilcoxon rank-sum (Mann-Whitney U) test per window set (W: tumor-normal pairs with the tumor higher);
##   _lungTAD_over_randomwin_ratio_lennorm_incl_self_normals_vs_tumors{_values.tsv,_wilcox.csv,_wilcox.pdf} (box-dot plot, dots by tissue)
## The per-tumor TSVs are looked up in outdir and then in further folders (4th argument); no NA window row for a chromosome without candidates
## run: FigS5B_job-script_comb_TADratio_200sqrd_20Mb.sh; args: gene outdir tumors (comma separated) [further folders with tumor TSVs, comma separated]
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')

library(GxG)
library(gGnome)
library(parallel)

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop('usage: Rscript <script> gene outdir tumors (comma separated) [tumor TSV folders, comma separated]', call. = FALSE)
xG <- args[1]
outdir <- sub('/*$', '/', args[2])
tumors <- strsplit(args[3], ',')[[1]]
tumor_dirs <- c(outdir, if (length(args) >= 4) sub('/*$', '/', strsplit(args[4], ',')[[1]]))
print(paste('gene', xG, '; output in', outdir, '; tumors', paste(tumors, collapse = ', '), '; tumor TSVs from', paste(tumor_dirs, collapse = ', ')))
model <- 'logistic'
norm_tag <- paste0('hiccov_', model, 'norm')
chrs <- c('chr7', 'chr8', 'chr11')
win_sets <- c('randomwin', 'randomwin_noampfilter')
CGC_gr_hg38 <- readRDS('../../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
Oncgr <- gr.chr(CGC_gr_hg38 %Q% (Gene_Symbol == xG))[1]
binsize = 5e3
Dmax = 2e7
cn_cutoff = 9
n_win = 200
# hg38 chromosome lengths: straw query = the whole chromosome
chr_len_hg38 <- c(chr7 = 159345973, chr8 = 145138636, chr11 = 135086622)

# normal lung TAD(s) overlapping xG; nested Arrowhead domains are merged to their outer span
Lung_norm_16f_TAD_1 <- fread(file.path(IN_DIR, 'Lung_norm_16f_TAD_1_ENCFF581SYY.bedpe'))
colnames(Lung_norm_16f_TAD_1)[1:3] = c('chr', 'start', 'end')
Lung_TAD_gr <- gr.chr(dt2gr(Lung_norm_16f_TAD_1))
TAD_oi <- reduce(Lung_TAD_gr[Lung_TAD_gr %^% Oncgr], ignore.strand = TRUE)
xchr <- as.character(seqnames(TAD_oi))
stopifnot(length(TAD_oi) == 1, xchr %in% chrs)
TAD_str <- paste0(xchr, ':', start(TAD_oi), '-', end(TAD_oi))
print(paste('TAD used (hg38):', TAD_str))

normHiCs <- c(breast = file.path(IN_DIR, 'ENCFF924WYB_BRCA_hg38_norm.hic'),
    esophagus = file.path(IN_DIR, 'ENCFF224HKR_squamousEpEsoph_hg38.hic'),
    lung = file.path(IN_DIR, 'ENCFF094KKX_lung16f_hg38.hic'))

add_ratios <- function(dt) {
    dt[, ratio := within / outside]
    dt[, ratio_incl_self := (within + self) / outside]
    # length normalized interaction count sums (self = 1 bin)
    dt[, within_lennorm := within / n_within]
    dt[, outside_lennorm := outside / n_outside]
    dt[, ratio_lennorm := within_lennorm / outside_lennorm]
    dt[, ratio_lennorm_incl_self := ((within + self) / (n_within + 1)) / outside_lennorm]
    dt
}
ratio_cols <- c('ratio', 'ratio_incl_self', 'ratio_lennorm', 'ratio_lennorm_incl_self')

# one normal: TAD bins (xchr) and both random window sets across the chromosomes
normal_ratios <- function(tis) {
    # whole-chromosome maps; straw bins are unordered: tile index along the chromosome for every map bin, and which tiles have a bin
    maps <- lapply(setNames(chrs, chrs), function(ch) {
        print(paste('reading in', tis, 'HiC data,', ch))
        xT_map = straw(normHiCs[[tis]], gr = GRanges(ch, IRanges(1, chr_len_hg38[[ch]])), res = binsize)
        gr_map_dt <- gr2dt(xT_map$gr)
        # straw appends the query ':start:end' to seqnames; strip whatever suffix it is
        gr_map_dt[, seqnames := sub(':.*$', '', seqnames), by = seqnames]
        xT_map$gr <- dt2gr(gr_map_dt)
        t_map <- start(xT_map$gr) %/% binsize + 1
        n_t <- max(t_map)
        list(map = xT_map, t_map = t_map, n_t = n_t, present_t = seq_len(n_t) %in% t_map)
    })
    # within/outside sums of map bin k of chromosome ch for the region in_region (logical per tile); partners <= Dmax away
    bin_sums <- function(k, ch, in_region) {
        m <- maps[[ch]]
        # virtual 4C of bin k
        xT_map_4C <- m$map[seq_along(m$map$gr) == k, drop = T]
        p <- start(xT_map_4C) %/% binsize + 1
        is_self <- p == m$t_map[k]
        is_near <- abs(p - m$t_map[k]) * binsize <= Dmax
        data.table(sample = paste(tis, 'normal'), tissue = tis, seqnames = ch, start = start(m$map$gr)[k], end = end(m$map$gr)[k],
            self = sum(xT_map_4C$value[is_self], na.rm = T),
            within = sum(xT_map_4C$value[in_region[p] & !is_self], na.rm = T),
            outside = sum(xT_map_4C$value[!in_region[p] & is_near], na.rm = T),
            n_within = sum(in_region[p] & !is_self), n_outside = sum(!in_region[p] & is_near))
    }
    # bins are assigned to the TAD by their midpoint so the 0-based bedpe start can't pull in the neighbouring bin
    mx <- maps[[xchr]]
    TAD_bins <- which(gr.mid(mx$map$gr) %^% TAD_oi)
    in_TAD_t <- seq_len(mx$n_t) %in% mx$t_map[TAD_bins]
    print(paste(tis, ':', length(TAD_bins), 'bins in TAD'))
    TAD_res <- rbindlist(lapply(TAD_bins, bin_sums, ch = xchr, in_region = in_TAD_t))

    # random windows of the TAD's size (same seed for every normal and set): no missing tile in the window and no TAD tile within window +-
    # Dmax (randomwin) or within the window (randomwin_noampfilter); candidates of the chromosomes pooled, drawn in random order, kept if not
    # overlapping a kept window of the same chromosome
    win_n <- length(TAD_bins)
    D_tiles <- Dmax / binsize
    rand_res <- list()
    for (set in win_sets) {
        cand <- rbindlist(lapply(chrs, function(ch) {
            m <- maps[[ch]]
            TAD_t <- if (ch == xchr) in_TAD_t else rep(FALSE, m$n_t)
            TAD_cs <- c(0, cumsum(TAD_t))
            miss_cs <- c(0, cumsum(!m$present_t))
            s <- (D_tiles + 1):(m$n_t - win_n + 1 - D_tiles)
            ok <- miss_cs[s + win_n] - miss_cs[s] == 0
            ok <- ok & if (set == 'randomwin') TAD_cs[s + win_n + D_tiles] - TAD_cs[s - D_tiles] == 0 else TAD_cs[s + win_n] - TAD_cs[s] == 0
            # no candidate on this chromosome: no row (data.table() would recycle chr and fill s with NA)
            if (!any(ok)) return(NULL)
            data.table(chr = ch, s = s[ok])
        }))
        set.seed(20261018)
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
        win_dt[, window := seq_len(.N)]
        print(paste(tis, ':', set, ':', nrow(cand), 'candidate window starts (', paste(chrs, sapply(chrs, function(ch) sum(cand$chr == ch)), collapse = ', '),
            '),', nrow(win_dt), 'random windows of', win_n, 'bins (', paste(chrs, sapply(chrs, function(ch) sum(win_dt$chr == ch)), collapse = ', '), ')'))
        rand_res[[set]] <- rbindlist(lapply(seq_len(nrow(win_dt)), function(w) {
            ch <- win_dt$chr[w]
            m <- maps[[ch]]
            win_t <- win_dt$s[w] + 0:(win_n - 1)
            win_bins <- match(win_t, m$t_map)
            win_str <- paste0(ch, ':', start(m$map$gr)[win_bins[1]], '-', end(m$map$gr)[win_bins[win_n]])
            data.table(window = win_dt$window[w], window_range = win_str, window_n_amp_tiles = 0L,
                rbindlist(lapply(win_bins, bin_sums, ch = ch, in_region = seq_len(m$n_t) %in% win_t)))
        }))
    }
    print(paste(tis, ': done'))
    list(TAD = TAD_res, rand = rand_res)
}
# normals: per-bin tables of a previous run of this gene in outdir are reused (same TAD, windows and seed), else computed and written
norm_files <- c(TAD = paste0(outdir, '20261024_', xG, '_normLungTAD_inout_ratio_normHiC_', Dmax / 1e6, 'Mb.tsv'),
    setNames(paste0(outdir, '20261024_', xG, '_', win_sets, '_inout_ratio_normHiC_', Dmax / 1e6, 'Mb.tsv'), win_sets))
if (all(file.exists(norm_files))) {
    print(paste('normals: reusing', paste(basename(norm_files), collapse = ', ')))
    norm_dt <- fread(norm_files[['TAD']])
    norm_rand <- lapply(setNames(win_sets, win_sets), function(set) fread(norm_files[[set]]))
} else {
    norm_res <- mclapply(names(normHiCs), normal_ratios, mc.cores = length(normHiCs))
    names(norm_res) <- names(normHiCs)
    failed <- names(norm_res)[sapply(norm_res, function(r) inherits(r, 'try-error') || is.null(r))]
    if (length(failed) > 0) stop(paste('normal(s) failed:', paste(failed, collapse = ', '), '-', paste(sapply(norm_res[failed], as.character), collapse = '; ')))
    norm_dt <- add_ratios(rbindlist(lapply(norm_res, `[[`, 'TAD')))
    norm_rand <- lapply(setNames(win_sets, win_sets), function(set) add_ratios(rbindlist(lapply(norm_res, function(r) r$rand[[set]]))))
    rm(norm_res)
    fwrite(norm_dt, norm_files[['TAD']], sep = '\t')
    for (set in win_sets) fwrite(norm_rand[[set]], norm_files[[set]], sep = '\t')
}
sample_lv <- c(paste(names(normHiCs), 'normal'), tumors)

print('reading in tumor ratios')
read_tumor <- function(sampleid, suffix) {
    f <- paste0(tumor_dirs, '20261024_', sampleid, '_', xG, suffix)
    f <- f[file.exists(f)]
    if (length(f) == 0) {
        warning(paste('missing', sampleid, suffix, 'in', paste(tumor_dirs, collapse = ', ')))
        return(NULL)
    }
    print(paste('tumor', sampleid, ':', f[1]))
    fread(f[1])
}
TAD_dt <- rbind(norm_dt, rbindlist(lapply(tumors, read_tumor, suffix = paste0('_lungTAD_hs1_inout_ratio_', norm_tag, '_CNgt', cn_cutoff, '.tsv'))), fill = TRUE)
TAD_dt[, sample := factor(sample, levels = sample_lv)]
fwrite(TAD_dt, paste0(outdir, '20261024_', xG, '_lungTAD_inout_ratio_normals', Dmax / 1e6, 'Mb_vs_tumors_', norm_tag, '_CNgt', cn_cutoff, '.tsv'), sep = '\t')
print('lung TAD, median per sample:')
print(TAD_dt[, c(.(n_bins = .N), lapply(.SD, function(x) median(x[is.finite(x)]))), keyby = sample, .SDcols = ratio_cols])

# one value per random window: median of its per-bin ratios; amplified window: a tile with JaBbA CN > cn_cutoff
rand_med <- lapply(setNames(win_sets, win_sets), function(set) {
    rand_dt <- rbind(norm_rand[[set]], rbindlist(lapply(tumors, read_tumor, suffix = paste0('_', set, '_hs1_inout_ratio_', norm_tag, '.tsv'))), fill = TRUE)
    rmd <- rand_dt[, c(.(n_bins = .N, n_amp_tiles = window_n_amp_tiles[1], max_jabba_cn = as.numeric(suppressWarnings(max(jabba_cn, na.rm = TRUE))),
        median_cn_model = median(cn_model)), lapply(.SD, function(x) median(x[is.finite(x)]))), by = .(sample, tissue, seqnames, window, window_range),
        .SDcols = ratio_cols]
    rmd[!is.finite(max_jabba_cn), max_jabba_cn := NA_real_]
    rmd[, amplified := n_amp_tiles > 0]
    rmd[, sample := factor(sample, levels = sample_lv)]
    fwrite(rmd, paste0(outdir, '20261024_', xG, '_', set, '_inout_ratio_median_per_window_', norm_tag, '.tsv'), sep = '\t')
    print(paste(set, ': random windows per sample and chromosome (amplified windows):'))
    print(rmd[, .(n_windows = .N, n_amplified = sum(amplified)), keyby = .(sample, seqnames)])
    print(paste(set, ': random windows, median of the window medians per sample:'))
    print(rmd[, c(.(n_windows = .N), lapply(.SD, function(x) median(x[is.finite(x)]))), keyby = sample, .SDcols = ratio_cols])
    rmd
})

# permutation test per sample, window set and in/out ratio: how often does a random window of the sample (median of its bins) reach the
# ratio of the lung TAD of the sample (median of its bins)? empirical p = (1 + #windows >= TAD) / (1 + #windows), one-sided
perm_dt <- rbindlist(lapply(win_sets, function(set) rbindlist(lapply(sample_lv[sample_lv %in% as.character(TAD_dt$sample)], function(s) rbindlist(lapply(ratio_cols, function(rc) {
    v_TAD <- TAD_dt[sample == s][[rc]]
    v_TAD <- v_TAD[is.finite(v_TAD)]
    rs <- rand_med[[set]][sample == s & is.finite(get(rc))]
    v_rand <- rs[[rc]]
    obs <- if (length(v_TAD) > 0) median(v_TAD) else NA_real_
    n_ge <- if (is.finite(obs)) sum(v_rand >= obs) else NA_integer_
    data.table(window_set = set, sample = s, tissue = TAD_dt[sample == s, tissue[1]], ratio_type = rc, n_TAD_bins = length(v_TAD), TAD_median_ratio = obs,
        n_random_windows = length(v_rand), random_windows_median_ratio = if (length(v_rand) > 0) median(v_rand) else NA_real_,
        n_random_windows_ge_TAD = n_ge, p_value = if (length(v_rand) > 0 && is.finite(obs)) (1 + n_ge) / (1 + length(v_rand)) else NA_real_,
        n_amplified_windows = sum(rs$amplified), n_amplified_windows_ge_TAD = if (is.finite(obs)) sum(rs$amplified & v_rand >= obs) else NA_integer_,
        amplified_windows_median_ratio = if (any(rs$amplified)) median(v_rand[rs$amplified]) else NA_real_)
}))))))
# Benjamini-Hochberg q-values over the samples of the gene, per window set and ratio type
perm_dt[, q_value := p.adjust(p_value, method = 'BH'), by = .(window_set, ratio_type)]
setcolorder(perm_dt, c(setdiff(names(perm_dt), c('n_amplified_windows', 'n_amplified_windows_ge_TAD', 'amplified_windows_median_ratio')),
    'n_amplified_windows', 'n_amplified_windows_ge_TAD', 'amplified_windows_median_ratio'))
fwrite(perm_dt, paste0(outdir, '20261024_', xG, '_lungTAD_inout_ratio_permtest_normals', Dmax / 1e6, 'Mb_vs_tumors_', norm_tag, '_CNgt', cn_cutoff, '.csv'))
print('permutation test, lung TAD vs random windows of the same sample:')
print(perm_dt[, .(window_set, sample, ratio_type, TAD_median_ratio = signif(TAD_median_ratio, 3), random_windows_median_ratio = signif(random_windows_median_ratio, 3),
    n_random_windows, n_random_windows_ge_TAD, p_value = signif(p_value, 3), q_value = signif(q_value, 3), n_amplified_windows, n_amplified_windows_ge_TAD)])

# red/blue as the BRCA/NSCLC tracks; normal TAD boxes filled, tumor TAD boxes open, random windows grey (amplified windows orange)
tis_cols <- c(breast = 'red', esophagus = 'blue', lung = 'blue')
plot_ratio <- function(TAD_dt, rmd, set, yfield, ylab) {
    samples_in <- sample_lv[sample_lv %in% c(as.character(TAD_dt$sample), as.character(rmd$sample))]
    vals <- list()
    amp_flags <- list()
    for (s in samples_in) {
        v_TAD <- TAD_dt[sample == s][[yfield]]
        rs <- rmd[sample == s & is.finite(get(yfield))]
        vals <- c(vals, list(v_TAD[is.finite(v_TAD)], rs[[yfield]]))
        amp_flags <- c(amp_flags, list(rep(FALSE, sum(is.finite(v_TAD))), rs$amplified))
    }
    tis_s <- sapply(samples_in, function(s) c(TAD_dt[sample == s, tissue], rmd[sample == s, tissue])[1])
    at <- as.vector(t(outer((seq_along(samples_in) - 1) * 3, 1:2, '+')))
    border <- as.vector(rbind(tis_cols[tis_s], 'grey40'))
    fill <- as.vector(rbind(ifelse(grepl('normal', samples_in), adjustcolor(tis_cols[tis_s], 0.15), 'white'), 'grey90'))
    # headroom for the significance stars
    y_rng <- range(unlist(vals))
    boxplot(vals, at = at, outline = F, col = fill, border = border, lwd = 1.5, xaxt = 'n', las = 1,
        xlim = c(0.5, max(at) + 0.5), ylim = y_rng + c(0, 0.08 * diff(y_rng)), xlab = '', ylab = ylab,
        main = paste0(xG, ' normal lung TAD ', TAD_str, ' (hg38; lifted to hs1 for tumors) vs random windows (', set, ') of the tested TAD\'s size on ',
            paste(chrs, collapse = ', '), ' (n per box above)\n', binsize / 1e3, 'kb bins, partners <= ', Dmax / 1e6, 'Mb; tumors: counts / CN effect of the ',
            model, ' HiC coverage CN model (chr7/8/11 fit up to ', Dmax / 1e6, 'Mb), TAD tiles with JaBbA CN > ', cn_cutoff), cex.main = 0.75)
    # one dot per TAD bin / per random window (orange: random window with a tile of JaBbA CN > cn_cutoff)
    set.seed(1)
    for (g in seq_along(vals)) {
        if (length(vals[[g]]) == 0) next
        dot_col <- if (border[g] == 'grey40') ifelse(amp_flags[[g]], 'darkorange', 'grey30') else border[g]
        points(at[g] + runif(length(vals[[g]]), -0.15, 0.15), vals[[g]], pch = 16, cex = 0.6, col = adjustcolor(dot_col, 0.6))
    }
    axis(1, at = (seq_along(samples_in) - 1) * 3 + 1.5, labels = samples_in, las = 2, tick = F)
    mtext(paste0('n=', sapply(vals, length)), side = 3, at = at, line = 0.2, cex = 0.5)
    abline(v = 3 * sum(grepl('normal', samples_in)), lty = 2, col = 'grey60')
    # permutation test of the sample (TAD vs its random windows): black star q < 0.1, grey star p < 0.05
    pv <- perm_dt[window_set == set & ratio_type == yfield]
    for (k in seq_along(samples_in)) {
        r <- pv[sample == samples_in[k]]
        if (nrow(r) == 0 || !is.finite(r$p_value)) next
        star_col <- if (r$q_value < 0.1) 'black' else if (r$p_value < 0.05) 'grey55' else NA
        if (!is.na(star_col)) points((k - 1) * 3 + 1.5, y_rng[2] + 0.05 * diff(y_rng), pch = 8, cex = 0.8, lwd = 1.5, col = star_col)
    }
    legend('topright', legend = c(paste(xG, 'TAD: one dot per 5kb bin'), 'random windows: one dot per window (median of its bins)',
        paste0('random window with a tile of JaBbA CN > ', cn_cutoff), 'TAD vs random windows: q < 0.1 (permutation, BH over the samples)',
        'TAD vs random windows: p < 0.05 (permutation)'), fill = c('white', 'grey90', NA, NA, NA), border = c('black', 'grey40', NA, NA, NA),
        pch = c(NA, NA, 16, 8, 8), col = c(NA, NA, 'darkorange', 'black', 'grey55'), cex = 0.7, bg = adjustcolor('white', 0.8), box.col = NA)
}

print('plotting')
pdf(file = paste0(outdir, '20261024_', xG, '_lungTAD_inout_ratio_normals', Dmax / 1e6, 'Mb_vs_tumors_', norm_tag, '_CNgt', cn_cutoff, '.pdf'), width = 11, height = 7)
par(mar = c(10, 6.5, 5, 1), mgp = c(4, 1, 0))
for (set in win_sets) {
    plot_ratio(TAD_dt, rand_med[[set]], set, 'ratio', 'within-TAD / outside-TAD contacts\n(self-contact excluded)')
    plot_ratio(TAD_dt, rand_med[[set]], set, 'ratio_incl_self', 'within-TAD / outside-TAD contacts\n(self-contact included)')
    plot_ratio(TAD_dt, rand_med[[set]], set, 'ratio_lennorm', 'within-TAD / outside-TAD contacts per bin\n(length normalized, self-contact excluded)')
    plot_ratio(TAD_dt, rand_med[[set]], set, 'ratio_lennorm_incl_self', 'within-TAD / outside-TAD contacts per bin\n(length normalized, self-contact included)')
}
dev.off()

# normals vs tumors: per sample, lung TAD median / median of the random window medians (ratio_lennorm_incl_self); three normal tissues vs all
# tumors, two-sided Wilcoxon rank-sum (Mann-Whitney U) test per window set
enr_rc <- 'ratio_lennorm_incl_self'
enr_prefix <- paste0(outdir, '20261024_', xG, '_lungTAD_over_randomwin_', enr_rc, '_normals_vs_tumors')
enr_dt <- perm_dt[ratio_type == enr_rc & is.finite(TAD_median_ratio) & is.finite(random_windows_median_ratio),
    .(window_set, sample, tissue, group = ifelse(grepl('normal', sample), 'normal', 'tumor'), TAD_median_ratio, random_windows_median_ratio,
    n_random_windows, TAD_over_random = TAD_median_ratio / random_windows_median_ratio)]
fwrite(enr_dt, paste0(enr_prefix, '_values.tsv'), sep = '\t')
wil_dt <- rbindlist(lapply(win_sets, function(set) {
    v_n <- enr_dt[window_set == set & group == 'normal', TAD_over_random]
    v_t <- enr_dt[window_set == set & group == 'tumor', TAD_over_random]
    wt <- wilcox.test(v_t, v_n, alternative = 'two.sided')
    data.table(gene = xG, window_set = set, ratio_type = enr_rc, value = 'lung TAD median / median of the random window medians, per sample',
        n_normal = length(v_n), n_tumor = length(v_t), median_normal = median(v_n), median_tumor = median(v_t), W_tumor_gt_normal = unname(wt$statistic),
        p_value = wt$p.value, test = paste(wt$method, '(two-sided)'))
}))
fwrite(wil_dt, paste0(enr_prefix, '_wilcox.csv'))
print('normals vs tumors, lung TAD / random windows (ratio_lennorm_incl_self):')
print(enr_dt[, .(window_set, sample, group, TAD_over_random = signif(TAD_over_random, 3))])
print(wil_dt)
pdf(file = paste0(enr_prefix, '_wilcox.pdf'), width = 5.5, height = 6.5)
par(mar = c(4, 6.5, 5, 1), mgp = c(4, 1, 0))
for (set in win_sets) {
    e <- enr_dt[window_set == set]
    w <- wil_dt[window_set == set]
    e[, xpos := ifelse(group == 'normal', 1, 2)]
    boxplot(TAD_over_random ~ factor(group, levels = c('normal', 'tumor')), data = e, outline = F, col = 'white', border = 'grey40', lwd = 1.5,
        xaxt = 'n', las = 1, ylim = range(c(1, e$TAD_over_random)), xlab = '', ylab = paste0('lung TAD / random windows (', set, ')\n(', enr_rc, ', medians)'),
        main = paste0(xG, ' normal lung TAD vs random windows (', set, '), partners <= ', Dmax / 1e6, 'Mb\nWilcoxon rank-sum, normal tissues vs tumors: W = ',
            w$W_tumor_gt_normal, ', p = ', signif(w$p_value, 3)), cex.main = 0.75)
    axis(1, at = 1:2, labels = c(paste0('normal tissues (n=', w$n_normal, ')'), paste0('tumors (n=', w$n_tumor, ')')), tick = F)
    abline(h = 1, lty = 2, col = 'grey60')
    set.seed(1)
    xj <- e$xpos + runif(nrow(e), -0.15, 0.15)
    points(xj, e$TAD_over_random, pch = 16, cex = 1.1, col = adjustcolor(tis_cols[e$tissue], 0.8))
    text(xj, e$TAD_over_random, sub(' normal', '', e$sample), pos = 4, cex = 0.55, col = 'grey30')
    legend('topright', legend = c('breast', 'lung / esophagus'), pch = 16, col = c('red', 'blue'), cex = 0.75, bg = adjustcolor('white', 0.8), box.col = NA)
}
dev.off()
print('Done')
