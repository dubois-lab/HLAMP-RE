## Fig S6: per-tumor plots of the NSD3 and EGFR HLAMP tumors (CN-normalized and raw HiC matrices with tracks) with DMR tracks:
## - plot windows: gene midpoint +-0.5/1/1.5Mb (hs1); TSS of the gene directly from the T2T/hs1 TSS rds (NSD3_TSSgr_T2T.rds, EGFR_TSS_T2T.rds;
##   no liftOver from hg38), first TSS
## - DMRs (hs1) with higher methylation in BRCA (dmrs_d_BRCAgt) for the BRCA tumors, with higher methylation in NSCLC (dmrs_d_NSCLCgt) for
##   the LUSC tumors
## - the CN-normalized v4C of the TSS saved as rds (list: v4C GRanges of counts / CN effect per tile, gt as plotted, y1, TSS, TAD, ...)
## Pages 1/2 window_oi, 3/4 window_wide, 5/6 window_wide2; CN-normalized matrix with its HiC coverage and the v4C of the gene TSS, then the raw
## count matrix and raw v4C; HiC colour scale and v4C y axis 0 to 4x the median v4C over the lifted lung TAD; counts / CN effect of the logistic
## 20Mb background model (200 x 200 cells)
## - contact map: only the tile pairs on the gene's chromosome with a tile in the plot region (widest plot window, lifted lung TAD and TSS) and
##   partners up to Dmax (raw 5kb map of the nf-core run or .hic file); HiC coverage, CN of the model and tile QC from the CN model
## - check against the TAD ratio run: TAD tiles, HiC coverage, CN of the model and kept TAD tiles vs
##   <IN_DIR>/TADratio/<gene>/20261024_<id>_<gene>_lungTAD_hs1_tiles_jabbaCN.tsv
## - output: <outdir>/20261028_<id>_<gene>_hiccov_logisticnorm_matrix_v4C_lungTAD_hs1_CNgt9_DMR.pdf and
##   <outdir>/20261028_<id>_<gene>_hiccov_logisticnorm_v4C_TSS_hs1.rds
## run: FigS6_job-script_NSD3_EGFR_tumor_plots.sh; args: sampleid ampChr tissue hic_dir gene outdir cn_model_file [map_id] (tissue: breast =
##   BRCA, lung = LUSC)
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')

library(GxG)
library(gGnome)
library(rtracklayer)

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop('usage: Rscript <script> sampleid ampChr tissue hic_dir gene outdir cn_model_file [map_id]', call. = FALSE)
sampleid <- args[1]
ampChr <- args[2]
tissue <- args[3]
hic_dir <- args[4]
xG <- args[5]
outdir <- sub('/*$', '/', args[6])
cn_model_file <- args[7]
map_id <- if (length(args) >= 8) args[8] else sampleid
stopifnot(xG %in% c('NSD3', 'EGFR'), tissue %in% c('breast', 'lung'))
model <- 'logistic'
print(paste0('Plots: ', sampleid, ' (', hic_dir, ', map ', map_id, '); gene ', xG, '; ', tissue, '; CN model ', cn_model_file))
amp <- c(chr8 = 'FGFR1', chr11 = 'CCND1')[[ampChr]]
norm_tag <- paste0('hiccov_', model, 'norm')
binsize = 5e3
Dmax = 2e7
cn_cutoff = 9
min_hiccov_frac = 0.1
clim_factor = 4
full_run_dir <- paste0(IN_DIR, '/TADratio/', xG, '/')

# white y.grid for all data tracks
fmt_gt <- function(gt) {
    gt$y.grid.col= "white"
    gt
}

hg38ToHs1_chain <- import.chain(file.path(IN_DIR, 'hg38ToHs1.over.chain'))
gencode_T2T_gt <- readRDS(file.path(IN_DIR, 'final_T2T_gencode_gt_USCS.rds'))
CGC_gr_T2T <- readRDS('../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')
CGC_gr_hg38 <- readRDS('../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
Oncgr <- gr.chr(CGC_gr_T2T[mcols(CGC_gr_T2T)$Gene_Symbol == xG])
xchr <- as.character(seqnames(Oncgr))[1]
# plot windows (hs1); TSS of the gene directly from the T2T/hs1 rds (no liftOver), first TSS
xG_mid <- round((min(start(Oncgr)) + max(end(Oncgr))) / 2)
window_oi <- GRanges(xchr, IRanges(xG_mid - 5e5, xG_mid + 5e5))
window_wide <- GRanges(xchr, IRanges(xG_mid - 1e6, xG_mid + 1e6))
window_wide2 <- GRanges(xchr, IRanges(xG_mid - 1.5e6, xG_mid + 1.5e6))
TSS_file <- c(NSD3 = file.path(IN_DIR, 'NSD3_TSSgr_T2T.rds'), EGFR = file.path(IN_DIR, 'EGFR_TSS_T2T.rds'))[[xG]]
xG_TSS_all <- readRDS(TSS_file)
if (!inherits(xG_TSS_all, 'GRanges')) xG_TSS_all <- dt2gr(as.data.table(xG_TSS_all))
xG_TSS_all <- gr.chr(xG_TSS_all)
print(paste0(length(xG_TSS_all), ' TSS(s) in ', basename(TSS_file), ': ', paste(gr.string(xG_TSS_all), collapse = ', ')))
xG_TSSgr <- xG_TSS_all[1]
stopifnot(as.character(seqnames(xG_TSSgr)) == xchr)

# DMR tracks (hs1): BRCA tumors the DMRs with higher methylation in BRCA, LUSC tumors those with
# higher methylation in NSCLC
xD <- c(NSD3 = 0.3, EGFR = 0.4)[[xG]]
dmrs_d <- fread(paste0(IN_DIR, "/", c(NSD3 = '20260615_NSD3_TM_dmrs.csv', EGFR = '20251003_EGFR_TM_dmrs_d_delt_0.4.csv')[[xG]]))
colnames(dmrs_d)[6:7] <- c("meanMethy_NSCLC", "meanMethy_BRCA")
dmrs_d_gr <- gr.chr(dt2gr(dmrs_d))
dmrs_d <- as.data.table(dmrs_d_gr)
# each DMR is shaded by the highest |diff.Methy| cutoff (xD / xD+0.1 / xD+0.2) it passes;
# categorical colormap as for the tissue-colored Fisher track (a continuous per-row grey() color column does not render in gTrack)
dmrs_d[, methyl_level := ifelse(abs(diff.Methy) > xD + 0.2, "dhi", ifelse(abs(diff.Methy) > xD + 0.1, "dmid", "dlo"))]
colormap_dmrs_d_grey <- list(methyl_level = c(dlo = grey(0.7), dmid = grey(0.4), dhi = grey(0.1)))
dmrs_d_NSCLCgt <- gTrack(dt2gr(dmrs_d[diff.Methy < -xD]), gr.colorfield = "color", colormaps = colormap_dmrs_d_grey)
dmrs_d_BRCAgt <- gTrack(dt2gr(dmrs_d[diff.Methy > xD]), gr.colorfield = "color", colormaps = colormap_dmrs_d_grey)
dmr_gt <- if (tissue == 'breast') dmrs_d_BRCAgt else dmrs_d_NSCLCgt
print(paste0('DMRs (|diff.Methy| > ', xD, '): ', nrow(dmrs_d[diff.Methy > xD]), ' higher in BRCA, ', nrow(dmrs_d[diff.Methy < -xD]), ' higher in NSCLC (',
    nrow(dmrs_d[seqnames == xchr & diff.Methy > xD]), ' / ', nrow(dmrs_d[seqnames == xchr & diff.Methy < -xD]), ' on ', xchr, '); ',
    if (tissue == 'breast') 'BRCA' else 'NSCLC', ' track for this tumor'))
rm(dmrs_d_gr)

# normal lung TADs lifted to hs1
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
Lung_norm_16f_TAD_1_hs1_gt <- gTrack(data = lift_hs1(Lung_TAD_gr[Lung_TAD_gr %^% gr.chr(CGC_gr_hg38)]), name = 'Lung_norm_16f_TAD_1_hs1', col = 'grey', height = 0.5)

# 5kb tiles of xchr (bins of the raw contact map; .hic: hs1 5kb bins of the nf-core run), tile k ends at k*binsize
is_hic <- grepl('\\.hic$', hic_dir)
bins_file <- if (is_hic) file.path(IN_DIR, 'hic/nfcore_results/contact_maps/bins/cooler_bins_5000.bed') else paste0(hic_dir, '/contact_maps/bins/cooler_bins_5000.bed')
cbins <- fread(bins_file, header = FALSE, select = 1:3, col.names = c('chr', 'start', 'end'))
idx <- which(cbins$chr == xchr)
stopifnot(length(idx) > 0, all(diff(idx) == 1))
xT_gr <- GRanges(xchr, IRanges(cbins$start[idx] + 1, cbins$end[idx]))
stopifnot(all(ceiling(end(xT_gr) / binsize) == seq_along(xT_gr)))
n_t <- length(xT_gr)
chr_len <- cbins$end[idx[n_t]]
rm(cbins)
D_tiles <- Dmax / binsize
tss_t <- which(xT_gr %^% xG_TSSgr)
print(paste0('TSS (hs1) ', gr.string(xG_TSSgr), ': tile(s) ', paste(tss_t, collapse = ', '), ' (', paste(gr.string(xT_gr[tss_t]), collapse = ', '), ')'))
stopifnot(length(tss_t) == 1)
# plot region: widest window, lifted TAD and TSS; its tiles and those of the widest window
region <- GRanges(xchr, IRanges(min(start(window_wide2), start(TAD_oi), start(xG_TSSgr)), max(end(window_wide2), end(TAD_oi), end(xG_TSSgr))))
reg_t <- which(xT_gr %^% region)
wide_t <- which(xT_gr %^% window_wide2)
print(paste0('plot region (hs1): ', gr.string(region), ' (', length(reg_t), ' tiles); tile pairs with a region tile and partners up to ', Dmax / 1e6, 'Mb'))

# tile pairs (local tile index, i <= j) with a region tile and j - i <= D_tiles
if (is_hic) {
    q_lo <- max(1, reg_t[1] - D_tiles)
    q_hi <- min(n_t, reg_t[length(reg_t)] + D_tiles)
    m <- straw(hic_dir, gr = GRanges(xchr, IRanges((q_lo - 1) * binsize + 1, min(chr_len, q_hi * binsize))), res = binsize, norm = 'NONE')
    t_map <- as.integer(start(m$gr) %/% binsize + 1)
    d <- m$dat
    it <- data.table(i = pmin(t_map[d$i], t_map[d$j]), j = pmax(t_map[d$i], t_map[d$j]), counts = as.integer(round(d$value)))
    rm(m, d)
} else {
    g0 <- idx[1] - 1L
    awk_prog <- sprintf('$1 > %d && $2 <= %d && $2 - $1 <= %d && (($1 > %d && $1 <= %d) || ($2 > %d && $2 <= %d)) { print $1 - %d, $2 - %d, $3 }',
        g0, g0 + n_t, D_tiles, g0 + reg_t[1] - 1L, g0 + reg_t[length(reg_t)], g0 + reg_t[1] - 1L, g0 + reg_t[length(reg_t)], g0, g0)
    it <- fread(cmd = paste0("awk '", awk_prog, "' ", hic_dir, '/contact_maps/txt/', map_id, '.5000_raw.txt'), col.names = c('i', 'j', 'counts'))
}
it <- it[j - i <= D_tiles & (i %in% reg_t | j %in% reg_t)]
stopifnot(it[, all(i <= j & i >= 1 & j <= n_t)])
print(paste(nrow(it), 'nonzero tile pairs,', sum(as.numeric(it$counts)), 'counts'))

# CN model: CN per tile (HiC coverage / median QC tile, matched by tile midpoint) and CN effect; HiC coverage and tile QC as in the TAD ratio script
cn_model <- readRDS(cn_model_file)
stopifnot(Dmax <= cn_model$Dmax)
th <- cn_model$fits[[model]]$par[cn_model$k + seq_len(cn_model$forms[[model]]$n)]
cn_effect <- function(ld, lx) exp(cn_model$forms[[model]]$eta(th, cn_model$design(ld, lx))$eta)
m_idx <- gr.match(gr.mid(xT_gr), cn_model$tile_gr)
total <- cn_model$hic_cov[m_idx]
total[is.na(total)] <- 0
cn_t <- (cn_model$hic_cov / cn_model$hiccov_median)[m_idx]
valid <- total > 0 & is.finite(cn_t) & cn_t > 0
good <- valid & total > min_hiccov_frac * median(total[total > 0])
jabba_res_gg = gG(jabba=paste0(IN_DIR, '/jabba/',amp,'/',tissue,'/',sampleid,'/jabba.gg.rds'))
jabba_nodes_gr <- gr.chr(gr.stripstrand(jabba_res_gg$nodes$gr))
tile_cn <- jabba_nodes_gr$cn[gr.match(gr.mid(xT_gr), jabba_nodes_gr)]
print(paste0('TSS tile: CN of the model ', round(cn_t[tss_t], 2), ', JaBbA CN ', tile_cn[tss_t]))

# normalized pairs (both tiles with a CN of the model) and raw counts of the widest window
raw_wide <- it[i %in% wide_t & j %in% wide_t, .(i, j, value = counts)]
pairs <- it[valid[i] & valid[j]]
pairs[, value := counts / cn_effect(log10(pmax(abs(j - i), 1) * binsize), log10(cn_t[i] * cn_t[j]))]
pairs <- pairs[, .(i, j, value)]
rm(it)

# lung TAD tiles kept (tile QC and JaBbA CN > cn_cutoff); check against the TAD ratio run
TAD_tiles <- which(gr.mid(xT_gr) %^% TAD_oi)
keep_tile <- good & !is.na(tile_cn) & tile_cn > cn_cutoff
TAD_kept <- TAD_tiles[keep_tile[TAD_tiles]]
print(paste(length(TAD_tiles), 'tiles in TAD,', length(TAD_kept), 'kept with JaBbA CN >', cn_cutoff, 'and passing tile QC'))
full_tsv <- paste0(full_run_dir, '20261024_', sampleid, '_', xG, '_lungTAD_hs1_tiles_jabbaCN.tsv')
if (file.exists(full_tsv)) {
    ft <- fread(full_tsv)
    same <- nrow(ft) == length(TAD_tiles) && all(ft$start == start(xT_gr)[TAD_tiles]) && isTRUE(all.equal(ft$hic_cov, total[TAD_tiles])) &&
        isTRUE(all.equal(ft$cn_model, cn_t[TAD_tiles])) && all(ft$kept == keep_tile[TAD_tiles])
    print(paste('check vs full run (', basename(full_tsv), '): TAD tiles, HiC coverage, CN of the model and kept tiles', if (same) 'identical' else 'DIFFER'))
    if (!same) warning(paste(sampleid, xG, ': TAD tiles differ from the full run'))
} else warning(paste('no full-run TAD tile table', full_tsv))

# diagnostic plots (with the DMR track; data tracks formatted with fmt_gt)
v4C_of <- function(gmat) gmat[gmat$gr %^% xG_TSSgr, drop = T]
clim_of <- function(v4C_gr, label) {
    clim_max <- clim_factor * median(v4C_gr$value[gr.mid(v4C_gr) %^% TAD_oi], na.rm = T)
    print(paste0('HiC colour scale max, ', label, ' (', clim_factor, 'x median ', xG, ' v4C over the TAD): ', clim_max))
    stopifnot(clim_max > 0)
    clim_max
}
gmat <- gM(gr = xT_gr, dat = pairs, fill = 0)
bedmeth_files <- paste0(IN_DIR, c('/bedmeth/', '/bedmethylRDS/'), sampleid, '.wf_mods.ungrouped.bedmethyl.rds')
mods.ungrouped.bedmethyl <- readRDS(bedmeth_files[file.exists(bedmeth_files)][1])
mods.ungrouped.bedmethyl_logmCe3_gt <- fmt_gt(gTrack(dt2gr(mods.ungrouped.bedmethyl[modified_base_code == 'm'& Nvalid_cov > 5]), y.field = 'neglog_mth', name = '-log5mC_%', bars = T, y1= 3))
rm(mods.ungrouped.bedmethyl)
# HiC coverage tracks over the plot region: raw (row sum over the chromosome) and normalized (row sum of counts / CN effect over the partners
# up to Dmax, self-contact once; tiles without a CN of the model 0)
cov_norm <- numeric(n_t)
s1 <- pairs[, .(t = sum(value)), by = i]
cov_norm[s1$i] <- s1$t
s2 <- pairs[i != j, .(t = sum(value)), by = j]
cov_norm[s2$j] <- cov_norm[s2$j] + s2$t
hic_cov_gr <- xT_gr[reg_t]
hic_cov_gr$hic_cov <- total[reg_t]
hic_cov_gr$hic_cov_norm <- cov_norm[reg_t]
hic_cov_gt <- fmt_gt(gTrack(data = hic_cov_gr, y.field = 'hic_cov', name = 'HiC_cov', bars = T, y0 = 0))
hic_cov_norm_gt <- fmt_gt(gTrack(data = hic_cov_gr, y.field = 'hic_cov_norm', name = paste0('HiC_cov_CNnorm_', Dmax / 1e6, 'Mb'), bars = T, y0 = 0))
TAD_gts <- Lung_norm_16f_TAD_1_hs1_gt
if (length(TAD_kept) > 0) TAD_gts <- c(TAD_gts, gTrack(data = reduce(xT_gr[TAD_kept]), name = paste0('TAD_CN>', cn_cutoff), col = 'black', height = 0.5))
# v4C y axis capped at the HiC colour scale max
v4C_xG_gr <- v4C_of(gmat)
clim_norm <- clim_of(v4C_xG_gr, 'normalized')
v4C_xG_gt <- fmt_gt(gTrack(data = v4C_xG_gr, y.field = 'value', bars = T, name = 'counts/CNeffect', y0 = 0, y1 = clim_norm, col = 'red'))
# CN-normalized v4C of the TSS (all tiles of xchr; tiles beyond Dmax of the TSS or without a CN of the model 0) for later use
v4C_rds <- paste0(outdir, '20261028_', sampleid, '_', xG, '_', norm_tag, '_v4C_TSS_hs1.rds')
saveRDS(list(v4C = v4C_xG_gr, gt = v4C_xG_gt, y1 = clim_norm, clim_factor = clim_factor, sample = sampleid, gene = xG, tissue = tissue,
    TSS_hs1 = xG_TSSgr, TSS_tile = xT_gr[tss_t], TAD_hs1 = TAD_oi, windows = c(window_oi, window_wide, window_wide2), Dmax = Dmax,
    cn_model_file = cn_model_file, hic = hic_dir, map_id = map_id), v4C_rds)
print(paste0('normalized v4C (', length(v4C_xG_gr), ' tiles, ', sum(v4C_xG_gr$value > 0), ' > 0) saved: ', v4C_rds))
gmat_raw <- gM(gr = xT_gr, dat = raw_wide, fill = 0)
v4C_raw_gr <- v4C_of(gmat_raw)
clim_raw <- clim_of(v4C_raw_gr, 'raw counts')
v4C_raw_gt <- fmt_gt(gTrack(data = v4C_raw_gr, y.field = 'value', bars = T, name = 'raw counts', y0 = 0, y1 = clim_raw, col = 'red'))
base_gts <- c(gencode_T2T_gt, dmr_gt, hic_cov_gt, fmt_gt(jabba_res_gg$gt), mods.ungrouped.bedmethyl_logmCe3_gt)
pdf(file = paste0(outdir, '20261028_', sampleid, '_', xG, '_', norm_tag, '_matrix_v4C_lungTAD_hs1_CNgt', cn_cutoff, '_DMR.pdf'), paper = 'a4', height = 24)
for (win in list(window_oi, window_wide, window_wide2)) {
    plot(c(base_gts, hic_cov_norm_gt, v4C_xG_gt, TAD_gts,
        gmat$gtrack(name = paste0('counts/(CN_i*CN_j)^gamma, ', model, ' HiCcov CN model (chr7/8/11)'), height = 90, clim = c(0, clim_norm))), win, cex.label = 0.1)
    plot(c(base_gts, v4C_raw_gt, TAD_gts,
        gmat_raw$gtrack(name = 'raw counts (raw 5kb contact map)', height = 90, clim = c(0, clim_raw))), win, cex.label = 0.1)
}
dev.off()
print(paste0('Done: ', sampleid, ' (', xG, ', ', tissue, ')'))
