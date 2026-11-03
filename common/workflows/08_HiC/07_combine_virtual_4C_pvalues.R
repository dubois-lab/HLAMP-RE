# Combines the virtual 4C p-values of the samples of one group per bin (Fisher's method), with q-values across the chromosome.
# Run from this folder: Rscript 07_combine_virtual_4C_pvalues.R <gene> <group> <sample> <sample> [<sample> ...]
# Writes comb_Pvals_v4C_<gene>_<group>_gr.rds

# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'results')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'results')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(gUtils)
library(poolr)
library(qvalue)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop('Usage: Rscript 07_combine_virtual_4C_pvalues.R <gene> <group> <sample> <sample> [<sample> ...]')
xG <- args[1]
xGroup <- args[2]
samples <- args[-(1:2)]

# The tracks of one gene share the same bins
v4C_dts <- lapply(samples, function(sid) gr2dt(readRDS(file.path(IN_DIR, paste0(sid, '_covCovariate_v4C_', xG, '_nlg_p.rds')))))
comb_pvals <- v4C_dts[[1]][, .(seqnames, start, end)]
p_cols <- paste0('p_', samples)
for (n in seq_along(samples)) {
  stopifnot(identical(v4C_dts[[n]]$start, comb_pvals$start))
  comb_pvals[, (p_cols[n]) := 10^-v4C_dts[[n]]$value]
}

comb_pvals[, combp := fisher(unlist(.BY))$p, by = p_cols]
comb_pvals[, comb_q := qvalue(combp)$qvalue]
comb_pvals[, nlg_combp := -log10(combp)]
comb_pvals[, nlg_comb_q := -log10(comb_q)]
saveRDS(dt2gr(comb_pvals), file.path(OUT_DIR, paste0('comb_Pvals_v4C_', xG, '_', xGroup, '_gr.rds')))
