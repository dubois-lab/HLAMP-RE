## =============================================================================
## FigS2A/B -- genome-wide amplification-fraction tracks
## CN-only vs JaBbA walks
## CCND1 / NSD3 / EGFR, BRCA + NSCLC(LUSC)
## -----------------------------------------------------------------------------
## Per gene, four stacked tracks (top to bottom): walks BRCA, walks NSCLC,
## CN-only BRCA, CN-only NSCLC. Bins overlapping the hs1 centromere/telomere
## blacklist are removed before plotting. chr_sizes/offset/genome_pos are not used by facet_grid(free_x),
## which plots each chromosome on its own local scale.
## =============================================================================
# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(GenomicRanges)
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(gGnome)

genes      <- c('CCND1', 'NSD3', 'EGFR')
date_tag   <- '20260910'
cohort_tag <- c(BRCA = 'BRCA', NSCLC = 'NSCLC_LUSC')

cn_dir    <- file.path(IN_DIR, 'CNdists')
walks_dir <- file.path(IN_DIR, 'jabbawlks')
out_dir   <- OUT_DIR

chr_levels <- paste0("chr", c(1:22, "X", "M"))

## hs1 chromosome sizes
chr_sizes <- fread('../../common/data/refGenome/hs1.chrom.sizes',
                    header = FALSE, col.names = c("seqnames", "chr_len"))
chr_sizes[, offset := shift(cumsum(as.numeric(chr_len)), fill = 0)]

## Factor here too: a character seqnames in any layer makes facet_grid fall back to lexical order.
chr_bounds <- chr_sizes[seqnames %in% chr_levels, .(seqnames = factor(seqnames, levels = chr_levels),
                                                    chr_len, xmin = 0, xmax = chr_len)]

## Centromere/telomere/acrocentric-p-arm mask: T2T satellite arrays have no usable
## depth, so JaBbA carries one segment's CN across them and they show up as tall
## artefactual bars genome-wide.
blacklist <- readRDS(file.path(cn_dir, 'blacklist/hs1_blacklist_centro_telo.rds'))

read_track <- function(dir, xG, cohort) {
  f <- file.path(dir, xG, sprintf('%s_%s_%s_CNMatm_combined_hs1.rds', date_tag, cohort_tag[[cohort]], xG))
  dt <- readRDS(f)[, .(seqnames, start, end, AmpPercent)]
  dt <- dt[!overlapsAny(GRanges(seqnames, IRanges(start, end)), blacklist)]
  dt <- chr_sizes[dt, on = "seqnames"]
  dt[, genome_pos := start + offset]
  dt[, seqnames := factor(seqnames, levels = chr_levels)]
  dt[]
}

make_track <- function(dt, colr) {
  ggplot(dt[AmpPercent > 0]) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = AmpPercent), color = colr) +
    geom_rect(data = chr_bounds, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = NA, color = NA) +
    scale_x_continuous(
      labels = function(x) paste0(round(x / 1e6), "Mb"),
      breaks = function(x) c(min(x), max(x)),
      expand = c(0, 0)
    ) +
    facet_grid(~seqnames, scales = "free_x", space = "free_x") +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    labs(x = "Chromosome", y = "Amplification fraction (% samples)", fill = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title = element_text(face = "bold"),
      panel.spacing.x = unit(0.5, "lines"),
      strip.background = element_blank(),
    )
}

plots <- lapply(genes, function(xG) {
  CN_B    <- read_track(cn_dir,    xG, "BRCA")
  CN_L    <- read_track(cn_dir,    xG, "NSCLC")
  walks_B <- read_track(walks_dir, xG, "BRCA")
  walks_L <- read_track(walks_dir, xG, "NSCLC")

  t1 <- make_track(CN_B,    "brown3")
  t2 <- make_track(CN_L,    "royalblue3")
  t3 <- make_track(walks_B, "brown3")
  t4 <- make_track(walks_L, "royalblue3")

  t_comb <- t3 / t4 / t1 / t2

  ggsave(
    filename = file.path(out_dir, sprintf('20260911_FigS2A_B_%s_genomeWide.pdf', xG)),
    plot = t_comb,
    width = 8,
    height = 12,
    units = "in"
  )
  cat(sprintf('[%s] wrote %s\n', xG,
              file.path(out_dir, sprintf('20260911_FigS2A_B_%s_genomeWide.pdf', xG))))

  t_comb
})
names(plots) <- genes
