# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)

CGC    <- readRDS('../../supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
CGC_dt <- as.data.table(CGC)

# Source directories
src_dir      <- IN_DIR
src_dir_fig1 <- IN_DIR

make_volcano <- function(xTT, xG, plot_title, highlight_genes = xG,
                         data_dir = src_dir) {
  resdt <- fread(file.path(data_dir,
                           paste0('20251231_DESeq2_ampVSnonamp_', xTT, '_', xG, '.csv')))

  resdt[padj < 0.05 & stat < 0, colr := 'blue']
  resdt[padj < 0.05 & stat > 0, colr := 'red']
  resdt[(padj < 0.05 & rn %in% CGC_dt$Gene_Symbol) | rn %in% highlight_genes | padj < 1e-7,
        labelid := rn]
  resdt[, nlgp := -log10(padj)]
  resdts <- resdt[abs(log2FoldChange) > 1]
  resdts[, opas := 0.2]
  resdts[rn %in% CGC_dt$Gene_Symbol, opas := 0.5]
  resdts[rn %in% highlight_genes, opas := 1]
  resdts[, baseMean_clamped := pmin(pmax(baseMean, 0.5e4), 3e4)]

  ggplot(resdts, aes(x = log2FoldChange, y = nlgp, size = baseMean_clamped)) +
    scale_size_continuous(limits = c(0.5e4, 3e4), range = c(1, 10),
                          breaks = c(5000, 15000, 30000),
                          name = "mean count") +
    scale_alpha_identity() +
    geom_point(aes(col = colr, size = baseMean_clamped, alpha = opas)) +
    geom_text_repel(aes(label = labelid), size = 6, force = 3, max.overlaps = 20,
                    family = "sans") +
    guides(colour = "none") +
    ylab(expression(-log[10](padj))) +
    ggtitle(plot_title) +
    theme_classic(base_family = "sans") +
    theme(text = element_text(size = 24, family = "sans"),
          plot.title = element_text(size = 24, family = "sans", hjust = 0.5))
}

# CCND1 panels
gg_CCND1_BRCA <- make_volcano('BRCA',  'CCND1', 'CCND1 - BRCA')
gg_CCND1_LUSC <- make_volcano('NSCLC', 'CCND1', 'CCND1 - LUSC')

# EGFR panels
gg_EGFR_BRCA  <- make_volcano('BRCA',  'EGFR',  'EGFR - BRCA')
gg_EGFR_LUSC  <- make_volcano('NSCLC', 'EGFR',  'EGFR - LUSC')

# FGFR1/NSD3 panels; highlight both co-amplified genes
gg_FGFR1_BRCA <- make_volcano('BRCA',  'FGFR1', 'FGFR1/NSD3 - BRCA',
                               highlight_genes = c('FGFR1', 'NSD3'),
                               data_dir = src_dir_fig1)
gg_FGFR1_LUSC <- make_volcano('NSCLC', 'FGFR1', 'FGFR1/NSD3 - LUSC',
                               highlight_genes = c('FGFR1', 'NSD3'),
                               data_dir = src_dir_fig1)

combined <- (gg_CCND1_BRCA | gg_CCND1_LUSC) /
            (gg_EGFR_BRCA  | gg_EGFR_LUSC)  /
            (gg_FGFR1_BRCA | gg_FGFR1_LUSC) +
            plot_layout(guides = 'collect')

pdf(file = file.path(OUT_DIR, '20260812_DESeq2_ampVSnonamp_CCND1_EGFR_FGFR1NSD3_BRCA_LUSC.pdf'),
    width = 16.54, height = 17.53)   # A3 landscape height for 3 rows
print(combined)
dev.off()
