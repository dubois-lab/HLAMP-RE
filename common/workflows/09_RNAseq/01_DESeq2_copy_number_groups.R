# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(DESeq2)
library(data.table)
library(ggplot2)
library(ggrepel)

# Usage: Rscript 01_DESeq2_copy_number_groups.R <gene> <tumor type> <AMP or GAIN>
# AMP compares amplified samples with samples of normal copy number
# GAIN also counts amplified samples as gained
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3, args[3] %in% c('AMP', 'GAIN'))
xG <- args[1]
xTT <- args[2]
xAlt <- args[3]

# One gene count table per batch
count_files <- list.files(IN_DIR, pattern = '^salmon\\.merged\\.gene_counts.*\\.tsv$', full.names = TRUE)
countmat_lng <- rbindlist(lapply(count_files, function(xf){
  counts_xb <- fread(xf)
  # Use the gene ID where a gene name occurs more than once
  counts_xb[!duplicated(counts_xb$gene_name), gene_nameDD := gene_name]
  counts_xb[duplicated(counts_xb$gene_name), gene_nameDD := gene_id]
  counts_xbm <- melt.data.table(data = counts_xb, id.vars = c('gene_name','gene_id', 'gene_nameDD'), variable.name = 'sample', value.name = 'count')
  return(counts_xbm)
}))

countmat_m <- dcast(data = countmat_lng, formula = gene_nameDD ~ sample, value.var = 'count')

# One row per sample, named as in the count tables, with the copy number status of each gene
# Leave the status empty for normal copy number
metaD <- fread('samplesheet_blank.csv', na.strings = '')
metaD[is.na(get(xG)), (xG) := 'normCN']
if (xAlt == 'GAIN') metaD[get(xG) == 'AMP', (xG) := 'GAIN']

CGC <- readRDS('../../../Figures/supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
CGC_dt <- as.data.table(CGC)

metaD_xTT <- metaD[group == xTT]

# make sure colnames and metadata order match
idxcco <- as.character(metaD_xTT$sample)
stopifnot(length(idxcco) > 0, idxcco %in% colnames(countmat_m))
countmat_m_ord <- countmat_m[, ..idxcco]

cts <- as.matrix(countmat_m_ord)
rownames(cts) <- countmat_m$gene_nameDD
ctsr <- round(cts)
coldata <- as.data.frame(metaD_xTT[,-1])
rownames(coldata) <- metaD_xTT$sample

# BRCA subtypes are a covariate
xDesign <- as.formula(if (xTT == 'BRCA') paste('~ subgroup +', xG) else paste('~', xG))
dds <- DESeqDataSetFromMatrix(countData = ctsr,
                              colData = coldata,
                              design = xDesign)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, contrast=c(xG,xAlt,"normCN"), parallel = TRUE)
resdt <- as.data.table(keep.rownames = T, res)
resdt[rn %in% CGC_dt$Gene_Symbol, fullname := CGC_dt[rn == Gene_Symbol, Name], by = rn]
resdt[rn %in% CGC_dt$Gene_Symbol, role := CGC_dt[rn == Gene_Symbol, Role_in_Cancer], by = rn]
resdt[rn %in% CGC_dt$Gene_Symbol, chr := CGC_dt[rn == Gene_Symbol, seqnames], by = rn]
resdt[rn %in% CGC_dt$Gene_Symbol, pos12 := CGC_dt[rn == Gene_Symbol, locC], by = rn]

# Fig. S1D reads the AMP tables under these names
out_name <- paste0('20251231_DESeq2_', if (xAlt == 'AMP') 'ampVSnonamp' else 'gainVSnormCN', '_', xTT, '_', xG)
write.csv(row.names = F, resdt, file.path(OUT_DIR, paste0(out_name, '.csv')))

#add colors
resdt[padj < 0.05 & stat < 0 ,colr := 'blue']
resdt[padj < 0.05 & stat > 0 ,colr := 'red']
resdt[padj < 0.01 | rn == xG ,labelid := rn]
resdt[, nlgp := -log10(padj)]
resdts <- resdt[abs(log2FoldChange) >1]
resdts[, opas := 1/11]
resdts[rn %in% CGC_dt$Gene_Symbol, opas := 1/10]

gg <- ggplot(resdts, aes(x=log2FoldChange, y=nlgp, size = baseMean )) + scale_size_continuous(range = c(1, 10)) +
  geom_point(aes(col=colr, size = baseMean , alpha = opas)) + geom_text_repel(aes(label=labelid), label.size = 0.75, force = 3, max.overlaps =20) + theme_classic() + theme(text = element_text(size=24))
pdf(file = file.path(OUT_DIR, paste0(out_name, '.pdf')), paper = 'a4r')
  print(gg)
dev.off()
