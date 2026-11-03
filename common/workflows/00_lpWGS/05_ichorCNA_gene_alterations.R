# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)

# Cancer Gene Census genes on each ichorCNA segment, written by 03_ichorCNA_amplified_oncogenes.R
cgc_segs <- readRDS(file.path(IN_DIR, '20240717_ichorCNA_cgc.seg.rds'))

idxC <- colnames(cgc_segs) %in% c('Gene_Symbol', 'ID', 'Corrected_Copy_Number', 'Corrected_Call')

cbio_varList <- cgc_segs[, ..idxC]
colnames(cbio_varList)[1:3] = c('Gene', 'Sample', 'Copy_Number')
table(cbio_varList$Corrected_Call)
cbio_varList[, Copy_Number := as.numeric(Copy_Number)]
cbio_varList[Copy_Number < 0.3, Alteration := 'HOMDEL']
cbio_varList[Copy_Number >= 0.3 & Copy_Number < 1.3, Alteration := 'HETLOSS', by = Copy_Number]
cbio_varList[ Copy_Number > 3.7, Alteration := 'GAIN']
cbio_varList[ Copy_Number > 10, Alteration := 'AMP']
cbio_varList[ , Type :="CNA"]
# format: Sample Gene Alteration Type
cbio_varList_ord <- cbio_varList[!is.na(Alteration), c(2,1,5,6,3,4)]

# Samples and genes shown in Fig. 1B and Fig. S1A
metaD_sequenced <- fread('../../../Figures/Figure_1/Fig1B/source_data/lpWGS_clinical_vars.tsv')
int_genes <- c('FGFR1', 'CCND1', 'MYC', 'EGFR', 'ERBB2', 'PIK3CA', 'GATA3', 'TP53', 'RB1', 'PTEN', 'BRCA1', 'KRAS', 'BRAF', 'MET', 'NKX2-1', 'TP63', 'FOXA1', 'SOX2', 'CDKN2A')
lpWGS_genomic <- cbio_varList_ord[Sample %in% metaD_sequenced$sid & Gene %in% int_genes, c(1:4)]
write.table(lpWGS_genomic, file.path(OUT_DIR, 'lpWGS_genomic.tsv'), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
