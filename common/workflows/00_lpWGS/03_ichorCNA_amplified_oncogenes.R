# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (xG in c('FGFR1', 'CCND1', 'MYC', 'EGFR', 'MDM2', 'ERBB2')) dir.create(file.path(OUT_DIR, xG), showWarnings = FALSE)

library(data.table)
library(parallel)
library(gUtils)

# ichorCNA results, one folder per sample (see 01_ichorCNA_snakemake.sh)
results_dir <- file.path(IN_DIR, 'ichorCNA/results')
fls <- list.files(results_dir)

IGVseg <- rbindlist(lapply(fls, function (xsid) {
  print(xsid)
  segX <- fread(paste0(results_dir, '/', xsid, '/', xsid, '.seg'))
  return(segX)
}))
# Add the segments of an earlier ichorCNA run
IGVseg_old <- fread(file.path(IN_DIR, '20240621_ichorCNA_IGV.seg'))
write.table(IGVseg, file.path(OUT_DIR, '20240717_ichorCNA_IGV.seg'),quote=FALSE, sep='\t', row.names=FALSE)
IGVseg_all <- rbind(IGVseg_old, IGVseg)
write.table(IGVseg_all, file.path(OUT_DIR, '20240717_ichorCNA_IGV_all.seg'),quote=FALSE, sep='\t', row.names=FALSE)

seg_all <- rbindlist(lapply(fls, function (xsid) {
  print(xsid)
  segX_txt <- fread(paste0(results_dir, '/', xsid, '/', xsid, '.seg.txt'))
  return(segX_txt)
}))
seg_all_old <- readRDS(file.path(IN_DIR, '20240621_ichorCNA.seg.rds'))
seg_alln <- rbind(seg_all_old, seg_all)
saveRDS(seg_alln, file.path(OUT_DIR, '20240717_ichorCNA.seg.rds'))

# Cancer Gene Census genes on each segment
cgc <- readRDS('../../../Figures/supp_Figure_3/FigS3A/source_data/CGC_hg38.rds')
cgc_hg38_dt <- gr2dt(cgc) 

seg_alln[, chrom := gsub(pattern = 'chr', replacement = '', x = chrom), by = chrom]

cgc_segs <- rbindlist(mclapply(unique(seg_alln$ID), function(xsid){
  print(xsid)
  seg_X <- seg_alln[ID == xsid]
  seg_X_gr <- dt2gr(seg_X)
  cgcX <- cgc %$% seg_X_gr
  cgcX_dt <- gr2dt(cgcX)
  return(cgcX_dt)
}, mc.cores = 10)) # adjust number of cores to fit machine
saveRDS(cgc_segs, file.path(OUT_DIR, '20240717_ichorCNA_cgc.seg.rds'))

# Amplified oncogene segments: corrected copy number above 7, or above 10 for high-level amplification
cgc_segs_ampCN <- cgc_segs[Corrected_Copy_Number > 7 & Role_in_Cancer != 'TSG']
sort(table(cgc_segs_ampCN$Gene_Symbol), decreasing = T)[1:25]
sort(table(cgc_segs_ampCN[Tier ==1, Gene_Symbol]), decreasing = T)[1:20]
cgc_segs_HLampCN <- cgc_segs[Corrected_Copy_Number > 10 & Role_in_Cancer != 'TSG']
sort(table(cgc_segs_HLampCN$Gene_Symbol), decreasing = T)[1:25]
sort(table(cgc_segs_HLampCN[Tier ==1, Gene_Symbol]), decreasing = T)[1:20]

cgc_segs_HLamplRcn <- cgc_segs[logR_Copy_Number >10 & Role_in_Cancer != 'TSG']
sort(table(cgc_segs_HLamplRcn[Tier ==1, Gene_Symbol]), decreasing = T) [1:20]
length(unique(seg_alln$ID))

int_IDs <- unique(c(cgc_segs_HLampCN[Gene_Symbol %in% c(names(sort(table(cgc_segs_HLampCN[Tier ==1, Gene_Symbol]), decreasing = T))[1:10], 'ERBB2', 'CCNE1', 'MDM2', 'MYC', 'FOXA1', 'NKX2-1', 'PPM1D'), ID], 
                  cgc_segs_HLamplRcn[Gene_Symbol %in% c(names(sort(table(cgc_segs_HLamplRcn[Tier ==1, Gene_Symbol]), decreasing = T))[1:10], 'ERBB2', 'CCNE1', 'MDM2', 'MYC', 'FOXA1'), ID]))
# write out IGV_segs 
IGVseg_oncAmps <- IGVseg_all[sample %in% int_IDs]
write.table(IGVseg_oncAmps, file.path(OUT_DIR, '20240717_ichorCNA_IGV_203oncamps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
write.csv2(cgc_segs[ID %in% int_IDs], file.path(OUT_DIR, '20240717_ichorCNA_203oncamps.csv'), row.names = F)

int_IDs2 <- unique(c(cgc_segs_HLampCN[Gene_Symbol %in% c('FGFR1','CCND1', 'MYC', 'EGFR', 'ERBB2', 'CCNE1', 'MDM2', 'MYC', 'FOXA1', 'NKX2-1', 'PPM1D'), ID], 
                  cgc_segs_HLamplRcn[Gene_Symbol %in% c('FGFR1','CCND1', 'MYC', 'EGFR', 'ERBB2', 'CCNE1', 'MDM2', 'MYC', 'FOXA1', 'NKX2-1', 'PPM1D'), ID]))
IGVseg_oncAmps2 <- IGVseg_all[sample %in% int_IDs2]
write.table(IGVseg_all[sample %in% int_IDs2], file.path(OUT_DIR, '20240717_ichorCNA_IGV_122oncamps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
write.csv2(cgc_segs[ID %in% int_IDs2], file.path(OUT_DIR, '20240717_ichorCNA_122oncamps.csv'), row.names = F)

# Sample annotation with columns sid, group and subgroup
metaD_sequenced <- fread('../../../Figures/Figure_1/Fig1B/source_data/lpWGS_clinical_vars.tsv')

metaD_sequenced_int2 <- metaD_sequenced[sid %in% int_IDs2]
table(metaD_sequenced_int2$group)
table(metaD_sequenced_int2$subgroup)

int_IDs3 <- unique(c(cgc_segs_HLampCN[Gene_Symbol %in% c('FGFR1','CCND1', 'MYC', 'EGFR', 'ERBB2',  'MDM2'), ID], 
                     cgc_segs_HLamplRcn[Gene_Symbol %in% c('FGFR1','CCND1', 'MYC', 'EGFR', 'ERBB2', 'MDM2'), ID]))
write.table(IGVseg_all[sample %in% int_IDs3], file.path(OUT_DIR, '20240717_ichorCNA_IGV110oncamps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
write.csv2(cgc_segs[ID %in% int_IDs3], file.path(OUT_DIR, '20240717_ichorCNA_110oncamps.csv'), row.names = F)

metaD_sequenced_int3 <- metaD_sequenced[sid %in% int_IDs3]
table(metaD_sequenced_int3$group)
table(metaD_sequenced_int3$subgroup)

FGFR1amps <- unique(c(cgc_segs_HLampCN[Gene_Symbol %in% c('FGFR1'), ID], cgc_segs_HLamplRcn[Gene_Symbol %in% c('FGFR1'), ID]))
table(metaD_sequenced_int2[sid %in% FGFR1amps, group])
table(metaD_sequenced_int2[sid %in% FGFR1amps, subgroup])
FGFR1amps_BRCA <- metaD_sequenced_int2[sid %in% FGFR1amps & group == 'BRCA', sid]
write.table(IGVseg_all[sample %in% FGFR1amps_BRCA], file.path(OUT_DIR, 'FGFR1/20240717_ichorCNA_IGV_BRCA_FGFR1amps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
FGFR1amps_LUSC <- metaD_sequenced_int2[sid %in% FGFR1amps & subgroup == 'LUSC', sid]
write.table(IGVseg_all[sample %in% FGFR1amps_LUSC], file.path(OUT_DIR, 'FGFR1/20240717_ichorCNA_IGV_LUSC_FGFR1amps.seg'),quote=FALSE, sep='\t', row.names=FALSE)

CCND1amps <- unique(c(cgc_segs_HLampCN[Gene_Symbol %in% c('CCND1'), ID], cgc_segs_HLamplRcn[Gene_Symbol %in% c('CCND1'), ID]))
table(metaD_sequenced_int2[sid %in% CCND1amps, group])
table(metaD_sequenced_int2[sid %in% CCND1amps, subgroup])
CCND1amps_BRCA <- metaD_sequenced_int2[sid %in% CCND1amps & group == 'BRCA', sid]
write.table(IGVseg_all[sample %in% CCND1amps_BRCA], file.path(OUT_DIR, 'CCND1/20240717_ichorCNA_IGV_BRCA_CCND1amps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
CCND1amps_LUSC <- metaD_sequenced_int2[sid %in% CCND1amps & subgroup == 'LUSC', sid]
write.table(IGVseg_all[sample %in% CCND1amps_LUSC], file.path(OUT_DIR, 'CCND1/20240717_ichorCNA_IGV_LUSC_CCND1amps.seg'),quote=FALSE, sep='\t', row.names=FALSE)

intIDsf <- unique(c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC))
c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC)[duplicated(c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC))]

MYCamps <- unique(c(cgc_segs_HLampCN[Gene_Symbol %in% c('MYC'), ID], cgc_segs_HLamplRcn[Gene_Symbol %in% c('MYC'), ID]))
table(metaD_sequenced_int2[sid %in% MYCamps, group])
table(metaD_sequenced_int2[sid %in% MYCamps, subgroup])
MYCamps_BRCA <- metaD_sequenced_int2[sid %in% MYCamps & group == 'BRCA', sid]
write.table(IGVseg_all[sample %in% MYCamps_BRCA], file.path(OUT_DIR, 'MYC/20240717_ichorCNA_IGV_BRCA_MYCamps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
MYCamps_LUSC <- metaD_sequenced_int2[sid %in% MYCamps & subgroup == 'LUSC', sid]
write.table(IGVseg_all[sample %in% MYCamps_LUSC], file.path(OUT_DIR, 'MYC/20240717_ichorCNA_IGV_LUSC_MYCamps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
MYCamps_LUAD <- metaD_sequenced_int2[sid %in% MYCamps & subgroup == 'LUAD', sid]
write.table(IGVseg_all[sample %in% MYCamps_LUAD], file.path(OUT_DIR, 'MYC/20240717_ichorCNA_IGV_LUAD_MYCamps.seg'),quote=FALSE, sep='\t', row.names=FALSE)

intIDsfm <- unique(c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, MYCamps))
c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, MYCamps)[duplicated(c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, MYCamps))]

EGFRamps <- unique(c(cgc_segs_HLampCN[Gene_Symbol %in% c('EGFR'), ID], cgc_segs_HLamplRcn[Gene_Symbol %in% c('EGFR'), ID]))
table(metaD_sequenced_int2[sid %in% EGFRamps, group])
table(metaD_sequenced_int2[sid %in% EGFRamps, subgroup])
EGFRamps_BRCA <- metaD_sequenced_int2[sid %in% EGFRamps & group == 'BRCA', sid]
write.table(IGVseg_all[sample %in% EGFRamps_BRCA], file.path(OUT_DIR, 'EGFR/20240717_ichorCNA_IGV_BRCA_EGFRamps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
EGFRamps_LUSC <- metaD_sequenced_int2[sid %in% EGFRamps & subgroup == 'LUSC', sid]
write.table(IGVseg_all[sample %in% EGFRamps_LUSC], file.path(OUT_DIR, 'EGFR/20240717_ichorCNA_IGV_LUSC_EGFRamps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
EGFRamps_LUAD <- metaD_sequenced_int2[sid %in% EGFRamps & subgroup == 'LUAD', sid]
write.table(IGVseg_all[sample %in% EGFRamps_LUAD], file.path(OUT_DIR, 'EGFR/20240717_ichorCNA_IGV_LUAD_EGFRamps.seg'),quote=FALSE, sep='\t', row.names=FALSE)

intIDsf_ME <- unique(c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, MYCamps, EGFRamps_BRCA, EGFRamps_LUSC))
c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, MYCamps, EGFRamps_BRCA, EGFRamps_LUSC)[duplicated(c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, MYCamps, EGFRamps_BRCA, EGFRamps_LUSC))]

MDM2amps <- unique(c(cgc_segs_HLampCN[Gene_Symbol %in% c('MDM2'), ID], cgc_segs_HLamplRcn[Gene_Symbol %in% c('MDM2'), ID]))
table(metaD_sequenced_int2[sid %in% MDM2amps, group])
table(metaD_sequenced_int2[sid %in% MDM2amps, subgroup])
MDM2amps_BRCA <- metaD_sequenced_int2[sid %in% MDM2amps & group == 'BRCA', sid]
write.table(IGVseg_all[sample %in% MDM2amps_BRCA], file.path(OUT_DIR, 'MDM2/20240717_ichorCNA_IGV_BRCA_MDM2amps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
MDM2amps_LUSC <- metaD_sequenced_int2[sid %in% MDM2amps & subgroup == 'LUSC', sid]
write.table(IGVseg_all[sample %in% MDM2amps_LUSC], file.path(OUT_DIR, 'MDM2/20240717_ichorCNA_IGV_LUSC_MDM2amps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
MDM2amps_LUAD <- metaD_sequenced_int2[sid %in% MDM2amps & subgroup == 'LUAD', sid]
write.table(IGVseg_all[sample %in% MDM2amps_LUAD], file.path(OUT_DIR, 'MDM2/20240717_ichorCNA_IGV_LUAD_MDM2amps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
intIDsfm <- unique(c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, MDM2amps))
c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, MDM2amps)[duplicated(c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, MDM2amps))]

ERBB2amps <- unique(c(cgc_segs_HLampCN[Gene_Symbol %in% c('ERBB2'), ID], cgc_segs_HLamplRcn[Gene_Symbol %in% c('ERBB2'), ID]))
table(metaD_sequenced_int2[sid %in% ERBB2amps, group])
table(metaD_sequenced_int2[sid %in% ERBB2amps, subgroup])
ERBB2amps_BRCA <- metaD_sequenced_int2[sid %in% ERBB2amps & group == 'BRCA', sid]
write.table(IGVseg_all[sample %in% ERBB2amps_BRCA], file.path(OUT_DIR, 'ERBB2/20240717_ichorCNA_IGV_BRCA_ERBB2amps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
ERBB2amps_LUSC <- metaD_sequenced_int2[sid %in% ERBB2amps & subgroup == 'LUSC', sid]
write.table(IGVseg_all[sample %in% ERBB2amps_LUSC], file.path(OUT_DIR, 'ERBB2/20240717_ichorCNA_IGV_LUSC_ERBB2amps.seg'),quote=FALSE, sep='\t', row.names=FALSE)
ERBB2amps_LUAD <- metaD_sequenced_int2[sid %in% ERBB2amps & subgroup == 'LUAD', sid]
write.table(IGVseg_all[sample %in% ERBB2amps_LUAD], file.path(OUT_DIR, 'ERBB2/20240717_ichorCNA_IGV_LUAD_ERBB2amps.seg'),quote=FALSE, sep='\t', row.names=FALSE)

intIDsfm <- unique(c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, ERBB2amps))
c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, ERBB2amps)[duplicated(c(FGFR1amps_BRCA, FGFR1amps_LUSC, CCND1amps_BRCA, CCND1amps_LUSC, ERBB2amps))]
