library(data.table)
library(parallel)

IDs <- fread('samplesheet_blank.csv',header=T)

mclapply(1:nrow(IDs), function(i) {
  sid <- IDs$sid[i]
  path <- IDs$bedmethylPath[i]
  message("sampleID: ", sid)

  tumor_bedmethyl_ungrouped <- fread(paste0(path,sid,'.bedmethyl.gz'))
  colnames(tumor_bedmethyl_ungrouped)[1:3] = c('seqnames', 'start', 'end')
  colnames(tumor_bedmethyl_ungrouped)[4:18] = c('modified_base_code', 'score', 'strand', 'startposition', 'endposition', 'color','Nvalid_cov', 'percent_modified', 'Nmod', 'Ncanonical', 'Nother_mod', 'Ndelete', 'Nfail', 'Ndiff', 'Nnocall')
  # seqnames formatting
  tumor_bedmethyl_ungrouped[, seqnames := gsub(pattern = 'chr', replacement = '', seqnames), by = seqnames]
  # add demethylation track for visualization
  tumor_bedmethyl_ungrouped[, neglog_mth := -log10(percent_modified/100)]

  saveRDS(tumor_bedmethyl_ungrouped,paste0('output/',sid,'.bedmethyl.rds'))
}, mc.cores = 8)
