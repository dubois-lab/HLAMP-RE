library(gUtils)
library(data.table)
library(parallel)

IDs = fread('samplesheet_tumor_blank.csv', header = T)

mclapply(1:nrow(IDs), function(i) {
  sid <- IDs$sid[i]
  seg_dt <- fread(paste0('../03_CNV_calling/outputs/',sid,'_tumor_v_hmmT_5kbBin_UCSC.call.cns'))
  seg_gr <- dt2gr(seg_dt)
  saveRDS(seg_gr, paste0('inputdata/segmentation/',sid,'_tumor_v_hmmT_5kbBin_UCSC.rds'))
}, mc.cores = 8)
