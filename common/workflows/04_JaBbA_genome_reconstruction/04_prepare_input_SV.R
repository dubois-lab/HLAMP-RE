library(StructuralVariantAnnotation)
library(data.table)
library(parallel)

sv_dir="../02_SNV_SV_calling/SV_PoN_filtering/output/"
IDs = fread('samplesheet_tumor_blank.csv', header = T)

mclapply(1:nrow(IDs), function(i) {
  sid <- IDs$sid[i]
  sid_sniffles_vcf <- VariantAnnotation::readVcf(paste0(sv_dir,sid,'.sniffles_sv.PoNfiltered.vcf.gz'))
  sid_sniffles_breakpointRanges <- breakpointRanges(sid_sniffles_vcf, inferMissingBreakends =T)
  sid_sniffles_breakendRanges <- breakendRanges(sid_sniffles_vcf)
  sid_bedpe_dt <- as.data.table(breakpointgr2bedpe(sid_sniffles_breakpointRanges))
  sid_bedpe_dt[, SVtype := unlist(strsplit(x = name, split = '[.]'))[2], by = name]
  print(table(sid_bedpe_dt$SVtype))
  sid_bedpe_dt[is.na(SVtype), SVtype := 'BND']
  sid_bedpe_dt[, SPAN := start2 - end1]
  sid_bedpe_dt[SVtype == 'DUP', SPAN := end2 - end1]
  write.table(sid_bedpe_dt, paste0('inputdata/SV_bedpe/',sid,'_sniffles_PoNfiltered.bedpe'), row.names = F, quote = F, sep = '\t', col.names = F)
}, mc.cores = 4, mc.preschedule = F)
