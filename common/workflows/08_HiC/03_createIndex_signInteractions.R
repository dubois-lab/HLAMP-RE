library(HiCDCPlus)
library(data.table)
library(parallel)
library(InteractionSet)

IDs <- fread('samplesheet_HiCDCPlus_blank.csv')

indexfile <- data.frame()
indexfile_10prct <- data.frame()

# create an Index of all significant interactions
for (i in 1:nrow(IDs)) {
  sid <- IDs$sid[i]
  print(paste0("Starting sample: ",sid))

  gi_list <- readRDS(paste0('results/',sid,'_covCovariate_wholeGenome_gi_list_obj.rds'))

  for (i in seq_along(gi_list)) {
    sig_indices <- which(mcols(gi_list[[i]])$qvalue <= 0.01)
    if(length(sig_indices) > 0 ) {
      dt_tmp <- as.data.frame(gi_list[[i]][sig_indices])
      df_tmp <- dt_tmp[, c('seqnames1', 'start1', 'start2')]
      indexfile <- rbind(indexfile, df_tmp)
      indexfile <- unique(indexfile)
    }
    sig_indices_10prct <- which(mcols(gi_list[[i]])$qvalue <= 0.1)
    if(length(sig_indices_10prct) > 0 ) {
      dt_tmp <- as.data.frame(gi_list[[i]][sig_indices_10prct])
      df_tmp <- dt_tmp[, c('seqnames1', 'start1', 'start2')]
      indexfile_10prct <- rbind(indexfile_10prct, df_tmp)
      indexfile_10prct <- unique(indexfile_10prct)
    }
  }
}

colnames(indexfile)<-c('chr','startI','startJ')
data.table::fwrite(indexfile, 'results/hicdcdiff/hicdcDiff_indexfile_1prctSig.txt.gz', sep='\t',row.names=FALSE,quote=FALSE)

colnames(indexfile_10prct)<-c('chr','startI','startJ')
data.table::fwrite(indexfile_10prct, 'results/hicdcdiff/hicdcDiff_indexfile_10prctSig.txt.gz', sep='\t',row.names=FALSE,quote=FALSE)
