library(HiCDCPlus)
library(data.table)
library(parallel)
library(GenomicRanges)
library(InteractionSet)

case <- commandArgs(trailingOnly = TRUE)
cores <- case[1]  # number of cores

IDs <- fread('samplesheet_HiCDCPlus_blank.csv')

hicPath <- "./results/hicpro/valid_pairs" # path to nf-core/hic hicpro validpairs results

for (i in 1:nrow(IDs)) {
  sid <- IDs$sid[i]
  covPath <- IDs$covPath[i]
  print(paste0("Starting sample: ",sid))

  gi_list <- generate_bintolen_gi_list(bintolen_path = paste0("source_data/T2T_LinkPrep_5kb_wholeGenome_noChrY_bintolen.txt.gz"), gen = "Hsapiens", gen_ver = "T2T.CHM13v2.0")
  gi_list <- add_hicpro_allvalidpairs_counts(gi_list, allvalidpairs_path = paste0(hicPath,"/",sid,".allValidPairs"), binned = TRUE)

  ## add coverage as covariate
  cov <- read.csv(covPath, header=F, skip=1, stringsAsFactors=FALSE)
  colnames(cov) <- c('chr', 'start', 'end', 'cov')
  cov <- cov[!grepl("^chrM", cov$chr), ]
  cov_df <- as.data.frame(cov)
  colnames(cov_df) <- c('chr', 'start', 'end', 'cov')
  gi_list <- add_1D_features(gi_list,cov_df)
  gi_list <- expand_1D_features(gi_list)

  ## run HiCDCPlus
  set.seed(1010)
  print("Starting HiCDCPlus_parallel...")
  gi_list <- HiCDCPlus_parallel(gi_list,ncore=cores)
  gi_list <- lapply(gi_list, function(x) {
    mcols(x)$normcounts <- mcols(x)$counts / mcols(x)$mu
    mcols(x)$zvalue <- (mcols(x)$counts - mcols(x)$mu) / mcols(x)$sdev
    return(x)
  })
  saveRDS(gi_list,paste0('results/',sid,'_covCovariate_wholeGenome_gi_list_obj.rds'))
  print(paste0("Done for sample: ",sid))
}
