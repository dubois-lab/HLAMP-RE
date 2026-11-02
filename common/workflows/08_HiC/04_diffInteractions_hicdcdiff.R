library(HiCDCPlus)
library(data.table)
library(parallel)
library(GenomicRanges)
library(InteractionSet)
library(DESeq2)

case <- commandArgs(trailingOnly = TRUE)
chromosomeL <- case[1] 

IDs <- fread('samplesheet_HiCDCPlus_blank.csv')

# create list of HiCDCPlus files for each group
hicplusPath <- './results/'
suffix <- '_covCovariate_wholeGenome_gi_list_obj.rds'
IDs[, input_path := file.path(hicplusPath,paste0(sid,suffix))]
input_paths <- split(IDs$input_path, IDs$tissue)

hicdcdiff(
  input_paths = input_paths,
  filter_file = 'results/hicdcdiff/hicdcDiff_indexfile_1prctSig.txt.gz',
  output_path = 'results/hicdcdiff/hicdcdiff_breastVSlung_fitTypeLocal_1prctSig',
  bin_type = 'Bins-uniform',
  binsize = 5000,
  chrs = chromosomeL,
  fitType = 'local',
  DESeq.save = FALSE,
  diagnostics=TRUE
)
hicdcdiff(
  input_paths = input_paths,
  filter_file = 'results/hicdcdiff/hicdcDiff_indexfile_10prctSig.txt.gz',
  output_path = 'results/hicdcdiff/hicdcdiff_breastVSlung_fitTypeLocal_10prctSig',
  bin_type = 'Bins-uniform',
  binsize = 5000,
  chrs = chromosomeL,
  fitType = 'local',
  DESeq.save = FALSE,
  diagnostics=TRUE
)
