# Writes the HiCDCPlus results of one sample for the chromosomes of interest as a table, read by 06_virtual_4C_tracks.R.
# Run from this folder: Rscript 05_interactions_table.R <sample> <chromosome> [<chromosome> ...]

# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'results')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'results')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(InteractionSet)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('Usage: Rscript 05_interactions_table.R <sample> <chromosome> [<chromosome> ...]')
sid <- args[1]
chroms <- args[-1]

gi_list <- readRDS(file.path(IN_DIR, paste0(sid, '_covCovariate_wholeGenome_gi_list_obj.rds')))
interactions_dt <- rbindlist(lapply(chroms, function(xChr) as.data.table(gi_list[[xChr]])))
write.csv(interactions_dt, file.path(OUT_DIR, paste0(sid, '_covCovariate_interactions_dt.csv')), row.names = FALSE)
