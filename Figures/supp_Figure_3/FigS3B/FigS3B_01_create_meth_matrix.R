library(data.table)
library(gGnome)
library(rtracklayer)

case <- commandArgs(trailingOnly = TRUE)
xG <- case[1]  # CCND1, NSD3, EGFR

IDs <- fread("samplesheet_blank.csv")

bedmethyls_dir <- "./source_data/" # replace with dir containing bedmethyls
b_l_80_spanning_window <- readRDS(paste0("source_data/",xG,"_80prct_consAmp_window.rds"))

gr_list <- lapply(IDs$sid, function(s) {
  f <- file.path(paste0(bedmethyls_dir,s,".wf_mods.ungrouped.bedmethyl.rds"))
  dt <- readRDS(f)
  gr <- gr.chr(dt2gr(dt[modified_base_code == "m" & Nvalid_cov > 5]))
  rm(dt)
  subsetByOverlaps(gr, b_l_80_spanning_window)
})

names(gr_list) <- IDs$sid

# create gemomic site IDs
site_id <- function(gr){
  paste0(as.character(seqnames(gr)),":",start(gr))
}
all_sites <- unique(unlist(lapply(gr_list,site_id)))

# create methylation matrix
meth_mat <- matrix(
              NA,
              nrow=length(gr_list),
              ncol=length(all_sites),
              dimnames=list(names(gr_list),all_sites))

# fill matrix
for (nm in names(gr_list)) {
  gr <- gr_list[[nm]]
  ids <- site_id(gr)
  meth_mat[nm,ids] <- gr$neglog_mth
}

saveRDS(meth_mat,paste0(xG,'_meth_matrix.rds'))
