library(uwot)
library(ggplot2)
library(stats)
library(data.table)

case <- commandArgs(trailingOnly = TRUE)
xG <- case[1]  # CCND1, NSD3, EGFR

IDs <- fread("samplesheet_blank.csv")

## prepare matrix
meth_mat <- readRDS(paste0(xG,'_meth_matrix.rds'))
meth_complete <- meth_mat[, colSums(is.na(meth_mat))==0]
# filter no variance (invariant) sites
vars <- apply(meth_complete,2,var)
keep_var <- is.finite(vars) & vars > 0
meth_var <- meth_complete[, keep_var]
# scale
meth_scaled <- scale(meth_var)

## run UMAP
umap_res <- umap(meth_scaled, n_neighbors=5, min_dist=0.1, metric="correlation")
#saveRDS(umap_res,paste0(xG,'_umap_res_n5.rds'))

## plot
umap_df <- data.frame(
  Sample = IDs$sid,
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  Group = IDs$tissue
)
ggplot(
  umap_df,
  aes(UMAP1, UMAP2,
  color = Group,
  label = Sample)
) +
  geom_point(size = 4) +
  geom_text(vjust = -0.7) +
  theme_classic()
