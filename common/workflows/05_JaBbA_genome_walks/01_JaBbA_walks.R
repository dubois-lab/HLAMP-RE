#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(gGnome)
  library(parallel)
})

gt.ge  <- readRDS('../../data/refGenome/T2T_gencode_gt_UCSC.rds')
CGC_gr <- readRDS('../../data/refGenome/cgc_T2T_gr_UCSC.rds')

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Error: No argument provided for sid.", call. = FALSE)
sid <- args[1]
print(sid)

OUT_DIR <- './walks'
GENES   <- c('CCND1', 'ERBB2', 'MYC', 'NSD3', 'WHSC1L1', 'EGFR')

jabbaX <- gG(jabba = paste0('../04_JaBbA_genome_reconstruction/results/',sid,'_JaBbA/jabba.gg.rds'))
CNmax <- max(na.omit(jabbaX$dt$cn))
CNxt  <- CNmax * 0.5

while (CNxt >= 8) {
  print(CNxt)
  highcopyX <- jabbaX[cn > CNxt]
  try({
    walks <- highcopyX$walks()
    saveRDS(walks, paste0(OUT_DIR, '/', sid, '_walks_cnLrg', CNxt, '.rds'))
    print('Walks done')

    walks_circ <- walks[circular == TRUE]
    if (length(walks_circ) > 0) {
      print('Found a circular amplicon')
      saveRDS(walks_circ, paste0(OUT_DIR, '/', sid, '_walksCirc_cnLrg', CNxt, '.rds'))
      walks_l <- walks_circ[walk.id %in% names(sort(walks_circ$lengths, decreasing = TRUE)[1:3])]
      pdf(file = paste0(OUT_DIR, '/', sid, '_walksCirc_cnLrg', CNxt, '.pdf'),
          paper = 'a4', height = 24)
      tryCatch(
        plot(c(gt.ge, jabbaX$gt, walks_l$gtrack(name = "Longest walks")),
             walks_l$footprint + 4e5, cex.label = 0.1),
        finally = dev.off()
      )
    }

    # Serial gene loop. Stop scanning further genes for this CNxt once a
    # circular oncogene has been written.
    found_circ_oncogene <- FALSE
    for (xG in GENES) {
      if (found_circ_oncogene) break
      tryCatch({
        print(xG)
        CGC_xGgr <- CGC_gr %Q% (Gene_Symbol == xG)
        if (length(CGC_xGgr) == 0) return(invisible())
        seqlevels(CGC_xGgr) <- gsub("^chr", "", seqlevels(CGC_xGgr))
        seqnames(CGC_xGgr)  <- gsub("^chr", "", seqnames(CGC_xGgr))

        walks_xGr <- walks %&% CGC_xGgr
        if (length(walks_xGr) > 0) {
          print(paste0('Found an ', xG, ' amplicon'))
          saveRDS(walks_xGr,
                  paste0(OUT_DIR, '/', sid, '_', xG, '_OncogeneWalks_cnLrg', CNxt, '.rds'))
          pdf(file = paste0(OUT_DIR, '/', sid, '_', xG, '_walks_cnLrg', CNxt, '.pdf'),
              paper = 'a4', height = 24)
          tryCatch(
            plot(c(gt.ge, jabbaX$gt, walks_xGr$gtrack(name = "walks")),
                 walks_xGr$footprint + 4e5, cex.label = 0.1),
            finally = dev.off()
          )
        }

        if (length(walks_circ) > 0) {
          walks_circ_xGr <- walks_circ %&% CGC_xGgr
          if (length(walks_circ_xGr) > 0) {
            walks_l_xGr <- walks_circ_xGr[walk.id %in% names(sort(walks_circ_xGr$lengths, decreasing = TRUE)[1:min(3, length(walks_circ_xGr))])]
            print(paste0('Found a circular ', xG, ' amplicon'))
            saveRDS(walks_circ_xGr,
                    paste0(OUT_DIR, '/', sid, '_', xG, '_OncogeneWalksCirc_cnLrg', CNxt, '.rds'))
            pdf(file = paste0(OUT_DIR, '/', sid, '_', xG, '_walksCirc_cnLrg', CNxt, '.pdf'),
                paper = 'a4', height = 24)
            tryCatch({
              plot(c(gt.ge, jabbaX$gt, walks_l_xGr$gtrack(name = "Longest walks")),
                   walks_l_xGr$footprint + 4e5, cex.label = 0.1)
              plot(c(gt.ge, jabbaX$gt, walks_circ_xGr$gtrack(name = "Oncogene walks", height = 40)),
                   walks_circ_xGr$footprint + 4e5, cex.label = 0.1)
            }, finally = dev.off())
            found_circ_oncogene <<- TRUE
          }
        }
      }, error = function(e) {
        message(sprintf("  [%s | CN=%s | %s] gene-block error: %s", sid, CNxt, xG, conditionMessage(e)))
      })
    }

    CNxt <- CNxt * 0.2
    if (CNxt < 10 && CNxt > 2) CNxt <- 8
    if (CNxt < 8)              print('done')
  })
}

