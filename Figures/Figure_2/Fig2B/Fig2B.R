library(gGnome)
library(gTrack)
library(plyranges)

gt.ge  <- readRDS('../../../common/data/refGenome/T2T_gencode_gt_UCSC.rds')
CGC_gr <- readRDS('../../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')

sid <- "C169"
jabba <- gG(jabba = paste0('source_data/',sid,'_jabba.gg.rds'))
circ_walk <- readRDS(paste0('source_data/',sid,''))

window <- GRanges(seqnames="11",ranges=IRanges(start=c(68990317,77759934),end=c(71152739,78874195)))

if (length(walks_xGr)>0) {  
  pdf(file='jabba_C169_cirWalk_CCND1.pdf', paper = 'a4', height = 24)
    plot(c(gt.ge, jabba$gt, circ_walk$gtrack(name = "walks",labels.suppress=TRUE)), window, cex.label = 0.1,
	 xaxis.interval=0.5e6, xaxis.unit=1e6, xaxis.suffix="Mb",
	 xaxis.cex.label=0.8,
	 yaxis.cex.label=0.6,
	 sep.bg.col="white",
	 yaxis.cex=0.6)
  dev.off()                                                                                                                                                                       
}                                                                                                                                                                      

