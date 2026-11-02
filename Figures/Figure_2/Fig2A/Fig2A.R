# Input folder: set IN_DIR in the environment or edit the default here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')

library(gGnome)
library(gTrack)
library(plyranges)

gt.ge  <- readRDS('../../../common/data/refGenome/T2T_gencode_gt_UCSC.rds')
CGC_gr <- readRDS('../../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')

sid <- "C96"
jabba <- gG(jabba = paste0(IN_DIR,'/',sid,'_jabba.gg.rds'))

SVs <- jJ(paste0(IN_DIR,'/',sid,'_sniffles_PoNfiltered.bedpe'))
cov <- dt2gr(fread(paste0(IN_DIR,'/',sid,'_tumor_5kbBin_USCS_cov.csv')))
cov300_gt <- gTrack(data = cov, y.field = 'ratio', name = 'coverage', y0=0, y1 = 300)

CNxt <- 8
highcopyX = jabba[cn>CNxt]
walks = highcopyX$walks()
xG <- "CCND1"
CGC_xGgr = CGC_gr %Q% (Gene_Symbol == xG) 
seqlevels(CGC_xGgr) <- gsub("^chr", "", seqlevels(CGC_xGgr))
seqnames(CGC_xGgr) <- gsub("^chr", "", seqnames(CGC_xGgr))
walks_xGr <- walks %&% CGC_xGgr                                                                                                                                         

window <- GRanges(seqnames="11",ranges=IRanges(start=c(61000000,68500000,73500000),end=c(62750000,71900000,74250000)))

if (length(walks_xGr)>0) {  
  pdf(file='jabba_C96_walks_CCND1.pdf', paper = 'a4', height = 24)
    plot(c(gt.ge,cov300_gt, jabba$gt, walks_xGr$gtrack(name = "walks",labels.suppress=TRUE)), window, cex.label = 0.1,
	 links=SVs$grl,
	 xaxis.interval=0.5e6, xaxis.unit=1e6, xaxis.suffix="Mb",
	 xaxis.cex.label=0.8,
	 yaxis.cex.label=0.6,
	 sep.bg.col="white",
	 yaxis.cex=0.6)
  dev.off()                                                                                                                                                                       
}                                                                                                                                                                      

