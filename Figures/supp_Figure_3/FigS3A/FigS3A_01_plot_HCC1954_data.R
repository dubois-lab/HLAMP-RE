library(gGnome)
library(gTrack)
library(rtracklayer)

gencode_hg38 = readRDS('source_data/gencode_hg38.rds')

# methylation data
HCC1954_bedmethyl_ungrouped <- readRDS('source_data/HCC1954.ungrouped.bedmethyl.rds')
HCC1954_bedmethyl_ungrouped_logmCe3_gt <- gTrack(dt2gr(HCC1954_bedmethyl_ungrouped[modified_base_code == 'm' & Nvalid_cov > 5]), y.field = 'neglog_mth',bars = T, y1= 3,yaxis.cex=0.7, name='neglog_meth')
# coverage
cnvkit_cov <- read.table("source_data/HCC1954_Tumor_ONT.GRCh38.sorted.targetcoverage.cnn", sep="\t",header=T,stringsAsFactors=F)
cnvkit_cov_gr <- dt2gr(cnvkit_cov) ## depth
cnvkit_cov_gr <- cnvkit_cov_gr[seqnames(cnvkit_cov_gr)=="chr11"]
cnvkit_cov_gt <- gTrack(data=cnvkit_cov_gr, y.field="depth", name="coverage", y0=0, y1=600, yaxis.cex=0.7)
# H3K27ac ChIP
bw_FE <- import('source_data/HCC1954_H3K27ac_ChIP_FE.bw')
chip_FE_gt <- gTrack(data = bw_FE, y.field = 'score', name = 'FE_score',bars = T,yaxis.cex=0.7,y1=8)
rm(bw_FE)

window_CCND1 <- GRanges(seqnames="chr11",ranges=IRanges(start=69000000,end=70500000))

pdf(file = paste0('FigS4A_HCC1954.pdf'),paper='a4r')
plot(c(gencode_hg38,
       chip_FE_gt,
       HCC1954_bedmethyl_ungrouped_logmCe3_gt,
       cnvkit_cov_gt
       ), window_CCND1, cex.label = 0.1,xaxis.interval=0.5e6, xaxis.unit=1e6, xaxis.suffix="Mb")
dev.off()
