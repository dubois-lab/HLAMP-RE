library(ichorCNA)
library(optparse)
option_list <- list(
  make_option(c("--id"), type = "character", default="test", help = "Patient ID. Default: [%default]"),
  make_option(c("--cores"), type = "numeric", default=24, help = "Cores. Default: 24"),
  make_option(c("--outDir"), type="character", default="./", help = "Output Directory. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

patientID <- opt$id
numcores <- opt$cores
outDir <- opt$outDir

print(patientID); 
# The read depth files are in 10 kb bins: set binSize to 10000 in the config file
run_ichorCNA(id = patientID, tumor_wig = paste0(outDir, "/readDepth/", patientID, ".bin10000.wig"), repTimeWig = '/ichorCNA/inst/extdata/Koren_repTiming_hg38_10kb.wig',
             gcWig = '/ichorCNA/inst/extdata/gc_hg38_10kb.wig', mapWig = '/ichorCNA/inst/extdata/map_hg38_10kb.wig', ploidy = 'c(2)', normal = 'c(0.25,0.5, 0.75)',
             maxCN = 10, chrs = 'c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")',
             chrTrain = 'c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")',
             genomeStyle = 'UCSC', genomeBuild = 'hg38', estimateScPrevalence = TRUE, maxFracGenomeSubclone = 0.5, maxFracCNASubclone = 0.7, scPenalty = 0.1, plotYLim = 'c(-2,10)',
             outDir = paste0(outDir, "/ichorCNA/", patientID,"/"), cores  = numcores , centromere = '/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt')
