library(DSS)
require(bsseq)
library(data.table)

bedmethyls_dir <- "./output/" ## see 01_prepare_bedmethyls_input.R
IDs <- fread('samplesheet_blank.csv',header=T)

# create input list of methylation files used for diff. meth. analysis
meth_list <- lapply(IDs$sid, function(s) {
  dt <- readRDS(paste0(bedmethyls_dir,s,'.bedmethyl.rds'))
  dt <- dt[modified_base_code == 'm', .(chr=seqnames,pos=start,N=Nvalid_cov,X=Nmod)]
  dt[, chr := sub("^chr", "", chr)]
  dt
})
names(meth_list) <- IDs$sid

#create BSseq object for DSS
BSobj = makeBSseqData(meth_list, IDs$sid)

# run differential methylation
g1 <- IDs$sid[IDs$tissue=="LUSC"]
g2 <- IDs$sid[IDs$tissue=="BRCA"]
dmlTest.sm = DMLtest(BSobj, group1=g1, group2=g2, smoothing=TRUE, ncores = 72) # adjust number of cores to fit machine

# dmls (diff. meth. sites)
dmls2 = callDML(dmlTest.sm, p.threshold=0.005, delta=0.2)

# dmrs (diff. meth. regions)
dmrs_d = callDMR(dmlTest.sm, p.threshold=0.01, minlen = 50, dis.merge = 500, pct.sig = 0.6, delta=0.4)
