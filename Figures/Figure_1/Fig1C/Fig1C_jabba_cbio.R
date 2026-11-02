# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')

library(data.table)
library(gGnome)
CGC <- readRDS(file.path(IN_DIR, '20190829cancergene_elements.rds'))
CGC_gr <- dt2gr(CGC[, c(1:3, 6)])
cbio_genedata_select_sr <- as.data.table(read.csv2(file.path(IN_DIR, '20241124_selected_genomic.csv')))
metaD_Zebanc <- fread(file.path(IN_DIR, "20240724_BRCA__Proben.txt"))

jabba_meta1 <- fread(file.path(IN_DIR, '20250225b_jabbawalks_res.csv'), fill = T)
jabba_meta2 <- fread(file.path(IN_DIR, '20250617b_jabbaMetaNSCLCs.csv'), fill = T)
jabba_meta_all <- rbind(jabba_meta1, jabba_meta2, use.names = F, fill = T)

BRCA_jabba_rds <- list.files(file.path(IN_DIR, 'final_jabba_calls'))
NSCLC_jabba_rds <- list.files(file.path(IN_DIR, 'jabba_rds'))
NSCLC_jabba_rds_ids <- gsub(pattern = '_jabba.simple.gg.rds', replacement = '', x = NSCLC_jabba_rds)

NSCLC_jabba_rds[!(NSCLC_jabba_rds_ids %in% jabba_meta_all$sid)]

xid =jabba_meta_all$sid[1]
xid ='KQ3CFX'
CGC_seg <- (lapply(c(jabba_meta_all$sid, 'KQ3CFX'), function(xid){
  print(xid)
  if (sum(grepl(pattern = xid, x = BRCA_jabba_rds)) >= 1) {
    Xid_jabba <- readRDS(paste0(IN_DIR, '/final_jabba_calls/', 
                                grep(pattern = xid, x = BRCA_jabba_rds, value = T)[1]))
    cn_dt <- Xid_jabba$dt
    CGC_xid_gr <- gr.nochr(CGC_gr) %$% dt2gr(cn_dt[, 1:6])
    CGC_xid_dt <- gr2dt(CGC_xid_gr)
    CGC_xid_dt[, sid := xid]
  }
  if (sum(grepl(pattern = xid, x = NSCLC_jabba_rds)) == 1) {
    Xid_jabba <- readRDS(paste0(IN_DIR, '/jabba_rds/',
                                grep(pattern = xid, x = NSCLC_jabba_rds, value = T)))
    cn_dt <- Xid_jabba$dt
    CGC_xid_gr <- gr.nochr(CGC_gr) %$% dt2gr(cn_dt[, 1:6])
    CGC_xid_dt <- gr2dt(CGC_xid_gr)
    CGC_xid_dt[, sid := xid]
  }
  return(CGC_xid_dt)
}))
# CGC_seg_dt <- rbindlist(CGC_seg)                                      

CGC_seg_dt <- fread(file.path(IN_DIR, '20250708jabba_cgc_seg_dt.csv'))

# cbio_gene_list <- CGC_seg_dt[, c(8, 6, 7)]
# cbio_gene_list[cn >= 5, Alteration := 'GAIN']
# cbio_gene_list[cn >= 10, Alteration := 'AMP']
# cbio_gene_list[cn <= 1, Alteration := 'HETLOSS']
# cbio_gene_list[cn <= 0.5, Alteration := 'HOMDEL']
# cbio_gene_list[, Type := 'CNA']

cbio_gene_list <- fread(file.path(IN_DIR, '20250708jabba_cbiogenelist.csv'))

# genelist: FGFR1,    CCND1  ,   EGFR   ,   PIK3CA   ,  SOX2, TP63, MYC, TP53, CDKN2A, PTEN, ATRX
cbio_meta <- as.data.table(read.csv2(
  file.path(IN_DIR, '20241124_all_clinical.csv')))
# cbio_meta_s <- cbio_meta[sid %in% cbio_gene_list$sid]
# cbio_meta_sm <- as.data.table(merge(cbio_meta_s, jabba_meta_all[, 1:3], by = 'sid', all.y =T))
# class(cbio_meta_sm)
# cbio_meta_sm[is.na(group), group := 'NSCLC']
# cbio_meta_sm[group == 'NSCLC', subgroup := 'LUSC']
cbio_meta_sm <- fread(file.path(IN_DIR, '20250708jabba_cbio_meta.tsv'))
cbio_meta_sm[sid == 'C23', subgroup := 'HER2_HR-']
cbio_meta_sm[sid == 'C23', PRgroup := 'Recurrent']
cbio_meta_sm[sid == 'Cb26', subgroup := 'HR+(lumB)']
cbio_meta_sm[sid == 'Cb26', neoT := 'NACT']
write.table(cbio_meta_sm, row.names = F, sep = '\t', quote = F, 
            file.path(IN_DIR, '20250829jabba_cbio_meta.tsv'))
cbio_meta_sm <- fread(file.path(IN_DIR, '20250829jabba_cbio_meta.tsv'))
write.table(cbio_meta_sm[, c(1:4)], row.names = F, sep = '\t', quote = F, 
            file.path(IN_DIR, '20250829jabba_cbio_meta_s.tsv'))



FGFR1amped_BRCA_metaD_mrg <- as.data.table(read.csv2(file.path(IN_DIR, '20250812_FGFR1amped_BRCA_metaD_mrgd.csv')))



write.table(cbio_meta_sm, row.names = F, sep = '\t', quote = F, 
            file.path(IN_DIR, '20250708jabba_cbio_meta.tsv'))






