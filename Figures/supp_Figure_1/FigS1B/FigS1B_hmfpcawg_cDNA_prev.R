# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)
hmf_scount <- fread(file.path(IN_DIR, '20220603hmf_counts.txt'))
pcawg_scount <- fread(file.path(IN_DIR, 'Proj_Samples_Percentages_merged.txt'))

pcawgMat <- fread(file.path(IN_DIR, '20210820T_VarMat_amp_Sum.txt'))
paca_pcawgAmpmat <- fread(file.path(IN_DIR, '20220317_cn>5_genecounttable.txt'))
paca_pcawg_ecAmpmat <- fread(file.path(IN_DIR, '20220317_cn>10_genecounttable.txt'))

idxp <- colnames(pcawgMat) %in% c("variant", "Ovary-AdenoCA", "Panc-Endocrine", "Breast-AdenoCA", "Lung-AdenoCA", "Panc-AdenoCA", "Lung-SCC", "Breast-LobularCA" )
pcawg_ecAmpMat <- pcawgMat[grepl(pattern = 'ExChr_amp', x = variant), ..idxp]
pcawg_AmpMat <- pcawgMat[!grepl(pattern = 'ExChr_amp', x = variant), ..idxp]
pcawg_ecAmpMat$variant <- gsub(pattern = ':ExChr_amp', replacement = '', x = pcawg_ecAmpMat$variant)
pcawg_AmpMat$variant <- gsub(pattern = ':amp', replacement = '', x = pcawg_AmpMat$variant)

# pcawg_ecAmpMat[, Breast := `Breast-AdenoCA` + `Breast-LobularCA`]
# pcawg_ecAmpMat[, Lung :=  `Lung-AdenoCA` + `Lung-SCC`]
# pcawg_ecAmpMat[, Ovary :=  `Ovary-AdenoCA`]
# pcawg_ecAmpMat[, Pancreas :=  `Panc-AdenoCA` + `Panc-Endocrine`]
# 
# pcawg_ecAmpMat[, Breast_preval := Breast / 209]
# pcawg_ecAmpMat[, Lung_preval := Lung / (48+36)]
# pcawg_ecAmpMat[, Ovary_preval := Ovary / 112]
# pcawg_ecAmpMat[, Pancreas_preval := Pancreas / 237]
# 
# pcawg_ecAmpMat[,pan_preval := Breast_preval+Lung_preval+Ovary_preval+Pancreas_preval]

hmf_ecAmpMat <- fread(file.path(IN_DIR, '20220528ecAmp_countsMat_unique.txt'))
colnames(hmf_ecAmpMat)[2:5] <- paste0(colnames(hmf_ecAmpMat)[2:5], '_hmf')
ecAmpMat <- merge(x = hmf_ecAmpMat[,-6:-7], y = pcawg_ecAmpMat, by.x = 'gene', by.y = 'variant', all = F)
ecAmpMat[, Breast := Breast_hmf + `Breast-AdenoCA` + `Breast-LobularCA`]
ecAmpMat[, Lung := Lung_hmf + `Lung-AdenoCA` + `Lung-SCC`]
ecAmpMat[, Ovary := Ovary_hmf + `Ovary-AdenoCA`]
ecAmpMat[, Pancreas := Pancreas_hmf + `Panc-AdenoCA` + `Panc-Endocrine`]

ecAmpMat[, Breast_preval := Breast / (209 +802)]
ecAmpMat[, Lung_preval := Lung / (48+36 +593)]
ecAmpMat[, Ovary_preval := Ovary / (112 +178)]
ecAmpMat[, Pancreas_preval := Pancreas / (237 +191)]


idxAmpMat_f <- colnames(ecAmpMat) %in% c("gene", "Breast", "Lung", "Ovary", "Pancreas", "Breast_preval", "Lung_preval", "Ovary_preval", "Pancreas_preval")
ecAmpMat_f <- ecAmpMat[,..idxAmpMat_f]
ecAmpMat_f[, Total_prev := Breast_preval + Lung_preval + Ovary_preval + Pancreas_preval]
ecAmpMat_f[, BrstLuTotal := Breast_preval + Lung_preval ]
ecAmpMat_f[, mean_prev := (Breast_preval + Lung_preval + Ovary_preval + Pancreas_preval)/4]
ecAmpMat_f[, BrstLuMean := (Breast_preval + Lung_preval)/2 ]
write.table(ecAmpMat_f, file.path(OUT_DIR, "20220604hmfpcawg_ecAmpCount_prev.txt"),sep = "\t", row.names = F, col.names = T, quote = F)
