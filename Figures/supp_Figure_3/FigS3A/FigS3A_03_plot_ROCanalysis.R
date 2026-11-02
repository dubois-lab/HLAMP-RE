# plot ROC curve
sens_high <- readRDS("source_data/HCC1954_ChIPvsMeth_AUC_sensitivity_covRange_high.rds")
spec_high <- readRDS("source_data/HCC1954_ChIPvsMeth_AUC_specificity_covRange_high.rds")
fpr_high <- 1 - spec_high

sens_medium <- readRDS("source_data/HCC1954_ChIPvsMeth_AUC_sensitivity_covRange_medium.rds")
spec_medium <- readRDS("source_data/HCC1954_ChIPvsMeth_AUC_specificity_covRange_medium.rds")
fpr_medium <- 1 - spec_medium

graphics::plot(fpr_medium, sens_medium, type = "l",
               xlab = "False Positive Rate (1 - specificity)",
               ylab = "True Positive Rate (sensitivity)",
               ylim=c(0,1),
               col="darkolivegreen4")
abline(0, 1, lty = 2, col = "gray")
graphics::lines(fpr_high,sens_high,col="royalblue")
graphics::legend("bottomright",
                 legend=c("cov := [125,250]","cov := [350,550]"),
                 col=c("darkolivegreen4","royalblue"),
                 lwd=2)

# compute AUC
o <- order(fpr_high)
x <- fpr_high[o]
y <- sens_high[o]

auc <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
