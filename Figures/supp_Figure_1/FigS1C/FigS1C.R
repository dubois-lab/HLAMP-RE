library(data.table)
library(dplyr)
library(ggplot2)

clinical <- fread('~/Downloads/cbio_ONT_clinical.tsv')
genomic <- fread('~/Downloads/CN_SNV_genomics.tsv')
colnames(genomic) <- c("sid","gene","amp","X")
n_NSCLC <- 30
n_BRCA <- 28

# get AMP counts
amp <- genomic[amp=="AMP"]
amp_clinical <- amp[clinical, group:=i.group, on="sid"]
amp_clinical[, .(total=.N,BRCA=sum(group=="BRCA"),NSCLC=sum(group=="NSCLC")),by=gene][order(-total)]

# plotting (min. sample size = 4)
df <- data.frame(
   sample = rep(c("FGFR1/NSD3 (8p11.23)","CCND1 (11q13.3)","EGFR (7p11.2)","PIK3CA (3q26.32)","SOX2/FXR1 (3q26.33)"), each = 2),
   type   = rep(c("BRCA", "NSCLC"), times = 5),
   value  = c(6/n_BRCA,10/n_NSCLC,6/n_BRCA,6/n_NSCLC,5/n_BRCA,4/n_NSCLC,0/n_BRCA,4/n_NSCLC,0/n_BRCA,4/n_NSCLC)
)
sample_order <- c("SOX2/FXR1 (3q26.33)","PIK3CA (3q26.32)","EGFR (7p11.2)","CCND1 (11q13.3)","FGFR1/NSD3 (8p11.23)")
df$sample <- factor(df$sample, levels = sample_order)
df$type <- factor(df$type, levels = c("BRCA", "NSCLC"))
df <- df %>% arrange(sample, type)
ggplot(df, aes(x = value, y = sample, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("NSCLC" = "royalblue3","BRCA" = "brown3"),
                    breaks=c("NSCLC", "BRCA"),
                    labels = c("NSCLC" = "NSCLC (n = 30)", "BRCA" = "BRCA (n = 28)")) +
  scale_x_continuous(expand=c(0,0)) +
  labs(x = "Samples [%]", y = "Amplified genes (cytoband)") +
  theme_minimal() + 
  theme(
    legend.title=element_blank(),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    axis.line=element_line(color="black"),
    axis.ticks=element_line(color="black"),
    axis.title=element_text(face="bold")
  )
