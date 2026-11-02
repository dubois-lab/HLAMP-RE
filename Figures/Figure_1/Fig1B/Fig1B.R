library(data.table)
library(dplyr)
library(ggplot2)

clinical <- fread('source_data/lpWGS_clinical.tsv')
genomic <- fread('source_data/lpWGS_genomic.tsv')
colnames(genomic) <- c("sid","gene","amp","X")

# get sample sizes
sid <- unique(genomic$sid)
sid_dt <- as.data.table(sid,use.row.names=F)
count <- sid_dt[clinical, group:=i.group, on="sid"]
sort(table(count[sid %in% sid]$group), decreasing = TRUE)
n_NSCLC <- 293
n_BRCA <- 166

# get AMP counts
amp <- genomic[amp=="AMP"]
amp_clinical <- amp[clinical, group:=i.group, on="sid"]  
amp_clinical[, .(total=.N,BRCA=sum(group=="BRCA"),NSCLC=sum(group=="NSCLC")),by=gene][order(-total)]

# plotting
df <- data.frame(
  sample = rep(c("FGFR1/NSD3 (8p11.23)","CCND1 (11q13.3)","SOX2 (3q26.33)","MYC (8q24.21)","PIK3CA (3q26.32)","EGFR (7p11.2)","ERBB2 (17q12)"), each = 2),
  type   = rep(c("BRCA", "NSCLC"), times = 7),
  value  = c(8/n_BRCA*100,21/n_NSCLC*100,10/n_BRCA*100,11/n_NSCLC*100,1/n_BRCA*100,17/n_NSCLC*100,12/n_BRCA*100,5/n_NSCLC*100,1/n_BRCA*100,15/n_NSCLC*100,5/n_BRCA*100,9/n_NSCLC*100,9/n_BRCA*100,2/n_NSCLC*100)
)
sample_order <- c("ERBB2 (17q12)","EGFR (7p11.2)","PIK3CA (3q26.32)","MYC (8q24.21)","SOX2 (3q26.33)","CCND1 (11q13.3)","FGFR1/NSD3 (8p11.23)")
df$sample <- factor(df$sample, levels = sample_order)
df$type <- factor(df$type, levels = c("BRCA", "NSCLC"))
df <- df %>% arrange(sample, type)
ggplot(df, aes(x = value, y = sample, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("NSCLC" = "royalblue3","BRCA" = "brown3"),
                    breaks=c("NSCLC", "BRCA"),
                    labels = c("NSCLC" = "NSCLC (n = 239)", "BRCA" = "BRCA (n = 166)")) +
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
