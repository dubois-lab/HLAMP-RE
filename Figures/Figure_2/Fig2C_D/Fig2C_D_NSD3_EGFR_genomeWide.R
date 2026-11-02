# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(gGnome)

xG <- "CCND1" # NSD3, EGFR

cn_consAmp_BRCA <- readRDS(paste0(IN_DIR, '/20260831_ONT_HMF_PCAWG_ampSum_cnBased_',xG,'_BRCA_binned10kb.rds'))
cn_consAmp_NSCLC <- readRDS(paste0(IN_DIR, '/20260831_ONT_HMF_PCAWG_ampSum_cnBased_',xG,'_NSCLC_binned10kb.rds'))
walks_consAmp_NSCLC <- readRDS(paste0(IN_DIR, '/20260825_NSCLC_LUSC_',xG,'_CNMatm_combined_hs1.rds'))
walks_consAmp_BRCA <- readRDS(paste0(IN_DIR, '/20260825_BRCA_',xG,'_CNMatm_combined_hs1.rds'))

CN_B <- as.data.frame(cn_consAmp_BRCA)
CN_L <- as.data.frame(cn_consAmp_NSCLC)
walks_B <- as.data.frame(walks_consAmp_BRCA)
walks_L <- as.data.frame(walks_consAmp_NSCLC)

chr_sizes <- fread('../../../common/data/refGenome/T2T_contig_lengths_UCSC.tsv')
colnames(chr_sizes) <- c("seqnames","chr_len")
chr_sizes <- chr_sizes %>%
  mutate(offset = lag(cumsum(as.numeric(chr_len)), default = 0))

CN_B_plot <- CN_B %>%
  left_join(chr_sizes, by="seqnames") %>%
  mutate(genome_pos=start+offset)

CN_L_plot <- CN_L %>%
  left_join(chr_sizes, by="seqnames") %>%
  mutate(genome_pos=start+offset)

walks_B_plot <- walks_B %>%
  left_join(chr_sizes, by="seqnames") %>%
  mutate(genome_pos=start+offset)

walks_L_plot <- walks_L %>%
  left_join(chr_sizes, by="seqnames") %>%
  mutate(genome_pos=start+offset)

CN_B_plot$seqnames <- factor(CN_B_plot$seqnames, levels=paste0("chr", c(1:22, "X", "M")))
CN_L_plot$seqnames <- factor(CN_L_plot$seqnames, levels=paste0("chr", c(1:22, "X", "M")))
walks_B_plot$seqnames <- factor(walks_B_plot$seqnames, levels=paste0("chr", c(1:22, "X", "M")))
walks_L_plot$seqnames <- factor(walks_L_plot$seqnames, levels=paste0("chr", c(1:22, "X", "M")))

chr_bounds <- CN_B_plot %>%
  distinct(seqnames,chr_len) %>%
  mutate(xmin=0,xmax=chr_len)


t1 <- ggplot(CN_B_plot %>% filter(AmpPercent > 0)) +
  geom_rect(aes(xmin=start,xmax=end,ymin = 0,ymax=AmpPercent),color = "brown3") +
  geom_rect(data=chr_bounds,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            inherit.aes=F,fill=NA,color=NA) +
  scale_x_continuous(
    labels = function(x) paste0(round(x/1e6), "Mb"),
    breaks = function(x) c(min(x),max(x)),
    expand = c(0,0)
  ) +
  facet_grid(~seqnames, scales="free_x", space="free_x") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Chromosome", y = "Amplification fraction (% samples)", fill = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle=45,hjust=1,size=8),
    axis.title=element_text(face="bold"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.background = element_blank(),
  )
t2 <- ggplot(CN_L_plot %>% filter(AmpPercent > 0)) +
  geom_rect(aes(xmin=start,xmax=end,ymin = 0,ymax=AmpPercent),color = "royalblue3") +
  geom_rect(data=chr_bounds,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            inherit.aes=F,fill=NA,color=NA) +
  scale_x_continuous(
    labels = function(x) paste0(round(x/1e6), "Mb"),
    breaks = function(x) c(min(x),max(x)),
    expand = c(0,0)
  ) +
  facet_grid(~seqnames, scales="free_x", space="free_x") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Chromosome", y = "Amplification fraction (% samples)", fill = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle=45,hjust=1,size=8),
    axis.title=element_text(face="bold"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.background = element_blank(),
  )
t3 <- ggplot(walks_B_plot %>% filter(AmpPercent > 0)) +
  geom_rect(aes(xmin=start,xmax=end,ymin = 0,ymax=AmpPercent),color = "brown3") +
  geom_rect(data=chr_bounds,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            inherit.aes=F,fill=NA,color=NA) +
  scale_x_continuous(
    labels = function(x) paste0(round(x/1e6), "Mb"),
    breaks = function(x) c(min(x),max(x)),
    expand = c(0,0)
  ) +
  facet_grid(~seqnames, scales="free_x", space="free_x") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Chromosome", y = "Amplification fraction (% samples)", fill = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle=45,hjust=1,size=8),
    axis.title=element_text(face="bold"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.background = element_blank(),
  )
t4 <- ggplot(walks_L_plot %>% filter(AmpPercent > 0)) +
  geom_rect(aes(xmin=start,xmax=end,ymin = 0,ymax=AmpPercent),color = "royalblue3") +
  geom_rect(data=chr_bounds,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            inherit.aes=F,fill=NA,color=NA) +
  scale_x_continuous(
    labels = function(x) paste0(round(x/1e6), "Mb"),
    breaks = function(x) c(min(x),max(x)),
    expand = c(0,0)
  ) +
  facet_grid(~seqnames, scales="free_x", space="free_x") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Chromosome", y = "Amplification fraction (% samples)", fill = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle=45,hjust=1,size=8),
    axis.title=element_text(face="bold"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.background = element_blank(),
  )

t_comb <- t3 / t4 / t1 / t2
t_comb
ggsave(
  filename = paste0(OUT_DIR, "/Fig2C_D_",xG,".pdf"),
  plot = t_comb,
  width = 8,
  height = 12,
  units = "in"
)

###### add second y-axis
sf_CN_B <- max(CN_B_plot$AmpFraction_smoothed,na.rm=T) / max(CN_B_plot$mean_support_smoothed, na.rm=T)
sf_CN_L <- max(CN_L_plot$AmpFraction_smoothed,na.rm=T) / max(CN_L_plot$mean_support_smoothed, na.rm=T)
sf_walks_B <- max(walks_B_plot$AmpFraction_smoothed,na.rm=T) / max(walks_B_plot$mean_support_smoothed, na.rm=T)
sf_walks_L <- max(walks_L_plot$AmpFraction_smoothed,na.rm=T) / max(walks_L_plot$mean_support_smoothed, na.rm=T)

t1 <- ggplot(CN_B_plot %>% filter(AmpFraction_smoothed > 0)) +
  geom_rect(aes(xmin=start,xmax=end,ymin = 0,ymax=AmpFraction_smoothed),color = "brown3") +
  geom_line(aes(x=(start+end)/2, y = mean_support_smoothed * sf_CN_B), color=NA, linewidth=60) +
  geom_rect(data=chr_bounds,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            inherit.aes=F,fill=NA,color=NA) +
  scale_x_continuous(
    labels = function(x) paste0(round(x/1e6), "Mb"),
    breaks = function(x) c(min(x),max(x)),
    expand = c(0,0)
  ) +
  scale_y_continuous(#trans="sqrt",
    sec.axis= sec_axis(~ . / sf_CN_B, name = "mean sample support")) +
  facet_grid(~seqnames, scales="free_x", space="free_x") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Chromosome", y = "Amplification fraction (% samples)", fill = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle=45,hjust=1,size=8),
    axis.title=element_text(face="bold"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.background = element_blank(),
  )
t2 <- ggplot(CN_L_plot %>% filter(AmpFraction_smoothed > 0)) +
  geom_rect(aes(xmin=start,xmax=end,ymin = 0,ymax=AmpFraction_smoothed),color = "royalblue3") +
  geom_line(aes(x=(start+end)/2, y = mean_support_smoothed * sf_CN_L), color=NA, linewidth=60) +
  geom_rect(data=chr_bounds,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            inherit.aes=F,fill=NA,color=NA) +
  scale_x_continuous(
    labels = function(x) paste0(round(x/1e6), "Mb"),
    breaks = function(x) c(min(x),max(x)),
    expand = c(0,0)
  ) +
  scale_y_continuous(#trans="sqrt",
    sec.axis= sec_axis(~ . / sf_CN_L, name = "mean sample support")) +
  facet_grid(~seqnames, scales="free_x", space="free_x") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Chromosome", y = "Amplification fraction (% samples)", fill = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle=45,hjust=1,size=8),
    axis.title=element_text(face="bold"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.background = element_blank(),
  )
sec_brakes <- sort(unique(c(scales::pretty_breaks(n=5)(range(walks_B_plot$mean_support_smoothed)),1)))
t3 <- ggplot(walks_B_plot %>% filter(AmpFraction_smoothed > 0)) +
  geom_rect(aes(xmin=start,xmax=end,ymin = 0,ymax=AmpFraction_smoothed),color = "brown3") +
  geom_line(aes(x=(start+end)/2, y = mean_support_smoothed * sf_walks_B), color=NA, linewidth=60) +
  geom_rect(data=chr_bounds,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            inherit.aes=F,fill=NA,color=NA) +
  scale_x_continuous(
    labels = function(x) paste0(round(x/1e6), "Mb"),
    breaks = function(x) c(min(x),max(x)),
    expand = c(0,0)
  ) +
  scale_y_continuous(#trans="sqrt",
    sec.axis= sec_axis(~ . / sf_walks_B, name = "mean sample support",
                       breaks = sec_brakes)) +
  facet_grid(~seqnames, scales="free_x", space="free_x") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Chromosome", y = "Amplification fraction (% samples)", fill = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle=45,hjust=1,size=8),
    axis.title=element_text(face="bold"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.background = element_blank(),
  )
sec_brakes <- sort(unique(c(scales::pretty_breaks(n=5)(range(walks_L_plot$mean_support_smoothed)),1)))
t4 <- ggplot(walks_L_plot %>% filter(AmpFraction_smoothed > 0)) +
  geom_rect(aes(xmin=start,xmax=end,ymin = 0,ymax=AmpFraction_smoothed),color = "royalblue3") +
  geom_line(aes(x=(start+end)/2, y = mean_support_smoothed * sf_walks_L), color=NA, linewidth=60) +
  geom_rect(data=chr_bounds,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            inherit.aes=F,fill=NA,color=NA) +
  scale_x_continuous(
    labels = function(x) paste0(round(x/1e6), "Mb"),
    breaks = function(x) c(min(x),max(x)),
    expand = c(0,0)
  ) +
  scale_y_continuous(#trans="sqrt",
    sec.axis= sec_axis(~ . / sf_walks_L, name = "mean sample support",
                       breaks = sec_brakes)) +
  facet_grid(~seqnames, scales="free_x", space="free_x") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Chromosome", y = "Amplification fraction (% samples)", fill = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle=45,hjust=1,size=8),
    axis.title=element_text(face="bold"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.background = element_blank(),
  )

t3 / t4 / t1 / t2

t_comb <- t3 / t4 / t1 / t2
t_comb
ggsave(
  filename = file.path(OUT_DIR, "Fig2C_D_NSD3_second_axis.pdf"),
  plot = t_comb,
  width = 8,
  height = 12,
  units = "in"
)
