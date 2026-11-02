# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggsignif)
library(pracma)

# ---------------------------------------------------------
# Step 1 — Prepare data
# ---------------------------------------------------------
df <- read.csv(
  file.path(IN_DIR, "20260206_ZR751_samplesheet_complete.csv")
)
names(df)[names(df) == "valuie"] <- "value"

df$group <- factor(df$group,
                   levels = c("gNT", "gCCND1_1", "gCCND1_2", "gCCND1_DepE1"))

df_144 <- df %>%
  filter(time <= 144, !is.na(value))

df_144$group <- factor(
  df_144$group,
  levels = c("gNT", "gCCND1_1", "gCCND1_2", "gCCND1_DepE1")
)

#take out g1
valid_groups <- c("gNT", "gCCND1_2", "gCCND1_DepE1")

df_144 <- df %>%
  filter(time <= 144,
         !is.na(value),
         group %in% valid_groups)

# ---------------------------------------------------------
# Step 2 — Fit model + Dunnett test
# ---------------------------------------------------------
model_dyn <- lmer(
  value ~ group * time +
    (1 + time | Biorep) +
    (1 | Biorep:well),
  data = df_144
)

emm_dyn <- emmeans(model_dyn, ~ group)

contr_dyn <- contrast(
  emm_dyn,
  method = "trt.vs.ctrl",
  ref = "gNT",
  adjust = "dunnett"
)

dyn_df <- as.data.frame(contr_dyn)

p_to_star <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

dyn_df$label <- sapply(dyn_df$p.value, p_to_star)


# ---------------------------------------------------------
# Step 3 — Compute AUC per well
# ---------------------------------------------------------
auc_well <- df_144 %>%
  group_by(Biorep, group, well) %>%
  arrange(time) %>%
  summarise(AUC = trapz(time, value), .groups = "drop")

# ---------------------------------------------------------
# Step 4 — Collapse to biorep-level AUC
# ---------------------------------------------------------
auc_biorep <- auc_well %>%
  group_by(Biorep, group) %>%
  summarise(AUC = mean(AUC), .groups = "drop") %>%
  mutate(group = factor(group, levels = levels(df_144$group)))

auc_biorep <- auc_biorep %>%
  mutate(Biorep = factor(Biorep))
# ---------------------------------------------------------
# Step 5 — Prepare significance labels for AUC plot
# ---------------------------------------------------------

tmp <- strsplit(dyn_df$contrast, " - ")

stat_df <- dyn_df %>%
  mutate(
    group2 = vapply(tmp, `[`, character(1), 1),
    group1 = vapply(tmp, `[`, character(1), 2),
    group1 = factor(group1, levels = levels(auc_biorep$group)),
    group2 = factor(group2, levels = levels(auc_biorep$group))
  )

# P value + significance stars
stat_df$label <- paste0(
  "P = ",
  format.pval(stat_df$p.value, digits = 3, eps = 0.001),
  " (", 
  stat_df$label,
  ")"
)

y_max <- max(auc_biorep$AUC)
stat_df$y_position <- seq(
  y_max * 1.05,
  y_max * 1.25,
  length.out = nrow(stat_df)
)

# ---------------------------------------------------------
# Step 6 — Plot: biorep-level AUC boxplot + dynamic-model p-values
# ---------------------------------------------------------
p_auc <- ggplot(
  auc_biorep,
  aes(x = group, y = AUC, fill = group)
) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.shape = NA) +
  geom_point(
    aes(shape = Biorep),
    size = 3,
    alpha = 0.9
  ) +
  scale_shape_manual(
    values = c("1" = 16, "2" = 17, "3" = 15)
  ) +
  scale_fill_manual(
    values = c(
      "gNT" = "steelblue",
      "gCCND1_1" = "red",
      "gCCND1_2" = "palegreen3",
      "gCCND1_DepE1" = "purple"
    )
  ) +
  geom_signif(
    data = stat_df,
    aes(
      xmin = group1,
      xmax = group2,
      annotations = label,
      y_position = y_position
    ),
    manual = TRUE,
    inherit.aes = FALSE,
    textsize = 4
  ) +
  labs(
    title = "ZR751 AUC per biological replicate",
    subtitle = paste(
      "P values from mixed model: value ~ group × time",
      "Visualization uses biorep-level AUC",
      sep = "\n"
    ),
    x = "Group",
    y = "AUC"
  ) +
  coord_cartesian(ylim = c(0, 500)) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    axis.ticks.length = unit(0.25, "cm")
  )


pdf(
  file = file.path(OUT_DIR, "20260909_ZR751_AUC_LMM_adjp.pdf"),
  paper = "a4",
  height = 70
)

print(p_auc)
dev.off()

