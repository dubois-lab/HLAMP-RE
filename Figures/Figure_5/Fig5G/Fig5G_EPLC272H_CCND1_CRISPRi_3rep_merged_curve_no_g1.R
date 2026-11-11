# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(readr)
library(dplyr)
library(ggplot2)
library(forcats)

samplesheet_path <- file.path(IN_DIR, "Fig.5GH_20261104_EPLC_samplesheet_complete.csv")

# -----------------------------
# Step 1: Read corrected tidy samplesheet
# -----------------------------
df <- read_csv(samplesheet_path) %>%
  mutate(
    Biorep = factor(Biorep),
    group  = fct_relevel(group, "gNT", "gCCND1_1", "gCCND1_2", "gCCND1_DepE1")
  )

# -----------------------------
# Step 2: Average technical replicates within group × biorep × time
# -----------------------------
df_biorep_avg <- df %>%
  filter(group != "gCCND1_1") %>%
  group_by(group, Biorep, time) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  filter(time <= 144)

# -----------------------------
# Step 3: Plot — points per biorep with group-specific loess curves
# -----------------------------
gg <- ggplot(df_biorep_avg, aes(x = time, y = value, color = group, shape = Biorep)) +
  geom_point(size = 2.5, alpha = 0.35) +
  geom_smooth(method = "loess", se = FALSE, aes(group = group), linewidth = 1) +
  scale_color_manual(
    values = c(
      "gNT" = "steelblue",
      "gCCND1_2" = "palegreen3",
      "gCCND1_DepE1" = "purple"
    )
  ) +
  labs(
    title = "EPLC Incucyte Scatterplot",
    subtitle = "Each point = mean across technical replicates per biological replicate",
    x = "Time (hours)",
    y = "Confluence (normalized)",
    caption = "Source: Fig.5GH_20261104_EPLC_samplesheet_complete.csv"
  ) +
  theme_bw()

# Save PDF plot
out_pdf <- file.path(OUT_DIR, "20261104_EPLC272H_CCND1_CRISPRi_3rep_merged_curve_no_g1.R.pdf")
ggsave(out_pdf, plot = gg, width = 7, height = 5)

# Also save PNG for quick viewing
out_png <- file.path(OUT_DIR, "20261104_EPLC272H_CCND1_CRISPRi_3rep_merged_curve_no_g1.png")
ggsave(out_png, plot = gg, width = 7, height = 5, dpi = 300)

message("Saved corrected merged curve PDF to: ", out_pdf)
message("Saved corrected merged curve PNG to: ", out_png)
