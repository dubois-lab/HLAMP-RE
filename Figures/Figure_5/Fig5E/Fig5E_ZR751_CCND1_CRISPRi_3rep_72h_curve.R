# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(readr)
library(dplyr)
library(ggplot2)
library(forcats)

# -----------------------------
# Step 1: Read tidy samplesheet
# -----------------------------
samplesheet_path <- file.path(IN_DIR, "Fig.5EF_20260206_ZR751_samplesheet_complete.csv")
df <- read_csv(samplesheet_path, show_col_types = FALSE) %>%
  mutate(
    Biorep = factor(Biorep),
    group  = fct_relevel(group, "gNT", "gCCND1_1", "gCCND1_2", "gCCND1_DepE1")
  )

# -----------------------------
# Step 2: Average technical replicates per biorep up to 72 h (exponential phase)
# -----------------------------
df_avg_72 <- df %>%
  filter(group != "gCCND1_1", time <= 72) %>%
  group_by(group, Biorep, time) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

# -----------------------------
# Step 3: Plot exponential proliferation curves (0 - 72 h)
# -----------------------------
gg_72 <- ggplot(df_avg_72, aes(x = time, y = value, color = group)) +
  geom_point(aes(shape = Biorep), size = 2.5, alpha = 0.35) +
  geom_smooth(method = "loess", se = FALSE, aes(group = group), linewidth = 1.2) +
  scale_color_manual(
    name = "Group",
    values = c(
      "gNT"          = "steelblue",
      "gCCND1_2"     = "palegreen3",
      "gCCND1_DepE1" = "purple"
    )
  ) +
  scale_shape_manual(
    name = "Biological Rep",
    values = c("1" = 16, "2" = 17, "3" = 15)
  ) +
  labs(
    title = "ZR751 Incucyte Proliferation Curve (0 - 72 h)",
    subtitle = "Exponential growth phase; points = biological replicates (mean across technical wells)",
    x = "Time (hours)",
    y = "Confluence (normalized to 0 h)"
  ) +
  theme_bw() +
  theme(
    legend.position   = "right",
    plot.title        = element_text(size = 15, face = "bold"),
    plot.subtitle     = element_text(size = 10.5),
    axis.title.x      = element_text(size = 13),
    axis.title.y      = element_text(size = 13),
    axis.text.x       = element_text(size = 12),
    axis.text.y       = element_text(size = 12)
  )

# Export PDF and PNG
out_pdf <- file.path(OUT_DIR, "20261106_ZR751_CCND1_CRISPRi_3rep_72h_curve.pdf")
out_png <- file.path(OUT_DIR, "20261106_ZR751_CCND1_CRISPRi_3rep_72h_curve.png")

ggsave(out_pdf, plot = gg_72, width = 7.5, height = 5.5)
ggsave(out_png, plot = gg_72, width = 7.5, height = 5.5, dpi = 300)

message("Saved 72h proliferation curve to:")
message("PDF: ", out_pdf)
message("PNG: ", out_png)
