# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------
# 1. Load packages
# --------------------------------------------------------------

library(tidyverse)
library(ggsignif)

# --------------------------------------------------------------
# 2. File paths
# --------------------------------------------------------------

csv_path <- file.path(IN_DIR, "20260909_EPLC_CCND1-E1_reporter_3biorep_quant.csv")

out_dir <- OUT_DIR

pdf_path <- file.path(
  out_dir,
  "20260909_CCND1_E1_EPLC_boxplot_log2_paired_ttest_ymax0.pdf"
)

# --------------------------------------------------------------
# 3. Read data
# --------------------------------------------------------------

df <- read_csv(
  csv_path,
  col_types = cols()
) %>%
  rename(
    replicate = Replicate,
    control = `Empty vector`,
    experimental = `CCND1-HRBC-E1`
  )

cat("\nRaw data:\n")
print(df)

# --------------------------------------------------------------
# 4. Log2 transformation
# --------------------------------------------------------------

df <- df %>%
  mutate(
    log2_control = log2(control),
    log2_experimental = log2(experimental)
  )

cat("\nLog2-transformed data:\n")
print(df)

# --------------------------------------------------------------
# 5. Paired t-test
# --------------------------------------------------------------
#
# Rep 1: control <-> experimental
# Rep 2: control <-> experimental
# Rep 3: control <-> experimental
#
# Therefore paired = TRUE.
# --------------------------------------------------------------

ttest_result <- t.test(
  df$log2_control,
  df$log2_experimental,
  paired = TRUE
)

p_value <- ttest_result$p.value

cat("\n========================================\n")
cat("Paired t-test on log2-transformed values\n")
cat("========================================\n")

print(ttest_result)

cat("\nP value =", p_value, "\n")
# --------------------------------------------------------------
# 6. Prepare data for plotting
# --------------------------------------------------------------

# Significance label
if (p_value <= 0.0001) {
  p_label <- "****"
} else if (p_value <= 0.001) {
  p_label <- "***"
} else if (p_value <= 0.01) {
  p_label <- "**"
} else if (p_value <= 0.05) {
  p_label <- "*"
} else {
  p_label <- "ns"
}

# Reshape data for ggplot
df_plot <- df %>%
  select(replicate, log2_control, log2_experimental) %>%
  pivot_longer(
    cols = c(log2_control, log2_experimental),
    names_to = "condition",
    values_to = "log2_value"
  ) %>%
  mutate(
    condition = factor(
      condition,
      levels = c("log2_control", "log2_experimental"),
      labels = c("Empty vector", "CCND1-HRBC-E1")
    )
  )

# Shapes for biological replicates
shape_map <- c(
  "Rep 1" = 16,
  "Rep 2" = 17,
  "Rep 3" = 15
)

# Colours for biological replicates
color_map <- c(
  "Rep 1" = "#E41A1C",
  "Rep 2" = "#377EB8",
  "Rep 3" = "#4DAF4A"
)

# Y-axis range
y_max <- 0
y_min <- min(df_plot$log2_value, na.rm = TRUE)

y_range <- y_max - y_min

if (y_range == 0) {
  y_range <- 1
}

p <- ggplot(
  df_plot,
  aes(x = condition, y = log2_value)
) +
  geom_boxplot(
    width = 0.5,
    fill = "#c6dbef",
    colour = "#3182bd",
    outlier.shape = NA
  ) +
  geom_point(
    aes(shape = replicate, colour = replicate),
    size = 3.5,
    stroke = 1.1
  ) +
  scale_shape_manual(values = shape_map) +
  scale_colour_manual(values = color_map) +
  
  # Paired t-test significance bracket
  geom_signif(
    comparisons = list(c("Empty vector", "CCND1-HRBC-E1")),
    annotations = paste0(
      "P = ",
      format.pval(p_value, digits = 3, eps = 0.001),
      " (", p_label, ")"
    ),
    y_position = -0.15 * y_range,
    tip_length = 0.03,
    textsize = 5
  ) +
  
  labs(
    x = NULL,
    y = expression(log[2]~"(Firefly:Renilla Expression Ratio)"),
    title = "CCND1-E1 reporter - EPLC",
    shape = "Biological replicate",
    colour = "Biological replicate"
  ) +
  coord_cartesian(
    ylim = c(NA, 0)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

print(p)

ggsave(
  filename = pdf_path,
  plot = p,
  device = "pdf",
  width = 6,
  height = 6,
  units = "in"
)