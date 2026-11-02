# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')

library(readr)
library(dplyr)
library(ggplot2)
library(forcats)

# -----------------------------
# Step 1: Read your tidy samplesheet
# -----------------------------
df <- read_csv(file.path(IN_DIR, "20251218_samplesheet_complete_v2.csv")) %>%
  mutate(
    Biorep = factor(Biorep),
    group  = fct_relevel(group, "gNT", "gCCND1_1", "gCCND1_2", "gCCND1_DepE1")
  )

# -----------------------------
# Step 2: Average technical replicates within group × biorep × time
# -----------------------------
# This collapses wells (e.g., A4–A8, B4–B8, ...) so each biorep-group-time has a single value
df_biorep_avg <- df %>%
  filter(group != "gCCND1_1") %>%
  group_by(group, Biorep, time) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
# -----------------------

# -----------------------------
# Step 3: Plot — 12 points per time (4 groups × 3 bioreps) with group-specific smooth curves
# -----------------------------
gg <- ggplot(df_biorep_avg, aes(x = time, y = value, color = group, shape = Biorep)) +
  geom_point(size = 2.5, alpha = 0.25) +
  geom_smooth(method = "loess", se = FALSE, aes(group = group), size = 1) +
  scale_color_manual(
    values = c(
      "gNT" = "steelblue",
      "gCCND1_2" = "palegreen3",
      "gCCND1_DepE1" = "purple"
    )
  ) +
  labs(
    title = "MFM223 Incucyte Scatterplot",
    subtitle = "Each point = mean across technical replicates per biological replicate",
    x = "Time (hours)",
    y = "Confluence (normalized)",
    caption = "Source: 20251218_samplesheet_complete_v2.csv"
  ) +
  theme_bw()
  theme(
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    axis.text.x  = element_text(size = 14), 
    axis.text.y  = element_text(size = 14),  
    axis.ticks.length = unit(0.25, "cm")    
)


print(gg)
