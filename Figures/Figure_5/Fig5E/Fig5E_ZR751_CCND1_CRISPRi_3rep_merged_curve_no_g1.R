# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')

library(readr)
library(dplyr)
library(ggplot2)
library(forcats)

df <- read_csv(file.path(IN_DIR, "20260206_ZR751_samplesheet_complete.csv")) %>%
  mutate(
    Biorep = factor(Biorep),
    group  = fct_relevel(group, "gNT", "gCCND1_1", "gCCND1_2", "gCCND1_DepE1")
  )

df_biorep_avg <- df %>%
  filter(group != "gCCND1_1") %>%
  group_by(group, Biorep, time) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  filter(time <= 144)

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
    title = "ZR751 Incucyte Scatterplot",
    subtitle = "Each point = mean across technical replicates per biological replicate",
    x = "Time (hours)",
    y = "Confluence (normalized)",
    caption = "Source: 20260206_samplesheet_complete.csv"
  ) +
  theme_bw()


print(gg)
