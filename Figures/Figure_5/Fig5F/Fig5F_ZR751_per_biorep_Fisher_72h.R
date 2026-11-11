# Input and output folders: set IN_DIR and OUT_DIR in the environment or edit the defaults here
IN_DIR <- Sys.getenv('IN_DIR', unset = 'input')
OUT_DIR <- Sys.getenv('OUT_DIR', unset = 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Fig. 5 - ZR751 CCND1 CRISPRi: each biological replicate analysed separately,
# 0-72 h, combined across replicates with Fisher's method
#
#   AUC       per-well trapezoidal AUC, Welch t-test vs gNT
#   emmeans   group difference at the mean time point (36 h), linear and log
#   emtrends  growth rate (slope over 0-72 h), linear and log
# emmeans and emtrends come from lmer(response ~ group * time_days +
# (1 + time_days | well)) with Kenward-Roger df (12).
#
# Fisher's method combines the two-sided unadjusted P values across the 3
# Bioreps, separately for each measure, when the difference from gNT has the
# same sign in all 3 Bioreps. Combined P values are Bonferroni-adjusted for the
# 2 guides.
# =============================================================================

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(pracma)
library(ggplot2)
library(ggsignif)

input_csv <- file.path(IN_DIR, "Fig.5EF_20260206_ZR751_samplesheet_complete.csv")
out_dir   <- OUT_DIR

time_cutoff <- 72

control <- "gNT"
guides  <- c("gCCND1_2", "gCCND1_DepE1")
valid_groups <- c(control, guides)
contrast_labels <- paste(guides, "-", control)
group_colors <- c(
  "gNT" = "steelblue",
  "gCCND1_2" = "palegreen3",
  "gCCND1_DepE1" = "purple"
)
guide_colors <- setNames(group_colors[guides], contrast_labels)

analysis_labels <- c(
  AUC = "AUC (confluence x h)",
  emmeans_linear = "emmeans, linear (confluence at 36 h)",
  emmeans_log = "emmeans, log (at 36 h)",
  emtrends_linear = "emtrends, linear (confluence / day)",
  emtrends_log = "emtrends, log (1 / day)"
)

# ---------------------------------------------------------
# Step 1 - Prepare data
# ---------------------------------------------------------
df_72 <- read.csv(input_csv) %>%
  filter(time <= time_cutoff, !is.na(value), group %in% valid_groups) %>%
  mutate(
    group = factor(group, levels = valid_groups),
    Biorep = factor(Biorep),
    well = factor(well),
    time_days = time / 24,
    log_value = log(value)
  )

# ---------------------------------------------------------
# Step 2 - Tests within each Biorep
# ---------------------------------------------------------
lmer_ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

contrast_rows <- function(grid, analysis, note) {
  raw <- as.data.frame(summary(
    contrast(grid, method = "trt.vs.ctrl", ref = control, adjust = "none"),
    infer = c(TRUE, TRUE)
  ))
  adj <- as.data.frame(
    contrast(grid, method = "trt.vs.ctrl", ref = control, adjust = "dunnettx")
  )
  data.frame(
    analysis = analysis, contrast = as.character(raw$contrast),
    estimate = raw$estimate, SE = raw$SE, df = raw$df,
    lower.CL = raw$lower.CL, upper.CL = raw$upper.CL,
    statistic = raw$t.ratio, p_raw = raw$p.value,
    p_adj = adj$p.value, adjustment = "Dunnett", note = note
  )
}

auc_contrasts <- function(aucs) {
  bind_rows(lapply(guides, function(g) {
    tt <- t.test(aucs$AUC[aucs$group == g], aucs$AUC[aucs$group == control])
    data.frame(
      analysis = "AUC", contrast = paste(g, "-", control),
      estimate = unname(tt$estimate[1] - tt$estimate[2]),
      SE = unname(tt$stderr), df = unname(tt$parameter),
      lower.CL = tt$conf.int[1], upper.CL = tt$conf.int[2],
      statistic = unname(tt$statistic), p_raw = tt$p.value
    )
  })) %>%
    mutate(p_adj = p.adjust(p_raw, method = "holm"), adjustment = "Holm",
           note = "per-well AUC, Welch")
}

# Fits with bobyqa first. If lme4 reports a convergence problem, the model is
# refitted with nloptwrap and the note column says so.
fit_lmm <- function(d, response) {
  model_formula <- as.formula(paste(response, "~ group * time_days + (1 + time_days | well)"))
  fit_with <- function(ctrl) {
    fit_warnings <- character()
    model <- withCallingHandlers(
      suppressMessages(lmer(model_formula, data = d, REML = TRUE, control = ctrl)),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(model = model, warnings = fit_warnings)
  }
  fit <- fit_with(lmer_ctrl)
  if (any(grepl("converge", fit$warnings))) {
    fit <- fit_with(lmerControl(optimizer = "nloptwrap"))
    fit$warnings <- c("refit with nloptwrap after a bobyqa convergence warning", fit$warnings)
  }
  fit
}

analyse_biorep <- function(b) {
  d <- filter(df_72, Biorep == b)

  aucs <- d %>%
    arrange(time) %>%
    group_by(group, well) %>%
    summarise(AUC = trapz(time, value), .groups = "drop")

  contrasts <- auc_contrasts(aucs)
  group_values <- aucs %>%
    group_by(group) %>%
    summarise(value = mean(AUC), .groups = "drop") %>%
    mutate(analysis = "AUC")

  for (scale_name in c("linear", "log")) {
    response <- if (scale_name == "linear") "value" else "log_value"
    fit <- fit_lmm(d, response)
    model <- fit$model
    note <- paste(c(if (isSingular(model)) "singular random effects", fit$warnings),
                  collapse = "; ")
    trends <- emtrends(model, ~ group, var = "time_days", lmer.df = "kenward-roger")

    contrasts <- bind_rows(
      contrasts,
      contrast_rows(suppressMessages(emmeans(model, ~ group, lmer.df = "kenward-roger")),
                    paste0("emmeans_", scale_name), note),
      contrast_rows(trends, paste0("emtrends_", scale_name), note)
    )
    group_values <- bind_rows(
      group_values,
      as.data.frame(trends) %>%
        transmute(group, value = time_days.trend, analysis = paste0("emtrends_", scale_name))
    )
  }

  list(contrasts = mutate(contrasts, Biorep = b),
       group_values = mutate(group_values, Biorep = b))
}

per_biorep <- lapply(levels(df_72$Biorep), analyse_biorep)

results <- bind_rows(lapply(per_biorep, `[[`, "contrasts")) %>%
  select(Biorep, analysis, contrast, estimate, SE, df, lower.CL, upper.CL,
         statistic, p_raw, p_adj, adjustment, note)

group_values <- bind_rows(lapply(per_biorep, `[[`, "group_values")) %>%
  mutate(Biorep = factor(Biorep)) %>%
  select(Biorep, analysis, group, value)

cat("\nPer-Biorep results, 0-72 h\n")
print(results %>% select(Biorep, analysis, contrast, estimate, SE, df, p_raw),
      digits = 3, row.names = FALSE)

write.csv(results,
          file.path(out_dir, "20261110_ZR751_per_biorep_72h_results.csv"),
          row.names = FALSE)
write.csv(group_values,
          file.path(out_dir, "20261110_ZR751_per_biorep_72h_group_values.csv"),
          row.names = FALSE)

# ---------------------------------------------------------
# Step 3 - Consistency across Bioreps and Fisher's method
# ---------------------------------------------------------
fisher_table <- results %>%
  group_by(analysis, contrast) %>%
  summarise(
    diff_B1 = estimate[Biorep == "1"],
    diff_B2 = estimate[Biorep == "2"],
    diff_B3 = estimate[Biorep == "3"],
    p_B1 = p_raw[Biorep == "1"],
    p_B2 = p_raw[Biorep == "2"],
    p_B3 = p_raw[Biorep == "3"],
    n_p_below_0.05 = sum(p_raw < 0.05),
    same_direction = n_distinct(sign(estimate)) == 1,
    .groups = "drop"
  ) %>%
  mutate(
    fisher_chisq = ifelse(same_direction, -2 * (log(p_B1) + log(p_B2) + log(p_B3)), NA),
    fisher_p = ifelse(same_direction, pchisq(fisher_chisq, df = 6, lower.tail = FALSE), NA)
  ) %>%
  group_by(analysis) %>%
  mutate(fisher_p_bonferroni = p.adjust(fisher_p, method = "bonferroni", n = 2)) %>%
  ungroup() %>%
  arrange(match(analysis, names(analysis_labels)), contrast)

cat("\nConsistency across Bioreps and Fisher's method\n")
print(as.data.frame(fisher_table), digits = 3, row.names = FALSE)

write.csv(fisher_table,
          file.path(out_dir, "20261110_ZR751_per_biorep_72h_fisher.csv"),
          row.names = FALSE)

# ---------------------------------------------------------
# Step 4 - Forest plot of the per-Biorep differences
# ---------------------------------------------------------
p_forest <- results %>%
  mutate(panel = factor(analysis_labels[analysis], levels = analysis_labels)) %>%
  ggplot(aes(x = estimate, y = Biorep, color = contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
  geom_point(size = 2.5) +
  facet_wrap(~ panel + contrast, ncol = 2, scales = "free_x") +
  scale_color_manual(values = guide_colors, guide = "none") +
  scale_y_discrete(limits = rev) +
  labs(title = "ZR751 per-Biorep differences vs gNT, 0-72 h",
       subtitle = "Each Biorep analysed separately, wells as replicates",
       x = "Difference vs gNT (95% CI, unadjusted)", y = "Biorep") +
  theme_bw() +
  theme(strip.text = element_text(size = 8))

ggsave(file.path(out_dir, "Fig.5_20261110_ZR751_per_biorep_72h_forest.pdf"),
       p_forest, width = 8, height = 9)

# ---------------------------------------------------------
# Step 5 - Box/dot plots with Fisher P values
# ---------------------------------------------------------
p_to_star <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

format_fisher <- function(p) {
  if (is.na(p)) return("not combined (direction differs)")
  paste0("Fisher P = ", formatC(p, format = "e", digits = 1), " (", p_to_star(p), ")")
}

box_plot <- function(analysis_name, title, y_label, test_note) {
  values <- filter(group_values, analysis == analysis_name)
  stats <- filter(fisher_table, analysis == analysis_name) %>%
    arrange(match(contrast, contrast_labels))
  labels <- vapply(stats$fisher_p_bonferroni, format_fisher, character(1))
  y_max <- max(values$value)

  ggplot(values, aes(x = group, y = value, fill = group)) +
    geom_boxplot(width = 0.6, alpha = 0.5, outlier.shape = NA) +
    geom_line(aes(group = Biorep), color = "grey60", linewidth = 0.4) +
    geom_point(aes(shape = Biorep), size = 3, alpha = 0.9) +
    scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15)) +
    scale_fill_manual(values = group_colors, guide = "none") +
    geom_signif(
      aes(x = group, y = value),
      inherit.aes = FALSE,
      comparisons = lapply(guides, function(g) c(control, g)),
      annotations = labels,
      y_position = c(y_max * 1.08, y_max * 1.24),
      textsize = 4
    ) +
    expand_limits(y = c(0, y_max * 1.32)) +
    labs(
      title = title,
      subtitle = paste(
        "P: Fisher's method across 3 Bioreps, Bonferroni-adjusted for 2 guides",
        test_note,
        "Points = one value per Biorep, lines join the same replicate",
        sep = "\n"
      ),
      x = "Group",
      y = y_label
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 11),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.ticks.length = unit(0.25, "cm")
    )
}

p_auc <- box_plot(
  "AUC", "ZR751 AUC per biological replicate, 0-72 h",
  "AUC (normalized confluence x hours)",
  "Per-Biorep test: Welch t-test on the AUC of 5 wells per group"
)
p_rate_linear <- box_plot(
  "emtrends_linear", "ZR751 growth rate per biological replicate, 0-72 h",
  "Growth rate (normalized confluence / day)",
  "Per-Biorep test: emtrends on lmer with 5 wells per group"
)
p_rate_log <- box_plot(
  "emtrends_log", "ZR751 specific growth rate per replicate, 0-72 h",
  "Specific growth rate (1/day)",
  "Per-Biorep test: emtrends on lmer of log confluence, 5 wells per group"
)

ggsave(file.path(out_dir, "Fig.5F_20261110_ZR751_72h_AUC_Fisher.pdf"),
       p_auc, width = 7, height = 6.5)
ggsave(file.path(out_dir, "Fig.5F_20261110_ZR751_72h_emtrends_linear_Fisher.pdf"),
       p_rate_linear, width = 7, height = 6.5)
ggsave(file.path(out_dir, "Fig.5F_20261110_ZR751_72h_emtrends_log_Fisher.pdf"),
       p_rate_log, width = 7, height = 6.5)
