suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

cat("\n")
cat("============================================================\n")
cat("AK1 ENSEMBLE QUANTITATIVE FIGURE (A/B/C)\n")
cat("============================================================\n\n")


# ============================================================
# PATHS
# ============================================================

base_dir <- file.path(
  "..",
  "results",
  "AK1_ensemble_sensitivity"
)

input_dir <- file.path(
  base_dir,
  "docking_analysis_v1"
)

out_dir <- file.path(
  base_dir,
  "quantitative_figure_v2"
)

fig_dir <- file.path(
  out_dir,
  "figures"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


summary_file <- file.path(
  input_dir,
  "AK1_7X7S_ensemble_summary.csv"
)

seed_file <- file.path(
  input_dir,
  "AK1_7X7S_ensemble_seed_results.csv"
)

ap5_file <- file.path(
  input_dir,
  "AK1_7X7S_AP5_region_metrics.csv"
)


# ============================================================
# CHECK FILES
# ============================================================

required_files <- c(
  summary_file,
  seed_file,
  ap5_file
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (
  length(
    missing_files
  ) > 0
) {

  stop(
    paste(
      "Missing required file(s):",
      paste(
        missing_files,
        collapse = "\n"
      ),
      sep = "\n"
    )
  )
}


# ============================================================
# READ DATA
# ============================================================

summary_df <- read_csv(
  summary_file,
  show_col_types = FALSE
)

seed_df <- read_csv(
  seed_file,
  show_col_types = FALSE
)

ap5_df <- read_csv(
  ap5_file,
  show_col_types = FALSE
)


# ============================================================
# ORDER OF RECEPTOR STATES
# ============================================================

state_order <- c(
  "1Z83 closed",
  "7X7S least-open",
  "7X7S medoid",
  "7X7S most-open"
)

display_labels <- c(
  "1Z83\nclosed",
  "7X7S\nleast-open",
  "7X7S\nmedoid",
  "7X7S\nmost-open"
)


# ============================================================
# STANDARDIZE SUMMARY DATA
# ============================================================

summary_plot <- summary_df %>%

  select(
    Receptor,
    Opening_metric_A,
    Mean_Vina_score,
    SD_Vina_score,
    Mean_ISO_AP5_centroid_distance_A,
    Mean_ISO_atoms_within_6A_of_AP5_percent
  ) %>%

  mutate(

    Receptor = factor(
      Receptor,
      levels = state_order,
      labels = display_labels
    )
  )


# ============================================================
# STANDARDIZE SEED-LEVEL SCORE DATA
# ============================================================

score_plot <- seed_df %>%

  select(
    Receptor,
    Seed,
    Best_Vina_score
  ) %>%

  mutate(

    Receptor = factor(
      Receptor,
      levels = state_order,
      labels = display_labels
    )
  )


# ============================================================
# STANDARDIZE AP5-REGION DATA
# ============================================================

ap5_plot <- ap5_df %>%

  select(
    Receptor,
    Seed,
    ISO_AP5_centroid_distance_A,
    ISO_atoms_within_6A_of_AP5_percent
  ) %>%

  mutate(

    Receptor = factor(
      Receptor,
      levels = state_order,
      labels = display_labels
    )
  )


# ============================================================
# CHECK GROUP SIZES
# ============================================================

cat(
  "Docking-score observations per receptor:\n"
)

print(
  table(
    score_plot$Receptor
  )
)

cat("\n")

cat(
  "AP5 localization observations per receptor:\n"
)

print(
  table(
    ap5_plot$Receptor
  )
)

cat("\n")


# ============================================================
# SCORE SUMMARY FROM THE 10 SEEDS
# ============================================================

score_stats <- score_plot %>%

  group_by(
    Receptor
  ) %>%

  summarise(

    n = n(),

    Mean =
      mean(
        Best_Vina_score
      ),

    SD =
      sd(
        Best_Vina_score
      ),

    Min =
      min(
        Best_Vina_score
      ),

    Max =
      max(
        Best_Vina_score
      ),

    .groups = "drop"
  )


# ============================================================
# AP5 LOCALIZATION SUMMARY
# ============================================================

ap5_stats <- ap5_plot %>%

  group_by(
    Receptor
  ) %>%

  summarise(

    n = n(),

    Mean_centroid =
      mean(
        ISO_AP5_centroid_distance_A
      ),

    SD_centroid =
      sd(
        ISO_AP5_centroid_distance_A
      ),

    Mean_within6 =
      mean(
        ISO_atoms_within_6A_of_AP5_percent
      ),

    .groups = "drop"
  )


# ============================================================
# MASTER SUMMARY TABLE
# ============================================================

master_summary <- summary_plot %>%

  select(
    Receptor,
    Opening_metric_A
  ) %>%

  left_join(
    score_stats,
    by = "Receptor"
  ) %>%

  left_join(
    ap5_stats,
    by = "Receptor"
  )


write_csv(
  master_summary,
  file.path(
    out_dir,
    "AK1_ensemble_quantitative_figure_summary.csv"
  )
)


# ============================================================
# COLORS
# ============================================================

receptor_colors <- c(
  "1Z83\nclosed" = "#666666",
  "7X7S\nleast-open" = "#91C788",
  "7X7S\nmedoid" = "#25BFC6",
  "7X7S\nmost-open" = "#D946EF"
)


# ============================================================
# COMMON THEME
# ============================================================

theme_pub <- theme_classic(
  base_size = 12
) +

  theme(

    axis.text.x = element_text(
      size = 10,
      face = "bold",
      color = "black",
      lineheight = 0.9
    ),

    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),

    axis.title = element_text(
      size = 11,
      face = "bold"
    ),

    plot.title = element_text(
      size = 14,
      face = "bold"
    ),

    plot.subtitle = element_text(
      size = 9
    ),

    plot.caption = element_text(
      size = 8,
      color = "gray30"
    ),

    legend.position = "none",

    plot.margin = margin(
      10,
      12,
      10,
      10
    )
  )


# ============================================================
# PANEL A
# RECEPTOR OPENING
# ============================================================

panel_A <- ggplot(
  summary_plot,
  aes(
    x = Receptor,
    y = Opening_metric_A,
    color = Receptor
  )
) +

  geom_segment(

    aes(
      xend = Receptor,
      y = 0,
      yend = Opening_metric_A
    ),

    color = "gray75",
    linewidth = 0.8
  ) +

  geom_point(
    size = 4.5
  ) +

  geom_text(

    aes(
      label = sprintf(
        "%.2f",
        Opening_metric_A
      )
    ),

    color = "black",
    vjust = -1,
    size = 3.7
  ) +

  scale_color_manual(
    values = receptor_colors
  ) +

  scale_y_continuous(

    limits = c(
      0,
      32
    ),

    breaks = seq(
      0,
      30,
      5
    ),

    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(

    title =
      "Receptor opening",

    subtitle =
      "Mean V47-R132 / V47-T145 C-alpha distance",

    x = NULL,

    y =
      "Opening metric (Å)"
  ) +

  theme_pub


# ============================================================
# PANEL B
# VINA SCORE DISTRIBUTIONS
# ============================================================

panel_B <- ggplot(
  score_plot,
  aes(
    x = Receptor,
    y = Best_Vina_score,
    color = Receptor
  )
) +

  geom_jitter(

    width = 0.10,
    height = 0,

    size = 2.8,
    alpha = 0.85
  ) +

  geom_errorbar(

    data = score_stats,

    aes(
      x = Receptor,
      ymin = Mean - SD,
      ymax = Mean + SD
    ),

    inherit.aes = FALSE,

    width = 0.16,
    linewidth = 0.8,
    color = "black"
  ) +

  geom_point(

    data = score_stats,

    aes(
      x = Receptor,
      y = Mean
    ),

    inherit.aes = FALSE,

    shape = 23,
    size = 4.2,
    stroke = 1,
    fill = "white",
    color = "black"
  ) +

  scale_color_manual(
    values = receptor_colors
  ) +

  labs(

    title =
      "Docking-score distribution",

    subtitle =
      "Ten independent docking seeds per receptor state",

    x = NULL,

    y =
      "Best Vina score (kcal/mol)",

    caption =
      "Diamond = mean; error bar = ±1 SD"
  ) +

  theme_pub


# ============================================================
# PANEL C
# AP5-REGION LOCALIZATION
# ============================================================

label_height <- (
  max(
    ap5_plot$ISO_AP5_centroid_distance_A
  )
  +
  0.65
)


panel_C <- ggplot(
  ap5_plot,
  aes(
    x = Receptor,
    y = ISO_AP5_centroid_distance_A,
    color = Receptor
  )
) +

  geom_jitter(

    width = 0.10,
    height = 0,

    size = 2.8,
    alpha = 0.85
  ) +

  geom_errorbar(

    data = ap5_stats,

    aes(
      x = Receptor,
      ymin = Mean_centroid - SD_centroid,
      ymax = Mean_centroid + SD_centroid
    ),

    inherit.aes = FALSE,

    width = 0.16,
    linewidth = 0.8,
    color = "black"
  ) +

  geom_point(

    data = ap5_stats,

    aes(
      x = Receptor,
      y = Mean_centroid
    ),

    inherit.aes = FALSE,

    shape = 23,
    size = 4.2,
    stroke = 1,
    fill = "white",
    color = "black"
  ) +

  geom_label(

    data = ap5_stats,

    aes(
      x = Receptor,
      y = label_height,
      label = paste0(
        "ISO atoms within 6 Å\n",
        sprintf(
          "%.1f%%",
          Mean_within6
        )
      )
    ),

    inherit.aes = FALSE,

    size = 3,
    color = "black",
    fill = "white",
    label.size = 0.25
  ) +

  scale_color_manual(
    values = receptor_colors
  ) +

  scale_y_continuous(

    limits = c(
      0,
      label_height + 0.65
    ),

    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(

    title =
      "Localization to the AP5 region",

    subtitle =
      "ISO centroid relative to crystallographic AP5",

    x = NULL,

    y =
      "ISO-AP5 centroid distance (Å)",

    caption =
      "Diamond = mean; error bar = ±1 SD"
  ) +

  theme_pub


# ============================================================
# SAVE INDIVIDUAL PANELS
# ============================================================

ggsave(

  file.path(
    fig_dir,
    "AK1_ensemble_panelA_opening_metric.png"
  ),

  panel_A,

  width = 6.5,
  height = 5.2,
  dpi = 600,
  bg = "white"
)


ggsave(

  file.path(
    fig_dir,
    "AK1_ensemble_panelA_opening_metric.pdf"
  ),

  panel_A,

  width = 6.5,
  height = 5.2,
  bg = "white"
)


ggsave(

  file.path(
    fig_dir,
    "AK1_ensemble_panelB_Vina_scores.png"
  ),

  panel_B,

  width = 6.5,
  height = 5.2,
  dpi = 600,
  bg = "white"
)


ggsave(

  file.path(
    fig_dir,
    "AK1_ensemble_panelB_Vina_scores.pdf"
  ),

  panel_B,

  width = 6.5,
  height = 5.2,
  bg = "white"
)


ggsave(

  file.path(
    fig_dir,
    "AK1_ensemble_panelC_AP5_localization.png"
  ),

  panel_C,

  width = 6.7,
  height = 5.2,
  dpi = 600,
  bg = "white"
)


ggsave(

  file.path(
    fig_dir,
    "AK1_ensemble_panelC_AP5_localization.pdf"
  ),

  panel_C,

  width = 6.7,
  height = 5.2,
  bg = "white"
)


# ============================================================
# COMBINED A/B/C FIGURE
# ============================================================

combined <- (
  panel_A
  |
  panel_B
  |
  panel_C
) +

  plot_annotation(

    title =
      "AK1 receptor-conformation sensitivity of isorhamnetin docking",

    tag_levels = "A",

    theme = theme(

      plot.title = element_text(
        size = 17,
        face = "bold",
        hjust = 0.5
      ),

      plot.tag = element_text(
        size = 16,
        face = "bold"
      )
    )
  )


ggsave(

  file.path(
    fig_dir,
    "AK1_ensemble_quantitative_ABC_combined.png"
  ),

  combined,

  width = 18.5,
  height = 6.1,
  dpi = 600,
  bg = "white"
)


ggsave(

  file.path(
    fig_dir,
    "AK1_ensemble_quantitative_ABC_combined.pdf"
  ),

  combined,

  width = 18.5,
  height = 6.1,
  bg = "white"
)


# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("AK1 ENSEMBLE QUANTITATIVE FIGURE COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Score summary:\n"
)

print(
  score_stats
)

cat("\n")

cat(
  "AP5 localization summary:\n"
)

print(
  ap5_stats
)

cat("\n")

cat(
  "Figures saved in:\n"
)

cat(
  fig_dir,
  "\n\n"
)

cat(
  "Combined figure:\n"
)

cat(
  file.path(
    fig_dir,
    "AK1_ensemble_quantitative_ABC_combined.png"
  ),
  "\n"
)

cat(
  file.path(
    fig_dir,
    "AK1_ensemble_quantitative_ABC_combined.pdf"
  ),
  "\n\n"
)

cat(
  "Analysis complete.\n"
)