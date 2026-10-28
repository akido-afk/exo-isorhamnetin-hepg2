# ============================================================
# AK1 WHOLE-PROTEIN BLIND-DOCKING VALIDATION FIGURE
#
# A = Mode-1 recovery RMSD across 10 blind-docking seeds
# B = Targeted vs blind Vina-score distributions
# C = Targeted-pose-family recovery frequencies
#
# Also creates a dedicated Excel workbook containing:
# - full 200-pose blind data
# - per-seed blind summary
# - overall blind summary
# - targeted 1Z83 scores
# - exact data used for panels A, B and C
# - method notes
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

needed_packages <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "patchwork",
  "openxlsx"
)


for (pkg in needed_packages) {

  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}


suppressPackageStartupMessages({

  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(openxlsx)

})


# ============================================================
# PATHS
# ============================================================

project_dir <- "C:/Users/akidi/Desktop/ISO_AK_Docking"


blind_dir <- file.path(
  project_dir,
  "results",
  "AK1_blind_docking",
  "full_10seed_analysis_v1"
)


targeted_dir <- file.path(
  project_dir,
  "results",
  "AK1_ensemble_sensitivity",
  "docking_analysis_v1"
)


output_dir <- file.path(
  blind_dir,
  "figures_v1"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# INPUT FILES
# ============================================================

blind_seed_file <- file.path(
  blind_dir,
  "AK1_BLIND_seed_summary.csv"
)


blind_all_file <- file.path(
  blind_dir,
  "AK1_BLIND_all_poses.csv"
)


blind_overall_file <- file.path(
  blind_dir,
  "AK1_BLIND_overall_summary.csv"
)


targeted_seed_file <- file.path(
  targeted_dir,
  "AK1_7X7S_ensemble_seed_results.csv"
)


required_files <- c(
  blind_seed_file,
  blind_all_file,
  blind_overall_file,
  targeted_seed_file
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
      "Missing required files:",
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

blind_seed <- read_csv(
  blind_seed_file,
  show_col_types = FALSE
)


blind_all <- read_csv(
  blind_all_file,
  show_col_types = FALSE
)


blind_overall <- read_csv(
  blind_overall_file,
  show_col_types = FALSE
)


targeted_all <- read_csv(
  targeted_seed_file,
  show_col_types = FALSE
)


# ============================================================
# TARGETED 1Z83 DATA ONLY
# ============================================================

targeted_closed <- targeted_all %>%

  filter(
    Receptor == "1Z83 closed"
  ) %>%

  select(
    Seed,
    Best_Vina_score
  ) %>%

  mutate(
    Search_strategy =
      "AP5-centered targeted"
  )


# ============================================================
# BLIND MODE-1 SCORES
# ============================================================

blind_scores <- blind_seed %>%

  select(
    Seed,
    Global_best_score
  ) %>%

  rename(
    Best_Vina_score =
      Global_best_score
  ) %>%

  mutate(
    Search_strategy =
      "Whole-protein blind"
  )


# ============================================================
# COMBINED SCORE DATA
# ============================================================

score_data <- bind_rows(
  targeted_closed,
  blind_scores
)


score_data$Search_strategy <- factor(

  score_data$Search_strategy,

  levels = c(
    "AP5-centered targeted",
    "Whole-protein blind"
  )
)


# ============================================================
# PANEL A DATA
# ============================================================

panel_A_data <- blind_seed %>%

  transmute(

    Seed = factor(
      Seed,
      levels = as.character(
        Seed
      )
    ),

    Seed_numeric =
      as.numeric(
        as.character(
          Seed
        )
      ),

    Mode1_RMSD_A =
      Mode1_nearest_family_RMSD_A,

    Nearest_family =
      Mode1_nearest_family,

    Recovered_2A =
      Mode1_targeted_like_2A,

    Recovered_3A =
      Mode1_targeted_like_3A
  )


mean_mode1_rmsd <- mean(
  panel_A_data$Mode1_RMSD_A
)


sd_mode1_rmsd <- sd(
  panel_A_data$Mode1_RMSD_A
)


max_mode1_rmsd <- max(
  panel_A_data$Mode1_RMSD_A
)


# ============================================================
# PANEL B SUMMARY
# ============================================================

score_summary <- score_data %>%

  group_by(
    Search_strategy
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

    Best =
      min(
        Best_Vina_score
      ),

    Worst =
      max(
        Best_Vina_score
      ),

    .groups = "drop"
  )


targeted_mean <- score_summary$Mean[
  score_summary$Search_strategy ==
    "AP5-centered targeted"
]


blind_mean <- score_summary$Mean[
  score_summary$Search_strategy ==
    "Whole-protein blind"
]


delta_mean <- (
  blind_mean -
  targeted_mean
)


# ============================================================
# PANEL C DATA
# ============================================================

n_seeds <- nrow(
  blind_seed
)


n_all_poses <- nrow(
  blind_all
)


mode1_2A_n <- sum(
  blind_seed$Mode1_targeted_like_2A
)


mode1_3A_n <- sum(
  blind_seed$Mode1_targeted_like_3A
)


all_2A_n <- sum(
  blind_all$Targeted_like_2A
)


all_3A_n <- sum(
  blind_all$Targeted_like_3A
)


panel_C_data <- data.frame(

  Population = c(
    "Global best\nmode 1",
    "Global best\nmode 1",
    "All blind\nmodes",
    "All blind\nmodes"
  ),

  Threshold = c(
    "<= 2 Å",
    "<= 3 Å",
    "<= 2 Å",
    "<= 3 Å"
  ),

  Recovered_n = c(
    mode1_2A_n,
    mode1_3A_n,
    all_2A_n,
    all_3A_n
  ),

  Total_n = c(
    n_seeds,
    n_seeds,
    n_all_poses,
    n_all_poses
  )
)


panel_C_data <- panel_C_data %>%

  mutate(

    Percent =
      100 *
      Recovered_n /
      Total_n,

    Label = paste0(
      Recovered_n,
      "/",
      Total_n,
      "\n",
      sprintf(
        "%.1f%%",
        Percent
      )
    ),

    Population = factor(

      Population,

      levels = c(
        "Global best\nmode 1",
        "All blind\nmodes"
      )
    ),

    Threshold = factor(

      Threshold,

      levels = c(
        "<= 2 Å",
        "<= 3 Å"
      )
    )
  )


# ============================================================
# PRESERVE EVERYTHING TO EXCEL
# ============================================================

excel_file <- file.path(
  output_dir,
  "AK1_BLIND_VALIDATION_FIGURE_DATA.xlsx"
)


wb <- createWorkbook()


# ------------------------------------------------------------
# Common styles
# ------------------------------------------------------------

header_style <- createStyle(

  fontColour = "#FFFFFF",
  fgFill = "#4472C4",

  halign = "center",
  valign = "center",

  textDecoration = "bold",

  border = "Bottom"
)


subheader_style <- createStyle(

  fontColour = "#FFFFFF",
  fgFill = "#5B9BD5",

  textDecoration = "bold"
)


number_style <- createStyle(
  numFmt = "0.000"
)


percent_style <- createStyle(
  numFmt = "0.0"
)


# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------

add_data_sheet <- function(
  workbook,
  sheet_name,
  data
) {

  addWorksheet(
    workbook,
    sheet_name
  )


  writeData(
    workbook,
    sheet_name,
    data,
    headerStyle = header_style
  )


  freezePane(
    workbook,
    sheet_name,
    firstRow = TRUE
  )


  setColWidths(
    workbook,
    sheet_name,
    cols = 1:ncol(data),
    widths = "auto"
  )


  if (
    ncol(data) > 0
  ) {

    addFilter(
      workbook,
      sheet_name,
      row = 1,
      cols = 1:ncol(data)
    )
  }
}


# ------------------------------------------------------------
# Full raw + summary sheets
# ------------------------------------------------------------

add_data_sheet(
  wb,
  "Blind_All_200_Poses",
  blind_all
)


add_data_sheet(
  wb,
  "Blind_Seed_Summary",
  blind_seed
)


add_data_sheet(
  wb,
  "Blind_Overall_Summary",
  blind_overall
)


add_data_sheet(
  wb,
  "Targeted_1Z83_Scores",
  targeted_closed
)


# ------------------------------------------------------------
# Exact figure data
# ------------------------------------------------------------

add_data_sheet(
  wb,
  "Figure_A_RMSD",
  panel_A_data
)


add_data_sheet(
  wb,
  "Figure_B_Scores",
  score_data
)


add_data_sheet(
  wb,
  "Figure_B_Summary",
  score_summary
)


add_data_sheet(
  wb,
  "Figure_C_Recovery",
  panel_C_data
)


# ------------------------------------------------------------
# Figure summary
# ------------------------------------------------------------

figure_summary <- data.frame(

  Metric = c(

    "Blind seeds",

    "Blind poses",

    "Mode-1 <=2 A recovery",

    "Mode-1 <=3 A recovery",

    "Mean mode-1 RMSD to targeted family (A)",

    "SD mode-1 RMSD (A)",

    "Targeted mean Vina score",

    "Blind mean Vina score",

    "Blind minus targeted mean score",

    "All poses <=2 A",

    "All poses <=3 A"
  ),

  Value = c(

    n_seeds,

    n_all_poses,

    paste0(
      mode1_2A_n,
      "/",
      n_seeds,
      " (",
      sprintf(
        "%.1f%%",
        100 * mode1_2A_n / n_seeds
      ),
      ")"
    ),

    paste0(
      mode1_3A_n,
      "/",
      n_seeds,
      " (",
      sprintf(
        "%.1f%%",
        100 * mode1_3A_n / n_seeds
      ),
      ")"
    ),

    sprintf(
      "%.3f",
      mean_mode1_rmsd
    ),

    sprintf(
      "%.3f",
      sd_mode1_rmsd
    ),

    sprintf(
      "%.3f",
      targeted_mean
    ),

    sprintf(
      "%.3f",
      blind_mean
    ),

    sprintf(
      "%.3f",
      delta_mean
    ),

    paste0(
      all_2A_n,
      "/",
      n_all_poses,
      " (",
      sprintf(
        "%.1f%%",
        100 * all_2A_n / n_all_poses
      ),
      ")"
    ),

    paste0(
      all_3A_n,
      "/",
      n_all_poses,
      " (",
      sprintf(
        "%.1f%%",
        100 * all_3A_n / n_all_poses
      ),
      ")"
    )
  )
)


add_data_sheet(
  wb,
  "Figure_Summary",
  figure_summary
)


# ------------------------------------------------------------
# Method notes
# ------------------------------------------------------------

method_notes <- data.frame(

  Item = c(

    "Blind docking",

    "Search space",

    "Sampling",

    "Targeted comparison",

    "Pose RMSD",

    "2 A threshold",

    "3 A threshold",

    "Score interpretation",

    "Main conclusion"
  ),

  Description = c(

    paste(
      "Whole-protein blind docking of neutral isorhamnetin",
      "against closed human AK1."
    ),

    paste(
      "AK1 whole-protein bounding box with 6 A padding",
      "per side. AP5 was not used to define the blind box."
    ),

    paste(
      "10 independent seeds, exhaustiveness 128,",
      "20 modes per seed, energy range 10 kcal/mol."
    ),

    paste(
      "Original AP5-centered neutral-ISO docking against",
      "the same AK1 receptor and using the same ligand preparation."
    ),

    paste(
      "Direct ligand heavy-atom positional RMSD in the",
      "shared receptor coordinate frame; no ligand superposition."
    ),

    "Stringent targeted-pose-family recovery criterion.",

    "Broader targeted-pose-family recovery criterion.",

    paste(
      "Vina values are docking scores and are not interpreted",
      "as experimental binding affinities or free energies."
    ),

    paste(
      "Whole-protein blind docking independently recovered",
      "the targeted dominant ISO pose as the global best",
      "solution in all 10 independent seeds."
    )
  )
)


add_data_sheet(
  wb,
  "Method_Notes",
  method_notes
)


# ------------------------------------------------------------
# Improve widths on text-heavy sheets
# ------------------------------------------------------------

setColWidths(
  wb,
  "Method_Notes",
  cols = 1,
  widths = 25
)


setColWidths(
  wb,
  "Method_Notes",
  cols = 2,
  widths = 95
)


setColWidths(
  wb,
  "Figure_Summary",
  cols = 1,
  widths = 45
)


setColWidths(
  wb,
  "Figure_Summary",
  cols = 2,
  widths = 28
)


saveWorkbook(
  wb,
  excel_file,
  overwrite = TRUE
)


# ============================================================
# FIGURE COLORS
# ============================================================

targeted_color <- "#F28E2B"
blind_color <- "#4E79A7"

threshold_2_color <- "#4E79A7"
threshold_3_color <- "#59A14F"


# ============================================================
# COMMON THEME
# ============================================================

theme_pub <- theme_classic(
  base_size = 12
) +

  theme(

    axis.text = element_text(
      color = "black"
    ),

    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),

    axis.text.y = element_text(
      size = 10
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
      size = 9.5
    ),

    plot.caption = element_text(
      size = 8,
      color = "gray30"
    ),

    legend.position = "bottom",

    legend.title = element_blank(),

    plot.margin = margin(
      10,
      12,
      10,
      10
    )
  )


# ============================================================
# PANEL A
# MODE-1 RMSD
# ============================================================

panel_A <- ggplot(

  panel_A_data,

  aes(
    x = Seed,
    y = Mode1_RMSD_A
  )
) +

  geom_hline(

    yintercept =
      mean_mode1_rmsd,

    linetype = "dashed",

    linewidth = 0.7,

    color = "gray40"
  ) +

  geom_segment(

    aes(
      xend = Seed,
      y = 0,
      yend = Mode1_RMSD_A
    ),

    linewidth = 0.7,

    color = "gray75"
  ) +

  geom_point(

    size = 3.8,

    color = blind_color
  ) +

  geom_text(

    aes(
      label = sprintf(
        "%.3f",
        Mode1_RMSD_A
      )
    ),

    vjust = -0.8,

    size = 3.1
  ) +

  annotate(

    "text",

    x = 1,

    y = mean_mode1_rmsd + 0.012,

    label = paste0(
      "Mean = ",
      sprintf(
        "%.3f Å",
        mean_mode1_rmsd
      )
    ),

    hjust = 0,

    size = 3.2,

    color = "gray25"
  ) +

  scale_y_continuous(

    limits = c(
      0,
      max_mode1_rmsd + 0.07
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
      "Global-best blind pose recovery",

    subtitle =
      "Mode 1 converged to the targeted dominant family in all 10 seeds",

    x =
      "Blind-docking seed",

    y =
      "RMSD to targeted pose family (Å)",

    caption =
      "All global-best blind poses were well below the 2 Å recovery criterion."
  ) +

  theme_pub +

  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        face = "plain",
        size = 9
      ),

    legend.position =
      "none"
  )


# ============================================================
# PANEL B
# TARGETED VS BLIND SCORES
# ============================================================

panel_B <- ggplot(

  score_data,

  aes(
    x = Search_strategy,
    y = Best_Vina_score,
    color = Search_strategy
  )
) +

  geom_jitter(

    width = 0.08,
    height = 0,

    size = 2.8,

    alpha = 0.9
  ) +

  geom_errorbar(

    data = score_summary,

    aes(
      x = Search_strategy,
      ymin = Mean - SD,
      ymax = Mean + SD
    ),

    inherit.aes = FALSE,

    width = 0.13,

    linewidth = 0.8,

    color = "black"
  ) +

  geom_point(

    data = score_summary,

    aes(
      x = Search_strategy,
      y = Mean
    ),

    inherit.aes = FALSE,

    shape = 23,

    size = 4.2,

    stroke = 1,

    fill = "white",

    color = "black"
  ) +

  geom_text(

    data = score_summary,

    aes(
      x = Search_strategy,
      y = Mean + 0.065,
      label = paste0(
        sprintf(
          "%.3f",
          Mean
        ),
        " ± ",
        sprintf(
          "%.3f",
          SD
        )
      )
    ),

    inherit.aes = FALSE,

    size = 3.2,

    color = "black"
  ) +

  scale_color_manual(

    values = c(

      "AP5-centered targeted" =
        targeted_color,

      "Whole-protein blind" =
        blind_color
    )
  ) +

  scale_x_discrete(

    labels = c(

      "AP5-centered targeted" =
        "Targeted",

      "Whole-protein blind" =
        "Blind"
    )
  ) +

  labs(

    title =
      "Docking-score reproducibility",

    subtitle = paste0(
      "Mean difference = ",
      sprintf(
        "%.3f kcal/mol",
        delta_mean
      )
    ),

    x =
      "Search strategy",

    y =
      "Best Vina score (kcal/mol)",

    caption =
      "Diamond = mean; error bar = ±1 SD. Same receptor, ligand and Vina scoring function."
  ) +

  theme_pub +

  theme(
    legend.position =
      "none"
  )


# ============================================================
# PANEL C
# RECOVERY FREQUENCY
# ============================================================

panel_C <- ggplot(

  panel_C_data,

  aes(
    x = Population,
    y = Percent,
    fill = Threshold
  )
) +

  geom_col(

    position = position_dodge(
      width = 0.72
    ),

    width = 0.62
  ) +

  geom_text(

    aes(
      label = Label
    ),

    position = position_dodge(
      width = 0.72
    ),

    vjust = -0.35,

    size = 3.5,

    lineheight = 0.9
  ) +

  scale_fill_manual(

    values = c(

      "<= 2 Å" =
        threshold_2_color,

      "<= 3 Å" =
        threshold_3_color
    ),

    labels = c(
      "<= 2 Å",
      "<= 3 Å"
    )
  ) +

  scale_y_continuous(

    limits = c(
      0,
      112
    ),

    breaks = c(
      0,
      20,
      40,
      60,
      80,
      100
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
      "Targeted-pose-family recovery",

    subtitle =
      "Recovery assessed after whole-protein docking",

    x =
      NULL,

    y =
      "Recovered poses (%)"
  ) +

  theme_pub


# ============================================================
# SAVE INDIVIDUAL PANELS
# ============================================================

ggsave(

  file.path(
    output_dir,
    "AK1_blind_panelA_mode1_RMSD.png"
  ),

  panel_A,

  width = 7.3,
  height = 5.3,

  dpi = 600,

  bg = "white"
)


ggsave(

  file.path(
    output_dir,
    "AK1_blind_panelA_mode1_RMSD.pdf"
  ),

  panel_A,

  width = 7.3,
  height = 5.3,

  bg = "white"
)


ggsave(

  file.path(
    output_dir,
    "AK1_blind_panelB_targeted_vs_blind_scores.png"
  ),

  panel_B,

  width = 6.6,
  height = 5.3,

  dpi = 600,

  bg = "white"
)


ggsave(

  file.path(
    output_dir,
    "AK1_blind_panelB_targeted_vs_blind_scores.pdf"
  ),

  panel_B,

  width = 6.6,
  height = 5.3,

  bg = "white"
)


ggsave(

  file.path(
    output_dir,
    "AK1_blind_panelC_recovery_frequency.png"
  ),

  panel_C,

  width = 6.5,
  height = 5.3,

  dpi = 600,

  bg = "white"
)


ggsave(

  file.path(
    output_dir,
    "AK1_blind_panelC_recovery_frequency.pdf"
  ),

  panel_C,

  width = 6.5,
  height = 5.3,

  bg = "white"
)


# ============================================================
# COMBINED FIGURE
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
      "Whole-protein blind docking independently recovers the targeted AK1 isorhamnetin pose",

    tag_levels =
      "A",

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
    output_dir,
    "AK1_BLIND_VALIDATION_ABC.png"
  ),

  combined,

  width = 19.0,
  height = 6.3,

  dpi = 600,

  bg = "white"
)


ggsave(

  file.path(
    output_dir,
    "AK1_BLIND_VALIDATION_ABC.pdf"
  ),

  combined,

  width = 19.0,
  height = 6.3,

  bg = "white"
)


# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("AK1 BLIND-DOCKING VALIDATION FIGURE COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Mode-1 recovery:\n"
)

cat(
  mode1_2A_n,
  "/",
  n_seeds,
  " seeds <=2 A\n",
  sep = ""
)

cat(
  "Mean mode-1 RMSD = ",
  sprintf(
    "%.3f",
    mean_mode1_rmsd
  ),
  " A\n\n",
  sep = ""
)


cat(
  "Targeted score = ",
  sprintf(
    "%.3f",
    targeted_mean
  ),
  " kcal/mol\n",
  sep = ""
)


cat(
  "Blind score    = ",
  sprintf(
    "%.3f",
    blind_mean
  ),
  " kcal/mol\n",
  sep = ""
)


cat(
  "Mean difference = ",
  sprintf(
    "%.3f",
    delta_mean
  ),
  " kcal/mol\n\n",
  sep = ""
)


cat(
  "All blind poses <=2 A: ",
  all_2A_n,
  "/",
  n_all_poses,
  " (",
  sprintf(
    "%.1f%%",
    100 * all_2A_n / n_all_poses
  ),
  ")\n",
  sep = ""
)


cat(
  "All blind poses <=3 A: ",
  all_3A_n,
  "/",
  n_all_poses,
  " (",
  sprintf(
    "%.1f%%",
    100 * all_3A_n / n_all_poses
  ),
  ")\n\n",
  sep = ""
)


cat(
  "Excel workbook:\n",
  excel_file,
  "\n\n"
)


cat(
  "Combined figure:\n",
  file.path(
    output_dir,
    "AK1_BLIND_VALIDATION_ABC.png"
  ),
  "\n"
)


cat(
  file.path(
    output_dir,
    "AK1_BLIND_VALIDATION_ABC.pdf"
  ),
  "\n\n"
)


cat(
  "All source data and figure data preserved.\n"
)