# ============================================================
# AK1 — ISO vs AP5 DOCKING-SCORE CONTEXT
# Publication version V3
#
# 10 independent Vina seeds per ligand.
# ISO dominant and alternative pose families distinguished.
#
# IMPORTANT:
# Vina scores are computational docking scores,
# NOT experimentally measured affinities/free energies.
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

packages <- c(
  "ggplot2",
  "dplyr",
  "openxlsx",
  "patchwork"
)

for (p in packages) {

  if (!requireNamespace(p, quietly = TRUE)) {

    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

library(ggplot2)
library(dplyr)
library(openxlsx)
library(patchwork)


# ============================================================
# PATHS
# ============================================================

docking_dir <- "../docking"

out_dir <- "../results/AK1_AP5_vs_ISO_scores"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# SEEDS
# ============================================================

seeds <- c(
  1001,
  2002,
  3003,
  4004,
  5005,
  6006,
  7007,
  8008,
  9009,
  10101
)


# ============================================================
# ISO CLUSTER MEMBERSHIP
# ============================================================

iso_dominant <- c(
  1001,
  2002,
  3003,
  4004,
  5005,
  7007,
  9009,
  10101
)

iso_alternative <- c(
  6006,
  8008
)


# ============================================================
# FUNCTION — READ BEST VINA SCORE
# ============================================================

read_best_vina_score <- function(file) {

  if (!file.exists(file)) {

    stop(
      paste(
        "File not found:",
        file
      )
    )
  }

  lines <- readLines(
    file,
    warn = FALSE
  )

  result_lines <- grep(
    "^REMARK VINA RESULT:",
    lines,
    value = TRUE
  )

  if (length(result_lines) == 0) {

    stop(
      paste(
        "No Vina result found in:",
        file
      )
    )
  }

  first_line <- result_lines[1]

  parts <- strsplit(
    trimws(first_line),
    "\\s+"
  )[[1]]

  as.numeric(
    parts[4]
  )
}


# ============================================================
# READ ISO SCORES
# ============================================================

iso_rows <- data.frame()

for (seed in seeds) {

  file <- file.path(
    docking_dir,
    paste0(
      "AK1_ISO_seed",
      seed,
      ".pdbqt"
    )
  )

  score <- read_best_vina_score(
    file
  )

  cluster <- ifelse(
    seed %in% iso_dominant,
    "ISO dominant",
    "ISO alternative"
  )

  iso_rows <- rbind(
    iso_rows,
    data.frame(
      Ligand = "ISO",
      Seed = seed,
      Score = score,
      Group = cluster
    )
  )
}


# ============================================================
# READ AP5 SCORES
# ============================================================

ap5_rows <- data.frame()

for (seed in seeds) {

  file <- file.path(
    docking_dir,
    paste0(
      "AK1_AP5_seed",
      seed,
      ".pdbqt"
    )
  )

  score <- read_best_vina_score(
    file
  )

  ap5_rows <- rbind(
    ap5_rows,
    data.frame(
      Ligand = "AP5",
      Seed = seed,
      Score = score,
      Group = "AP5 reference"
    )
  )
}


# ============================================================
# COMBINE
# ============================================================

all_scores <- bind_rows(
  iso_rows,
  ap5_rows
)

all_scores <- all_scores %>%
  mutate(

    Ligand = factor(
      Ligand,
      levels = c(
        "ISO",
        "AP5"
      )
    ),

    Group = factor(
      Group,
      levels = c(
        "ISO dominant",
        "ISO alternative",
        "AP5 reference"
      )
    )
  )


# ============================================================
# SUMMARY STATISTICS
# ============================================================

iso_mean <- mean(
  iso_rows$Score
)

iso_sd <- sd(
  iso_rows$Score
)

iso_median <- median(
  iso_rows$Score
)


iso_dom_scores <- iso_rows$Score[
  iso_rows$Group == "ISO dominant"
]

iso_alt_scores <- iso_rows$Score[
  iso_rows$Group == "ISO alternative"
]


iso_dom_mean <- mean(
  iso_dom_scores
)

iso_dom_sd <- sd(
  iso_dom_scores
)


iso_alt_mean <- mean(
  iso_alt_scores
)

iso_alt_sd <- sd(
  iso_alt_scores
)


ap5_mean <- mean(
  ap5_rows$Score
)

ap5_sd <- sd(
  ap5_rows$Score
)

ap5_median <- median(
  ap5_rows$Score
)


mean_score_separation <- abs(
  ap5_mean - iso_mean
)


# ============================================================
# SUMMARY TABLE
# ============================================================

summary_complete <- data.frame(

  Group = c(
    "ISO all 10",
    "ISO dominant 8",
    "ISO alternative 2",
    "AP5 all 10"
  ),

  n = c(
    10,
    8,
    2,
    10
  ),

  Mean = c(
    iso_mean,
    iso_dom_mean,
    iso_alt_mean,
    ap5_mean
  ),

  Median = c(
    iso_median,
    median(iso_dom_scores),
    median(iso_alt_scores),
    ap5_median
  ),

  SD = c(
    iso_sd,
    iso_dom_sd,
    iso_alt_sd,
    ap5_sd
  ),

  Best = c(
    min(iso_rows$Score),
    min(iso_dom_scores),
    min(iso_alt_scores),
    min(ap5_rows$Score)
  ),

  Worst = c(
    max(iso_rows$Score),
    max(iso_dom_scores),
    max(iso_alt_scores),
    max(ap5_rows$Score)
  ),

  Range = c(
    max(iso_rows$Score) -
      min(iso_rows$Score),

    max(iso_dom_scores) -
      min(iso_dom_scores),

    max(iso_alt_scores) -
      min(iso_alt_scores),

    max(ap5_rows$Score) -
      min(ap5_rows$Score)
  )
)


# ============================================================
# SEEDWISE COMPARISON
# ============================================================

paired <- iso_rows %>%
  select(
    Seed,
    ISO_score = Score,
    ISO_cluster = Group
  ) %>%
  left_join(
    ap5_rows %>%
      select(
        Seed,
        AP5_score = Score
      ),
    by = "Seed"
  ) %>%
  mutate(
    AP5_minus_ISO =
      AP5_score - ISO_score
  )


# ============================================================
# SAVE CSV FILES
# ============================================================

scores_csv <- file.path(
  out_dir,
  "AK1_AP5_vs_ISO_all_seed_scores.csv"
)

summary_csv <- file.path(
  out_dir,
  "AK1_AP5_vs_ISO_summary_statistics.csv"
)

paired_csv <- file.path(
  out_dir,
  "AK1_AP5_vs_ISO_seedwise_comparison.csv"
)


write.csv(
  all_scores,
  scores_csv,
  row.names = FALSE
)

write.csv(
  summary_complete,
  summary_csv,
  row.names = FALSE
)

write.csv(
  paired,
  paired_csv,
  row.names = FALSE
)


# ============================================================
# SAVE EXCEL
# ============================================================

xlsx_file <- file.path(
  out_dir,
  "AK1_AP5_vs_ISO_10seed_score_comparison.xlsx"
)

wb <- createWorkbook()


addWorksheet(
  wb,
  "All_Scores"
)

writeData(
  wb,
  "All_Scores",
  all_scores
)


addWorksheet(
  wb,
  "Summary"
)

writeData(
  wb,
  "Summary",
  summary_complete
)


addWorksheet(
  wb,
  "Seedwise_Comparison"
)

writeData(
  wb,
  "Seedwise_Comparison",
  paired
)


addWorksheet(
  wb,
  "Notes"
)


notes <- data.frame(

  Item = c(
    "Purpose",
    "ISO",
    "AP5",
    "AP5 redocking validation",
    "Score interpretation",
    "Comparison limitation",
    "Next analysis"
  ),

  Description = c(

    paste0(
      "Docking-score context for ISO relative to the ",
      "crystallographic AK1 reference ligand AP5."
    ),

    paste0(
      "Ten independent Zn-free ISO docking searches. ",
      "Eight formed the dominant pose family and two ",
      "formed the alternative family."
    ),

    paste0(
      "Ten AP5 searches performed using the same AK1 receptor, ",
      "box, exhaustiveness, number of modes, and seed series."
    ),

    paste0(
      "Previous AP5 redocking RMSD relative to crystallographic ",
      "AP5 = 1.776 A."
    ),

    paste0(
      "Vina values are computational docking scores and should ",
      "not be interpreted as experimentally measured affinities."
    ),

    paste0(
      "ISO and AP5 differ markedly in molecular size, charge, ",
      "and chemistry. Raw score differences therefore provide ",
      "scoring context rather than quantitative experimental ",
      "binding-affinity ratios."
    ),

    paste0(
      "The ten AP5 best poses should next be compared against ",
      "crystallographic AP5 because the AP5 score distribution ",
      "appears to contain more than one score band."
    )
  )
)


writeData(
  wb,
  "Notes",
  notes
)


setColWidths(
  wb,
  "All_Scores",
  cols = 1:ncol(all_scores),
  widths = "auto"
)

setColWidths(
  wb,
  "Summary",
  cols = 1:ncol(summary_complete),
  widths = "auto"
)

setColWidths(
  wb,
  "Seedwise_Comparison",
  cols = 1:ncol(paired),
  widths = "auto"
)

setColWidths(
  wb,
  "Notes",
  cols = 1,
  widths = 28
)

setColWidths(
  wb,
  "Notes",
  cols = 2,
  widths = 105
)


saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


# ============================================================
# MEAN ± SD DATA
# ============================================================

mean_sd_df <- all_scores %>%
  group_by(
    Ligand
  ) %>%
  summarise(
    Mean = mean(Score),
    SD = sd(Score),
    .groups = "drop"
  )


# ============================================================
# FIXED JITTER POSITIONS
#
# Slightly wider than before so the 10 ISO values can be seen.
# ============================================================

set.seed(1001)

all_scores <- all_scores %>%
  mutate(

    x_numeric = case_when(
      Ligand == "ISO" ~ 1,
      Ligand == "AP5" ~ 2
    ),

    x_jitter =
      x_numeric +
      runif(
        n(),
        min = -0.11,
        max = 0.11
      )
  )


# ============================================================
# MAIN FIGURE
# ============================================================

p <- ggplot() +

  # ----------------------------------------------------------
  # ISO DOMINANT
  # ----------------------------------------------------------

  geom_point(
    data = subset(
      all_scores,
      Group == "ISO dominant"
    ),
    aes(
      x = x_jitter,
      y = Score
    ),
    shape = 21,
    size = 3.8,
    stroke = 0.8,
    fill = "#E69F00",
    color = "black"
  ) +

  # ----------------------------------------------------------
  # ISO ALTERNATIVE
  # ----------------------------------------------------------

  geom_point(
    data = subset(
      all_scores,
      Group == "ISO alternative"
    ),
    aes(
      x = x_jitter,
      y = Score
    ),
    shape = 24,
    size = 4.0,
    stroke = 0.8,
    fill = "#56B4E9",
    color = "black"
  ) +

  # ----------------------------------------------------------
  # AP5
  # ----------------------------------------------------------

  geom_point(
    data = subset(
      all_scores,
      Group == "AP5 reference"
    ),
    aes(
      x = x_jitter,
      y = Score
    ),
    shape = 22,
    size = 3.8,
    stroke = 0.8,
    fill = "#00BFC4",
    color = "black"
  ) +

  # ----------------------------------------------------------
  # MEAN ± SD
  # ----------------------------------------------------------

  geom_errorbar(
    data = mean_sd_df,
    aes(
      x = as.numeric(Ligand),
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    width = 0.19,
    linewidth = 1.0,
    color = "black"
  ) +

  geom_point(
    data = mean_sd_df,
    aes(
      x = as.numeric(Ligand),
      y = Mean
    ),
    shape = 23,
    size = 5.4,
    stroke = 1.15,
    fill = "white",
    color = "black"
  ) +

  # ==========================================================
  # ISO SUMMARY
  #
  # MOVED WELL BELOW THE ISO POINT CLOUD.
  # ==========================================================

  annotate(
    "label",
    x = 1.00,
    y = -9.18,
    label = paste0(
      "ISO\n",
      "all 10: ",
      sprintf("%.3f", iso_mean),
      " ± ",
      sprintf("%.3f", iso_sd),
      "\n",
      "dominant 8: ",
      sprintf("%.3f", iso_dom_mean),
      " ± ",
      sprintf("%.3f", iso_dom_sd)
    ),
    size = 3.65,
    lineheight = 1.05,
    fill = "white",
    label.size = 0.3
  ) +

  # ==========================================================
  # AP5 SUMMARY
  #
  # MOVED BELOW THE AP5 DISTRIBUTION.
  # ==========================================================

  annotate(
    "label",
    x = 2.00,
    y = -11.66,
    label = paste0(
      "AP5\n",
      sprintf("%.3f", ap5_mean),
      " ± ",
      sprintf("%.3f", ap5_sd)
    ),
    size = 3.65,
    lineheight = 1.05,
    fill = "white",
    label.size = 0.3
  ) +

  # ==========================================================
  # MEAN SCORE SEPARATION
  # ==========================================================

  annotate(
    "segment",
    x = 1.18,
    xend = 1.82,
    y = -10.02,
    yend = -10.02,
    linewidth = 0.65,
    color = "grey35"
  ) +

  annotate(
    "segment",
    x = 1.18,
    xend = 1.18,
    y = -9.96,
    yend = -10.08,
    linewidth = 0.65,
    color = "grey35"
  ) +

  annotate(
    "segment",
    x = 1.82,
    xend = 1.82,
    y = -9.96,
    yend = -10.08,
    linewidth = 0.65,
    color = "grey35"
  ) +

  annotate(
    "label",
    x = 1.50,
    y = -10.02,
    label = paste0(
      "\u0394 mean score = ",
      sprintf(
        "%.2f",
        mean_score_separation
      ),
      " kcal/mol"
    ),
    size = 3.55,
    fill = "white",
    label.size = 0
  ) +

  # ==========================================================
  # AXES
  # ==========================================================

  scale_x_continuous(
    limits = c(
      0.55,
      2.45
    ),
    breaks = c(
      1,
      2
    ),
    labels = c(
      "ISO",
      "AP5"
    ),
    expand = c(
      0,
      0
    )
  ) +

  scale_y_continuous(
    limits = c(
      -11.90,
      -8.40
    ),
    breaks = seq(
      -11.5,
      -8.5,
      by = 0.5
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(

    title =
      "AK1 docking-score context",

    subtitle =
      paste0(
        "ISO and crystallographic AP5 reference across ",
        "10 independent Vina searches"
      ),

    x = NULL,

    y =
      "Best Vina docking score (kcal/mol)"
  ) +

  theme_classic(
    base_size = 15
  ) +

  theme(

    axis.text.x =
      element_text(
        size = 17,
        face = "bold",
        color = "black",
        margin = margin(
          t = 8
        )
      ),

    axis.text.y =
      element_text(
        size = 13,
        color = "black"
      ),

    axis.title.y =
      element_text(
        size = 16,
        face = "bold",
        margin = margin(
          r = 10
        )
      ),

    axis.line =
      element_line(
        linewidth = 0.8,
        color = "black"
      ),

    axis.ticks =
      element_line(
        linewidth = 0.7,
        color = "black"
      ),

    plot.title =
      element_text(
        size = 19,
        face = "bold",
        hjust = 0,
        margin = margin(
          b = 4
        )
      ),

    plot.subtitle =
      element_text(
        size = 11.5,
        hjust = 0,
        color = "grey25",
        margin = margin(
          b = 15
        )
      ),

    plot.margin =
      margin(
        t = 20,
        r = 25,
        b = 16,
        l = 22
      ),

    plot.title.position =
      "plot"
  )


# ============================================================
# COMPACT BOTTOM LEGEND
# ============================================================

legend_df <- data.frame(

  x = c(
    1.0,
    2.25,
    3.65,
    5.10
  ),

  y = c(
    1,
    1,
    1,
    1
  ),

  Label = c(
    "ISO dominant",
    "ISO alternative",
    "AP5 reference",
    "Mean ± SD"
  )
)


legend_plot <- ggplot() +

  # ISO dominant
  geom_point(
    aes(
      x = 1.0,
      y = 1
    ),
    shape = 21,
    size = 3.8,
    fill = "#E69F00",
    color = "black"
  ) +

  # ISO alternative
  geom_point(
    aes(
      x = 2.25,
      y = 1
    ),
    shape = 24,
    size = 4.0,
    fill = "#56B4E9",
    color = "black"
  ) +

  # AP5
  geom_point(
    aes(
      x = 3.65,
      y = 1
    ),
    shape = 22,
    size = 3.8,
    fill = "#00BFC4",
    color = "black"
  ) +

  # Mean
  geom_point(
    aes(
      x = 5.10,
      y = 1
    ),
    shape = 23,
    size = 4.5,
    fill = "white",
    color = "black"
  ) +

  geom_text(
    aes(
      x = 1.12,
      y = 1,
      label = "ISO dominant"
    ),
    hjust = 0,
    size = 4.0
  ) +

  geom_text(
    aes(
      x = 2.37,
      y = 1,
      label = "ISO alternative"
    ),
    hjust = 0,
    size = 4.0
  ) +

  geom_text(
    aes(
      x = 3.77,
      y = 1,
      label = "AP5 reference"
    ),
    hjust = 0,
    size = 4.0
  ) +

  geom_text(
    aes(
      x = 5.22,
      y = 1,
      label = "Mean ± SD"
    ),
    hjust = 0,
    size = 4.0
  ) +

  xlim(
    0.8,
    6.15
  ) +

  ylim(
    0.75,
    1.25
  ) +

  theme_void()


# ============================================================
# COMBINE
# ============================================================

final_plot <- p / legend_plot +

  plot_layout(
    heights = c(
      14,
      1
    )
  )


# ============================================================
# SAVE FIGURE V3
# ============================================================

png_file <- file.path(
  out_dir,
  "AK1_AP5_vs_ISO_10seed_score_comparison_V3.png"
)

pdf_file <- file.path(
  out_dir,
  "AK1_AP5_vs_ISO_10seed_score_comparison_V3.pdf"
)


ggsave(
  filename = png_file,
  plot = final_plot,
  width = 8.8,
  height = 8.2,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = pdf_file,
  plot = final_plot,
  width = 8.8,
  height = 8.2,
  units = "in",
  bg = "white"
)


# ============================================================
# FINAL CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("AK1 ISO vs AP5 SCORE COMPARISON — V3\n")
cat("============================================================\n\n")

cat(
  sprintf(
    "ISO all 10:      %.3f ± %.3f kcal/mol\n",
    iso_mean,
    iso_sd
  )
)

cat(
  sprintf(
    "ISO dominant 8:  %.3f ± %.3f kcal/mol\n",
    iso_dom_mean,
    iso_dom_sd
  )
)

cat(
  sprintf(
    "AP5 all 10:      %.3f ± %.3f kcal/mol\n",
    ap5_mean,
    ap5_sd
  )
)

cat(
  sprintf(
    "Mean-score separation: %.3f kcal/mol\n",
    mean_score_separation
  )
)

cat("\nSaved:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat(xlsx_file, "\n")
cat(scores_csv, "\n")
cat(summary_csv, "\n")
cat(paired_csv, "\n")