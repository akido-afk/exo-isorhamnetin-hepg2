# ============================================================
# FIGURE 8 — PANEL D
# AK1 multi-seed docking reproducibility
# Axis scale
# ============================================================

packages <- c("readxl", "ggplot2", "dplyr")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

library(readxl)
library(ggplot2)
library(dplyr)


# ============================================================
# PATHS
# ============================================================

project <- "C:/Users/akidi/Desktop/ISO_AK_Docking"

input_file <- file.path(
  project,
  "results",
  "AK1_ISO_COMPREHENSIVE_DOCKING_RESULTS.xlsx"
)

output_dir <- file.path(
  project,
  "results",
  "AK1_panelD"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# READ REAL DATA
# ============================================================

dat <- read_excel(
  input_file,
  sheet = "02_Seed_Results",
  skip = 2
)

names(dat)[1:3] <- c(
  "Seed",
  "Cluster",
  "Score"
)

dat <- dat %>%
  select(Seed, Cluster, Score) %>%
  filter(
    !is.na(Seed),
    !is.na(Score)
  )

dat$Seed <- factor(
  dat$Seed,
  levels = dat$Seed
)

dat$Cluster <- factor(
  dat$Cluster,
  levels = c(
    "Dominant",
    "Alternative"
  )
)


# ============================================================
# STATISTICS
# ============================================================

mean_score <- mean(dat$Score)
median_score <- median(dat$Score)
sd_score <- sd(dat$Score)

best_score <- min(dat$Score)
worst_score <- max(dat$Score)
score_range <- worst_score - best_score

dominant_n <- sum(
  dat$Cluster == "Dominant"
)

total_n <- nrow(dat)


cat("\n")
cat("============================================\n")
cat("AK1 MULTI-SEED REPRODUCIBILITY\n")
cat("============================================\n")

cat("n =", total_n, "\n")
cat("Mean =", sprintf("%.3f", mean_score), "kcal/mol\n")
cat("Median =", sprintf("%.3f", median_score), "kcal/mol\n")
cat("SD =", sprintf("%.3f", sd_score), "kcal/mol\n")
cat("Range =", sprintf("%.3f", score_range), "kcal/mol\n")
cat(
  "Dominant pose cluster =",
  dominant_n,
  "/",
  total_n,
  "\n"
)


# ============================================================
# PLOT
#
# IMPORTANT:
# Fixed broader Y scale prevents exaggeration of
# the very small score differences.
# ============================================================

p <- ggplot(
  dat,
  aes(
    x = Seed,
    y = Score,
    fill = Cluster
  )
) +

  # Mean ± SD band
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = mean_score - sd_score,
    ymax = mean_score + sd_score,
    fill = "grey70",
    alpha = 0.22
  ) +

  # Mean line
  geom_hline(
    yintercept = mean_score,
    linewidth = 0.7,
    color = "grey35"
  ) +

  # Median line
  geom_hline(
    yintercept = median_score,
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +

  # Seed values
  geom_point(
    shape = 21,
    size = 5,
    stroke = 0.8,
    color = "black"
  ) +

  # Colors represent POSE CLUSTER, not energetic groups
  scale_fill_manual(
    values = c(
      "Dominant" = "#E69F00",
      "Alternative" = "#56B4E9"
    ),
    labels = c(
      "Dominant pose cluster (8/10)",
      "Alternative pose cluster (2/10)"
    )
  ) +

  # BROADER FIXED Y SCALE
  coord_cartesian(
    ylim = c(-9.00, -8.50)
  ) +

  scale_y_continuous(
    breaks = seq(
      -9.0,
      -8.5,
      by = 0.1
    )
  ) +

  labs(
    x = "Independent Vina seed",
    y = "Best Vina score (kcal/mol)",
    fill = NULL
  ) +

  # Statistical summary
  annotate(
    "label",
    x = 6.1,
    y = -8.54,
    label = paste0(
      "Mean = ",
      sprintf("%.3f", mean_score),
      " kcal/mol\n",
      "Median = ",
      sprintf("%.3f", median_score),
      " kcal/mol\n",
      "SD = ",
      sprintf("%.3f", sd_score),
      " kcal/mol\n",
      "Total range = ",
      sprintf("%.3f", score_range),
      " kcal/mol"
    ),
    hjust = 0,
    vjust = 1,
    size = 4.0,
    label.size = 0.3,
    fill = "white"
  ) +

  # Key interpretation
  annotate(
    "text",
    x = 5.5,
    y = -8.94,
    label = "Highly consistent scores across 10 independent searches",
    fontface = "bold",
    size = 4.3
  ) +

  theme_classic(
    base_size = 15
  ) +

  theme(
    axis.text = element_text(
      color = "black"
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),

    axis.title = element_text(
      face = "bold",
      color = "black"
    ),

    legend.position = "top",

    legend.text = element_text(
      size = 11
    ),

    plot.margin = margin(
      12,
      18,
      12,
      12
    )
  )


# ============================================================
# SAVE AS NEW FILES
# ============================================================

png_file <- file.path(
  output_dir,
  "AK1_panelD_multiseed_reproducibility_V2.png"
)

ggsave(
  filename = png_file,
  plot = p,
  width = 7.2,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


pdf_file <- file.path(
  output_dir,
  "AK1_panelD_multiseed_reproducibility_V2.pdf"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 7.2,
  height = 5.5,
  units = "in"
)


cat("\nSaved PNG:\n")
cat(png_file, "\n")

cat("\nSaved PDF:\n")
cat(pdf_file, "\n")