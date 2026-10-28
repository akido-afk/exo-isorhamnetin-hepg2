# ============================================================
# FIGURE 8 — PANEL D
# AK1 multi-seed docking reproducibility
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
# READ THE REAL 10-SEED DATA
# ============================================================

dat <- read_excel(
  input_file,
  sheet = "02_Seed_Results",
  skip = 2
)

# First three columns are:
# Seed / Cluster / Best Vina score
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
# SUMMARY STATISTICS
# ============================================================

mean_score <- mean(dat$Score)
median_score <- median(dat$Score)
sd_score <- sd(dat$Score)

dominant_n <- sum(
  dat$Cluster == "Dominant"
)

total_n <- nrow(dat)


cat("\n")
cat("============================================\n")
cat("AK1 MULTI-SEED REPRODUCIBILITY\n")
cat("============================================\n")

cat(
  "n =", total_n, "\n"
)

cat(
  "Mean =",
  sprintf("%.3f", mean_score),
  "kcal/mol\n"
)

cat(
  "Median =",
  sprintf("%.3f", median_score),
  "kcal/mol\n"
)

cat(
  "SD =",
  sprintf("%.3f", sd_score),
  "kcal/mol\n"
)

cat(
  "Dominant cluster =",
  dominant_n,
  "/",
  total_n,
  "\n"
)


# ============================================================
# PLOT
# ============================================================

p <- ggplot(
  dat,
  aes(
    x = Seed,
    y = Score,
    fill = Cluster
  )
) +

  # Median reference line
  geom_hline(
    yintercept = median_score,
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +

  # Individual independent docking runs
  geom_point(
    shape = 21,
    size = 5,
    stroke = 0.8,
    color = "black"
  ) +

  # Cluster colors
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

  # Keep scale tight enough to show the small variation
  scale_y_continuous(
    breaks = seq(
      floor(min(dat$Score) * 100) / 100,
      ceiling(max(dat$Score) * 100) / 100,
      by = 0.01
    ),
    expand = expansion(
      mult = c(0.12, 0.35)
    )
  ) +

  labs(
    x = "Independent Vina seed",
    y = "Best Vina score (kcal/mol)",
    fill = NULL
  ) +

  # Summary annotation
  annotate(
    "label",
    x = 7.5,
    y = max(dat$Score) + 0.012,
    label = paste0(
      "Median = ",
      sprintf("%.3f", median_score),
      " kcal/mol\n",
      "Mean = ",
      sprintf("%.3f", mean_score),
      " kcal/mol\n",
      "SD = ",
      sprintf("%.3f", sd_score),
      " kcal/mol\n",
      "Dominant cluster = ",
      dominant_n,
      "/",
      total_n
    ),
    hjust = 0,
    vjust = 0,
    size = 4.2,
    label.size = 0.35,
    fill = "white"
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
      color = "black",
      face = "bold"
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
# SAVE HIGH-RESOLUTION PNG
# ============================================================

png_file <- file.path(
  output_dir,
  "AK1_panelD_multiseed_reproducibility.png"
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


# ============================================================
# SAVE VECTOR PDF
# ============================================================

pdf_file <- file.path(
  output_dir,
  "AK1_panelD_multiseed_reproducibility.pdf"
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