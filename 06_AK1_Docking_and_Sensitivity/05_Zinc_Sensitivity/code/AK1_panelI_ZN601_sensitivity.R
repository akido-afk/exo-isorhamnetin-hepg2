# ============================================================
# FIGURE 8 — PANEL I
# AK1 receptor-state sensitivity to crystallographic ZN601
# ============================================================

packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
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
library(tidyr)
library(patchwork)


# ============================================================
# PATHS
# ============================================================

project <- "C:/Users/akidi/Desktop/ISO_AK_Docking"

base_dir <- file.path(
  project,
  "results",
  "AK1_ZN_sensitivity",
  "ZN601_10seed"
)

score_file <- file.path(
  base_dir,
  "AK1_ZN601_scores.csv"
)

rmsd_file <- file.path(
  base_dir,
  "AK1_ZN601_pose_RMSD_matrix.csv"
)

family_file <- file.path(
  base_dir,
  "AK1_ZN601_vs_Znfree_families.csv"
)

zn_file <- file.path(
  base_dir,
  "AK1_ZN601_Zn_distances.csv"
)

output_dir <- file.path(
  project,
  "results",
  "AK1_panelI"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# READ DATA
# ============================================================

scores <- read.csv(
  score_file,
  check.names = FALSE
)

families <- read.csv(
  family_file,
  check.names = FALSE
)

zn <- read.csv(
  zn_file,
  check.names = FALSE
)

rmsd <- read.csv(
  rmsd_file,
  check.names = FALSE,
  row.names = 1
)

rmsd <- as.matrix(rmsd)

mode(rmsd) <- "numeric"


# ============================================================
# SUMMARY VALUES
# ============================================================

znfree_mean <- mean(
  scores$`Zn-free score`
)

zn_mean <- mean(
  scores$`ZN601-retained score`
)

paired_shift <- mean(
  scores$`Score difference ZN601-minus-Znfree`
)

zn_dist_mean <- mean(
  zn$`Closest polar distance to ZN601 (A)`
)

zn_dist_min <- min(
  zn$`Closest polar distance to ZN601 (A)`
)

zn_dist_max <- max(
  zn$`Closest polar distance to ZN601 (A)`
)


# ============================================================
# I1 — PAIRED SCORE COMPARISON
# ============================================================

score_long <- scores %>%
  select(
    Seed,
    `Zn-free score`,
    `ZN601-retained score`
  ) %>%
  pivot_longer(
    cols = -Seed,
    names_to = "State",
    values_to = "Score"
  ) %>%
  mutate(
    State = recode(
      State,
      `Zn-free score` = "Zn-free",
      `ZN601-retained score` = "ZN601 retained"
    )
  )


p1 <- ggplot(
  score_long,
  aes(
    x = State,
    y = Score,
    group = Seed
  )
) +

  geom_line(
    color = "grey65",
    linewidth = 0.7
  ) +

  geom_point(
    aes(fill = State),
    shape = 21,
    size = 3.8,
    color = "black",
    stroke = 0.7
  ) +

  scale_fill_manual(
    values = c(
      "Zn-free" = "#E69F00",
      "ZN601 retained" = "#56B4E9"
    ),
    guide = "none"
  ) +

  coord_cartesian(
    ylim = c(
      -8.86,
      -8.72
    )
  ) +

  labs(
    x = NULL,
    y = "Best Vina score (kcal/mol)"
  ) +

  annotate(
    "label",
    x = 1.5,
    y = -8.735,
    label = paste0(
      "Mean shift = ",
      sprintf("%.3f", paired_shift),
      " kcal/mol"
    ),
    size = 3.8,
    fill = "white"
  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(
    axis.text = element_text(
      color = "black"
    ),

    axis.title = element_text(
      face = "bold"
    )
  )


# ============================================================
# I2 — ZN601 10×10 RMSD HEATMAP
# ============================================================

rmsd_df <- as.data.frame(rmsd)

rmsd_df$Seed_y <- rownames(rmsd_df)

rmsd_long <- rmsd_df %>%
  pivot_longer(
    cols = -Seed_y,
    names_to = "Seed_x",
    values_to = "RMSD"
  )


seed_order <- colnames(rmsd)

rmsd_long$Seed_x <- factor(
  rmsd_long$Seed_x,
  levels = seed_order
)

rmsd_long$Seed_y <- factor(
  rmsd_long$Seed_y,
  levels = rev(seed_order)
)


p2 <- ggplot(
  rmsd_long,
  aes(
    x = Seed_x,
    y = Seed_y,
    fill = RMSD
  )
) +

  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +

  geom_text(
    aes(
      label = sprintf("%.2f", RMSD)
    ),
    size = 2.5
  ) +

  scale_fill_gradient(
    low = "#2166AC",
    high = "#F7F7F7",
    name = "RMSD\n(Å)"
  ) +

  labs(
    x = "ZN601 Vina seed",
    y = "ZN601 Vina seed"
  ) +

  coord_fixed() +

  theme_minimal(
    base_size = 11
  ) +

  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),

    axis.text.y = element_text(
      color = "black"
    ),

    axis.title = element_text(
      face = "bold"
    )
  )


# ============================================================
# I3 — DISTANCE TO ORIGINAL Zn-FREE POSE FAMILIES
# ============================================================

family_long <- families %>%
  select(
    `ZN601 seed`,
    `Mean RMSD to original dominant (A)`,
    `Mean RMSD to original alternative (A)`
  ) %>%
  pivot_longer(
    cols = -`ZN601 seed`,
    names_to = "Original_family",
    values_to = "RMSD"
  ) %>%
  mutate(
    Original_family = recode(
      Original_family,
      `Mean RMSD to original dominant (A)` =
        "Original dominant\nZn-free family",
      `Mean RMSD to original alternative (A)` =
        "Original alternative\nZn-free family"
    )
  )


p3 <- ggplot(
  family_long,
  aes(
    x = Original_family,
    y = RMSD
  )
) +

  geom_jitter(
    width = 0.08,
    size = 3,
    alpha = 0.65,
    color = "grey30"
  ) +

  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 5,
    fill = "#56B4E9",
    color = "black"
  ) +

  labs(
    x = NULL,
    y = "RMSD to Zn-free pose family (Å)"
  ) +

  annotate(
    "label",
    x = 1.5,
    y = max(family_long$RMSD) + 0.15,
    label = paste0(
      "ZN601–ISO polar distance:\n",
      sprintf("%.3f", zn_dist_mean),
      " Å mean (",
      sprintf("%.3f", zn_dist_min),
      "–",
      sprintf("%.3f", zn_dist_max),
      " Å)"
    ),
    size = 3.5,
    fill = "white"
  ) +

  coord_cartesian(
    ylim = c(
      5.5,
      8.25
    )
  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(
    axis.text = element_text(
      color = "black"
    ),

    axis.title = element_text(
      face = "bold"
    )
  )


# ============================================================
# COMBINE
# ============================================================

combined <- (
  p1 |
  p2 |
  p3
) +

  plot_annotation(
    tag_levels = "a"
  ) &

  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    )
  )


# ============================================================
# SAVE
# ============================================================

png_file <- file.path(
  output_dir,
  "AK1_panelI_ZN601_receptor_state_sensitivity.png"
)

pdf_file <- file.path(
  output_dir,
  "AK1_panelI_ZN601_receptor_state_sensitivity.pdf"
)


ggsave(
  filename = png_file,
  plot = combined,
  width = 15,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = combined,
  width = 15,
  height = 5.5,
  units = "in"
)


cat("\nSaved:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")