# ============================================================
# FIGURE 8 — PANEL E
# Pairwise RMSD heatmap of AK1–ISO best poses
# ============================================================

packages <- c("ggplot2", "dplyr", "tidyr")

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


# ============================================================
# PATHS
# ============================================================

project <- "C:/Users/akidi/Desktop/ISO_AK_Docking"

input_file <- file.path(
  project,
  "results",
  "AK1_ISO_pose_RMSD_matrix.csv"
)

output_dir <- file.path(
  project,
  "results",
  "AK1_panelE"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# READ REAL RMSD MATRIX
# ============================================================

rmsd <- read.csv(
  input_file,
  check.names = FALSE,
  row.names = 1
)

rmsd <- as.matrix(rmsd)

mode(rmsd) <- "numeric"

# make sure row and column names are identical
rownames(rmsd) <- as.character(rownames(rmsd))
colnames(rmsd) <- as.character(colnames(rmsd))


# ============================================================
# HIERARCHICAL CLUSTERING
# ============================================================

hc <- hclust(
  as.dist(rmsd),
  method = "average"
)

cluster_order <- hc$labels[hc$order]

cat("\n")
cat("============================================\n")
cat("AK1 POSE RMSD CLUSTER ORDER\n")
cat("============================================\n")
cat(cluster_order, sep = "  ")
cat("\n\n")


# ============================================================
# KNOWN CLUSTER MEMBERSHIP FROM THE ACTUAL ANALYSIS
# ============================================================

dominant <- c(
  "1001",
  "2002",
  "3003",
  "4004",
  "5005",
  "7007",
  "9009",
  "10101"
)

alternative <- c(
  "6006",
  "8008"
)


# ============================================================
# ORDER SEEDS TO MAKE THE TWO CLUSTERS VISUALLY CLEAR
#
# Within each cluster, preserve hierarchical ordering.
# ============================================================

dominant_order <- cluster_order[
  cluster_order %in% dominant
]

alternative_order <- cluster_order[
  cluster_order %in% alternative
]

display_order <- c(
  dominant_order,
  alternative_order
)


# ============================================================
# CONVERT MATRIX TO LONG FORMAT
# ============================================================

df <- as.data.frame(rmsd)

df$Seed_y <- rownames(df)

long <- df %>%
  pivot_longer(
    cols = -Seed_y,
    names_to = "Seed_x",
    values_to = "RMSD"
  )


# ============================================================
# FACTOR ORDER
#
# Reverse Y order so top-left reads naturally.
# ============================================================

long$Seed_x <- factor(
  long$Seed_x,
  levels = display_order
)

long$Seed_y <- factor(
  long$Seed_y,
  levels = rev(display_order)
)


# ============================================================
# AXIS LABELS
# Mark representative seed 9009 with a star
# ============================================================

x_labels <- setNames(
  display_order,
  display_order
)

x_labels["9009"] <- "9009 *"

y_order <- rev(display_order)

y_labels <- setNames(
  y_order,
  y_order
)

y_labels["9009"] <- "9009 *"


# ============================================================
# CLUSTER BOUNDARY
# ============================================================

n_dom <- length(dominant_order)
n_total <- length(display_order)

# Because Y axis is reversed, dominant block occupies
# upper-left / lower-right symmetrically in matrix.
# We'll draw separation lines after the 8th pose.


# ============================================================
# PLOT
# ============================================================

p <- ggplot(
  long,
  aes(
    x = Seed_x,
    y = Seed_y,
    fill = RMSD
  )
) +

  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +

  # numeric RMSD values inside cells
  geom_text(
    aes(
      label = sprintf("%.2f", RMSD)
    ),
    size = 3.2,
    color = "black"
  ) +

  # RMSD color scale
  scale_fill_gradientn(
    colours = c(
      "#2166AC",
      "#67A9CF",
      "#D1E5F0",
      "#F7F7F7",
      "#FDDBC7",
      "#EF8A62",
      "#B2182B"
    ),
    values = scales::rescale(
      c(
        0,
        0.5,
        1,
        2,
        3,
        5,
        max(long$RMSD)
      )
    ),
    limits = c(
      0,
      max(long$RMSD)
    ),
    name = "Pairwise RMSD\n(Å)"
  ) +

  # Separation between dominant 8 and alternative 2
  geom_vline(
    xintercept = n_dom + 0.5,
    linewidth = 1.2,
    color = "black"
  ) +

  # Y separation needs opposite coordinate because Y is reversed
  geom_hline(
    yintercept = length(alternative_order) + 0.5,
    linewidth = 1.2,
    color = "black"
  ) +

  scale_x_discrete(
    labels = x_labels,
    drop = FALSE
  ) +

  scale_y_discrete(
    labels = y_labels,
    drop = FALSE
  ) +

  labs(
    x = "Independent Vina seed",
    y = "Independent Vina seed"
  ) +

  coord_fixed() +

  theme_minimal(
    base_size = 14
  ) +

  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black",
      size = 11
    ),

    axis.text.y = element_text(
      color = "black",
      size = 11
    ),

    axis.title = element_text(
      face = "bold",
      color = "black"
    ),

    legend.title = element_text(
      face = "bold"
    ),

    legend.position = "right",

    plot.margin = margin(
      12,
      15,
      12,
      12
    )
  )


# ============================================================
# SAVE BASE HEATMAP
# ============================================================

png_file <- file.path(
  output_dir,
  "AK1_panelE_pose_RMSD_heatmap.png"
)

ggsave(
  filename = png_file,
  plot = p,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


pdf_file <- file.path(
  output_dir,
  "AK1_panelE_pose_RMSD_heatmap.pdf"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 7.5,
  height = 6.5,
  units = "in"
)


# ============================================================
# PRINT SUMMARY
# ============================================================

pairwise_values <- rmsd[
  upper.tri(rmsd)
]

cat("Mean pairwise RMSD =",
    sprintf("%.3f", mean(pairwise_values)),
    "A\n")

cat("Minimum pairwise RMSD =",
    sprintf("%.3f", min(pairwise_values)),
    "A\n")

cat("Maximum pairwise RMSD =",
    sprintf("%.3f", max(pairwise_values)),
    "A\n")

cat("\nDominant cluster (8/10):\n")
cat(dominant, sep = " ")
cat("\n")

cat("\nAlternative cluster (2/10):\n")
cat(alternative, sep = " ")
cat("\n")

cat("\nRepresentative dominant pose: seed 9009\n")
cat("Mean RMSD within dominant cluster: 0.039 A\n")

cat("\nSaved PNG:\n")
cat(png_file, "\n")

cat("\nSaved PDF:\n")
cat(pdf_file, "\n")