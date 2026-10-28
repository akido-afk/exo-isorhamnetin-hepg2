# ============================================================
# FIGURE 8 — PANEL E, EXPLANATORY VERSION
# Pairwise RMSD among the top-ranked AK1–ISO poses
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
# READ RMSD MATRIX
# ============================================================

rmsd <- read.csv(
  input_file,
  check.names = FALSE,
  row.names = 1
)

rmsd <- as.matrix(rmsd)
mode(rmsd) <- "numeric"

rownames(rmsd) <- as.character(rownames(rmsd))
colnames(rmsd) <- as.character(colnames(rmsd))


# ============================================================
# CLUSTERS
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
# HIERARCHICAL ORDER WITHIN EACH CLUSTER
# ============================================================

hc <- hclust(
  as.dist(rmsd),
  method = "average"
)

cluster_order <- hc$labels[hc$order]

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
# CALCULATE THE IMPORTANT RMSD RANGES AUTOMATICALLY
# ============================================================

dom_mat <- rmsd[
  dominant,
  dominant,
  drop = FALSE
]

dom_values <- dom_mat[
  upper.tri(dom_mat)
]

alt_mat <- rmsd[
  alternative,
  alternative,
  drop = FALSE
]

alt_values <- alt_mat[
  upper.tri(alt_mat)
]

between_values <- rmsd[
  dominant,
  alternative,
  drop = FALSE
]

dom_min <- min(dom_values)
dom_max <- max(dom_values)

alt_min <- min(alt_values)
alt_max <- max(alt_values)

between_min <- min(between_values)
between_max <- max(between_values)


# ============================================================
# CONVERT TO LONG FORMAT WITH NUMERIC POSITIONS
# ============================================================

plot_matrix <- rmsd[
  display_order,
  display_order
]

long <- as.data.frame(plot_matrix)

long$Seed_y <- rownames(long)

long <- long %>%
  pivot_longer(
    cols = -Seed_y,
    names_to = "Seed_x",
    values_to = "RMSD"
  )


# Numeric positions
x_lookup <- setNames(
  seq_along(display_order),
  display_order
)

# Reverse Y position so first seed is at top
y_lookup <- setNames(
  rev(seq_along(display_order)),
  display_order
)

long$x <- x_lookup[long$Seed_x]
long$y <- y_lookup[long$Seed_y]


# ============================================================
# AXIS LABELS
# ============================================================

x_axis_labels <- display_order
x_axis_labels[
  x_axis_labels == "9009"
] <- "9009 *"

y_display <- rev(display_order)

y_axis_labels <- y_display
y_axis_labels[
  y_axis_labels == "9009"
] <- "9009 *"


# ============================================================
# COLORS
# Same cluster colors as Panel D
# ============================================================

dominant_color <- "#E69F00"
alternative_color <- "#56B4E9"


# ============================================================
# PLOT
# ============================================================

p <- ggplot(
  long,
  aes(
    x = x,
    y = y,
    fill = RMSD
  )
) +

  # Heatmap cells
  geom_tile(
    color = "white",
    linewidth = 0.6,
    width = 1,
    height = 1
  ) +

  # RMSD values
  geom_text(
    aes(
      label = sprintf("%.2f", RMSD)
    ),
    size = 3.0,
    color = "black"
  ) +

  # RMSD scale
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
    name = "Pairwise pose RMSD (Å)"
  ) +

  # ========================================================
  # CLUSTER DIVIDERS
  # ========================================================

  geom_vline(
    xintercept = 8.5,
    linewidth = 1.2,
    color = "black"
  ) +

  geom_hline(
    yintercept = 2.5,
    linewidth = 1.2,
    color = "black"
  ) +

  # ========================================================
  # CLUSTER COLOR STRIPS — TOP
  # ========================================================

  annotate(
    "rect",
    xmin = 0.5,
    xmax = 8.5,
    ymin = 10.62,
    ymax = 10.84,
    fill = dominant_color
  ) +

  annotate(
    "rect",
    xmin = 8.5,
    xmax = 10.5,
    ymin = 10.62,
    ymax = 10.84,
    fill = alternative_color
  ) +

  # ========================================================
  # CLUSTER LABELS ABOVE HEATMAP
  # ========================================================

  annotate(
    "text",
    x = 4.5,
    y = 11.15,
    label = "Dominant pose cluster (8/10)",
    fontface = "bold",
    size = 4.2
  ) +

  annotate(
    "text",
    x = 9.5,
    y = 11.15,
    label = "Alternative\ncluster (2/10)",
    fontface = "bold",
    size = 3.8
  ) +

  # ========================================================
  # EXPLANATORY CALLOUTS
  # ========================================================

  annotate(
    "label",
    x = 11.75,
    y = 8.4,
    label = paste0(
      "Within dominant cluster\n",
      sprintf("%.2f", dom_min),
      "–",
      sprintf("%.2f", dom_max),
      " Å"
    ),
    hjust = 0,
    size = 4.0,
    fill = "white"
  ) +

  annotate(
    "label",
    x = 11.75,
    y = 6.0,
    label = paste0(
      "Between clusters\n",
      sprintf("%.2f", between_min),
      "–",
      sprintf("%.2f", between_max),
      " Å"
    ),
    hjust = 0,
    size = 4.0,
    fill = "white"
  ) +

  annotate(
    "label",
    x = 11.75,
    y = 3.6,
    label = paste0(
      "Within alternative cluster\n",
      sprintf("%.2f", alt_min),
      "–",
      sprintf("%.2f", alt_max),
      " Å"
    ),
    hjust = 0,
    size = 4.0,
    fill = "white"
  ) +

  # ========================================================
  # EXPLANATORY FOOTNOTE
  # ========================================================

  annotate(
    "text",
    x = 5.5,
    y = -0.35,
    label = paste0(
      "Each cell compares the top-ranked ISO pose from the row seed ",
      "with the top-ranked pose from the column seed.\n",
      "Diagonal = self-comparison (0 Å); * seed 9009 = representative dominant pose."
    ),
    size = 3.6,
    hjust = 0.5
  ) +

  # ========================================================
  # AXES
  # ========================================================

  scale_x_continuous(
    breaks = 1:10,
    labels = x_axis_labels,
    limits = c(
      0.5,
      15.2
    ),
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    breaks = 1:10,
    labels = y_axis_labels,
    limits = c(
      -0.8,
      11.5
    ),
    expand = c(0, 0)
  ) +

  labs(
    x = "Best pose from Vina seed",
    y = "Best pose from Vina seed"
  ) +

  coord_fixed(
    ratio = 1,
    clip = "off"
  ) +

  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5
    )
  ) +

  theme_minimal(
    base_size = 14
  ) +

  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black",
      size = 10.5
    ),

    axis.text.y = element_text(
      color = "black",
      size = 10.5
    ),

    axis.title = element_text(
      face = "bold",
      color = "black"
    ),

    legend.title = element_text(
      face = "bold"
    ),

    plot.margin = margin(
      25,
      25,
      45,
      15
    )
  )


# ============================================================
# SAVE AS NEW FILES
# Original Panel E remains untouched
# ============================================================

png_file <- file.path(
  output_dir,
  "AK1_panelE_pose_RMSD_heatmap_EXPLAINED.png"
)

ggsave(
  filename = png_file,
  plot = p,
  width = 11.0,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white"
)


pdf_file <- file.path(
  output_dir,
  "AK1_panelE_pose_RMSD_heatmap_EXPLAINED.pdf"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 11.0,
  height = 7.2,
  units = "in"
)


# ============================================================
# PRINT SUMMARY
# ============================================================

cat("\n")
cat("Dominant-cluster RMSD range:",
    sprintf("%.2f", dom_min),
    "-",
    sprintf("%.2f", dom_max),
    "A\n")

cat("Alternative-cluster RMSD range:",
    sprintf("%.2f", alt_min),
    "-",
    sprintf("%.2f", alt_max),
    "A\n")

cat("Between-cluster RMSD range:",
    sprintf("%.2f", between_min),
    "-",
    sprintf("%.2f", between_max),
    "A\n")

cat("\nSaved:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")