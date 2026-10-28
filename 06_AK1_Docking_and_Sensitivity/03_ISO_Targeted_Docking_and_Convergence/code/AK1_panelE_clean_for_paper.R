# ============================================================
# FIGURE 8 — PANEL E (CLEAN PAPER VERSION)
# Pairwise RMSD heatmap of best AK1–ISO poses across seeds
# ============================================================

packages <- c("ggplot2", "dplyr", "tidyr", "scales")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

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

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# READ MATRIX
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

dominant <- c("1001", "2002", "3003", "4004", "5005", "7007", "9009", "10101")
alternative <- c("6006", "8008")

# ============================================================
# ORDER
# ============================================================

hc <- hclust(as.dist(rmsd), method = "average")
cluster_order <- hc$labels[hc$order]

dominant_order <- cluster_order[cluster_order %in% dominant]
alternative_order <- cluster_order[cluster_order %in% alternative]

display_order <- c(dominant_order, alternative_order)

plot_matrix <- rmsd[display_order, display_order]

# ============================================================
# SUMMARY VALUES
# ============================================================

dom_mat <- rmsd[dominant, dominant, drop = FALSE]
dom_vals <- dom_mat[upper.tri(dom_mat)]

alt_mat <- rmsd[alternative, alternative, drop = FALSE]
alt_vals <- alt_mat[upper.tri(alt_mat)]

between_vals <- rmsd[dominant, alternative, drop = FALSE]

dom_min <- min(dom_vals)
dom_max <- max(dom_vals)

alt_min <- min(alt_vals)
alt_max <- max(alt_vals)

between_min <- min(between_vals)
between_max <- max(between_vals)

# ============================================================
# LONG FORMAT
# ============================================================

long <- as.data.frame(plot_matrix)
long$Seed_y <- rownames(long)

long <- long %>%
  pivot_longer(
    cols = -Seed_y,
    names_to = "Seed_x",
    values_to = "RMSD"
  )

x_lookup <- setNames(seq_along(display_order), display_order)
y_lookup <- setNames(rev(seq_along(display_order)), display_order)

long$x <- x_lookup[long$Seed_x]
long$y <- y_lookup[long$Seed_y]

# Axis labels
x_axis_labels <- display_order
x_axis_labels[x_axis_labels == "9009"] <- "9009*"

y_display <- rev(display_order)
y_axis_labels <- y_display
y_axis_labels[y_axis_labels == "9009"] <- "9009*"

# Colors
dominant_color <- "#E69F00"
alternative_color <- "#56B4E9"

# ============================================================
# PLOT
# ============================================================

p <- ggplot(long, aes(x = x, y = y, fill = RMSD)) +

  geom_tile(
    color = "white",
    linewidth = 0.7,
    width = 1,
    height = 1
  ) +

  geom_text(
    aes(label = sprintf("%.2f", RMSD)),
    size = 3.2,
    color = "black"
  ) +

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
    values = rescale(c(0, 0.5, 1, 2, 3, 5, max(long$RMSD))),
    limits = c(0, max(long$RMSD)),
    name = "Pairwise RMSD (Å)"
  ) +

  # Cluster dividing lines
  geom_vline(xintercept = 8.5, linewidth = 1.4, color = "black") +
  geom_hline(yintercept = 2.5, linewidth = 1.4, color = "black") +

  # Top cluster bars
  annotate("rect", xmin = 0.5, xmax = 8.5, ymin = 10.15, ymax = 10.35,
           fill = dominant_color, color = NA) +
  annotate("rect", xmin = 8.5, xmax = 10.5, ymin = 10.15, ymax = 10.35,
           fill = alternative_color, color = NA) +

  # Cluster labels
  annotate("text", x = 4.5, y = 10.62,
           label = "Dominant pose cluster (8/10)",
           size = 4.4, fontface = "bold") +

  annotate("text", x = 9.5, y = 10.62,
           label = "Alternative\ncluster (2/10)",
           size = 3.8, fontface = "bold", lineheight = 0.95) +

  # Right-side explanatory boxes
  annotate(
    "label",
    x = 12.0,
    y = 8.5,
    label = paste0(
      "Within dominant cluster\n",
      sprintf("%.2f", dom_min), "–", sprintf("%.2f", dom_max), " Å"
    ),
    hjust = 0,
    size = 4.0,
    label.size = 0.5,
    fill = "white"
  ) +

  annotate(
    "label",
    x = 12.0,
    y = 6.0,
    label = paste0(
      "Between clusters\n",
      sprintf("%.2f", between_min), "–", sprintf("%.2f", between_max), " Å"
    ),
    hjust = 0,
    size = 4.0,
    label.size = 0.5,
    fill = "white"
  ) +

  annotate(
    "label",
    x = 12.0,
    y = 3.5,
    label = paste0(
      "Within alternative cluster\n",
      sprintf("%.2f", alt_min), "–", sprintf("%.2f", alt_max), " Å"
    ),
    hjust = 0,
    size = 4.0,
    label.size = 0.5,
    fill = "white"
  ) +

  scale_x_continuous(
    breaks = 1:10,
    labels = x_axis_labels,
    limits = c(0.5, 15.4),
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    breaks = 1:10,
    labels = y_axis_labels,
    limits = c(0.5, 10.9),
    expand = c(0, 0)
  ) +

  labs(
    x = "Best pose from Vina seed",
    y = "Best pose from Vina seed"
  ) +

  coord_fixed(clip = "off") +

  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(70, "pt"),
      barwidth  = unit(12, "pt")
    )
  ) +

  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 11
    ),

    axis.text.y = element_text(
      color = "black",
      size = 11
    ),

    axis.title = element_text(
      face = "bold",
      color = "black",
      size = 16
    ),

    legend.title = element_text(
      face = "bold",
      size = 14
    ),

    legend.text = element_text(size = 11),

    plot.margin = margin(25, 120, 25, 20)
  )

# ============================================================
# SAVE
# ============================================================

png_file <- file.path(output_dir, "AK1_panelE_pose_RMSD_heatmap_CLEAN.png")
pdf_file <- file.path(output_dir, "AK1_panelE_pose_RMSD_heatmap_CLEAN.pdf")

ggsave(
  filename = png_file,
  plot = p,
  width = 12.5,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 12.5,
  height = 7.2,
  units = "in",
  bg = "white"
)

cat("\nSaved:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")