# ============================================================
# AK1 Panel E (clean paper version)
# Pairwise RMSD heatmap of best ISO poses across 10 Vina seeds
# ============================================================

# ----------------------------
# Packages
# ----------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")

library(ggplot2)
library(grid)

# ----------------------------
# Input / output
# ----------------------------
INPUT_CSV  <- "..\\results\\AK1_ISO_pose_RMSD_matrix.csv"
OUTPUT_PNG <- "..\\results\\AK1_panelE_clean_v2.png"

# ----------------------------
# Read matrix
# ----------------------------
mat <- as.matrix(read.csv(INPUT_CSV, row.names = 1, check.names = FALSE))

# Force numeric
mode(mat) <- "numeric"

# ----------------------------
# Desired plotting order
# Dominant cluster first, alternative cluster last
# 9009 is marked as representative dominant pose
# ----------------------------
seed_order <- c("10101", "1001", "2002", "5005", "3003", "9009", "4004", "7007", "6006", "8008")
axis_labels <- c("10101", "1001", "2002", "5005", "3003", "9009*", "4004", "7007", "6006", "8008")

mat <- mat[seed_order, seed_order, drop = FALSE]

# ----------------------------
# Build plotting dataframe
# y is reversed so first seed appears at top
# ----------------------------
plot_df <- expand.grid(
  x = 1:ncol(mat),
  y_index = 1:nrow(mat)
)

plot_df$seed_x <- rep(colnames(mat), each = nrow(mat))
plot_df$seed_y <- rep(rownames(mat), times = ncol(mat))

# top-to-bottom layout
plot_df$y <- 11 - plot_df$y_index

# matrix values
plot_df$value <- as.vector(t(mat))

# label formatting
plot_df$label <- sprintf("%.2f", plot_df$value)

# ----------------------------
# Cluster summary values
# ----------------------------
dominant_idx    <- 1:8
alternative_idx <- 9:10

dominant_block <- mat[dominant_idx, dominant_idx, drop = FALSE]
alt_block      <- mat[alternative_idx, alternative_idx, drop = FALSE]
between_block  <- mat[dominant_idx, alternative_idx, drop = FALSE]

dominant_vals <- dominant_block[upper.tri(dominant_block)]
alt_vals      <- alt_block[upper.tri(alt_block)]
between_vals  <- as.vector(between_block)

dom_min <- min(dominant_vals)
dom_max <- max(dominant_vals)

alt_min <- min(alt_vals)
alt_max <- max(alt_vals)

bet_min <- min(between_vals)
bet_max <- max(between_vals)

# ----------------------------
# Colors
# ----------------------------
dominant_color    <- "#E69F00"   # orange-gold
alternative_color <- "#56B4E9"   # light blue

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(plot_df, aes(x = x, y = y, fill = value)) +

  # heatmap tiles
  geom_tile(color = "white", linewidth = 0.8) +

  # numbers inside tiles
  geom_text(aes(label = label), size = 4.2) +

  # divide the two clusters
  geom_vline(xintercept = 8.5, linewidth = 1.6, color = "black") +
  geom_hline(yintercept = 2.5, linewidth = 1.6, color = "black") +

  # top cluster bars
  annotate("rect",
           xmin = 0.5, xmax = 8.5,
           ymin = 10.58, ymax = 10.76,
           fill = dominant_color, color = NA) +

  annotate("rect",
           xmin = 8.5, xmax = 10.5,
           ymin = 10.58, ymax = 10.76,
           fill = alternative_color, color = NA) +

  # top cluster titles
  annotate("text",
           x = 4.5, y = 10.98,
           label = "Dominant pose cluster (8/10)",
           size = 5.2, fontface = "bold") +

  annotate("text",
           x = 9.5, y = 10.98,
           label = "Alternative\ncluster (2/10)",
           size = 4.2, fontface = "bold", lineheight = 0.92) +

  # explanatory annotation boxes on right
  annotate("label",
           x = 13.0, y = 8.2,
           label = paste0("Within dominant cluster\n",
                          sprintf("%.2f", dom_min), "\u2013", sprintf("%.2f", dom_max), " \u00C5"),
           size = 4.9, hjust = 0, vjust = 0.5,
           label.size = 0.5, fill = "white") +

  annotate("label",
           x = 13.0, y = 5.8,
           label = paste0("Between clusters\n",
                          sprintf("%.2f", bet_min), "\u2013", sprintf("%.2f", bet_max), " \u00C5"),
           size = 4.9, hjust = 0, vjust = 0.5,
           label.size = 0.5, fill = "white") +

  annotate("label",
           x = 13.0, y = 3.5,
           label = paste0("Within alternative cluster\n",
                          sprintf("%.2f", alt_min), "\u2013", sprintf("%.2f", alt_max), " \u00C5"),
           size = 4.9, hjust = 0, vjust = 0.5,
           label.size = 0.5, fill = "white") +

  # manual centered x-axis title under heatmap only
  annotate("text",
           x = 5.5, y = -0.85,
           label = "Best pose from Vina seed",
           size = 8.0, fontface = "bold") +

  # footnote
  annotate("text",
           x = 5.5, y = -0.38,
           label = "Each cell compares the top-ranked ISO pose from the row seed with the top-ranked pose from the column seed.\nDiagonal = self-comparison (0 Å); * seed 9009 = representative dominant pose.",
           size = 4.5, lineheight = 1.05) +

  # fill scale
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 2.0,
    limits = c(0, 7.4),
    name = "Pairwise pose RMSD (\u00C5)"
  ) +

  # axes
  scale_x_continuous(
    breaks = 1:10,
    labels = axis_labels,
    limits = c(0.5, 16.8),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = 10:1,
    labels = axis_labels,
    limits = c(-1.1, 11.15),
    expand = c(0, 0)
  ) +

  labs(
    x = NULL,
    y = "Best pose from Vina seed"
  ) +

  coord_fixed(ratio = 1) +

  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = c(0.89, 0.56),
    legend.justification = c(0.5, 0.5),
    plot.margin = margin(30, 150, 55, 20)
  ) +

  guides(fill = guide_colorbar(
    barheight = unit(3.5, "cm"),
    barwidth  = unit(0.7, "cm"),
    title.position = "top"
  ))

# ----------------------------
# Save
# ----------------------------
ggsave(
  filename = OUTPUT_PNG,
  plot = p,
  width = 14,
  height = 8.5,
  dpi = 300,
  bg = "white"
)

cat("Saved cleaned panel to:\n", normalizePath(OUTPUT_PNG), "\n")