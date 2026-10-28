library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)

# =========================
# INPUT MATRIX
# =========================
seed_order <- c("10101","1001","2002","5005","3003","9009*","4004","7007","6006","8008")

mat <- matrix(c(
  0.00,0.07,0.04,0.05,0.06,0.05,0.10,0.09,7.35,7.35,
  0.07,0.00,0.05,0.05,0.04,0.04,0.06,0.08,7.35,7.35,
  0.04,0.05,0.00,0.03,0.03,0.02,0.07,0.06,7.35,7.35,
  0.05,0.05,0.03,0.00,0.03,0.03,0.06,0.04,7.35,7.34,
  0.06,0.04,0.03,0.03,0.00,0.02,0.05,0.06,7.36,7.35,
  0.05,0.04,0.02,0.03,0.02,0.00,0.05,0.06,7.35,7.35,
  0.10,0.06,0.07,0.06,0.05,0.05,0.00,0.06,7.36,7.35,
  0.09,0.08,0.06,0.04,0.06,0.06,0.06,0.00,7.34,7.33,
  7.35,7.35,7.35,7.35,7.36,7.35,7.36,7.34,0.00,0.06,
  7.35,7.35,7.35,7.34,7.35,7.35,7.35,7.33,0.06,0.00
), nrow = 10, byrow = TRUE)

rownames(mat) <- seed_order
colnames(mat) <- seed_order

df <- as.data.frame(mat)
df$RowSeed <- factor(rownames(df), levels = rev(seed_order))

df_long <- df %>%
  pivot_longer(-RowSeed, names_to = "ColSeed", values_to = "RMSD") %>%
  mutate(
    ColSeed = factor(ColSeed, levels = seed_order),
    RowNum = as.numeric(RowSeed),
    ColNum = as.numeric(ColSeed)
  )

# Convert to plotting positions
df_long$RowPlot <- factor(df_long$RowSeed, levels = rev(seed_order))
df_long$ColPlot <- factor(df_long$ColSeed, levels = seed_order)

# =========================
# BASIC NUMBERS
# =========================
n_total <- 10
n_dom <- 8
n_alt <- 2

# Centers of top labels
dom_center <- (1 + n_dom) / 2        # 4.5
alt_center <- (n_dom + 1 + n_total) / 2   # 9.5

# Top annotation positions
bar_y <- 10.35
text_y <- 10.72

# =========================
# PLOT
# =========================
p <- ggplot(df_long, aes(x = ColPlot, y = RowPlot, fill = RMSD)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", RMSD)), size = 4.8) +

  scale_fill_gradient2(
    low = "#2E6EB5",
    mid = "#F4E6D7",
    high = "#C5162F",
    midpoint = 2,
    name = "Pairwise pose RMSD (Å)"
  ) +

  # Thick black separator lines between dominant and alternative clusters
  geom_vline(xintercept = 8.5, linewidth = 1.6, color = "black") +
  geom_hline(yintercept = 2.5, linewidth = 1.6, color = "black") +

  # Top colored bars
  annotate("segment", x = 0.5, xend = 8.5, y = bar_y, yend = bar_y,
           color = "#E69F00", linewidth = 6, lineend = "butt") +
  annotate("segment", x = 8.5, xend = 10.5, y = bar_y, yend = bar_y,
           color = "#56B4E9", linewidth = 6, lineend = "butt") +

  # Top cluster labels - MOVED UPWARD
  annotate("text", x = dom_center, y = text_y,
           label = "Dominant pose cluster (8/10)",
           size = 7.2, fontface = "bold") +
  annotate("text", x = alt_center, y = text_y,
           label = "Alternative\ncluster (2/10)",
           size = 6.1, fontface = "bold", lineheight = 0.95) +

  # Right-side explanatory boxes
  annotate("label", x = 12.2, y = 8.0,
           label = "Within dominant cluster\n0.02–0.10 Å",
           hjust = 0, size = 6.2, label.size = 0.6, fill = "white") +
  annotate("label", x = 12.2, y = 5.7,
           label = "Between clusters\n7.33–7.36 Å",
           hjust = 0, size = 6.2, label.size = 0.6, fill = "white") +
  annotate("label", x = 12.2, y = 3.5,
           label = "Within alternative cluster\n0.06–0.06 Å",
           hjust = 0, size = 6.2, label.size = 0.6, fill = "white") +

  labs(
    x = "Best pose from Vina seed",
    y = "Best pose from Vina seed"
  ) +

  coord_cartesian(
    xlim = c(0.5, 15.5),
    ylim = c(0.5, 11.0),
    clip = "off"
  ) +

  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(
      size = 20, face = "bold",
      hjust = 0.5,                # center the bottom title
      margin = margin(t = 18)     # give it a bit of space from axis text
    ),
    axis.title.y = element_text(
      size = 20, face = "bold",
      margin = margin(r = 14)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 45, hjust = 1, vjust = 1
    ),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = c(0.90, 0.53),
    plot.margin = margin(35, 120, 25, 25)  # extra top and right space
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(100, "pt"),
    barwidth = unit(18, "pt")
  ))

# Save
ggsave(
  filename = "../results/AK1_panelE_pose_RMSD_heatmap_clean.png",
  plot = p,
  width = 14,
  height = 9,
  dpi = 300,
  bg = "white"
)