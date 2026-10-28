# ============================================================
# FIGURE 8 — PANEL F
# Stability of recurrent AK1–ISO polar contacts
# across the dominant 8-pose cluster
# ============================================================

packages <- c("ggplot2", "dplyr")

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


# ============================================================
# PATHS
# ============================================================

project <- "C:/Users/akidi/Desktop/ISO_AK_Docking"

input_file <- file.path(
  project,
  "results",
  "AK1_dominant8_contact_frequency",
  "AK1_dominant8_recurrent_polar_distances_long.csv"
)

output_dir <- file.path(
  project,
  "results",
  "AK1_panelF"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# READ REAL DATA
# ============================================================

dat <- read.csv(
  input_file,
  check.names = FALSE
)

names(dat)[3] <- "Distance"


# ============================================================
# ORDER RESIDUES BY MEAN DISTANCE
# ============================================================

summary_dat <- dat %>%
  group_by(Residue) %>%
  summarise(
    Mean = mean(Distance),
    SD = sd(Distance),
    Min = min(Distance),
    Max = max(Distance),
    n = n()
  ) %>%
  arrange(Mean)


dat$Residue <- factor(
  dat$Residue,
  levels = summary_dat$Residue
)

summary_dat$Residue <- factor(
  summary_dat$Residue,
  levels = summary_dat$Residue
)


# ============================================================
# PRINT SUMMARY
# ============================================================

print(summary_dat)


# ============================================================
# PLOT
# ============================================================

p <- ggplot(
  dat,
  aes(
    x = Residue,
    y = Distance
  )
) +

  # 3.6 A screening cutoff
  geom_hline(
    yintercept = 3.6,
    linetype = "dashed",
    linewidth = 0.8,
    color = "grey35"
  ) +

  # individual dominant-pose measurements
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 3.2,
    alpha = 0.65,
    color = "grey35"
  ) +

  # mean +/- SD
  geom_errorbar(
    data = summary_dat,
    mapping = aes(
      x = Residue,
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    width = 0.18,
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +

  geom_point(
    data = summary_dat,
    aes(
      x = Residue,
      y = Mean
    ),
    shape = 23,
    size = 5,
    fill = "#E69F00",
    color = "black",
    inherit.aes = FALSE
  ) +

  # n=8/8 annotation
  geom_text(
    data = summary_dat,
    aes(
      x = Residue,
      y = 3.72,
      label = "8/8"
    ),
    size = 4.2,
    fontface = "bold",
    inherit.aes = FALSE
  ) +

  annotate(
    "text",
    x = Inf,
    y = 3.60,
    label = "3.6 Å polar-contact screening cutoff",
    hjust = 1.05,
    vjust = -0.6,
    size = 4.0,
    color = "grey25"
  ) +

  labs(
    x = "AK1 residue",
    y = "Minimum polar heavy-atom distance (Å)"
  ) +

  coord_cartesian(
    ylim = c(
      min(dat$Distance) - 0.10,
      3.82
    )
  ) +

  theme_classic(
    base_size = 15
  ) +

  theme(
    axis.text = element_text(
      color = "black"
    ),

    axis.text.x = element_text(
      face = "bold"
    ),

    axis.title = element_text(
      face = "bold",
      color = "black"
    ),

    plot.margin = margin(
      15,
      20,
      15,
      15
    )
  )


# ============================================================
# SAVE
# ============================================================

png_file <- file.path(
  output_dir,
  "AK1_panelF_recurrent_polar_contact_stability.png"
)

pdf_file <- file.path(
  output_dir,
  "AK1_panelF_recurrent_polar_contact_stability.pdf"
)

ggsave(
  filename = png_file,
  plot = p,
  width = 7.5,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 7.5,
  height = 5.5,
  units = "in"
)


cat("\nSaved PNG:\n")
cat(png_file, "\n")

cat("\nSaved PDF:\n")
cat(pdf_file, "\n")