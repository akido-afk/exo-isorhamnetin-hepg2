# ============================================================
# FIGURE 8 — PANEL G
# Native AP5 vs isorhamnetin polar-contact fingerprint
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
  "AK1_AP5_ISO_fingerprint",
  "AK1_AP5_vs_ISO_residue_fingerprint.csv"
)

output_dir <- file.path(
  project,
  "results",
  "AK1_panelG"
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


# ============================================================
# KEEP ONLY POLAR-CONTACT FINGERPRINT
#
# AP5: crystallographic polar contact
# ISO: consensus polar contact >= 6/8 dominant poses
# ============================================================

polar <- dat %>%
  mutate(
    AP5_present =
      `AP5 polar contact` == 1,

    ISO_present =
      `ISO polar contact n/8` >= 6
  ) %>%

  filter(
    AP5_present |
    ISO_present
  )


# ============================================================
# CLASSIFY RESIDUES
# ============================================================

polar <- polar %>%
  mutate(
    Class = case_when(
      AP5_present & ISO_present ~ "Shared",
      !AP5_present & ISO_present ~ "ISO only",
      AP5_present & !ISO_present ~ "AP5 only",
      TRUE ~ "Other"
    )
  )


# ============================================================
# EXTRACT RESIDUE NUMBER FOR ORDERING
# ============================================================

polar <- polar %>%
  mutate(
    Residue_number =
      as.numeric(
        gsub(
          "[^0-9]",
          "",
          Residue
        )
      )
  )


# ============================================================
# ORDER:
# Shared first
# ISO-only second
# AP5-only third
# ============================================================

class_order <- c(
  "Shared",
  "ISO only",
  "AP5 only"
)

polar$Class <- factor(
  polar$Class,
  levels = class_order
)

polar <- polar %>%
  arrange(
    Class,
    Residue_number
  )


# ============================================================
# SUMMARY
# ============================================================

shared_n <- sum(
  polar$Class == "Shared"
)

iso_only_n <- sum(
  polar$Class == "ISO only"
)

ap5_only_n <- sum(
  polar$Class == "AP5 only"
)

ap5_set <- polar$Residue[
  polar$AP5_present
]

iso_set <- polar$Residue[
  polar$ISO_present
]

jaccard <- length(
  intersect(
    ap5_set,
    iso_set
  )
) /
length(
  union(
    ap5_set,
    iso_set
  )
)


cat("\n")
cat("============================================\n")
cat("AK1 AP5 vs ISO POLAR FINGERPRINT\n")
cat("============================================\n")

cat(
  "Shared residues:",
  paste(
    polar$Residue[
      polar$Class == "Shared"
    ],
    collapse = ", "
  ),
  "\n"
)

cat(
  "ISO-only residues:",
  paste(
    polar$Residue[
      polar$Class == "ISO only"
    ],
    collapse = ", "
  ),
  "\n"
)

cat(
  "AP5-only residues:",
  paste(
    polar$Residue[
      polar$Class == "AP5 only"
    ],
    collapse = ", "
  ),
  "\n"
)

cat(
  "Jaccard overlap =",
  sprintf("%.2f", jaccard),
  "\n"
)


# ============================================================
# LONG FORMAT
# ============================================================

plot_dat <- polar %>%
  select(
    Residue,
    Class,
    AP5_present,
    ISO_present
  ) %>%

  pivot_longer(
    cols = c(
      AP5_present,
      ISO_present
    ),
    names_to = "Ligand",
    values_to = "Present"
  ) %>%

  mutate(
    Ligand = recode(
      Ligand,
      AP5_present =
        "Native AP5",
      ISO_present =
        "ISO dominant\ncluster"
    )
  )


# ============================================================
# TILE COLOR
# ============================================================

plot_dat <- plot_dat %>%
  mutate(
    Tile = case_when(
      !Present ~ "Absent",

      Present &
        Ligand == "Native AP5"
        ~ "AP5",

      Present &
        Ligand == "ISO dominant\ncluster"
        ~ "ISO"
    )
  )


# ============================================================
# KEEP RESIDUE ORDER WITHIN FACETS
# ============================================================

plot_dat$Residue <- factor(
  plot_dat$Residue,
  levels = rev(
    polar$Residue
  )
)


# ============================================================
# PLOT
# ============================================================

p <- ggplot(
  plot_dat,
  aes(
    x = Ligand,
    y = Residue,
    fill = Tile
  )
) +

  geom_tile(
    color = "grey45",
    linewidth = 0.6,
    width = 0.88,
    height = 0.88
  ) +

  # mark present contacts
  geom_text(
    data = subset(
      plot_dat,
      Present
    ),
    aes(
      label = ifelse(
        Ligand == "ISO dominant\ncluster",
        "8/8",
        "✓"
      )
    ),
    color = "black",
    fontface = "bold",
    size = 4.2
  ) +

  scale_fill_manual(
    values = c(
      "AP5" = "#00BFC4",
      "ISO" = "#E69F00",
      "Absent" = "white"
    ),
    guide = "none"
  ) +

  facet_grid(
    Class ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +

  labs(
    x = NULL,
    y = NULL,
    title = NULL,
    subtitle = paste0(
      "Polar-contact fingerprint overlap: ",
      shared_n,
      " shared, ",
      iso_only_n,
      " ISO-only, ",
      ap5_only_n,
      " AP5-only residues; Jaccard = ",
      sprintf("%.2f", jaccard)
    )
  ) +

  theme_minimal(
    base_size = 14
  ) +

  theme(
    panel.grid = element_blank(),

    panel.spacing.y =
      unit(
        0.6,
        "lines"
      ),

    axis.text.x =
      element_text(
        face = "bold",
        size = 13,
        color = "black",
        margin = margin(
          t = 8
        )
      ),

    axis.text.y =
      element_text(
        face = "bold",
        size = 11,
        color = "black"
      ),

    strip.placement = "outside",

    strip.text.y.left =
      element_text(
        face = "bold",
        size = 12,
        angle = 0,
        color = "black"
      ),

    strip.background =
      element_rect(
        fill = "grey92",
        color = "grey55"
      ),

    plot.subtitle =
      element_text(
        size = 12.5,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 12
        )
      ),

    plot.margin =
      margin(
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
  "AK1_panelG_AP5_ISO_polar_fingerprint.png"
)

pdf_file <- file.path(
  output_dir,
  "AK1_panelG_AP5_ISO_polar_fingerprint.pdf"
)


ggsave(
  filename = png_file,
  plot = p,
  width = 6.5,
  height = 9.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 6.5,
  height = 9.0,
  units = "in"
)


cat("\n")
cat("FILES SAVED\n")
cat("===========\n")

cat(png_file, "\n")
cat(pdf_file, "\n")