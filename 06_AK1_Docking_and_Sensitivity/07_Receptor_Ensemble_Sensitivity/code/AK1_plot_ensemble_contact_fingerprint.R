# ============================================================
# AK1 RECEPTOR-CONFORMATION CONTACT FINGERPRINT FIGURES
# ============================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}

library(ggplot2)


# ============================================================
# PATHS
# ============================================================

project_dir <- "C:/Users/akidi/Desktop/ISO_AK_Docking"

input_dir <- file.path(
  project_dir,
  "results",
  "AK1_ensemble_sensitivity",
  "mapped_contact_fingerprint_v1"
)

figure_dir <- file.path(
  input_dir,
  "figures"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


recurrence_file <- file.path(
  input_dir,
  "AK1_ensemble_contact_recurrence_mapped_to_1Z83.csv"
)

overlap_file <- file.path(
  input_dir,
  "AK1_ensemble_contact_overlap_vs_closed.csv"
)


# ============================================================
# READ DATA
# ============================================================

recurrence <- read.csv(
  recurrence_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

overlap <- read.csv(
  overlap_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# GROUP ORDER / DISPLAY LABELS
# ============================================================

group_order <- c(
  "1Z83 closed dominant",
  "7X7S least-open",
  "7X7S medoid",
  "7X7S most-open"
)

group_labels <- c(
  "1Z83 closed\ndominant",
  "7X7S\nleast-open",
  "7X7S\nmedoid",
  "7X7S\nmost-open"
)


# ============================================================
# HELPER: RESIDUE NUMBER
# ============================================================

get_resnum <- function(x) {

  as.numeric(
    gsub(
      "[^0-9]",
      "",
      x
    )
  )
}


# ============================================================
# HEATMAP FUNCTION
# ============================================================

make_contact_heatmap <- function(
  data,
  contact_type,
  cutoff_text,
  output_prefix
) {

  if (contact_type == "POLAR") {

    percent_col <- "Polar_contact_percent"
    count_col <- "Polar_contact_n"

  } else if (contact_type == "ANY") {

    percent_col <- "Any_contact_percent"
    count_col <- "Any_contact_n"

  } else {

    stop(
      "contact_type must be POLAR or ANY"
    )
  }


  # ----------------------------------------------------------
  # Keep residues recurrent >=80% in at least one state
  # ----------------------------------------------------------

  max_by_residue <- aggregate(
    data[[percent_col]],
    by = list(
      Reference_residue = data$Reference_residue
    ),
    FUN = max,
    na.rm = TRUE
  )

  names(max_by_residue)[2] <- "Max_percent"


  keep_residues <- max_by_residue$Reference_residue[
    max_by_residue$Max_percent >= 80
  ]


  plot_data <- data[
    data$Reference_residue %in% keep_residues,
    ,
    drop = FALSE
  ]


  # ----------------------------------------------------------
  # Order receptors
  # ----------------------------------------------------------

  plot_data$Group <- factor(
    plot_data$Group,
    levels = group_order,
    labels = group_labels
  )


  # ----------------------------------------------------------
  # Order residues numerically
  # ----------------------------------------------------------

  residue_table <- unique(
    data.frame(
      residue = plot_data$Reference_residue,
      number = get_resnum(
        plot_data$Reference_residue
      ),
      stringsAsFactors = FALSE
    )
  )


  residue_table <- residue_table[
    order(
      residue_table$number,
      residue_table$residue
    ),
    ,
    drop = FALSE
  ]


  residue_order <- residue_table$residue


  plot_data$Reference_residue <- factor(
    plot_data$Reference_residue,
    levels = rev(
      residue_order
    )
  )


  # ----------------------------------------------------------
  # Values
  # ----------------------------------------------------------

  plot_data$Percent <- plot_data[[percent_col]]

  plot_data$Count <- plot_data[[count_col]]


  plot_data$Label <- ifelse(
    plot_data$Count > 0,
    paste0(
      plot_data$Count,
      "/",
      plot_data$n_poses
    ),
    ""
  )


  # ----------------------------------------------------------
  # SAVE EXACT PLOTTED DATA
  # ----------------------------------------------------------

  write.csv(
    plot_data,
    file.path(
      figure_dir,
      paste0(
        output_prefix,
        "_plot_data.csv"
      )
    ),
    row.names = FALSE
  )


  # ----------------------------------------------------------
  # PLOT
  # ----------------------------------------------------------

  p <- ggplot(
    plot_data,
    aes(
      x = Group,
      y = Reference_residue,
      fill = Percent
    )
  ) +

    geom_tile(
      linewidth = 0.45,
      color = "white"
    ) +

    geom_text(
      aes(
        label = Label
      ),
      size = 3.5,
      fontface = "bold"
    ) +

    scale_fill_gradient(
      low = "white",
      high = "#4472C4",
      limits = c(
        0,
        100
      ),
      breaks = c(
        0,
        20,
        40,
        60,
        80,
        100
      ),
      name = "Recurrence\n(%)"
    ) +

    labs(
      x = NULL,
      y = "AK1 residue (1Z83 numbering)",
      title = ifelse(
        contact_type == "POLAR",
        "AK1 polar-contact recurrence across receptor conformations",
        "AK1 residue-contact recurrence across receptor conformations"
      ),
      subtitle = paste0(
        cutoff_text,
        "; cells show number of poses with contact / total poses"
      )
    ) +

    theme_classic(
      base_size = 12
    ) +

    theme(
      axis.text.x = element_text(
        size = 11,
        face = "bold",
        lineheight = 0.9
      ),

      axis.text.y = element_text(
        size = 10
      ),

      axis.title.y = element_text(
        size = 11
      ),

      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0
      ),

      plot.subtitle = element_text(
        size = 10,
        hjust = 0
      ),

      legend.title = element_text(
        size = 10,
        face = "bold"
      ),

      legend.text = element_text(
        size = 9
      ),

      panel.grid = element_blank(),

      plot.margin = margin(
        12,
        16,
        12,
        12
      )
    )


  # ----------------------------------------------------------
  # Dynamic figure height
  # ----------------------------------------------------------

  n_residues <- length(
    residue_order
  )

  fig_height <- max(
    5.5,
    2.5 + n_residues * 0.30
  )


  # ----------------------------------------------------------
  # PNG
  # ----------------------------------------------------------

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        output_prefix,
        ".png"
      )
    ),
    plot = p,
    width = 8.2,
    height = fig_height,
    units = "in",
    dpi = 600,
    bg = "white"
  )


  # ----------------------------------------------------------
  # PDF
  # ----------------------------------------------------------

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        output_prefix,
        ".pdf"
      )
    ),
    plot = p,
    width = 8.2,
    height = fig_height,
    units = "in",
    bg = "white"
  )


  return(
    invisible(
      p
    )
  )
}


# ============================================================
# POLAR HEATMAP
# ============================================================

make_contact_heatmap(
  data = recurrence,
  contact_type = "POLAR",
  cutoff_text = "Polar O/N/S proximity <= 3.6 A",
  output_prefix = "AK1_ensemble_POLAR_contact_heatmap"
)


# ============================================================
# ANY-CONTACT HEATMAP
# ============================================================

make_contact_heatmap(
  data = recurrence,
  contact_type = "ANY",
  cutoff_text = "Heavy-atom contact <= 4.0 A",
  output_prefix = "AK1_ensemble_ANY_contact_heatmap"
)


# ============================================================
# CLOSED-FINGERPRINT RETENTION PLOT
# ============================================================

retention <- overlap[
  ,
  c(
    "Comparison",
    "Contact_type",
    "Closed_reference_retained_percent"
  )
]


retention$Receptor <- retention$Comparison


retention$Receptor[
  grepl(
    "least-open",
    retention$Comparison
  )
] <- "Least-open"


retention$Receptor[
  grepl(
    "medoid",
    retention$Comparison
  )
] <- "Medoid"


retention$Receptor[
  grepl(
    "most-open",
    retention$Comparison
  )
] <- "Most-open"


retention$Receptor <- factor(
  retention$Receptor,
  levels = c(
    "Least-open",
    "Medoid",
    "Most-open"
  )
)


retention$Contact_type <- factor(
  retention$Contact_type,
  levels = c(
    "Any",
    "Polar"
  ),
  labels = c(
    "Any contact",
    "Polar proximity"
  )
)


write.csv(
  retention,
  file.path(
    figure_dir,
    "AK1_ensemble_closed_fingerprint_retention_plot_data.csv"
  ),
  row.names = FALSE
)


p_retention <- ggplot(
  retention,
  aes(
    x = Receptor,
    y = Closed_reference_retained_percent,
    fill = Contact_type
  )
) +

  geom_col(
    position = position_dodge(
      width = 0.72
    ),
    width = 0.64
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        Closed_reference_retained_percent
      )
    ),
    position = position_dodge(
      width = 0.72
    ),
    vjust = -0.35,
    size = 3.6
  ) +

  scale_y_continuous(
    limits = c(
      0,
      80
    ),
    breaks = seq(
      0,
      80,
      20
    ),
    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +

  labs(
    x = "7X7S receptor conformation",
    y = "Closed-state fingerprint retained (%)",
    fill = NULL,
    title = "Retention of the closed-state ISO contact fingerprint"
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    axis.text.x = element_text(
      size = 11,
      face = "bold"
    ),

    axis.text.y = element_text(
      size = 10
    ),

    axis.title = element_text(
      size = 11
    ),

    plot.title = element_text(
      size = 15,
      face = "bold"
    ),

    legend.position = "bottom",

    legend.text = element_text(
      size = 10
    ),

    plot.margin = margin(
      12,
      16,
      12,
      12
    )
  )


ggsave(
  filename = file.path(
    figure_dir,
    "AK1_ensemble_closed_fingerprint_retention.png"
  ),
  plot = p_retention,
  width = 7.3,
  height = 5.3,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = file.path(
    figure_dir,
    "AK1_ensemble_closed_fingerprint_retention.pdf"
  ),
  plot = p_retention,
  width = 7.3,
  height = 5.3,
  units = "in",
  bg = "white"
)


# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("AK1 ENSEMBLE CONTACT FIGURES COMPLETE\n")
cat("============================================================\n")
cat("\n")

cat(
  "Figure folder:\n",
  figure_dir,
  "\n\n"
)

cat(
  "1. AK1_ensemble_POLAR_contact_heatmap.png/pdf\n"
)

cat(
  "2. AK1_ensemble_ANY_contact_heatmap.png/pdf\n"
)

cat(
  "3. AK1_ensemble_closed_fingerprint_retention.png/pdf\n"
)

cat("\n")

cat(
  "Plot-data CSV files were also saved.\n"
)