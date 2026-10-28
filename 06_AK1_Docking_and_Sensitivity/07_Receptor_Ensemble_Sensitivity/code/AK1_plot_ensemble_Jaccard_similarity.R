suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

# ============================================================
# PATHS
# ============================================================

input_dir <- file.path(
  "..",
  "results",
  "AK1_ensemble_sensitivity",
  "mapped_contact_fingerprint_v1"
)

input_file <- file.path(
  input_dir,
  "AK1_ensemble_contact_overlap_vs_closed.csv"
)

output_dir <- file.path(
  input_dir,
  "figures"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# READ DATA
# ============================================================

df <- read_csv(
  input_file,
  show_col_types = FALSE
)


# ============================================================
# CLEAN RECEPTOR LABELS
# ============================================================

df <- df %>%
  mutate(

    Receptor = case_when(

      grepl(
        "least-open",
        Comparison,
        ignore.case = TRUE
      ) ~ "Least-open",

      grepl(
        "medoid",
        Comparison,
        ignore.case = TRUE
      ) ~ "Medoid",

      grepl(
        "most-open",
        Comparison,
        ignore.case = TRUE
      ) ~ "Most-open",

      TRUE ~ Comparison
    ),

    Receptor = factor(
      Receptor,
      levels = c(
        "Least-open",
        "Medoid",
        "Most-open"
      )
    ),

    Contact_type = factor(
      Contact_type,
      levels = c(
        "Any",
        "Polar"
      ),
      labels = c(
        "Any contact",
        "Polar proximity"
      )
    )
  )


# ============================================================
# SAVE EXACT PLOT DATA
# ============================================================

write_csv(
  df,
  file.path(
    output_dir,
    "AK1_ensemble_Jaccard_similarity_plot_data.csv"
  )
)


# ============================================================
# PLOT
# ============================================================

p <- ggplot(
  df,
  aes(
    x = Receptor,
    y = Jaccard,
    fill = Contact_type
  )
) +

  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.66
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        Jaccard
      )
    ),
    position = position_dodge(
      width = 0.75
    ),
    vjust = -0.45,
    size = 4.2
  ) +

  scale_y_continuous(
    limits = c(
      0,
      0.45
    ),
    breaks = seq(
      0,
      0.4,
      0.1
    ),
    expand = expansion(
      mult = c(
        0,
        0.04
      )
    )
  ) +

  labs(
    title =
      "Similarity of ISO contact fingerprints to the closed AK1 state",

    subtitle =
      "Jaccard similarity accounts for both retained and newly acquired contacts",

    x =
      "7X7S receptor conformation",

    y =
      "Jaccard similarity",

    fill =
      NULL
  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(

    axis.text.x = element_text(
      size = 12,
      face = "bold"
    ),

    axis.text.y = element_text(
      size = 11
    ),

    axis.title = element_text(
      size = 12
    ),

    plot.title = element_text(
      size = 16,
      face = "bold"
    ),

    plot.subtitle = element_text(
      size = 10
    ),

    legend.position =
      "bottom",

    legend.text = element_text(
      size = 11
    ),

    plot.margin = margin(
      12,
      16,
      12,
      12
    )
  )


# ============================================================
# SAVE
# ============================================================

ggsave(
  filename = file.path(
    output_dir,
    "AK1_ensemble_Jaccard_similarity.png"
  ),
  plot = p,
  width = 7.5,
  height = 5.6,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = file.path(
    output_dir,
    "AK1_ensemble_Jaccard_similarity.pdf"
  ),
  plot = p,
  width = 7.5,
  height = 5.6,
  units = "in",
  bg = "white"
)


cat("\n")
cat("============================================================\n")
cat("AK1 JACCARD SIMILARITY FIGURE COMPLETE\n")
cat("============================================================\n\n")

print(
  df %>%
    select(
      Receptor,
      Contact_type,
      Jaccard,
      Shared_n,
      Closed_reference_n,
      Open_state_n
    )
)

cat("\nSaved:\n")

cat(
  file.path(
    output_dir,
    "AK1_ensemble_Jaccard_similarity.png"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "AK1_ensemble_Jaccard_similarity.pdf"
  ),
  "\n"
)