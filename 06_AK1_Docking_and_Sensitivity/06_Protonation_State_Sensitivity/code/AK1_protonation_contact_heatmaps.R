# ============================================================
# AK1 PROTONATION-STATE CONTACT FINGERPRINT HEATMAPS
#
# E1 = Polar-contact recurrence
# E2 = Any-contact recurrence
#
# Four groups:
#   Neutral dominant
#   Neutral alternative
#   7-O- monoanion
#   7-O-/4'-O- dianion
#
# Figures saved separately.
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "openxlsx"
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
library(openxlsx)


# ============================================================
# PATHS
# ============================================================

input_dir <- paste0(
  "../results/",
  "AK1_ligand_state_sensitivity/",
  "contact_fingerprint_v1"
)

output_dir <- paste0(
  "../results/",
  "AK1_ligand_state_sensitivity/",
  "contact_fingerprint_figures"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


input_file <- file.path(
  input_dir,
  "AK1_protonation_state_contact_recurrence.csv"
)


if (!file.exists(input_file)) {

  stop(
    paste(
      "Input file not found:",
      input_file
    )
  )
}


# ============================================================
# READ DATA
# ============================================================

dat <- read.csv(
  input_file,
  check.names = FALSE
)


# ============================================================
# EXPECTED GROUP ORDER
# ============================================================

group_levels <- c(
  "Neutral dominant",
  "Neutral alternative",
  "7-O- monoanion",
  "7-O-/4'-O- dianion"
)


dat$Group <- factor(
  dat$Group,
  levels = group_levels
)


# ============================================================
# VALIDATION
# ============================================================

required_cols <- c(
  "Group",
  "Residue",
  "n_poses",
  "Any_contact_n",
  "Any_contact_percent",
  "Polar_contact_n",
  "Polar_contact_percent"
)


if (!all(
  required_cols %in% colnames(dat)
)) {

  stop(
    paste(
      "Missing required columns. Found:",
      paste(
        colnames(dat),
        collapse = ", "
      )
    )
  )
}


# ============================================================
# POLAR CONTACT DATA
#
# Keep residues with at least one polar contact
# in any of the four groups.
# ============================================================

polar_residues <- dat %>%

  group_by(
    Residue
  ) %>%

  summarise(
    MaxPolar = max(
      Polar_contact_n,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%

  filter(
    MaxPolar > 0
  ) %>%

  pull(
    Residue
  )


# ============================================================
# MANUAL POLAR RESIDUE ORDER
#
# Ordered to make the family patterns easy to see.
# Any additional residue discovered by the input data is added
# automatically at the bottom.
# ============================================================

preferred_polar_order <- c(
  "LYS21",
  "ARG97",
  "ARG149",
  "ARG132",
  "ARG138",
  "MET61",
  "ASP93",
  "SER38",
  "GLY16",
  "GLY22",
  "PRO17"
)


extra_polar <- setdiff(
  polar_residues,
  preferred_polar_order
)


polar_order <- c(
  preferred_polar_order[
    preferred_polar_order %in% polar_residues
  ],
  sort(
    extra_polar
  )
)


# ============================================================
# PREPARE POLAR HEATMAP
# ============================================================

polar_plot_data <- dat %>%

  filter(
    Residue %in% polar_residues
  ) %>%

  mutate(

    Cell_label = ifelse(
      Polar_contact_n > 0,
      paste0(
        Polar_contact_n,
        "/",
        n_poses
      ),
      ""
    ),

    Residue = factor(
      Residue,
      levels = rev(
        polar_order
      )
    )
  )


# ============================================================
# POLAR HEATMAP
# ============================================================

p_polar <- ggplot(
  polar_plot_data,
  aes(
    x = Group,
    y = Residue,
    fill = Polar_contact_percent
  )
) +

  geom_tile(
    color = "grey40",
    linewidth = 0.7,
    width = 0.94,
    height = 0.90
  ) +

  geom_text(
    aes(
      label = Cell_label
    ),
    size = 4.6,
    fontface = "bold",
    color = "black"
  ) +

  scale_fill_gradient(
    low = "white",
    high = "#E69F00",
    limits = c(
      0,
      100
    ),
    breaks = c(
      0,
      25,
      50,
      75,
      100
    ),
    labels = function(x) {
      paste0(
        x,
        "%"
      )
    },
    name = "Recurrence"
  ) +

  labs(

    title =
      "AK1 polar-contact fingerprint across ISO protonation states",

    subtitle =
      paste0(
        "Cell labels show contact recurrence (n/N); ",
        "color intensity represents recurrence percentage"
      ),

    x = NULL,

    y = NULL,

    caption =
      paste0(
        "Polar contact = O/N/S heavy-atom proximity \u2264 3.6 \u00C5. ",
        "These proximities are not automatically interpreted as hydrogen bonds."
      )
  ) +

  theme_minimal(
    base_size = 14
  ) +

  theme(

    panel.grid =
      element_blank(),

    axis.text.x =
      element_text(
        size = 12,
        face = "bold",
        color = "black",
        angle = 18,
        hjust = 1,
        vjust = 1,
        margin = margin(
          t = 8
        )
      ),

    axis.text.y =
      element_text(
        size = 12.5,
        face = "bold",
        color = "black"
      ),

    legend.title =
      element_text(
        size = 11.5,
        face = "bold"
      ),

    legend.text =
      element_text(
        size = 10.5
      ),

    plot.title =
      element_text(
        size = 18,
        face = "bold",
        hjust = 0
      ),

    plot.subtitle =
      element_text(
        size = 11,
        color = "grey25",
        hjust = 0,
        margin = margin(
          b = 14
        )
      ),

    plot.caption =
      element_text(
        size = 9.5,
        color = "grey30",
        hjust = 0,
        margin = margin(
          t = 12
        )
      ),

    plot.margin =
      margin(
        25,
        30,
        25,
        30
      ),

    plot.title.position =
      "plot",

    plot.caption.position =
      "plot"
  )


# ============================================================
# SAVE POLAR HEATMAP
# ============================================================

polar_png <- file.path(
  output_dir,
  "AK1_protonation_state_POLAR_contact_heatmap.png"
)

polar_pdf <- file.path(
  output_dir,
  "AK1_protonation_state_POLAR_contact_heatmap.pdf"
)


ggsave(
  filename = polar_png,
  plot = p_polar,
  width = 11.5,
  height = 8.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = polar_pdf,
  plot = p_polar,
  width = 11.5,
  height = 8.5,
  units = "in",
  bg = "white"
)


# ============================================================
# ANY-CONTACT DATA
#
# Keep residues with at least one any-contact occurrence
# in any group.
# ============================================================

any_residues <- dat %>%

  group_by(
    Residue
  ) %>%

  summarise(
    MaxAny = max(
      Any_contact_n,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%

  filter(
    MaxAny > 0
  ) %>%

  pull(
    Residue
  )


# ============================================================
# ORDER ANY-CONTACT RESIDUES
#
# Prioritize biologically recurring / family-informative
# residues, then append any additional residues automatically.
# ============================================================

preferred_any_order <- c(
  "LYS21",
  "ARG97",
  "ARG149",
  "ARG132",
  "ARG138",
  "MET61",
  "ASP93",
  "ASP140",
  "GLY15",
  "GLY16",
  "GLY18",
  "GLY22",
  "GLY64",
  "GLY94",
  "PRO17",
  "SER38",
  "TYR153"
)


extra_any <- setdiff(
  any_residues,
  preferred_any_order
)


any_order <- c(
  preferred_any_order[
    preferred_any_order %in% any_residues
  ],
  sort(
    extra_any
  )
)


# ============================================================
# PREPARE ANY-CONTACT HEATMAP
# ============================================================

any_plot_data <- dat %>%

  filter(
    Residue %in% any_residues
  ) %>%

  mutate(

    Cell_label = ifelse(
      Any_contact_n > 0,
      paste0(
        Any_contact_n,
        "/",
        n_poses
      ),
      ""
    ),

    Residue = factor(
      Residue,
      levels = rev(
        any_order
      )
    )
  )


# ============================================================
# ANY-CONTACT HEATMAP
# ============================================================

p_any <- ggplot(
  any_plot_data,
  aes(
    x = Group,
    y = Residue,
    fill = Any_contact_percent
  )
) +

  geom_tile(
    color = "grey40",
    linewidth = 0.7,
    width = 0.94,
    height = 0.90
  ) +

  geom_text(
    aes(
      label = Cell_label
    ),
    size = 4.3,
    fontface = "bold",
    color = "black"
  ) +

  scale_fill_gradient(
    low = "white",
    high = "#56B4E9",
    limits = c(
      0,
      100
    ),
    breaks = c(
      0,
      25,
      50,
      75,
      100
    ),
    labels = function(x) {
      paste0(
        x,
        "%"
      )
    },
    name = "Recurrence"
  ) +

  labs(

    title =
      "AK1 residue-contact fingerprint across ISO protonation states",

    subtitle =
      paste0(
        "Cell labels show contact recurrence (n/N); ",
        "color intensity represents recurrence percentage"
      ),

    x = NULL,

    y = NULL,

    caption =
      "Any contact = minimum ligand\u2013residue heavy-atom distance \u2264 4.0 \u00C5."
  ) +

  theme_minimal(
    base_size = 14
  ) +

  theme(

    panel.grid =
      element_blank(),

    axis.text.x =
      element_text(
        size = 12,
        face = "bold",
        color = "black",
        angle = 18,
        hjust = 1,
        vjust = 1,
        margin = margin(
          t = 8
        )
      ),

    axis.text.y =
      element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),

    legend.title =
      element_text(
        size = 11.5,
        face = "bold"
      ),

    legend.text =
      element_text(
        size = 10.5
      ),

    plot.title =
      element_text(
        size = 18,
        face = "bold",
        hjust = 0
      ),

    plot.subtitle =
      element_text(
        size = 11,
        color = "grey25",
        hjust = 0,
        margin = margin(
          b = 14
        )
      ),

    plot.caption =
      element_text(
        size = 9.5,
        color = "grey30",
        hjust = 0,
        margin = margin(
          t = 12
        )
      ),

    plot.margin =
      margin(
        25,
        30,
        25,
        30
      ),

    plot.title.position =
      "plot",

    plot.caption.position =
      "plot"
  )


# ============================================================
# SAVE ANY-CONTACT HEATMAP
# ============================================================

any_png <- file.path(
  output_dir,
  "AK1_protonation_state_ANY_contact_heatmap.png"
)

any_pdf <- file.path(
  output_dir,
  "AK1_protonation_state_ANY_contact_heatmap.pdf"
)


ggsave(
  filename = any_png,
  plot = p_any,
  width = 11.5,
  height = 10.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = any_pdf,
  plot = p_any,
  width = 11.5,
  height = 10.5,
  units = "in",
  bg = "white"
)


# ============================================================
# QUANTITATIVE OVERLAP SUMMARY
# ============================================================

get_recurrent_set <- function(
  group_name,
  type = "polar",
  threshold = 80
) {

  group_data <- dat %>%
    filter(
      Group == group_name
    )

  if (type == "polar") {

    return(
      group_data %>%
        filter(
          Polar_contact_percent >= threshold
        ) %>%
        pull(
          Residue
        )
    )

  } else {

    return(
      group_data %>%
        filter(
          Any_contact_percent >= threshold
        ) %>%
        pull(
          Residue
        )
    )
  }
}


ND_any <- get_recurrent_set(
  "Neutral dominant",
  "any"
)

NA_any <- get_recurrent_set(
  "Neutral alternative",
  "any"
)

MONO_any <- get_recurrent_set(
  "7-O- monoanion",
  "any"
)

DI_any <- get_recurrent_set(
  "7-O-/4'-O- dianion",
  "any"
)


ND_polar <- get_recurrent_set(
  "Neutral dominant",
  "polar"
)

NA_polar <- get_recurrent_set(
  "Neutral alternative",
  "polar"
)

MONO_polar <- get_recurrent_set(
  "7-O- monoanion",
  "polar"
)

DI_polar <- get_recurrent_set(
  "7-O-/4'-O- dianion",
  "polar"
)


overlap_summary <- data.frame(

  Comparison = c(
    "Monoanion vs neutral alternative",
    "Monoanion vs neutral alternative",
    "Dianion vs neutral dominant",
    "Dianion vs neutral dominant"
  ),

  Contact_type = c(
    "Any",
    "Polar",
    "Any",
    "Polar"
  ),

  Reference_n = c(
    length(NA_any),
    length(NA_polar),
    length(ND_any),
    length(ND_polar)
  ),

  Test_n = c(
    length(MONO_any),
    length(MONO_polar),
    length(DI_any),
    length(DI_polar)
  ),

  Shared_n = c(
    length(
      intersect(
        NA_any,
        MONO_any
      )
    ),

    length(
      intersect(
        NA_polar,
        MONO_polar
      )
    ),

    length(
      intersect(
        ND_any,
        DI_any
      )
    ),

    length(
      intersect(
        ND_polar,
        DI_polar
      )
    )
  )
)


overlap_summary <- overlap_summary %>%

  mutate(

    Reference_retained_percent =
      100 *
      Shared_n /
      Reference_n,

    Jaccard =
      c(
        length(
          intersect(
            NA_any,
            MONO_any
          )
        ) /
          length(
            union(
              NA_any,
              MONO_any
            )
          ),

        length(
          intersect(
            NA_polar,
            MONO_polar
          )
        ) /
          length(
            union(
              NA_polar,
              MONO_polar
            )
          ),

        length(
          intersect(
            ND_any,
            DI_any
          )
        ) /
          length(
            union(
              ND_any,
              DI_any
            )
          ),

        length(
          intersect(
            ND_polar,
            DI_polar
          )
        ) /
          length(
            union(
              ND_polar,
              DI_polar
            )
          )
      )
  )


# ============================================================
# SAVE CLEAN FIGURE DATA
# ============================================================

polar_data_csv <- file.path(
  output_dir,
  "AK1_protonation_POLAR_heatmap_data.csv"
)

any_data_csv <- file.path(
  output_dir,
  "AK1_protonation_ANY_heatmap_data.csv"
)

overlap_csv <- file.path(
  output_dir,
  "AK1_protonation_contact_overlap_summary.csv"
)


write.csv(
  polar_plot_data,
  polar_data_csv,
  row.names = FALSE
)

write.csv(
  any_plot_data,
  any_data_csv,
  row.names = FALSE
)

write.csv(
  overlap_summary,
  overlap_csv,
  row.names = FALSE
)


# ============================================================
# EXCEL
# ============================================================

xlsx_file <- file.path(
  output_dir,
  "AK1_protonation_contact_heatmap_data.xlsx"
)


wb <- createWorkbook()


addWorksheet(
  wb,
  "Polar_Heatmap_Data"
)

writeData(
  wb,
  "Polar_Heatmap_Data",
  polar_plot_data
)


addWorksheet(
  wb,
  "Any_Heatmap_Data"
)

writeData(
  wb,
  "Any_Heatmap_Data",
  any_plot_data
)


addWorksheet(
  wb,
  "Overlap_Summary"
)

writeData(
  wb,
  "Overlap_Summary",
  overlap_summary
)


notes <- data.frame(

  Item = c(
    "Monoanion structural result",
    "Monoanion contact result",
    "Dianion structural result",
    "Dianion any-contact result",
    "Dianion polar-contact result",
    "Important terminology"
  ),

  Description = c(

    paste0(
      "The 7-O- monoanion matched the neutral alternative ",
      "pose family with mean nearest positional RMSD of 0.049 A."
    ),

    paste0(
      "The monoanion recurrent any-contact and polar-contact ",
      "fingerprints were identical to those of the neutral ",
      "alternative reference set."
    ),

    paste0(
      "The dianion matched the neutral dominant pose family ",
      "with mean nearest positional RMSD of 0.392 A."
    ),

    paste0(
      "The dianion retained 12 of 13 recurrent any-contact ",
      "residues from the neutral dominant family (92.3%)."
    ),

    paste0(
      "The dianion retained 5 of 6 recurrent polar-contact ",
      "residues from the neutral dominant family (83.3%). ",
      "ARG149 lost recurrent polar proximity, whereas GLY16, ",
      "GLY22 and PRO17 became recurrent polar-proximity residues."
    ),

    paste0(
      "Polar contact means O/N/S heavy-atom proximity within ",
      "3.6 A and should not automatically be described as a ",
      "hydrogen bond. For glycine and proline residues, the ",
      "identified polar atom can be part of the peptide backbone."
    )
  )
)


addWorksheet(
  wb,
  "Method_Notes"
)

writeData(
  wb,
  "Method_Notes",
  notes
)


setColWidths(
  wb,
  "Polar_Heatmap_Data",
  cols = 1:ncol(polar_plot_data),
  widths = "auto"
)

setColWidths(
  wb,
  "Any_Heatmap_Data",
  cols = 1:ncol(any_plot_data),
  widths = "auto"
)

setColWidths(
  wb,
  "Overlap_Summary",
  cols = 1:ncol(overlap_summary),
  widths = "auto"
)

setColWidths(
  wb,
  "Method_Notes",
  cols = 1,
  widths = 30
)

setColWidths(
  wb,
  "Method_Notes",
  cols = 2,
  widths = 110
)


saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("AK1 PROTONATION CONTACT HEATMAPS COMPLETE\n")
cat("============================================================\n\n")


cat("Overlap summary:\n")
print(
  overlap_summary
)


cat("\nPolar heatmap:\n")
cat(polar_png, "\n")


cat("\nAny-contact heatmap:\n")
cat(any_png, "\n")


cat("\nExcel:\n")
cat(xlsx_file, "\n")