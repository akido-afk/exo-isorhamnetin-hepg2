# ============================================================
# AK1 PANEL L
# Native-ligand interaction fingerprint comparison
#
# OUTPUTS ARE SAVED AS TWO SEPARATE FIGURES:
#
# L-A = Polar-contact residue fingerprint
# L-B = Quantitative fingerprint overlap
#
# ISO dominant ensemble vs crystallographic AP5
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
# OUTPUT DIRECTORY
# ============================================================

out_dir <- "../results/AK1_panelL_native_ligand_fingerprint"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# VERIFIED RESIDUE SETS
# ============================================================

# ------------------------------------------------------------
# ISO recurrent ANY contacts across dominant ensemble
# ------------------------------------------------------------

iso_any <- c(
  "ARG132",
  "ARG138",
  "ARG149",
  "ARG97",
  "ASP140",
  "ASP93",
  "GLY15",
  "GLY16",
  "GLY18",
  "GLY64",
  "LYS21",
  "MET61",
  "PRO17"
)


# ------------------------------------------------------------
# Native AP5 ANY-contact fingerprint
# ------------------------------------------------------------

ap5_any <- c(
  "ARG128",
  "ARG132",
  "ARG138",
  "ARG149",
  "ARG44",
  "ARG97",
  "ASP140",
  "GLN101",
  "GLN65",
  "GLY16",
  "GLY177",
  "GLY18",
  "GLY20",
  "GLY22",
  "GLY40",
  "GLY94",
  "LEU129",
  "LEU43",
  "LEU66",
  "LYS21",
  "MET61",
  "PRO17",
  "SER178",
  "SER19",
  "THR23",
  "THR39",
  "TYR95",
  "VAL179",
  "VAL67",
  "VAL72"
)


# ------------------------------------------------------------
# ISO recurrent POLAR contacts
# ------------------------------------------------------------

iso_polar <- c(
  "LYS21",
  "ARG132",
  "ARG138",
  "ARG149",
  "MET61",
  "ARG97"
)


# ------------------------------------------------------------
# Native AP5 POLAR-contact fingerprint
# ------------------------------------------------------------

ap5_polar <- c(
  "ARG128",
  "ARG132",
  "ARG138",
  "ARG149",
  "ARG44",
  "GLN101",
  "GLN65",
  "GLY177",
  "GLY18",
  "GLY20",
  "GLY22",
  "GLY94",
  "LYS21",
  "SER19",
  "THR23",
  "THR39",
  "VAL179",
  "VAL67"
)


# ============================================================
# CALCULATE OVERLAPS
# ============================================================

shared_any <- intersect(
  iso_any,
  ap5_any
)

iso_only_any <- setdiff(
  iso_any,
  ap5_any
)

ap5_only_any <- setdiff(
  ap5_any,
  iso_any
)


shared_polar <- intersect(
  iso_polar,
  ap5_polar
)

iso_only_polar <- setdiff(
  iso_polar,
  ap5_polar
)

ap5_only_polar <- setdiff(
  ap5_polar,
  iso_polar
)


# ============================================================
# SUMMARY METRICS
# ============================================================

any_iso_coverage <-
  length(shared_any) /
  length(iso_any) *
  100


polar_iso_coverage <-
  length(shared_polar) /
  length(iso_polar) *
  100


any_jaccard <-
  length(shared_any) /
  length(
    union(
      iso_any,
      ap5_any
    )
  )


polar_jaccard <-
  length(shared_polar) /
  length(
    union(
      iso_polar,
      ap5_polar
    )
  )


# ============================================================
# VALIDATION
# ============================================================

stopifnot(
  length(iso_any) == 13
)

stopifnot(
  length(ap5_any) == 30
)

stopifnot(
  length(shared_any) == 10
)

stopifnot(
  length(iso_only_any) == 3
)

stopifnot(
  length(ap5_only_any) == 20
)

stopifnot(
  length(iso_polar) == 6
)

stopifnot(
  length(ap5_polar) == 18
)

stopifnot(
  length(shared_polar) == 4
)

stopifnot(
  length(iso_only_polar) == 2
)

stopifnot(
  length(ap5_only_polar) == 14
)


# ============================================================
# PRINT NUMERICAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("AK1 NATIVE-LIGAND FINGERPRINT COMPARISON\n")
cat("============================================================\n\n")

cat(
  sprintf(
    "ANY CONTACTS: %d/%d ISO contacts shared with AP5 = %.1f%%\n",
    length(shared_any),
    length(iso_any),
    any_iso_coverage
  )
)

cat(
  sprintf(
    "ANY-contact Jaccard = %.3f\n\n",
    any_jaccard
  )
)

cat(
  sprintf(
    "POLAR CONTACTS: %d/%d ISO polar contacts shared with AP5 = %.1f%%\n",
    length(shared_polar),
    length(iso_polar),
    polar_iso_coverage
  )
)

cat(
  sprintf(
    "Polar-contact Jaccard = %.3f\n",
    polar_jaccard
  )
)


# ============================================================
# PANEL L-A
# POLAR-CONTACT RESIDUE FINGERPRINT
# ============================================================


# ------------------------------------------------------------
# RESIDUE ORDER
# ------------------------------------------------------------

shared_order <- c(
  "LYS21",
  "ARG132",
  "ARG138",
  "ARG149"
)

iso_only_order <- c(
  "MET61",
  "ARG97"
)

ap5_only_order <- c(
  "GLY18",
  "SER19",
  "GLY20",
  "GLY22",
  "THR23",
  "THR39",
  "ARG44",
  "GLN65",
  "VAL67",
  "GLY94",
  "GLN101",
  "ARG128",
  "GLY177",
  "VAL179"
)


polar_order <- c(
  shared_order,
  iso_only_order,
  ap5_only_order
)


polar_class <- data.frame(

  Residue = polar_order,

  Category = c(

    rep(
      "Shared",
      length(shared_order)
    ),

    rep(
      "ISO only",
      length(iso_only_order)
    ),

    rep(
      "AP5 only",
      length(ap5_only_order)
    )
  ),

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# LONG FORMAT
# ------------------------------------------------------------

polar_matrix <- polar_class %>%

  mutate(

    ISO_present =
      Residue %in% iso_polar,

    AP5_present =
      Residue %in% ap5_polar
  ) %>%

  pivot_longer(

    cols = c(
      ISO_present,
      AP5_present
    ),

    names_to =
      "Ligand",

    values_to =
      "Present"
  ) %>%

  mutate(

    Ligand = recode(

      Ligand,

      "ISO_present" =
        "ISO dominant\ncluster",

      "AP5_present" =
        "Native AP5"
    ),

    Category = factor(

      Category,

      levels = c(
        "Shared",
        "ISO only",
        "AP5 only"
      )
    ),

    Residue = factor(

      Residue,

      levels = rev(
        polar_order
      )
    ),

    Ligand = factor(

      Ligand,

      levels = c(
        "ISO dominant\ncluster",
        "Native AP5"
      )
    ),

    Cell_type = case_when(

      Present &
        Ligand ==
        "ISO dominant\ncluster" ~
        "ISO",

      Present &
        Ligand ==
        "Native AP5" ~
        "AP5",

      TRUE ~
        "Absent"
    ),

    Cell_label = case_when(

      Present &
        Ligand ==
        "ISO dominant\ncluster" ~
        "8/8",

      Present &
        Ligand ==
        "Native AP5" ~
        "\u2713",

      TRUE ~
        ""
    )
  )


# ------------------------------------------------------------
# BUILD PANEL L-A
# ------------------------------------------------------------

pA <- ggplot(

  polar_matrix,

  aes(
    x = Ligand,
    y = Residue
  )

) +

  geom_tile(

    aes(
      fill = Cell_type
    ),

    color = "grey50",

    linewidth = 0.7,

    width = 0.92,

    height = 0.90
  ) +

  geom_text(

    aes(
      label = Cell_label
    ),

    size = 5.1,

    fontface = "bold"
  ) +

  facet_grid(

    Category ~ .,

    scales = "free_y",

    space = "free_y",

    switch = "y"
  ) +

  scale_fill_manual(

    values = c(

      "ISO" =
        "#E69F00",

      "AP5" =
        "#00BFC4",

      "Absent" =
        "#F2F2F2"
    ),

    guide = "none"
  ) +

  labs(

    title =
      "AK1 polar-contact fingerprint",

    subtitle =
      "ISO dominant ensemble compared with crystallographic AP5",

    x = NULL,

    y = NULL,

    caption =
      paste0(
        "Shared polar contacts: LYS21, ARG132, ARG138 and ARG149. ",
        "ISO-specific polar contacts: MET61 and ARG97."
      )
  ) +

  theme_minimal(
    base_size = 15
  ) +

  theme(

    panel.grid =
      element_blank(),

    axis.text.x =
      element_text(
        size = 15,
        face = "bold",
        color = "black",
        margin = margin(
          t = 8
        )
      ),

    axis.text.y =
      element_text(
        size = 13.5,
        face = "bold",
        color = "black"
      ),

    plot.title =
      element_text(
        size = 20,
        face = "bold",
        hjust = 0,
        margin = margin(
          b = 5
        )
      ),

    plot.subtitle =
      element_text(
        size = 12,
        color = "grey25",
        hjust = 0,
        margin = margin(
          b = 14
        )
      ),

    plot.caption =
      element_text(
        size = 10,
        hjust = 0,
        color = "grey25",
        margin = margin(
          t = 12
        )
      ),

    strip.text.y.left =
      element_text(
        angle = 0,
        size = 12.5,
        face = "bold",
        color = "black"
      ),

    strip.background =
      element_rect(
        fill = "#E8E8E8",
        color = "grey60",
        linewidth = 0.5
      ),

    strip.placement =
      "outside",

    panel.spacing.y =
      unit(
        0.45,
        "lines"
      ),

    plot.margin =
      margin(
        25,
        30,
        30,
        40
      ),

    plot.title.position =
      "plot",

    plot.caption.position =
      "plot"
  )


# ============================================================
# SAVE PANEL L-A
# ============================================================

panelA_png <- file.path(
  out_dir,
  "AK1_panelL_A_polar_contact_fingerprint.png"
)

panelA_pdf <- file.path(
  out_dir,
  "AK1_panelL_A_polar_contact_fingerprint.pdf"
)


ggsave(

  filename =
    panelA_png,

  plot =
    pA,

  width =
    9.5,

  height =
    11,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


ggsave(

  filename =
    panelA_pdf,

  plot =
    pA,

  width =
    9.5,

  height =
    11,

  units =
    "in",

  bg =
    "white"
)


# ============================================================
# PANEL L-B DATA
# QUANTITATIVE OVERLAP
# ============================================================

overlap_summary <- data.frame(

  Contact_type = c(
    "Any contact",
    "Polar contact"
  ),

  ISO_total = c(
    length(iso_any),
    length(iso_polar)
  ),

  Shared_with_AP5 = c(
    length(shared_any),
    length(shared_polar)
  ),

  ISO_coverage_percent = c(
    any_iso_coverage,
    polar_iso_coverage
  ),

  Jaccard = c(
    any_jaccard,
    polar_jaccard
  ),

  stringsAsFactors = FALSE

) %>%

  mutate(

    Contact_type = factor(

      Contact_type,

      levels = c(
        "Any contact",
        "Polar contact"
      )
    ),

    Fraction_label =
      paste0(
        Shared_with_AP5,
        "/",
        ISO_total
      ),

    Percent_label =
      paste0(
        sprintf(
          "%.1f",
          ISO_coverage_percent
        ),
        "%"
      ),

    Main_label =
      paste0(
        Fraction_label,
        "  (",
        Percent_label,
        ")"
      ),

    Jaccard_label =
      paste0(
        "Jaccard = ",
        sprintf(
          "%.3f",
          Jaccard
        )
      )
  )


# ============================================================
# PANEL L-B
# ============================================================

pB <- ggplot(

  overlap_summary,

  aes(
    x = Contact_type,
    y = ISO_coverage_percent
  )

) +

  geom_col(

    width = 0.56,

    fill = "#E69F00",

    color = "black",

    linewidth = 0.85
  ) +

  geom_text(

    aes(
      label = Main_label
    ),

    vjust = -0.65,

    size = 5.0,

    fontface = "bold"
  ) +

  geom_label(

    aes(
      y = 12,
      label = Jaccard_label
    ),

    size = 4.1,

    label.size = 0.3,

    fill = "white"
  ) +

  scale_y_continuous(

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

    expand = c(
      0,
      0
    )
  ) +

  labs(

    title =
      "AK1 fingerprint overlap with native AP5",

    subtitle =
      paste0(
        "Fraction of recurrent ISO contacts that are also present\n",
        "in the crystallographic AP5 contact environment"
      ),

    x = NULL,

    y =
      "Recurrent ISO contacts shared with AP5",

    caption =
      paste0(
        "Any contact: 10/13 ISO residues overlap with AP5.  ",
        "Polar contact: 4/6 ISO residues overlap with AP5.\n",
        "Jaccard values describe symmetric overlap of the complete ",
        "ISO and AP5 residue sets."
      )
  ) +

  theme_classic(
    base_size = 16
  ) +

  theme(

    axis.text.x =
      element_text(
        size = 16,
        face = "bold",
        color = "black",
        margin = margin(
          t = 10
        )
      ),

    axis.text.y =
      element_text(
        size = 13,
        color = "black"
      ),

    axis.title.y =
      element_text(
        size = 15,
        face = "bold",
        margin = margin(
          r = 10
        )
      ),

    plot.title =
      element_text(
        size = 20,
        face = "bold",
        hjust = 0,
        margin = margin(
          b = 5
        )
      ),

    plot.subtitle =
      element_text(
        size = 12,
        color = "grey25",
        hjust = 0,
        margin = margin(
          b = 18
        )
      ),

    plot.caption =
      element_text(
        size = 10.5,
        hjust = 0,
        color = "grey25",
        margin = margin(
          t = 16
        )
      ),

    plot.margin =
      margin(
        30,
        40,
        35,
        35
      ),

    plot.title.position =
      "plot",

    plot.caption.position =
      "plot"
  )


# ============================================================
# SAVE PANEL L-B
# ============================================================

panelB_png <- file.path(
  out_dir,
  "AK1_panelL_B_fingerprint_overlap.png"
)

panelB_pdf <- file.path(
  out_dir,
  "AK1_panelL_B_fingerprint_overlap.pdf"
)


ggsave(

  filename =
    panelB_png,

  plot =
    pB,

  width =
    9,

  height =
    7.5,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


ggsave(

  filename =
    panelB_pdf,

  plot =
    pB,

  width =
    9,

  height =
    7.5,

  units =
    "in",

  bg =
    "white"
)


# ============================================================
# RESIDUE-LEVEL TABLES
# ============================================================

all_any_residues <- sort(
  union(
    iso_any,
    ap5_any
  )
)


any_table <- data.frame(

  Residue =
    all_any_residues,

  ISO =
    all_any_residues %in% iso_any,

  AP5 =
    all_any_residues %in% ap5_any,

  stringsAsFactors =
    FALSE

) %>%

  mutate(

    Classification = case_when(

      ISO & AP5 ~
        "Shared",

      ISO & !AP5 ~
        "ISO only",

      !ISO & AP5 ~
        "AP5 only"
    )
  )


all_polar_residues <- sort(
  union(
    iso_polar,
    ap5_polar
  )
)


polar_table <- data.frame(

  Residue =
    all_polar_residues,

  ISO =
    all_polar_residues %in% iso_polar,

  AP5 =
    all_polar_residues %in% ap5_polar,

  stringsAsFactors =
    FALSE

) %>%

  mutate(

    Classification = case_when(

      ISO & AP5 ~
        "Shared",

      ISO & !AP5 ~
        "ISO only",

      !ISO & AP5 ~
        "AP5 only"
    )
  )


# ============================================================
# SAVE CSV
# ============================================================

summary_csv <- file.path(
  out_dir,
  "AK1_panelL_fingerprint_overlap_summary.csv"
)

any_csv <- file.path(
  out_dir,
  "AK1_panelL_any_contact_fingerprint.csv"
)

polar_csv <- file.path(
  out_dir,
  "AK1_panelL_polar_contact_fingerprint.csv"
)


write.csv(
  overlap_summary,
  summary_csv,
  row.names = FALSE
)

write.csv(
  any_table,
  any_csv,
  row.names = FALSE
)

write.csv(
  polar_table,
  polar_csv,
  row.names = FALSE
)


# ============================================================
# SAVE EXCEL WORKBOOK
# ============================================================

xlsx_file <- file.path(
  out_dir,
  "AK1_panelL_native_ligand_fingerprint_data.xlsx"
)


wb <- createWorkbook()


# ------------------------------------------------------------
# OVERLAP SUMMARY
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Overlap_Summary"
)

writeData(
  wb,
  "Overlap_Summary",
  overlap_summary
)


# ------------------------------------------------------------
# ANY CONTACT
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Any_Contact_Fingerprint"
)

writeData(
  wb,
  "Any_Contact_Fingerprint",
  any_table
)


# ------------------------------------------------------------
# POLAR CONTACT
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Polar_Contact_Fingerprint"
)

writeData(
  wb,
  "Polar_Contact_Fingerprint",
  polar_table
)


# ------------------------------------------------------------
# RESIDUE SETS
# ------------------------------------------------------------

key_sets <- data.frame(

  Category = c(
    "Any shared",
    "Any ISO only",
    "Any AP5 only",
    "Polar shared",
    "Polar ISO only",
    "Polar AP5 only"
  ),

  Residues = c(

    paste(
      shared_any,
      collapse = ", "
    ),

    paste(
      iso_only_any,
      collapse = ", "
    ),

    paste(
      ap5_only_any,
      collapse = ", "
    ),

    paste(
      shared_polar,
      collapse = ", "
    ),

    paste(
      iso_only_polar,
      collapse = ", "
    ),

    paste(
      ap5_only_polar,
      collapse = ", "
    )
  )
)


addWorksheet(
  wb,
  "Residue_Sets"
)

writeData(
  wb,
  "Residue_Sets",
  key_sets
)


# ------------------------------------------------------------
# NOTES
# ------------------------------------------------------------

notes <- data.frame(

  Item = c(
    "Main interpretation",
    "Any-contact overlap",
    "Polar-contact overlap",
    "Jaccard",
    "Important limitation"
  ),

  Description = c(

    paste0(
      "Most recurrent ISO contacts fall within the native AP5 ",
      "contact environment, supporting localization of ISO ",
      "within the experimentally defined AK1 nucleotide-binding region."
    ),

    paste0(
      length(shared_any),
      " of ",
      length(iso_any),
      " recurrent ISO contact residues (",
      sprintf(
        "%.1f",
        any_iso_coverage
      ),
      "%) were also contacted by AP5."
    ),

    paste0(
      length(shared_polar),
      " of ",
      length(iso_polar),
      " recurrent ISO polar-contact residues (",
      sprintf(
        "%.1f",
        polar_iso_coverage
      ),
      "%) were also contacted by AP5."
    ),

    paste0(
      "Jaccard values are ",
      sprintf(
        "%.3f",
        any_jaccard
      ),
      " for any contacts and ",
      sprintf(
        "%.3f",
        polar_jaccard
      ),
      " for polar contacts."
    ),

    paste0(
      "Fingerprint overlap supports shared pocket occupancy but ",
      "does not demonstrate identical binding mode, mechanism, ",
      "or experimental affinity."
    )
  )
)


addWorksheet(
  wb,
  "Notes"
)

writeData(
  wb,
  "Notes",
  notes
)


# ------------------------------------------------------------
# COLUMN WIDTHS
# ------------------------------------------------------------

setColWidths(
  wb,
  "Overlap_Summary",
  cols = 1:ncol(overlap_summary),
  widths = "auto"
)

setColWidths(
  wb,
  "Any_Contact_Fingerprint",
  cols = 1:ncol(any_table),
  widths = "auto"
)

setColWidths(
  wb,
  "Polar_Contact_Fingerprint",
  cols = 1:ncol(polar_table),
  widths = "auto"
)

setColWidths(
  wb,
  "Residue_Sets",
  cols = 1,
  widths = 22
)

setColWidths(
  wb,
  "Residue_Sets",
  cols = 2,
  widths = 110
)

setColWidths(
  wb,
  "Notes",
  cols = 1,
  widths = 28
)

setColWidths(
  wb,
  "Notes",
  cols = 2,
  widths = 110
)


saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PANEL L COMPLETE — TWO SEPARATE FIGURES\n")
cat("============================================================\n\n")

cat(
  sprintf(
    "ANY CONTACTS: %d/%d = %.1f%% shared with AP5\n",
    length(shared_any),
    length(iso_any),
    any_iso_coverage
  )
)

cat(
  sprintf(
    "ANY-contact Jaccard = %.3f\n\n",
    any_jaccard
  )
)

cat(
  sprintf(
    "POLAR CONTACTS: %d/%d = %.1f%% shared with AP5\n",
    length(shared_polar),
    length(iso_polar),
    polar_iso_coverage
  )
)

cat(
  sprintf(
    "Polar-contact Jaccard = %.3f\n",
    polar_jaccard
  )
)

cat("\nPanel L-A:\n")
cat(panelA_png, "\n")

cat("\nPanel L-B:\n")
cat(panelB_png, "\n")

cat("\nExcel:\n")
cat(xlsx_file, "\n")