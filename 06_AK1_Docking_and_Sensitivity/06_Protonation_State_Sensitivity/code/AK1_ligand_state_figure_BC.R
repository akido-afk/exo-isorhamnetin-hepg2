# ============================================================
# AK1 LIGAND-STATE SENSITIVITY
#
# PANEL B:
#   Pose-family occupancy across protonation states
#
# PANEL C:
#   Per-seed RMSD to neutral dominant and alternative
#   reference pose families
#
# Figures are saved SEPARATELY.
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

analysis_dir <- paste0(
  "../results/",
  "AK1_ligand_state_sensitivity/",
  "analysis_v1"
)

figure_dir <- paste0(
  "../results/",
  "AK1_ligand_state_sensitivity/",
  "figures"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# INPUT FILES
# ============================================================

neutral_matrix_file <- file.path(
  analysis_dir,
  "AK1_ISO_Neutral_pairwise_RMSD.csv"
)

mono_file <- file.path(
  analysis_dir,
  "AK1_monoanion_vs_neutral_pose_families_per_seed.csv"
)

dianion_file <- file.path(
  analysis_dir,
  "AK1_dianion_vs_neutral_pose_families_per_seed.csv"
)


# ============================================================
# CHECK INPUT FILES
# ============================================================

for (f in c(
  neutral_matrix_file,
  mono_file,
  dianion_file
)) {

  if (!file.exists(f)) {

    stop(
      paste(
        "Missing input file:",
        f
      )
    )
  }
}


# ============================================================
# SEEDS / NEUTRAL REFERENCE FAMILIES
# ============================================================

all_seeds <- c(
  1001,
  2002,
  3003,
  4004,
  5005,
  6006,
  7007,
  8008,
  9009,
  10101
)


dominant_seeds <- c(
  1001,
  2002,
  3003,
  4004,
  5005,
  7007,
  9009,
  10101
)


alternative_seeds <- c(
  6006,
  8008
)


# ============================================================
# READ NEUTRAL PAIRWISE RMSD MATRIX
# ============================================================

neutral_raw <- read.csv(
  neutral_matrix_file,
  check.names = FALSE
)


rownames(
  neutral_raw
) <- as.character(
  neutral_raw$Seed
)


neutral_matrix <- as.matrix(
  neutral_raw[
    ,
    -1,
    drop = FALSE
  ]
)


mode(
  neutral_matrix
) <- "numeric"


# ============================================================
# VALIDATE NEUTRAL RMSD MATRIX
# ============================================================

if (!all(
  as.character(all_seeds)
  %in%
  rownames(neutral_matrix)
)) {

  stop(
    "Not all expected neutral seeds are present as matrix rows."
  )
}


if (!all(
  as.character(all_seeds)
  %in%
  colnames(neutral_matrix)
)) {

  stop(
    "Not all expected neutral seeds are present as matrix columns."
  )
}


# ============================================================
# HELPER:
# FIND NEAREST MEMBER OF A REFERENCE FAMILY
#
# drop = FALSE is critical here.
# It prevents R from turning a one-column matrix into
# an unnamed scalar when only one candidate remains.
# ============================================================

get_nearest_reference <- function(
  query_seed,
  candidate_seeds
) {

  if (length(candidate_seeds) == 0) {

    stop(
      paste(
        "No reference candidates available for seed",
        query_seed
      )
    )
  }


  submatrix <- neutral_matrix[
    as.character(query_seed),
    as.character(candidate_seeds),
    drop = FALSE
  ]


  values <- as.numeric(
    submatrix[1, ]
  )


  names(values) <- colnames(
    submatrix
  )


  min_index <- which.min(
    values
  )


  list(

    seed = as.integer(
      names(values)[
        min_index
      ]
    ),

    rmsd = as.numeric(
      values[
        min_index
      ]
    )
  )
}


# ============================================================
# NEUTRAL STATE:
# DISTANCE TO EACH REFERENCE FAMILY
#
# IMPORTANT:
# If a neutral pose belongs to the family being tested,
# that exact seed itself is excluded.
#
# This prevents a trivial self-RMSD of 0 A.
# ============================================================

neutral_rows <- data.frame()


for (seed in all_seeds) {


  # ----------------------------------------------------------
  # DOMINANT REFERENCE FAMILY
  # ----------------------------------------------------------

  dom_candidates <- dominant_seeds


  if (seed %in% dominant_seeds) {

    dom_candidates <- setdiff(
      dominant_seeds,
      seed
    )
  }


  nearest_dom <- get_nearest_reference(
    query_seed = seed,
    candidate_seeds = dom_candidates
  )


  # ----------------------------------------------------------
  # ALTERNATIVE REFERENCE FAMILY
  # ----------------------------------------------------------

  alt_candidates <- alternative_seeds


  if (seed %in% alternative_seeds) {

    alt_candidates <- setdiff(
      alternative_seeds,
      seed
    )
  }


  nearest_alt <- get_nearest_reference(
    query_seed = seed,
    candidate_seeds = alt_candidates
  )


  # ----------------------------------------------------------
  # ORIGINAL NEUTRAL FAMILY ASSIGNMENT
  # ----------------------------------------------------------

  observed_family <- ifelse(
    seed %in% dominant_seeds,
    "Dominant",
    "Alternative"
  )


  # ----------------------------------------------------------
  # ADD ROW
  # ----------------------------------------------------------

  neutral_rows <- rbind(

    neutral_rows,

    data.frame(

      State =
        "Neutral",

      Seed =
        seed,

      Observed_family =
        observed_family,

      Nearest_dominant_seed =
        nearest_dom$seed,

      Nearest_dominant_RMSD_A =
        nearest_dom$rmsd,

      Nearest_alternative_seed =
        nearest_alt$seed,

      Nearest_alternative_RMSD_A =
        nearest_alt$rmsd,

      stringsAsFactors =
        FALSE
    )
  )
}


# ============================================================
# PRINT NEUTRAL CALCULATIONS FOR VALIDATION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("NEUTRAL REFERENCE-FAMILY CALCULATIONS\n")
cat("============================================================\n\n")

print(
  neutral_rows
)


# ============================================================
# READ MONOANION RESULTS
# ============================================================

mono_raw <- read.csv(
  mono_file,
  check.names = FALSE
)


required_mono_cols <- c(
  "Monoanion_seed",
  "Nearest_dominant_seed",
  "Nearest_dominant_RMSD_A",
  "Nearest_alternative_seed",
  "Nearest_alternative_RMSD_A"
)


if (!all(
  required_mono_cols
  %in%
  colnames(mono_raw)
)) {

  stop(
    paste(
      "Monoanion CSV does not contain all required columns.\n",
      "Found:",
      paste(
        colnames(mono_raw),
        collapse = ", "
      )
    )
  )
}


mono_rows <- data.frame(

  State =
    "7-O\u207B monoanion",

  Seed =
    mono_raw$Monoanion_seed,

  Observed_family =
    "Alternative",

  Nearest_dominant_seed =
    mono_raw$Nearest_dominant_seed,

  Nearest_dominant_RMSD_A =
    mono_raw$Nearest_dominant_RMSD_A,

  Nearest_alternative_seed =
    mono_raw$Nearest_alternative_seed,

  Nearest_alternative_RMSD_A =
    mono_raw$Nearest_alternative_RMSD_A,

  stringsAsFactors =
    FALSE
)


# ============================================================
# READ DIANION RESULTS
# ============================================================

dianion_raw <- read.csv(
  dianion_file,
  check.names = FALSE
)


required_dianion_cols <- c(
  "Dianion_seed",
  "Nearest_dominant_seed",
  "Nearest_dominant_RMSD_A",
  "Nearest_alternative_seed",
  "Nearest_alternative_RMSD_A"
)


if (!all(
  required_dianion_cols
  %in%
  colnames(dianion_raw)
)) {

  stop(
    paste(
      "Dianion CSV does not contain all required columns.\n",
      "Found:",
      paste(
        colnames(dianion_raw),
        collapse = ", "
      )
    )
  )
}


dianion_rows <- data.frame(

  State =
    "7-O\u207B/4\u2032-O\u207B dianion",

  Seed =
    dianion_raw$Dianion_seed,

  Observed_family =
    "Dominant",

  Nearest_dominant_seed =
    dianion_raw$Nearest_dominant_seed,

  Nearest_dominant_RMSD_A =
    dianion_raw$Nearest_dominant_RMSD_A,

  Nearest_alternative_seed =
    dianion_raw$Nearest_alternative_seed,

  Nearest_alternative_RMSD_A =
    dianion_raw$Nearest_alternative_RMSD_A,

  stringsAsFactors =
    FALSE
)


# ============================================================
# COMBINE ALL THREE STATES
# ============================================================

all_state_rows <- bind_rows(
  neutral_rows,
  mono_rows,
  dianion_rows
)


state_levels <- c(
  "Neutral",
  "7-O\u207B monoanion",
  "7-O\u207B/4\u2032-O\u207B dianion"
)


all_state_rows$State <- factor(
  all_state_rows$State,
  levels = state_levels
)


# ============================================================
# VALIDATE 10 SEEDS PER STATE
# ============================================================

state_counts <- all_state_rows %>%
  count(
    State
  )


if (!all(
  state_counts$n == 10
)) {

  stop(
    "Expected exactly 10 seed results for each ligand state."
  )
}


# ============================================================
# PANEL B DATA:
# POSE-FAMILY OCCUPANCY
# ============================================================

occupancy <- all_state_rows %>%

  count(
    State,
    Observed_family,
    name = "n"
  ) %>%

  complete(

    State,

    Observed_family = c(
      "Dominant",
      "Alternative"
    ),

    fill = list(
      n = 0
    )
  ) %>%

  group_by(
    State
  ) %>%

  mutate(

    Percent =
      100 *
      n /
      sum(n),

    Label =
      ifelse(
        n > 0,
        paste0(
          n,
          "/10\n",
          sprintf(
            "%.0f%%",
            Percent
          )
        ),
        ""
      )
  ) %>%

  ungroup()


occupancy$Observed_family <- factor(

  occupancy$Observed_family,

  levels = c(
    "Dominant",
    "Alternative"
  )
)


# ============================================================
# PANEL B:
# POSE-FAMILY OCCUPANCY
# ============================================================

pB <- ggplot(

  occupancy,

  aes(
    x = State,
    y = Percent,
    fill = Observed_family
  )

) +

  geom_col(

    width = 0.64,

    color = "black",

    linewidth = 0.75
  ) +

  geom_text(

    aes(
      label = Label
    ),

    position = position_stack(
      vjust = 0.5
    ),

    size = 5.0,

    fontface = "bold"
  ) +

  scale_fill_manual(

    values = c(

      "Dominant" =
        "#E69F00",

      "Alternative" =
        "#56B4E9"
    ),

    name =
      "Neutral reference\npose family"
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
      "AK1 pose-family occupancy across ISO protonation states",

    subtitle =
      paste0(
        "Neutral ISO sampled two pose families, whereas each charged ",
        "state converged exclusively to one pre-existing family"
      ),

    x =
      NULL,

    y =
      "Best poses assigned to neutral reference family"
  ) +

  theme_classic(
    base_size = 15
  ) +

  theme(

    axis.text.x =
      element_text(
        size = 14,
        face = "bold",
        color = "black",
        margin = margin(
          t = 9
        )
      ),

    axis.text.y =
      element_text(
        size = 12,
        color = "black"
      ),

    axis.title.y =
      element_text(
        size = 14,
        face = "bold",
        margin = margin(
          r = 10
        )
      ),

    legend.title =
      element_text(
        size = 12,
        face = "bold"
      ),

    legend.text =
      element_text(
        size = 12
      ),

    legend.position =
      "right",

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
          b = 15
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
      "plot"
  )


# ============================================================
# SAVE PANEL B
# ============================================================

panelB_png <- file.path(
  figure_dir,
  "AK1_ligand_state_panelB_pose_family_occupancy.png"
)

panelB_pdf <- file.path(
  figure_dir,
  "AK1_ligand_state_panelB_pose_family_occupancy.pdf"
)


ggsave(
  filename = panelB_png,
  plot = pB,
  width = 10,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = panelB_pdf,
  plot = pB,
  width = 10,
  height = 7.5,
  units = "in",
  bg = "white"
)


# ============================================================
# PANEL C DATA:
# PER-SEED RMSD TO BOTH NEUTRAL REFERENCE FAMILIES
# ============================================================

rmsd_long <- all_state_rows %>%

  select(
    State,
    Seed,
    Observed_family,
    Nearest_dominant_RMSD_A,
    Nearest_alternative_RMSD_A
  ) %>%

  pivot_longer(

    cols = c(
      Nearest_dominant_RMSD_A,
      Nearest_alternative_RMSD_A
    ),

    names_to =
      "Reference_family",

    values_to =
      "Nearest_RMSD_A"
  ) %>%

  mutate(

    Reference_family = recode(

      Reference_family,

      "Nearest_dominant_RMSD_A" =
        "Dominant family",

      "Nearest_alternative_RMSD_A" =
        "Alternative family"
    ),

    Reference_family = factor(

      Reference_family,

      levels = c(
        "Dominant family",
        "Alternative family"
      )
    )
  )


# ============================================================
# <= 2 A COUNTS
# ============================================================

within2_summary <- rmsd_long %>%

  group_by(
    State,
    Reference_family
  ) %>%

  summarise(

    n_within_2A =
      sum(
        Nearest_RMSD_A <= 2.0
      ),

    .groups =
      "drop"
  ) %>%

  mutate(

    Label =
      paste0(
        "\u22642 \u00C5: ",
        n_within_2A,
        "/10"
      )
  )


# ============================================================
# RMSD SUMMARY STATISTICS
# ============================================================

rmsd_summary <- rmsd_long %>%

  group_by(
    State,
    Reference_family
  ) %>%

  summarise(

    Mean_RMSD_A =
      mean(
        Nearest_RMSD_A
      ),

    Median_RMSD_A =
      median(
        Nearest_RMSD_A
      ),

    SD_RMSD_A =
      sd(
        Nearest_RMSD_A
      ),

    Minimum_RMSD_A =
      min(
        Nearest_RMSD_A
      ),

    Maximum_RMSD_A =
      max(
        Nearest_RMSD_A
      ),

    .groups =
      "drop"
  )


# ============================================================
# PANEL C:
# PER-SEED FAMILY RMSD
# ============================================================

set.seed(
  1001
)


pC <- ggplot(

  rmsd_long,

  aes(
    x = Reference_family,
    y = Nearest_RMSD_A
  )

) +

  geom_hline(

    yintercept = 2.0,

    linetype = "dashed",

    linewidth = 0.8,

    color = "grey35"
  ) +

  geom_hline(

    yintercept = 3.0,

    linetype = "dotted",

    linewidth = 0.65,

    color = "grey60"
  ) +

  geom_jitter(

    aes(
      fill = Reference_family
    ),

    shape = 21,

    width = 0.10,

    height = 0,

    size = 4.2,

    stroke = 0.8,

    color = "black"
  ) +

  stat_summary(

    fun = median,

    geom = "crossbar",

    width = 0.42,

    linewidth = 0.9,

    color = "black",

    fatten = 0
  ) +

  geom_label(

    data =
      within2_summary,

    aes(
      x = Reference_family,
      y = 7.85,
      label = Label
    ),

    inherit.aes = FALSE,

    size = 3.7,

    fontface = "bold",

    fill = "white",

    label.size = 0.25
  ) +

  facet_wrap(
    ~ State,
    nrow = 1
  ) +

  scale_fill_manual(

    values = c(

      "Dominant family" =
        "#E69F00",

      "Alternative family" =
        "#56B4E9"
    ),

    guide =
      "none"
  ) +

  scale_y_continuous(

    limits = c(
      0,
      8.25
    ),

    breaks = seq(
      0,
      8,
      by = 1
    ),

    expand = c(
      0,
      0
    )
  ) +

  labs(

    title =
      "Per-seed structural similarity to the neutral AK1 pose families",

    subtitle =
      paste0(
        "Nearest heavy-atom positional RMSD to each neutral reference family; ",
        "exact neutral self-comparisons were excluded"
      ),

    x =
      NULL,

    y =
      "Nearest positional RMSD (\u00C5)",

    caption =
      paste0(
        "Dashed line: 2 \u00C5 family-assignment threshold; ",
        "dotted line: 3 \u00C5 sensitivity threshold. ",
        "Horizontal black bars indicate medians."
      )
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(

    strip.background =
      element_rect(
        fill = "#EFEFEF",
        color = "grey55",
        linewidth = 0.5
      ),

    strip.text =
      element_text(
        size = 12.5,
        face = "bold",
        color = "black",
        margin = margin(
          8,
          5,
          8,
          5
        )
      ),

    axis.text.x =
      element_text(
        size = 11.5,
        face = "bold",
        color = "black",
        margin = margin(
          t = 7
        )
      ),

    axis.text.y =
      element_text(
        size = 11.5,
        color = "black"
      ),

    axis.title.y =
      element_text(
        size = 13.5,
        face = "bold",
        margin = margin(
          r = 10
        )
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

    panel.spacing =
      unit(
        1.2,
        "lines"
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
# SAVE PANEL C
# ============================================================

panelC_png <- file.path(
  figure_dir,
  "AK1_ligand_state_panelC_family_RMSD.png"
)

panelC_pdf <- file.path(
  figure_dir,
  "AK1_ligand_state_panelC_family_RMSD.pdf"
)


ggsave(
  filename = panelC_png,
  plot = pC,
  width = 14,
  height = 7.8,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = panelC_pdf,
  plot = pC,
  width = 14,
  height = 7.8,
  units = "in",
  bg = "white"
)


# ============================================================
# SAVE CLEANED CSV DATA
# ============================================================

occupancy_csv <- file.path(
  figure_dir,
  "AK1_ligand_state_pose_family_occupancy.csv"
)

rmsd_csv <- file.path(
  figure_dir,
  "AK1_ligand_state_family_RMSD_all_seeds.csv"
)

within2_csv <- file.path(
  figure_dir,
  "AK1_ligand_state_family_RMSD_within2A_summary.csv"
)

rmsd_summary_csv <- file.path(
  figure_dir,
  "AK1_ligand_state_family_RMSD_summary_statistics.csv"
)


write.csv(
  occupancy,
  occupancy_csv,
  row.names = FALSE
)

write.csv(
  rmsd_long,
  rmsd_csv,
  row.names = FALSE
)

write.csv(
  within2_summary,
  within2_csv,
  row.names = FALSE
)

write.csv(
  rmsd_summary,
  rmsd_summary_csv,
  row.names = FALSE
)


# ============================================================
# EXCEL WORKBOOK
# ============================================================

xlsx_file <- file.path(
  figure_dir,
  "AK1_ligand_state_figure_BC_data.xlsx"
)


wb <- createWorkbook()


addWorksheet(
  wb,
  "Pose_Family_Occupancy"
)

writeData(
  wb,
  "Pose_Family_Occupancy",
  occupancy
)


addWorksheet(
  wb,
  "Per_Seed_Family_RMSD"
)

writeData(
  wb,
  "Per_Seed_Family_RMSD",
  rmsd_long
)


addWorksheet(
  wb,
  "RMSD_Summary"
)

writeData(
  wb,
  "RMSD_Summary",
  rmsd_summary
)


addWorksheet(
  wb,
  "Within_2A_Summary"
)

writeData(
  wb,
  "Within_2A_Summary",
  within2_summary
)


addWorksheet(
  wb,
  "Neutral_Calculations"
)

writeData(
  wb,
  "Neutral_Calculations",
  neutral_rows
)


notes <- data.frame(

  Item = c(
    "Neutral state",
    "7-O- monoanion",
    "7-O-/4'-O- dianion",
    "Neutral self-comparison",
    "RMSD convention",
    "Interpretation"
  ),

  Description = c(

    paste0(
      "Neutral ISO produced an 8/10 dominant pose family and ",
      "2/10 alternative pose family."
    ),

    paste0(
      "All 10 monoanion best poses corresponded to the neutral ",
      "alternative family. Mean nearest RMSD to that family was ",
      "0.049 A; mean nearest RMSD to the dominant family was ",
      "7.338 A."
    ),

    paste0(
      "All 10 dianion best poses corresponded to the neutral ",
      "dominant family. Mean nearest RMSD to that family was ",
      "0.392 A; mean nearest RMSD to the alternative family was ",
      "7.249 A."
    ),

    paste0(
      "For neutral poses, the exact same seed was excluded when ",
      "calculating the nearest RMSD to the pose family containing ",
      "that seed, preventing trivial self-RMSD values of 0 A."
    ),

    paste0(
      "Direct heavy-atom positional RMSD was calculated in the ",
      "common AK1 receptor coordinate frame without post-docking ",
      "structural superposition."
    ),

    paste0(
      "The results indicate protonation-state-dependent selection ",
      "between two pose families already sampled by neutral ISO. ",
      "They do not establish experimental state populations or ",
      "thermodynamic stabilization."
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
  "Pose_Family_Occupancy",
  cols = 1:ncol(occupancy),
  widths = "auto"
)

setColWidths(
  wb,
  "Per_Seed_Family_RMSD",
  cols = 1:ncol(rmsd_long),
  widths = "auto"
)

setColWidths(
  wb,
  "RMSD_Summary",
  cols = 1:ncol(rmsd_summary),
  widths = "auto"
)

setColWidths(
  wb,
  "Within_2A_Summary",
  cols = 1:ncol(within2_summary),
  widths = "auto"
)

setColWidths(
  wb,
  "Neutral_Calculations",
  cols = 1:ncol(neutral_rows),
  widths = "auto"
)

setColWidths(
  wb,
  "Method_Notes",
  cols = 1,
  widths = 28
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
cat("AK1 LIGAND-STATE FIGURE B/C COMPLETE\n")
cat("============================================================\n\n")


cat("POSE-FAMILY OCCUPANCY\n")
cat("---------------------\n")

print(
  occupancy
)


cat("\n")
cat("RMSD SUMMARY\n")
cat("------------\n")

print(
  rmsd_summary
)


cat("\n")
cat("<= 2 A FAMILY MATCHES\n")
cat("---------------------\n")

print(
  within2_summary
)


cat("\nSaved Panel B:\n")
cat(panelB_png, "\n")
cat(panelB_pdf, "\n")


cat("\nSaved Panel C:\n")
cat(panelC_png, "\n")
cat(panelC_pdf, "\n")


cat("\nSaved workbook:\n")
cat(xlsx_file, "\n")


cat("\nAnalysis and figure generation complete.\n")