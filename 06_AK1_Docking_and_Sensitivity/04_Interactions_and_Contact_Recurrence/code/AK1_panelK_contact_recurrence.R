# ============================================================
# AK1 PANEL K
# Interaction-frequency analysis across the dominant 8 ISO poses
# ============================================================

packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "openxlsx"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)

# ============================================================
# PATHS
# ============================================================

input_csv <- "../results/AK1_dominant8_contact_frequency/AK1_dominant8_contact_frequency.csv"
out_dir   <- "../results/AK1_panelK_contact_recurrence"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# READ DATA
# ============================================================

raw <- read.csv(
  input_csv,
  check.names = TRUE,
  stringsAsFactors = FALSE
)

cat("\nDetected columns:\n")
print(colnames(raw))

# ============================================================
# CLEAN TABLE
# ============================================================

freq <- raw %>%
  transmute(
    Residue                = as.character(Residue),
    Any_n                  = as.numeric(Any.contact.n),
    Any_percent            = as.numeric(Any.contact..),
    Any_mean_distance      = as.numeric(Any.mean.min.distance),
    Any_shortest_distance  = as.numeric(Any.shortest.distance),
    Polar_n                = as.numeric(Polar.contact.n),
    Polar_percent          = as.numeric(Polar.contact..),
    Polar_mean_distance    = as.numeric(Polar.mean.min.distance),
    Polar_shortest_distance= as.numeric(Polar.shortest.distance)
  )

# ============================================================
# CLASSIFICATION
# ============================================================

freq <- freq %>%
  mutate(
    Fully_recurrent_polar = Polar_n == 8,
    Fully_recurrent_any   = Any_n == 8,
    Contact_class = case_when(
      Polar_n == 8 ~ "Recurrent polar",
      Any_n == 8   ~ "Recurrent non-polar/proximity",
      TRUE         ~ "Occasional contact"
    )
  )

# ============================================================
# ORDER RESIDUES
# ============================================================

freq <- freq %>%
  arrange(
    desc(Polar_n),
    desc(Any_n),
    Polar_mean_distance,
    Residue
  )

residue_order <- freq$Residue

# ============================================================
# LONG FORMAT FOR PLOTTING
# ============================================================

plot_df <- bind_rows(
  freq %>%
    transmute(
      Residue,
      Contact_type = "Any contact",
      Count = Any_n
    ),
  freq %>%
    transmute(
      Residue,
      Contact_type = "Polar contact",
      Count = Polar_n
    )
) %>%
  mutate(
    Contact_type = factor(
      Contact_type,
      levels = c("Any contact", "Polar contact")
    ),
    Residue = factor(
      Residue,
      levels = rev(residue_order)
    ),
    Label = paste0(Count, "/8")
  )

# ============================================================
# CORE RECURRENT POLAR RESIDUES
# ============================================================

core <- freq %>%
  filter(Polar_n == 8)

core_names <- paste(core$Residue, collapse = ", ")

# ============================================================
# PRINT CORE RESULTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FULLY RECURRENT POLAR-CONTACT RESIDUES\n")
cat("============================================================\n")
print(
  core %>%
    select(
      Residue,
      Polar_n,
      Polar_mean_distance,
      Polar_shortest_distance
    )
)

# ============================================================
# FIGURE
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = Contact_type,
    y = Residue,
    fill = Count
  )
) +
  geom_tile(
    color = "grey65",
    linewidth = 0.7,
    width = 0.90,
    height = 0.88
  ) +
  geom_tile(
    data = subset(
      plot_df,
      Contact_type == "Polar contact" & Count == 8
    ),
    fill = NA,
    color = "black",
    linewidth = 1.4,
    width = 0.90,
    height = 0.88
  ) +
  geom_text(
    aes(label = Label),
    size = 5.2,
    fontface = "bold"
  ) +
  scale_fill_gradientn(
    colours = c("#F2F2F2", "#F7D98A", "#E69F00"),
    values  = c(0, 0.5, 1),
    limits  = c(0, 8),
    breaks  = c(0, 2, 4, 6, 8),
    name    = "Number of\nposes"
  ) +
  labs(
    title = "Interaction recurrence across the dominant AK1–ISO ensemble",
    subtitle = paste0(
      "Six residues showed polar-contact recurrence in all 8 dominant poses: ",
      core_names
    ),
    x = NULL,
    y = "AK1 residue",
    caption = paste0(
      "Values indicate the number of dominant poses containing each contact (n = 8). ",
      "Black outlines mark polar contacts present in all 8 poses."
    )
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 17,
      face = "bold"
    ),
    plot.title = element_text(
      size = 17,        # smaller than before
      face = "bold",
      hjust = 0,
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 11.5,      # smaller subtitle
      hjust = 0,
      margin = margin(b = 14)
    ),
    plot.caption = element_text(
      size = 10,
      hjust = 0,
      color = "grey25",
      margin = margin(t = 12)
    ),
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    legend.text = element_text(
      size = 12
    ),
    legend.position = "right",
    plot.margin = margin(
      t = 20,
      r = 38,
      b = 28,
      l = 24
    ),
    plot.title.position = "plot",
    plot.caption.position = "plot"
  )

# ============================================================
# SAVE FIGURE
# ============================================================

png_file <- file.path(
  out_dir,
  "AK1_panelK_dominant8_contact_recurrence.png"
)

pdf_file <- file.path(
  out_dir,
  "AK1_panelK_dominant8_contact_recurrence.pdf"
)

ggsave(
  filename = png_file,
  plot = p,
  width = 14.5,     # slightly wider
  height = 9.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 14.5,
  height = 9.5,
  units = "in",
  bg = "white"
)

# ============================================================
# SAVE CLEAN DATA — CSV
# ============================================================

csv_file <- file.path(
  out_dir,
  "AK1_panelK_dominant8_contact_recurrence_data.csv"
)

write.csv(
  freq,
  csv_file,
  row.names = FALSE
)

# ============================================================
# SAVE CLEAN DATA — EXCEL
# ============================================================

xlsx_file <- file.path(
  out_dir,
  "AK1_panelK_dominant8_contact_recurrence_data.xlsx"
)

wb <- createWorkbook()

addWorksheet(wb, "Contact_Recurrence")
writeData(wb, "Contact_Recurrence", freq)

addWorksheet(wb, "Core_Polar_8of8")
writeData(wb, "Core_Polar_8of8", core)

addWorksheet(wb, "Notes")

notes <- data.frame(
  Item = c(
    "Main result",
    "Dominant pose ensemble",
    "Any contact",
    "Polar contact",
    "Interpretation limitation"
  ),
  Description = c(
    paste0(
      "Six residues showed polar-contact recurrence in all 8 dominant ISO poses: ",
      core_names, "."
    ),
    "Eight Zn-free ISO poses belonging to the previously identified dominant docking cluster.",
    "Any qualifying ligand-residue heavy-atom proximity from the automated contact analysis.",
    "Distance/type-based candidate polar interaction from the automated contact analysis.",
    paste0(
      "The 8/8 recurrence demonstrates stability within the dominant docking pose family. ",
      "Because the eight dominant poses are themselves geometrically very similar, ",
      "this should not be interpreted as eight independent experimental confirmations."
    )
  )
)

writeData(wb, "Notes", notes)

setColWidths(wb, "Contact_Recurrence", cols = 1:ncol(freq), widths = "auto")
setColWidths(wb, "Core_Polar_8of8", cols = 1:ncol(core), widths = "auto")
setColWidths(wb, "Notes", cols = 1, widths = 28)
setColWidths(wb, "Notes", cols = 2, widths = 100)

saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)

# ============================================================
# FINAL CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("AK1 PANEL K COMPLETE\n")
cat("============================================================\n")

cat("\nCore 8/8 polar-contact residues:\n")
for (i in seq_len(nrow(core))) {
  cat(
    sprintf(
      "%-8s  %d/8   mean minimum distance = %.3f A   shortest = %.3f A\n",
      core$Residue[i],
      core$Polar_n[i],
      core$Polar_mean_distance[i],
      core$Polar_shortest_distance[i]
    )
  )
}

cat("\nSaved files:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat(csv_file, "\n")
cat(xlsx_file, "\n")