library(ggplot2)
library(dplyr)
library(tidyr)

# =========================================================
# INPUT / OUTPUT
# =========================================================
input_csv <- "../results/AK1_dominant8_contact_frequency/AK1_dominant8_contact_frequency.csv"
out_dir   <- "../results/AK1_panelG"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# HELPERS
# =========================================================
parse_ratio_numerator <- function(x) {
  x <- trimws(as.character(x))
  as.numeric(sub("/.*", "", x))
}

clean_distance <- function(x) {
  x <- as.character(x)
  x[x == "-" | x == "" | is.na(x)] <- NA
  suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
}

# =========================================================
# READ TABLE
# Expected CSV columns:
# Residue, Any, Polar, Polar., Min.polar.dist.
# =========================================================
raw <- read.csv(input_csv, check.names = TRUE, stringsAsFactors = FALSE)

# Safety check for expected columns
print("Detected columns:")
print(colnames(raw))

freq <- raw %>%
  transmute(
    Residue = Residue,
    Any = Any,
    Polar = Polar,
    PolarPct = if ("Polar.." %in% colnames(raw)) `Polar..` else if ("Polar." %in% colnames(raw)) `Polar.` else NA,
    MinPolarDist = if ("Min.polar.dist." %in% colnames(raw)) `Min.polar.dist.` else if ("Min.polar.dist" %in% colnames(raw)) `Min.polar.dist` else NA,
    Any_n = parse_ratio_numerator(Any),
    Polar_n = parse_ratio_numerator(Polar)
  ) %>%
  mutate(
    MinPolarDist_num = clean_distance(MinPolarDist),
    RecurrentPolar = Polar_n == 8
  )

# =========================================================
# ORDER RESIDUES
# Main ordering:
# 1) highest polar recurrence
# 2) highest any-contact recurrence
# 3) shortest min polar distance
# =========================================================
freq <- freq %>%
  arrange(desc(Polar_n), desc(Any_n), MinPolarDist_num)

freq <- freq %>%
  mutate(
    ResidueLabel = ifelse(RecurrentPolar, paste0(Residue, " *"), Residue)
  )

residue_order <- freq$ResidueLabel

# =========================================================
# LONG FORMAT FOR TILE PLOT
# =========================================================
plot_df <- bind_rows(
  freq %>%
    transmute(
      ResidueLabel,
      ContactType = "Any contact",
      Count = Any_n
    ),
  freq %>%
    transmute(
      ResidueLabel,
      ContactType = "Polar contact",
      Count = Polar_n
    )
) %>%
  mutate(
    ContactType = factor(ContactType, levels = c("Any contact", "Polar contact")),
    ResidueLabel = factor(ResidueLabel, levels = rev(residue_order)),
    Label = paste0(Count, "/8")
  )

# recurrent polar residues for subtitle
recurrent_residues <- freq %>%
  filter(RecurrentPolar) %>%
  pull(Residue)

subtitle_text <- paste0(
  length(recurrent_residues),
  " residues showed polar-contact recurrence in 8/8 dominant poses: ",
  paste(recurrent_residues, collapse = ", ")
)

# =========================================================
# PLOT
# =========================================================
p <- ggplot(plot_df, aes(x = ContactType, y = ResidueLabel, fill = Count)) +
  geom_tile(color = "white", linewidth = 1.0, width = 0.92, height = 0.92) +

  # outline fully recurrent polar residues
  geom_tile(
    data = subset(plot_df, ContactType == "Polar contact" & Count == 8),
    fill = NA, color = "black", linewidth = 1.2, width = 0.92, height = 0.92
  ) +

  geom_text(aes(label = Label), size = 5.3, fontface = "bold") +

  scale_fill_gradientn(
    colours = c("#F2F2F2", "#F7D98A", "#E69F00"),
    values = c(0, 0.5, 1),
    limits = c(0, 8),
    breaks = c(0, 2, 4, 6, 8),
    name = "Recurrence\nacross 8 poses"
  ) +

  labs(
    title = "AK1 dominant-cluster contact recurrence",
    subtitle = subtitle_text,
    x = NULL,
    y = "AK1 residue",
    caption = paste0(
      "* indicates residues with polar contact in all 8 dominant poses.\n",
      "Any contact = any heavy-atom contact; Polar contact = candidate H-bond/polar proximity."
    )
  ) +

  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 12)),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),
    plot.caption = element_text(size = 11, hjust = 0),
    plot.margin = margin(20, 20, 20, 20)
  )

# =========================================================
# SAVE FIGURES
# =========================================================
png_file <- file.path(out_dir, "AK1_panelG_contact_recurrence.png")
pdf_file <- file.path(out_dir, "AK1_panelG_contact_recurrence.pdf")

ggsave(
  filename = png_file,
  plot = p,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 10,
  height = 8,
  bg = "white"
)

# =========================================================
# CONSOLE SUMMARY
# =========================================================
cat("============================================================\n")
cat("AK1 PANEL G — DOMINANT 8-POSE CONTACT RECURRENCE\n")
cat("============================================================\n\n")

print(freq)

cat("\nSaved files:\n")
cat(normalizePath(png_file), "\n")
cat(normalizePath(pdf_file), "\n")