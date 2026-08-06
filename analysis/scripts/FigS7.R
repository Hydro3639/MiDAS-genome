############################################################
## Figure S7 — Trait conservation and species-level resolution
##
## Description:
##   This script visualizes trait conservation across core genera
##   and evaluates whether genus-level mixed traits are resolved
##   at the species level.
##
## Panels:
##   a) Genus-by-pathway trait conservation heatmap
##   b) Fraction of core genera classified as mixed for each pathway
##   c) Comparison of mixed fractions at genus level and among
##      species represented by at least three or five HQ MAGs
##
## Inputs:
##   - sources/FigS7a.tsv
##       Rank-resolved pathway conservation table
##
##   - sources/FigS7b.tsv
##       Rank-resolved trait-status summary table
##
##   - sources/FigS7c.tsv
##       Species-level follow-up for genus-mixed pathways,
##       restricted to species represented by at least three MAGs
##
##   - sources/FigS7d.tsv
##       Species-level follow-up for genus-mixed pathways,
##       restricted to species represented by at least five MAGs
##
## Outputs:
##   - figures/FigS7.pdf
##   - figures/FigS7.png
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ----------------------------- #
# Input and output paths
# ----------------------------- #

conservation_file <- "sources/FigS7a.tsv"
status_file       <- "sources/FigS7b.tsv"
species_min3_file <- "sources/FigS7c.tsv"
species_min5_file <- "sources/FigS7d.tsv"

output_pdf <- "figures/FigS7.pdf"

if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}

required_files <- c(
  conservation_file,
  status_file,
  species_min3_file,
  species_min5_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input file(s): ",
    paste(missing_files, collapse = ", ")
  )
}

# ----------------------------- #
# Read input data
# ----------------------------- #

conservation <- readr::read_tsv(
  conservation_file,
  show_col_types = FALSE
)

status_summary <- readr::read_tsv(
  status_file,
  show_col_types = FALSE
)

species_min3 <- readr::read_tsv(
  species_min3_file,
  show_col_types = FALSE
)

species_min5 <- readr::read_tsv(
  species_min5_file,
  show_col_types = FALSE
)

# ----------------------------- #
# Validate input columns
# ----------------------------- #

required_conservation_columns <- c(
  "Tax_rank",
  "Taxon",
  "Pathway",
  "Status"
)

required_status_columns <- c(
  "Tax_rank",
  "Pathway",
  "Status",
  "n_taxa",
  "fraction"
)

required_species_columns <- c(
  "Pathway",
  "Status"
)

missing_conservation_columns <- setdiff(
  required_conservation_columns,
  colnames(conservation)
)

missing_status_columns <- setdiff(
  required_status_columns,
  colnames(status_summary)
)

missing_species_min3_columns <- setdiff(
  required_species_columns,
  colnames(species_min3)
)

missing_species_min5_columns <- setdiff(
  required_species_columns,
  colnames(species_min5)
)

if (length(missing_conservation_columns) > 0) {
  stop(
    "Missing required columns in ",
    conservation_file,
    ": ",
    paste(missing_conservation_columns, collapse = ", ")
  )
}

if (length(missing_status_columns) > 0) {
  stop(
    "Missing required columns in ",
    status_file,
    ": ",
    paste(missing_status_columns, collapse = ", ")
  )
}

if (length(missing_species_min3_columns) > 0) {
  stop(
    "Missing required columns in ",
    species_min3_file,
    ": ",
    paste(missing_species_min3_columns, collapse = ", ")
  )
}

if (length(missing_species_min5_columns) > 0) {
  stop(
    "Missing required columns in ",
    species_min5_file,
    ": ",
    paste(missing_species_min5_columns, collapse = ", ")
  )
}

# ----------------------------- #
# Shared plot settings
# ----------------------------- #

status_levels <- c(
  "Conserved present",
  "Mixed",
  "Conserved absent"
)

status_colors <- c(
  "Conserved present" = "#2166AC",
  "Mixed"             = "#F4A582",
  "Conserved absent"  = "#D9D9D9"
)

pathway_label_map <- c(
  "Glycolysis" = "Glycolysis",
  "Glyconeogenesis" = "Gluconeogenesis",
  "Pentose_phosphate_pathway" = "Pentose phosphate pathway",
  "Pyruvate_oxidation" = "Pyruvate oxidation",
  "TCA_cycle" = "TCA cycle",
  "Glyoxylate_cycle" = "Glyoxylate cycle",
  "Fatty_acid_beta_oxidation" = "Fatty acid beta-oxidation",
  "Fatty_acid_beta-oxidation" = "Fatty acid beta-oxidation",
  "Oxidative_phosphorylation" = "Oxidative phosphorylation",
  "Glycogen_biosynthesis" = "Glycogen biosynthesis",
  "Glycogen_degradation" = "Glycogen degradation",
  "PHA_accumulation" = "PHA accumulation",
  "PHA_degradation" = "PHA degradation",
  "Calvin_cycle" = "Calvin cycle",
  "rTCA_cycle" = "rTCA cycle",
  "Ammonia_oxidation" = "Ammonia oxidation",
  "Comammox" = "Comammox",
  "Nitrite_oxidation" = "Nitrite oxidation",
  "Denitrification" = "Complete denitrification",
  "DNRA" = "DNRA",
  "Assimilatory_sulfate_reduction" =
    "Assimilatory sulfate reduction",
  "SOX" = "Sox",
  "Sox" = "Sox",
  "Oxidative_dsr" = "Oxidative Dsr",
  "Oxidative dsr" = "Oxidative Dsr"
)

preferred_pathway_order <- c(
  "Calvin cycle",
  "rTCA cycle",
  "Pyruvate oxidation",
  "TCA cycle",
  "Glycolysis",
  "Gluconeogenesis",
  "Pentose phosphate pathway",
  "Glyoxylate cycle",
  "Fatty acid beta-oxidation",
  "Oxidative phosphorylation",
  "Glycogen biosynthesis",
  "Glycogen degradation",
  "PHA accumulation",
  "PHA degradation",
  "Ammonia oxidation",
  "Nitrite oxidation",
  "Comammox",
  "Complete denitrification",
  "DNRA",
  "Assimilatory sulfate reduction",
  "Sox",
  "Oxidative Dsr"
)

clean_pathway_labels <- function(pathway) {
  dplyr::recode(
    pathway,
    !!!pathway_label_map,
    .default = pathway
  )
}

axis_theme <- theme(
  axis.title = element_text(
    color = "black",
    size = 12
  ),
  axis.text = element_text(
    color = "black",
    size = 10.5
  ),
  legend.title = element_text(
    color = "black",
    size = 11
  ),
  legend.text = element_text(
    color = "black",
    size = 10
  ),
  plot.tag = element_text(
    color = "black",
    size = 14,
    face = "bold"
  )
)

# ============================================================
# Panel a: Genus-level trait conservation heatmap
# ============================================================

figs7a_data <- conservation %>%
  filter(Tax_rank == "Genus") %>%
  mutate(
    Status = factor(
      Status,
      levels = status_levels
    ),
    Pathway_label = clean_pathway_labels(Pathway)
  )

if (nrow(figs7a_data) == 0) {
  stop(
    "No genus-level rows were found in ",
    conservation_file,
    "."
  )
}

# Retain the preferred pathway order for pathways present in the data.
pathway_levels <- preferred_pathway_order[
  preferred_pathway_order %in%
    unique(figs7a_data$Pathway_label)
]

# Append pathways not included in the predefined order.
pathway_levels <- c(
  pathway_levels,
  setdiff(
    unique(figs7a_data$Pathway_label),
    pathway_levels
  )
)

# Order genera by the number of conserved-present traits,
# followed by the number of mixed traits.
genus_order <- figs7a_data %>%
  group_by(Taxon) %>%
  summarise(
    n_present = sum(
      Status == "Conserved present",
      na.rm = TRUE
    ),
    n_mixed = sum(
      Status == "Mixed",
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_present),
    desc(n_mixed),
    Taxon
  ) %>%
  pull(Taxon)

figs7a_data <- figs7a_data %>%
  mutate(
    Pathway_label = factor(
      Pathway_label,
      levels = pathway_levels
    ),
    Taxon = factor(
      Taxon,
      levels = rev(genus_order)
    )
  )

p_a <- ggplot(
  figs7a_data,
  aes(
    x = Pathway_label,
    y = Taxon,
    fill = Status
  )
) +
  geom_tile(color = NA) +
  scale_fill_manual(
    values = status_colors,
    breaks = status_levels,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Trait status"
  ) +
  theme_bw(base_size = 8) +
  axis_theme +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 10.5
    ),
    legend.position = "top",
    plot.margin = margin(5, 5, 5, 5)
  )

p_a

# ============================================================
# Optional alternative Panel a: Stacked genus fractions
# this figure was not included in this version
# ============================================================

p_a_alt <- status_summary %>%
  filter(Tax_rank == "Genus") %>%
  mutate(
    Status = factor(
      Status,
      levels = status_levels
    ),
    Pathway_label = clean_pathway_labels(Pathway),
    Pathway_label = factor(
      Pathway_label,
      levels = pathway_levels
    )
  ) %>%
  ggplot(
    aes(
      x = Pathway_label,
      y = fraction,
      fill = Status
    )
  ) +
  geom_col(width = 0.75) +
  scale_fill_manual(
    values = status_colors,
    breaks = status_levels,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      scale = 100,
      accuracy = 1,
      suffix = "%"
    )
  ) +
  labs(
    x = NULL,
    y = "Fraction of 277 core genera",
    fill = "Trait status"
  ) +
  theme_bw(base_size = 10) +
  axis_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 10.5
    ),
    axis.text.y = element_text(
      color = "black",
      size = 10.5
    ),
    axis.title.y = element_text(
      color = "black",
      size = 12
    ),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

p_a_alt

# ============================================================
# Panel b: Genus-level mixed fraction by pathway
# ============================================================

genus_mixed_rows <- status_summary %>%
  filter(
    Tax_rank == "Genus",
    Status == "Mixed"
  ) %>%
  mutate(
    Pathway_label = clean_pathway_labels(Pathway)
  )

# Include pathways with no mixed genera as zero-valued rows.
all_genus_pathways <- status_summary %>%
  filter(Tax_rank == "Genus") %>%
  distinct(Pathway) %>%
  mutate(
    Pathway_label = clean_pathway_labels(Pathway)
  )

figs7b_data <- all_genus_pathways %>%
  left_join(
    genus_mixed_rows %>%
      select(
        Pathway,
        mixed_fraction = fraction,
        n_mixed = n_taxa
      ),
    by = "Pathway"
  ) %>%
  mutate(
    mixed_fraction = replace_na(
      mixed_fraction,
      0
    ),
    n_mixed = replace_na(
      n_mixed,
      0
    )
  ) %>%
  arrange(mixed_fraction) %>%
  mutate(
    Pathway_label = factor(
      Pathway_label,
      levels = unique(Pathway_label)
    )
  )

# Save the shared pathway order for Panels b and c.
pathway_order_bc <- levels(
  figs7b_data$Pathway_label
)

bar_color <- "#F4A582"

level_fill_colors <- c(
  "Genus"              = "#3C5488FF",
  "Species ≥3 HQ MAGs" = "#00A087FF",
  "Species ≥5 HQ MAGs" = "#E64B35FF"
)

p_b <- ggplot(
  figs7b_data,
  aes(
    x = mixed_fraction,
    y = Pathway_label
  )
) +
  geom_col(
    width = 0.75,
    fill = bar_color
  ) +
  geom_text(
    aes(label = n_mixed),
    hjust = -0.5,
    size = 3.2,
    color = "black"
  ) +
  scale_x_continuous(
    labels = scales::label_number(
      scale = 100,
      accuracy = 1,
      suffix = "%"
    ),
    expand = expansion(
      mult = c(0, 0.15)
    )
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "Fraction of mixed genera",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  axis_theme +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

p_b

# ============================================================
# Panel c: Species-level follow-up of genus-mixed pathways
# ============================================================

calculate_species_mixed_fraction <- function(
    data,
    level_label
) {
  data %>%
    group_by(
      Pathway,
      Status
    ) %>%
    summarise(
      n_taxa = n(),
      .groups = "drop"
    ) %>%
    group_by(Pathway) %>%
    mutate(
      total_taxa = sum(n_taxa),
      fraction = n_taxa / total_taxa
    ) %>%
    ungroup() %>%
    filter(Status == "Mixed") %>%
    transmute(
      Pathway,
      Level = level_label,
      mixed_fraction = fraction,
      n_mixed = n_taxa,
      total_taxa
    )
}

genus_mixed <- status_summary %>%
  filter(
    Tax_rank == "Genus",
    Status == "Mixed"
  ) %>%
  transmute(
    Pathway,
    Level = "Genus",
    mixed_fraction = fraction,
    n_mixed = n_taxa,
    total_taxa
  )

species_mixed_min3 <- calculate_species_mixed_fraction(
  species_min3,
  "Species ≥3 HQ MAGs"
)

species_mixed_min5 <- calculate_species_mixed_fraction(
  species_min5,
  "Species ≥5 HQ MAGs"
)

figs7c_data <- bind_rows(
  genus_mixed,
  species_mixed_min3,
  species_mixed_min5
) %>%
  mutate(
    Pathway_label = clean_pathway_labels(Pathway),
    Level = factor(
      Level,
      levels = c(
        "Genus",
        "Species ≥3 HQ MAGs",
        "Species ≥5 HQ MAGs"
      )
    )
  )

# Retain only pathways classified as mixed at genus level.
pathways_mixed_at_genus <- genus_mixed %>%
  filter(mixed_fraction > 0) %>%
  pull(Pathway)

complete_species_grid <- tidyr::expand_grid(
  Pathway = pathways_mixed_at_genus,
  Level = c(
    "Genus",
    "Species ≥3 HQ MAGs",
    "Species ≥5 HQ MAGs"
  )
)

figs7c_data <- complete_species_grid %>%
  left_join(
    figs7c_data,
    by = c(
      "Pathway",
      "Level"
    )
  ) %>%
  mutate(
    mixed_fraction = replace_na(
      mixed_fraction,
      0
    ),
    n_mixed = replace_na(
      n_mixed,
      0
    ),
    Pathway_label = clean_pathway_labels(Pathway),
    Level = factor(
      Level,
      levels = c(
        "Genus",
        "Species ≥3 HQ MAGs",
        "Species ≥5 HQ MAGs"
      )
    ),
    Pathway_label = factor(
      Pathway_label,
      levels = pathway_order_bc
    )
  )

p_c <- ggplot(
  figs7c_data,
  aes(
    x = mixed_fraction,
    y = Pathway_label
  )
) +
  geom_point(
    aes(
      shape = Level,
      fill = Level
    ),
    size = 2.8,
    stroke = 0.45,
    color = "grey15"
  ) +
  scale_shape_manual(
    values = c(
      "Genus"              = 21,
      "Species ≥3 HQ MAGs" = 24,
      "Species ≥5 HQ MAGs" = 22
    )
  ) +
  scale_fill_manual(
    values = level_fill_colors
  ) +
  scale_x_continuous(
    labels = scales::label_number(
      scale = 100,
      accuracy = 1,
      suffix = "%"
    ),
    expand = expansion(
      mult = c(0, 0.1)
    )
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "Mixed fraction",
    y = NULL,
    shape = NULL,
    fill = NULL
  ) +
  guides(
    shape = guide_legend(nrow = 1),
    fill = guide_legend(nrow = 1)
  ) +
  theme_bw(base_size = 9) +
  axis_theme +
  theme(
    panel.grid.major.y = element_line(
      linewidth = 0.2,
      color = "grey90"
    ),
    panel.grid.major.x = element_line(
      linewidth = 0.2,
      color = "grey90"
    ),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    plot.margin = margin(5, 5, 5, 5)
  )

p_c

# ============================================================
# Combine panels
# ============================================================

# Hide the pathway labels in Panel c while preserving its
# discrete y-axis structure and alignment with Panel b.
p_c_without_y_labels <- p_c +
  theme(
    axis.title.x = element_text(
      color = "black",
      size = 12
    ),
    axis.text.x = element_text(
      color = "black",
      size = 10.5
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(
      color = "black",
      size = 10
    ),
    legend.title = element_text(
      color = "black",
      size = 11
    ),
    plot.margin = margin(5, 5, 5, 5)
  )

p_b_aligned <- p_b +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

figs7 <- (
  patchwork::free(
    p_a,
    side = "l"
  ) /
    (
      p_b_aligned |
        p_c_without_y_labels
    )
) +
  plot_layout(
    heights = c(1, 1.2),
    widths = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "a"
  )

figs7

# ----------------------------- #
# Export
# ----------------------------- #

ggsave(
  filename = output_pdf,
  plot = figs7,
  width = 12,
  height = 12
)

# ----------------------------- #
# Record the R environment
# ----------------------------- #

if (!dir.exists("../envs")) {
  dir.create(
    "../envs",
    recursive = TRUE
  )
}

writeLines(
  capture.output(sessionInfo()),
  "../envs/FigS7.sessionInfo.txt"
)
