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
##   c) Corrected species mixed fractions versus within-genus null means,
##      at three or five HQ MAGs per species (specific mixed parent pairs)
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
##   - sources/FigS7_null_by_pathway.tsv (completed null summaries)
##
## Outputs:
##   - figures/FigS7.pdf
##   # - figures/FigS7.png (600 dpi; export is commented out)
##   - ../envs/FigS7.sessionInfo.txt
## This script reads prepared sources and generates the figure only
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

# the input files
conservation_file <- "sources/FigS7a.tsv"
status_file       <- "sources/FigS7b.tsv"
species_min3_file <- "sources/FigS7c.tsv"
species_min5_file <- "sources/FigS7d.tsv"
null_file         <- "sources/FigS7_null_by_pathway.tsv"

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("../envs", recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  conservation_file,
  status_file,
  species_min3_file,
  species_min5_file,
  null_file
)


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
species_null <- readr::read_tsv(null_file, show_col_types = FALSE)

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
  "Taxon", "Genus", "Pathway", "Status", "n_MAGs"
)

# Validate that the prepared species and null tables use the same eligible groups.
parent_pairs <- conservation %>% filter(Tax_rank == "Genus", Status == "Mixed") %>%
  transmute(Genus = Taxon, Pathway) %>% distinct()
for (cutoff in c(3L, 5L)) {
  x <- if (cutoff == 3L) species_min3 else species_min5
  if (nrow(anti_join(x, parent_pairs, by = c("Genus", "Pathway"))) > 0L) {
    stop("Species follow-up contains rows outside specific mixed genus-pathway pairs.")
  }
  stopifnot(all(x$n_MAGs >= cutoff),
    !anyDuplicated(x[c("Taxon", "Genus", "Pathway")]),
    all(x$Status %in% c("Conserved present", "Mixed", "Conserved absent")))
}
needed_null <- c("minimum_MAGs", "Pathway", "n_combinations", "n_mixed", "observed", "null_mean")

# ----------------------------- #
# Shared plot settings
# ----------------------------- #

figure_width <- 12
figure_height <- 12
#png_dpi <- 600

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
  "Gluconeogenesis" = "Gluconeogenesis",
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
    face = "plain"
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

if (interactive()) print(p_a)

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
  scale_y_discrete(drop = FALSE, expand = expansion(add = 0.6)) +
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

if (interactive()) print(p_b)

# ============================================================
# Panel c: Corrected species follow-up versus matched null expectations
# ============================================================
# Compare observed and null species fractions at the same sample sizes.
# Four pathways have no eligible species follow-up. Leave them blank, not at 0%.

figs7c_summary <- species_null %>% mutate(
  Pathway_label = factor(clean_pathway_labels(Pathway), levels = pathway_order_bc),
  Level = if_else(minimum_MAGs == 3L, "Species ≥3 HQ MAGs", "Species ≥5 HQ MAGs"),
  y = as.numeric(Pathway_label) + if_else(minimum_MAGs == 3L, 0.13, -0.13))
figs7c_data <- bind_rows(
  transmute(figs7c_summary, Pathway, Pathway_label, Level, y,
    Comparison = "Observed", mixed_fraction = observed),
  transmute(figs7c_summary, Pathway, Pathway_label, Level, y,
    Comparison = "Randomized null", mixed_fraction = null_mean)
) %>% mutate(Key = paste(Comparison, Level, sep = ": "))
key_levels <- c("Observed: Species ≥3 HQ MAGs", "Randomized null: Species ≥3 HQ MAGs",
                "Observed: Species ≥5 HQ MAGs", "Randomized null: Species ≥5 HQ MAGs")
key_labels <- c("Observed ≥3 HQ MAGs", "Null ≥3 HQ MAGs", "Observed ≥5 HQ MAGs", "Null ≥5 HQ MAGs")
figs7c_data$Key <- factor(figs7c_data$Key, levels = key_levels)
key_colors <- setNames(c("#00A087FF", "#00A087FF", "#E64B35FF", "#E64B35FF"), key_levels)
key_fills <- setNames(c("#00A087FF", "white", "#E64B35FF", "white"), key_levels)
key_shapes <- setNames(c(24, 24, 22, 22), key_levels)
p_c <- ggplot(figs7c_data, aes(x = mixed_fraction, y = y)) +
  geom_segment(data = figs7c_summary,
    aes(x = observed, xend = null_mean, y = y, yend = y, colour = Level),
    inherit.aes = FALSE, linewidth = 0.35, alpha = 0.6, show.legend = FALSE) +
  geom_point(aes(shape = Key, fill = Key), size = 2.8, stroke = 0.5, colour = "grey15") +
  scale_shape_manual(values = key_shapes, breaks = key_levels, labels = key_labels) +
  scale_fill_manual(values = key_fills, breaks = key_levels, labels = key_labels) +
  scale_colour_manual(values = level_fill_colors, guide = "none") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1),
    limits = c(-0.02, 1.04), breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq_along(pathway_order_bc), labels = pathway_order_bc,
    limits = c(0.4, length(pathway_order_bc) + 0.6), expand = expansion(mult = 0)) +
  labs(x = "Mixed species–pathway fraction", y = NULL, shape = NULL, fill = NULL) +
  guides(shape = guide_legend(nrow = 2, byrow = TRUE),
         fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw(base_size = 9) + axis_theme +
  theme(panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
    panel.grid.major.x = element_line(linewidth = 0.2, color = "grey90"),
    panel.grid.minor = element_blank(), legend.position = "top",
    legend.text = element_text(size = 8), legend.box = "horizontal",
    plot.margin = margin(5, 5, 5, 5))
if (interactive()) print(p_c)

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

if (interactive()) print(figs7)

# ----------------------------- #
# Export: PDF by default; optional PNG is commented out
# ----------------------------- #

ggsave("figures/FigS7.pdf", figs7,
  width = figure_width, height = figure_height, units = "in", bg = "white")

#ggsave("figures/FigS7.png", figs7, device = grDevices::png,
#  width = figure_width, height = figure_height, units = "in",
#  dpi = png_dpi, type = "cairo", bg = "white")

writeLines(capture.output(sessionInfo()), "../envs/FigS7.sessionInfo.txt")
