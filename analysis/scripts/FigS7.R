# ============================================================
# Fig. 4 | Genus-level taxonomy provides a pathway-dependent
# framework for functional interpretation
#
# a: Genus × trait conservation heatmap
# b: Mixed fraction across taxonomic ranks
# c: Species-level follow-up for genus-mixed pathways
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)
library(scales)
library(patchwork)

# -----------------------------
# 1. Input files
# -----------------------------

conservation_file <- "../rank_resolved_pathway_conservation.tsv"
status_file <- "../rank_resolved_status_summary.tsv"
species_min3_file <- "../species_followup_for_genus_mixed_pathways_min3MAGs.tsv"
species_min5_file <- "../species_followup_for_genus_mixed_pathways_min5MAGs.tsv"

conservation <- read_tsv(conservation_file, show_col_types = FALSE)
status_summary <- read_tsv(status_file, show_col_types = FALSE)
species_min3 <- read_tsv(species_min3_file, show_col_types = FALSE)
species_min5 <- read_tsv(species_min5_file, show_col_types = FALSE)

# -----------------------------
# 2. Settings
# -----------------------------

rank_levels <- c("Genus", "Family", "Order", "Class", "Phylum")

status_levels <- c(
  "Conserved present",
  "Mixed",
  "Conserved absent"
)

status_cols <- c(
  "Conserved present" = "#2166AC",
  "Mixed" = "#F4A582",
  "Conserved absent" = "#D9D9D9"
)

# Adjust labels if your pathway names differ
pathway_label_map <- c(
  "Glycolysis" = "Glycolysis",
  "Glyconeogenesis" = "Glyconeogenesis",
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
  "Assimilatory_sulfate_reduction" = "Assimilatory sulfate reduction",
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
  "Glyconeogenesis",
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


clean_pathway_labels <- function(x) {
  recode(x, !!!pathway_label_map, .default = x)
}
## axis titles
axis_black_big_theme <- theme(
  axis.title = element_text(color = "black", size = 12),
  axis.text  = element_text(color = "black", size = 10.5),
  legend.title = element_text(color = "black", size = 11),
  legend.text  = element_text(color = "black", size = 10),
  plot.tag = element_text(color = "black", size = 14, face = "bold")
)


# ============================================================
# Fig. Sx: Overall mixed fraction across taxonomic ranks
# ============================================================

# ============================================================
# Rank-resolved mixed fraction distribution by pathway
# Violin + boxplot + points + number of mixed combinations
# ============================================================

library(ggrepel)

fig4sx_dist_df <- status_summary %>%
  filter(Tax_rank %in% rank_levels) %>%
  mutate(
    Tax_rank = factor(Tax_rank, levels = rank_levels),
    Pathway_label = clean_pathway_labels(Pathway)
  ) %>%
  group_by(Tax_rank, Pathway, Pathway_label) %>%
  summarise(
    n_mixed = sum(n_taxa[Status == "Mixed"], na.rm = TRUE),
    total_combinations = sum(n_taxa, na.rm = TRUE),
    mixed_fraction = n_mixed / total_combinations,
    .groups = "drop"
  )

# Overall mixed fraction per taxonomic rank
fig4sx_overall_df <- status_summary %>%
  filter(Tax_rank %in% rank_levels) %>%
  group_by(Tax_rank) %>%
  summarise(
    mixed_combinations = sum(n_taxa[Status == "Mixed"], na.rm = TRUE),
    total_combinations = sum(n_taxa, na.rm = TRUE),
    mixed_fraction = mixed_combinations / total_combinations,
    .groups = "drop"
  ) %>%
  mutate(
    Tax_rank = factor(Tax_rank, levels = rank_levels),
    label = paste0(
      scales::number(mixed_fraction * 100, accuracy = 0.1),
      "\n",
      "n = ", mixed_combinations
    )
  )

p_Sx_violin_line <- ggplot(
  fig4sx_dist_df,
  aes(x = Tax_rank, y = mixed_fraction)
) +
  geom_violin(
    width = 0.85,
    trim = FALSE,
    fill = "#4DBBD5FF",
    alpha = 0.35,
    color = NA
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.4
  ) +
  geom_point(
    aes(size = n_mixed),
    shape = 21,
    fill = "#3C5488FF",
    color = "black",
    stroke = 0.35,
    alpha = 0.75,
    position = position_jitter(
      width = 0.08,
      height = 0,
      seed = 1
    )
  ) +
  geom_line(
    data = fig4sx_overall_df,
    aes(x = Tax_rank, y = mixed_fraction, group = 1),
    inherit.aes = FALSE,
    linewidth = 0.6,
    color = "black"
  ) +
  geom_point(
    data = fig4sx_overall_df,
    aes(x = Tax_rank, y = mixed_fraction),
    inherit.aes = FALSE,
    size = 2.8,
    color = "black"
  ) +
  geom_text(
    data = fig4sx_overall_df,
    aes(
      x = Tax_rank,
      y = mixed_fraction,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = -0.9,
    size = 3.2,
    color = "black",
    lineheight = 0.85
  ) +
  scale_size_continuous(
    name = "Mixed combinations",
    range = c(1.5, 4.5)
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      scale = 100,
      accuracy = 1
    ),
    limits = c(
      0,
      max(fig4sx_dist_df$mixed_fraction, fig4sx_overall_df$mixed_fraction, na.rm = TRUE) * 1.25
    ),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    x = "Taxonomic rank",
    y = "Mixed fraction across pathways (%)"
  ) +
  theme_bw(base_size = 11) +
  axis_black_big_theme +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black",
      size = 10.5
    ),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

p_Sx_violin_line


ggsave(
  "FigSx_Overall_mixed_fraction_across_taxonomic_ranks.pdf",
  p_Sx,
  width = 12,
  height = 12
)

ggsave(
  "FigSx_Overall_mixed_fraction_across_taxonomic_ranks.png",
  p_Sx,
  width = 12,
  height = 12,
  dpi = 400
)
