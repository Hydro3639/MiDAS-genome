
############################################################
## Fig. 1d — rRNA and BGC recovery in MiDAS and GWASMC
##
## Input:
##   - sources/Fig1d_rRNA_whole_catalog.tsv
##     Required columns:
##       Domain_subset, Group, Feature, Positive, Total, Percentage
##   - sources/Fig1d_BGC_whole_catalog_genomes.tsv
##     Required columns:
##       genome, Catalog, Quality, Domain, region_count,
##       complete_region_count
##   - sources/Fig1d_BGC_whole_catalog_summary.tsv
##     Required columns:
##       Domain_subset, Group, n_MAGs, total_regions,
##       total_complete_regions, region_count_median,
##       region_count_q25, region_count_q75, complete_count_median,
##       complete_count_q25, complete_count_q75
##
## Output:
##   - figures/Fig1d.pdf
##   ##- figures/Fig1d.png
##   ##- figures/Fig1d.rRNA_plot_data.tsv
##   ##- figures/Fig1d.BGC_summary.tsv
##   - ../envs/Fig1d.sessionInfo.txt
##
## Notes:
##   - HQ and MQ groups follow each study's original definitions.
##   - so, for GWASMC HQ MAGs, RNA genes were not required.
##   - Row order: MiDAS-HQ, MiDAS-MQ, GWASMC-HQ, GWASMC-MQ.
##   - Complete predicted BGC regions have contig_edge=False.
############################################################

library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)

# ----------------------------- #
# Input checks and output folders
# ----------------------------- #

input_files <- c(
  "sources/Fig1d_rRNA_whole_catalog.tsv",
  "sources/Fig1d_BGC_whole_catalog_genomes.tsv",
  "sources/Fig1d_BGC_whole_catalog_summary.tsv"
)
## check if we are in the right path
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Missing input files: ", paste(missing_files, collapse = ", "))
}
dir.create("figures", showWarnings = FALSE)
dir.create("../envs", showWarnings = FALSE)

# ----------------------------- #
# Shared settings
# ----------------------------- #

domain_subset <- "All"  # Change to "Bacteria" for bacteria only
panel_tag <- "d"        # One label for the complete five-plot row

group_levels <- c(
  "MiDAS-HQ", "MiDAS-MQ",
  "GWASMC-HQ", "GWASMC-MQ"
)

group_colors <- c(
  "MiDAS-HQ"   = "#ed665dFF",
  "GWASMC-HQ"  = "#59a14f",
  "MiDAS-MQ"   = "#f4b3ad",
  "GWASMC-MQ"  = "#b7d9a8"
)

figure_width <- 11.5
figure_height <- 3
base_text_size <- 12
value_text_size <- 3.5
median_y_offset <- 0.47  # Increase the number to move BGC median labels higher.
panel_widths <- c(1, 1, 1, 1.4, 1.4)

theme_base <- function(show_y_text = FALSE) {
  theme_classic(base_size = base_text_size, base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, margin = margin(b = 1)),
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      axis.line.x = element_line(linewidth = 0.4, colour = "black"),
      axis.line.y = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks.y = element_blank(),
      axis.ticks.length = unit(0.12, "cm"),
      axis.text = element_text(colour = "black"),
      axis.text.y = if (show_y_text) element_text(colour = "black") else element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(margin = margin(t = 6)),
      plot.margin = margin(7, 8, 5, 5),
      plot.tag = element_text(size = 11, face = "bold"),
      plot.tag.position = "topleft"
    )
}

# The same scale and expansion align category rows in all five plots.
group_axis <- function() {
  scale_y_discrete(
    limits = rev(group_levels),
    labels = identity,
    expand = expansion(add = c(0.85, 0.60)),
    name = NULL
  )
}

# ----------------------------- #
# Load and validate rRNA data
# ----------------------------- #

rrna_raw <- read.delim("sources/Fig1d_rRNA_whole_catalog.tsv", stringsAsFactors = FALSE, check.names = FALSE)
rrna_required <- c("Domain_subset", "Group", "Feature", "Positive", "Total", "Percentage")
feature_levels <- c("16S rRNA", "23S rRNA", "Intact operon")

rrna_data <- rrna_raw %>%
  filter(
    Domain_subset == domain_subset,
    Group %in% group_levels,
    Feature %in% feature_levels
  ) %>%
  mutate(
    Group = factor(Group, levels = rev(group_levels)),
    Feature = factor(Feature, levels = feature_levels),
    percentage_label = sprintf("%.1f%%", Percentage),
    label_x = ifelse(Percentage >= 85, Percentage - 1.6, Percentage + 1.5),
    label_hjust = ifelse(Percentage >= 85, 1, 0)
  )


rrna_totals <- rrna_data %>%
  group_by(Group) %>%
  summarise(Total = first(Total), n_totals = n_distinct(Total), .groups = "drop")

# ----------------------------- #
# Load and validate BGC data
# ----------------------------- #

bgc_raw <- read.delim(
  "sources/Fig1d_BGC_whole_catalog_genomes.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
)
summary_raw <- read.delim(
  "sources/Fig1d_BGC_whole_catalog_summary.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
)
bgc_required <- c("genome", "Catalog", "Quality", "Domain", "region_count", "complete_region_count")
summary_required <- c(
  "Domain_subset", "Group", "n_MAGs", "total_regions", "total_complete_regions",
  "region_count_median", "region_count_q25", "region_count_q75",
  "complete_count_median", "complete_count_q25", "complete_count_q75"
)

bgc_selected <- if (domain_subset == "All") bgc_raw else filter(bgc_raw, Domain == domain_subset)

# Keep zero-BGC MAGs. Each MAG contributes once, to its original quality group.
# Pooled HQ+MQ rows are not used here to make it similar with Fig1c
bgc_data <- bgc_selected %>%
  mutate(Group = paste(Catalog, Quality, sep = "-")) %>%
  filter(Group %in% group_levels) %>%
  mutate(Group = factor(Group, levels = rev(group_levels)))

bgc_summary <- bgc_data %>%
  group_by(Group) %>%
  summarise(
    n_MAGs = n(),
    total_regions = sum(region_count),
    total_complete_regions = sum(complete_region_count),
    region_count_median = median(region_count),
    region_count_q25 = quantile(region_count, 0.25, names = FALSE),
    region_count_q75 = quantile(region_count, 0.75, names = FALSE),
    complete_count_median = median(complete_region_count),
    complete_count_q25 = quantile(complete_region_count, 0.25, names = FALSE),
    complete_count_q75 = quantile(complete_region_count, 0.75, names = FALSE),
    .groups = "drop"
  )

# Recompute BGC summaries in R and compare with the supplied summaries.
reference <- filter(summary_raw, Domain_subset == domain_subset)
for (i in seq_len(nrow(bgc_summary))) {
  row <- bgc_summary[i, ]
  expected <- reference[reference$Group == as.character(row$Group), ]
  if (nrow(expected) != 1) stop("Missing/duplicate BGC reference summary row.")
  for (key in setdiff(names(bgc_summary), "Group")) {
    if (!isTRUE(all.equal(as.numeric(row[[key]]), as.numeric(expected[[key]]), tolerance = 1e-9))) {
      stop("R/source summary disagreement: ", as.character(row$Group), ": ", key)
    }
  }
}

# Both analyses must use the same numbers of MAGs in every comparison group.
denominator_check <- left_join(bgc_summary, rrna_totals, by = "Group")
if (anyNA(denominator_check$Total) || any(denominator_check$n_MAGs != denominator_check$Total)) {
  stop("rRNA and BGC comparison-group denominators differ.")
}

# ----------------------------- #
# Three rRNA bar plots
# ----------------------------- #

make_rrna_panel <- function(feature, title, show_y_text = FALSE) {
  ggplot(filter(rrna_data, Feature == feature), aes(x = Percentage, y = Group, fill = Group)) +
    geom_col(width = 0.66, colour = NA, orientation = "y") +
    geom_text(
      aes(x = label_x, label = percentage_label, hjust = label_hjust),
      size = value_text_size, colour = "black"
    ) +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_x_continuous(
      limits = c(0, 105), breaks = c(0, 25, 50, 75, 100),
      expand = expansion(mult = 0), name = "MAGs (%)"
    ) +
    group_axis() +
    coord_cartesian(clip = "off") +  # Show complete percentage labels beyond the panel edge.
    labs(title = title) +
    theme_base(show_y_text) +
    theme(plot.margin = margin(7, 16, 5, 5))
}

# ----------------------------- #
# Two BGC violin/box plots
# ----------------------------- #

# Complete BGCs: contig_edge=False, irrespective of the MAG's HQ/MQ status.
# Both BGC plots use the same square-root count axis.
count_limit <- max(bgc_data$region_count)
count_breaks <- c(0, 1, 5, 10, 20, 40, 80, 160, 320)
count_breaks <- count_breaks[count_breaks <= count_limit]

make_bgc_panel <- function(metric, median_column, title) {
  count_labels <- bgc_summary %>%
    mutate(
      median_value = .data[[median_column]],
      label = format(median_value, trim = TRUE, scientific = FALSE)
    )

  ggplot(bgc_data, aes(x = .data[[metric]], y = Group, fill = Group)) +
    geom_violin(
      width = 0.72, trim = TRUE, scale = "width", orientation = "y", linewidth = 0.3
    ) +
    geom_boxplot(
      width = 0.12, outlier.shape = NA, orientation = "y",
      fill = "white", colour = "black", linewidth = 0.3
    ) +
    geom_text(
      data = count_labels,
      aes(x = median_value, y = Group, label = label), inherit.aes = FALSE,
      position = position_nudge(y = median_y_offset), size = value_text_size, colour = "black"
    ) +
    scale_fill_manual(values = group_colors, guide = "none") +
    group_axis() +
    scale_x_sqrt(
      limits = c(0, count_limit), breaks = count_breaks,
      expand = expansion(mult = c(0.06, 0.06)), name = "BGC regions per MAG"
    ) +
    labs(title = title) +
    theme_base()
}

# ----------------------------- #
# Assemble the single row (panel d)
# ----------------------------- #

#p16S <- make_rrna_panel("16S rRNA", "16S rRNA", show_y_text = TRUE) + labs(tag = panel_tag)
p16S <- make_rrna_panel("16S rRNA", "16S rRNA", show_y_text = TRUE)
p23S <- make_rrna_panel("23S rRNA", "23S rRNA")
p_operon <- make_rrna_panel("Intact operon", "Intact operon")
p_total <- make_bgc_panel("region_count", "region_count_median", "Total BGCs")
p_complete <- make_bgc_panel("complete_region_count", "complete_count_median", "Complete BGCs")

p_combined <- (p16S | p23S | p_operon | p_total | p_complete) +
  plot_layout(nrow = 1, widths = panel_widths)

if (interactive()) print(p_combined)

# ----------------------------- #
# Export
# ----------------------------- #

ggsave(
  p_combined, filename = "figures/Fig1d.pdf", width = figure_width, height = figure_height,
  units = "in", bg = "white"
)
#ggsave(
#  p_combined, filename = "figures/Fig1d.png", width = figure_width, height = figure_height,
#  units = "in", dpi = 600, bg = "white"
#)
#write.table(
#  rrna_data, "figures/Fig1d.rRNA_plot_data.tsv",
#  sep = "\t", quote = FALSE, row.names = FALSE
#)
#write.table(
#  bgc_summary, "figures/Fig1d.BGC_summary.tsv",
#  sep = "\t", quote = FALSE, row.names = FALSE
#)
writeLines(capture.output(sessionInfo()), "../envs/Fig1d.sessionInfo.txt")
