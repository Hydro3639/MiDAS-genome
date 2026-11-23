############################################################
## Fig. 1c: Distribution of genome quality metrics
##
## Input:
##   - sources/Fig1c-I.tsv
##     Required columns:
##       bin, checkmx_completeness, checkmx_contamination,
##       Contig_N50, contigs, quality, source
##   - sources/Fig1c-II.tsv
##     Required columns:
##       bin, Cov, quality, source
##   - sources/Fig1c-gwasmc_rep.checkm1.checkm2.stats.tsv
##     Required columns:
##       bin, checkm2_completeness, checkm2_contamination,
##       checkm1_completeness, checkm1_contamination,
##       n_contigs, ctg_L50
##
## Output:
##   - figures/Fig1c.pdf
##   ##- figures/Fig1c.png
##   ##- figures/Fig1c.plot_data.tsv
##   ##- figures/Fig1c.GWASMC_selected_scores.tsv
##   - ../envs/Fig1c.sessionInfo.txt
############################################################

library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)

# ----------------------------- #
# Shared settings
# ----------------------------- #

random_seed <- 42
points_per_group <- 1000
figure_width <- 28
figure_height <- 8
base_text_size <- 30

# Top-to-bottom order; MiDAS is between the two reference catalogues.
group_levels <- c(
  "GTDB-HQ", "GTDB-MQ",
  "MiDAS-HQ", "MiDAS-MQ",
  "GWASMC-HQ", "GWASMC-MQ"
)
group_positions <- setNames(c(1.9, 1.6, 1.2, 0.9, 0.5, 0.2), group_levels)

group_colors <- c(
  "GTDB-HQ"   = "#5a93b2",
  "GTDB-MQ"   = "#a6cee3FF",
  "MiDAS-HQ"  = "#ed665dFF",
  "MiDAS-MQ"  = "#f4b3ad",
  "GWASMC-HQ" = "#59a14f",
  "GWASMC-MQ" = "#b7d9a8"
)

dir.create("figures", showWarnings = FALSE)
dir.create("../envs", showWarnings = FALSE)

theme_base <- function(show_y_text = TRUE) {
  theme_classic(base_size = base_text_size) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.line.x = element_line(linewidth = 0.8, colour = "black"),
      axis.line.y = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks.y = element_blank(),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text = element_text(colour = "black"),
      axis.text.y = if (show_y_text) element_text(colour = "black") else element_blank(),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_blank(),
      plot.margin = margin(10, 15, 10, 10)
    )
}

# The same y scale also places coverage beside the MiDAS rows.
group_axis <- function() {
  scale_y_continuous(
    breaks = rev(unname(group_positions)),
    labels = rev(group_levels),
    limits = c(0, 2.1),
    expand = expansion(add = 0),
    name = NULL
  )
}

add_box_and_points <- function(p) {
  p +
    geom_boxplot(width = 0.25, outlier.shape = NA, orientation = "y") +
    geom_point(
      position = position_jitter(width = 0, height = 0.08, seed = random_seed),
      size = 1, alpha = 0.1, colour = "black"
    ) +
    group_axis() +
    scale_fill_manual(values = group_colors, guide = "none")
}

# ----------------------------- #
# Load the original MiDAS and GTDB plotting data
# the data were prepared on the server using shuf command
# ----------------------------- #

original_data <- read.delim(
  "sources/Fig1c-I.tsv",
  stringsAsFactors = FALSE,
  na.strings = c("NA", "None")
) %>%
  dplyr::rename(
    checkmx_completeness  = checkm2_completeness, # for the combination, not  real checkm2 matrix
    checkmx_contamination = checkm2_contamination
  )
coverage_data <- read.delim("sources/Fig1c-II.tsv", stringsAsFactors = FALSE)

# ----------------------------- #
# Select paired CheckM/CheckM2 scores for GWASMC
# The Same rule applied for GTDB/MiDAS to have a fair comparison
# ----------------------------- #

# Prioritise an HQ assessment, then an MQ assessment. If both tools
# qualify at the same tier, maximise completeness - 5 * contamination.
gwasmc_scores <- read.delim(
  "sources/Fig1c-gwasmc_rep.checkm1.checkm2.stats.tsv",
  stringsAsFactors = FALSE
) %>%
  mutate(
    hq2 = checkm2_completeness >= 90 & checkm2_contamination <= 5,
    hq1 = checkm1_completeness >= 90 & checkm1_contamination <= 5,
    mq2 = checkm2_completeness >= 50 & checkm2_contamination <= 10,
    mq1 = checkm1_completeness >= 50 & checkm1_contamination <= 10,
    use_checkm2 = case_when(
      hq2 & !hq1 ~ TRUE,
      hq1 & !hq2 ~ FALSE,
      mq2 & !mq1 ~ TRUE,
      mq1 & !mq2 ~ FALSE,
      TRUE ~ checkm2_completeness - 5 * checkm2_contamination >=
        checkm1_completeness - 5 * checkm1_contamination
    ),
    checkmx_completeness = if_else(
      use_checkm2, checkm2_completeness, checkm1_completeness
    ),
    checkmx_contamination = if_else(
      use_checkm2, checkm2_contamination, checkm1_contamination
    ),
    selected_assessment = if_else(use_checkm2, "CheckM2", "CheckM1"),
    quality = if_else(hq1 | hq2, "HQ", "MQ"),
    source = "GWASMC",
    Contig_N50 = ctg_L50,
    contigs = n_contigs
  ) %>%
  filter(mq1 | mq2) %>%
  select(
    bin, checkmx_completeness, checkmx_contamination, Contig_N50,
    contigs, quality, source, selected_assessment,
    checkm2_completeness, checkm2_contamination,
    checkm1_completeness, checkm1_contamination
  )

# Original MiDAS / GTDB samples are prepared on the server,
# Only GWASMC is subsampled here for the ms. revision
set.seed(random_seed)
gwasmc_plot <- gwasmc_scores %>%
  group_by(quality) %>%
  slice_sample(n = points_per_group) %>%
  ungroup() %>%
  select(all_of(names(original_data)))

plot_data <- bind_rows(original_data, gwasmc_plot) %>%
  mutate(
    group = factor(paste(source, quality, sep = "-"), levels = group_levels),
    y_position = unname(group_positions[as.character(group)]),
    contigs = as.numeric(contigs),
    Contig_N50 = as.numeric(Contig_N50)
  ) %>%
  arrange(group)

coverage_data <- coverage_data %>%
  mutate(
    group = factor(paste(source, quality, sep = "-"), levels = group_levels),
    y_position = unname(group_positions[as.character(group)])
  )

# ----------------------------- #
# Panels 1-2: Completeness and contamination
# ----------------------------- #

p1 <- ggplot(plot_data, aes(x = checkmx_completeness, y = y_position, fill = group)) %>%
  add_box_and_points() +
  scale_x_continuous(limits = c(50, 100), name = "Completeness (%)") +
  theme_base(show_y_text = TRUE)

p2 <- ggplot(plot_data, aes(x = checkmx_contamination, y = y_position, fill = group)) %>%
  add_box_and_points() +
  scale_x_continuous(limits = c(0, 10), name = "Contamination (%)") +
  theme_base(show_y_text = FALSE)

# ----------------------------- #
# Panels 3-4: Contig count and N50
# ----------------------------- #

# Leave log-axis limits unset so all observations enter the boxplots.
p3 <- ggplot(plot_data, aes(x = contigs, y = y_position, fill = group)) %>%
  add_box_and_points() +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c("1", "10", "100", "1000", "10000"),
    name = "Number of contigs"
  ) +
  theme_base(show_y_text = FALSE)

p4 <- ggplot(plot_data, aes(x = Contig_N50 / 1e6, y = y_position, fill = group)) %>%
  add_box_and_points() +
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.1, 1, 10),
    labels = c("0.001", "0.01", "0.1", "1.0", "10.0"),
    name = "N50 (Mbp)"
  ) +
  theme_base(show_y_text = FALSE)

# ----------------------------- #
# Panel 5: Sequencing coverage for MiDAS
# We only have seq. cov info for MiDAS
# ----------------------------- #

p5 <- ggplot(coverage_data, aes(x = Cov, y = y_position, fill = group)) %>%
  add_box_and_points() +
  scale_x_log10(
    breaks = c(5, 10, 20, 100, 1000),
    labels = c("5", "10", "20", "100", "1000"),
    name = "Sequencing coverage"
  ) +
  theme_base(show_y_text = FALSE)

# ----------------------------- #
# Assemble and export
# ----------------------------- #

p_combined <- (p1 | p2 | p3 | p4 | p5) +
  plot_layout(widths = c(1.5, 1, 1, 1, 1))

ggsave(
  "figures/Fig1c.pdf", plot = p_combined,
  width = figure_width, height = figure_height, device = cairo_pdf
)
#ggsave(
#  "figures/Fig1c.png", plot = p_combined,
#  width = figure_width, height = figure_height, dpi = 300,
#  bg = "white"
#)

#write.table(
#  plot_data, "figures/Fig1c.plot_data.tsv",
#  sep = "\t", quote = FALSE, row.names = FALSE
#)
#write.table(
#  gwasmc_scores, "figures/Fig1c.GWASMC_selected_scores.tsv",
#  sep = "\t", quote = FALSE, row.names = FALSE
#)
writeLines(capture.output(sessionInfo()), "../envs/Fig1c.sessionInfo.txt")
