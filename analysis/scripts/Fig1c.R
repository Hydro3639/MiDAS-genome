############################################################
## Fig. 1c: Distribution of MAG quality metrics
##
## This script generates Figure 1c (multi-panel) including:
##  - CheckM2 completeness
##  - CheckM2 contamination
##  - Number of contigs (log10)
##  - Contig N50 in Mbp (log10; computed as Contig_N50/1e6)
##  - Sequencing coverage (MiDAS only; log10 transform)
##
## Inputs:
##   - sources/Fig1c-I.tsv
##   - sources/Fig1c-II.tsv
##
## Output:
##   - figures/Figure1c.pdf
############################################################

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

# ----------------------------- #
# Shared settings
# ----------------------------- #

group_colors <- c(
  "MiDAS-HQ" = "#ed665dFF",
  "GTDB-HQ"  = "#5a93b2",
  "MiDAS-MQ" = "#f4b3ad",
  "GTDB-MQ"  = "#a6cee3FF"
)

# Base theme used across panels
theme_base <- function(show_y_text = TRUE) {
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    axis.line.x = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks.x = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(color = "black"),
    text = element_text(size = 30),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.ticks.y = element_blank(),
    axis.text.y = if (show_y_text) element_text(color = "black") else element_blank(),
    axis.title.y = element_blank()
  )
}

# Common layers: horizontal boxplot + jittered points
add_box_and_points <- function(p) {
  p +
    geom_boxplot(
      width = 0.25,
      outlier.shape = NA,
      orientation = "y"
    ) +
    geom_point(
      position = position_jitter(width = 0, height = 0.08),
      size = 1,
      alpha = 0.3,
      color = "black"
    )
}

# ----------------------------- #
# Load and preprocess: Fig1c-I
# ----------------------------- #

data <- read.delim("sources/Fig1c-I.tsv")

data <- data %>%
  mutate(
    group = paste(source, quality, sep = "-"),
    y_position = case_when(
      group == "MiDAS-MQ" ~ 0.2,
      group == "MiDAS-HQ" ~ 0.5,
      group == "GTDB-MQ"  ~ 0.9,
      group == "GTDB-HQ"  ~ 1.2
    ),
    group = factor(group, levels = c("MiDAS-MQ", "MiDAS-HQ", "GTDB-HQ", "GTDB-MQ"))
  )

# ----------------------------- #
# Panel 1: CheckM2 completeness
# (keeps original y labels visible)
# ----------------------------- #

p1 <- ggplot(data, aes(x = checkm2_completeness, y = y_position, fill = group)) %>%
  add_box_and_points() +
  scale_x_continuous(
    limits = c(50, 100),
    name = "Completeness (%)"
  ) +
  scale_y_continuous(
    breaks = c(0.2, 0.5, 0.9, 1.2),
    labels = c("MiDAS-MQ", "MiDAS-HQ", "GTDB-MQ", "GTDB-HQ"),
    name = NULL,
    expand = expansion(add = 0.2)
  ) +
  scale_fill_manual(values = group_colors, guide = "none") +
  theme_base(show_y_text = TRUE)

# ----------------------------- #
# Panel 2: CheckM2 contamination
# (keeps original behavior: remove y text)
# ----------------------------- #

p2 <- ggplot(data, aes(x = checkm2_contamination, y = y_position, fill = group)) %>%
  add_box_and_points() +
  scale_x_continuous(
    limits = c(0, 10),
    name = "Contamination (%)"
  ) +
  scale_y_continuous(
    breaks = c(0.2, 0.5, 0.9, 1.2),
    labels = c("MiDAS-MQ", "MiDAS-HQ", "GTDB-MQ", "GTDB-HQ"),
    name = NULL,
    expand = expansion(add = 0.2)
  ) +
  scale_fill_manual(values = group_colors, guide = "none") +
  theme_base(show_y_text = FALSE)

# ----------------------------- #
# Panel 3: Contig count per MAG
# (keeps original conversion + x log10)
# ----------------------------- #

data$contigs <- as.numeric(as.character(data$contigs))

p3 <- ggplot(data, aes(x = contigs, y = y_position, fill = group)) %>%
  add_box_and_points() +
  scale_x_log10(
    limits = c(1, 1024),
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000"),
    name = "Number of contigs"
  ) +
  scale_y_continuous(
    breaks = c(0.2, 0.5, 0.9, 1.2),
    labels = c("MiDAS-MQ", "MiDAS-HQ", "GTDB-MQ", "GTDB-HQ"),
    name = NULL,
    expand = expansion(add = 0.2)
  ) +
  scale_fill_manual(values = group_colors, guide = "none") +
  theme_base(show_y_text = FALSE)

# ----------------------------- #
# Panel 4: Contig N50 (Mbp)
# (keeps original: x = Contig_N50/1e6 + log10 scale)
# ----------------------------- #

p4 <- ggplot(data, aes(x = Contig_N50 / 1000000, y = y_position, fill = group)) %>%
  add_box_and_points() +
  scale_x_log10(
    limits = c(0.005, 15),
    breaks = c(0.01, 0.1, 1, 10),
    labels = c("0.01", "0.1", "1.0", "10.0"),
    name = "N50 (Mbp)"
  ) +
  scale_y_continuous(
    breaks = c(0.2, 0.5, 0.9, 1.2),
    labels = c("MiDAS-MQ", "MiDAS-HQ", "GTDB-MQ", "GTDB-HQ"),
    name = NULL,
    expand = expansion(add = 0.2)
  ) +
  scale_fill_manual(values = group_colors, guide = "none") +
  theme_base(show_y_text = FALSE)

# ----------------------------- #
# Load and preprocess: Fig1c-II (coverage; MiDAS only)
# ----------------------------- #

data2 <- read.delim("sources/Fig1c-II.tsv")

data2 <- data2 %>%
  mutate(
    group = paste(source, quality, sep = "-"),
    y_position = case_when(
      group == "MiDAS-MQ" ~ 0.2,
      group == "MiDAS-HQ" ~ 0.5
    ),
    group = factor(group, levels = c("MiDAS-MQ", "MiDAS-HQ"))
  )

# ----------------------------- #
# Panel 5: Sequencing coverage
# (keeps original log10 transform and axis settings)
# ----------------------------- #

p5 <- ggplot(data2, aes(x = Cov, y = y_position, fill = group)) %>%
  add_box_and_points() +
  scale_x_continuous(
    trans = log10_trans(),
    limits = c(3, 1300),
    breaks = c(5, 10, 20, 100, 1000),
    labels = c("5", "10", "20", "100", "1000"),
    name = "Sequencing coverage"
  ) +
  scale_y_continuous(
    breaks = c(0.2, 0.5),
    labels = c("MiDAS-MQ", "MiDAS-HQ"),
    name = NULL,
    expand = expansion(add = 0.2)
  ) +
  scale_fill_manual(
    values = c(
      "MiDAS-HQ" = "#ed665dFF",
      "MiDAS-MQ" = "#f4b3ad"
    ),
    guide = "none"
  ) +
  theme_base(show_y_text = FALSE)

# ----------------------------- #
# Assemble and export
# (keeps original patchwork structure)
# ----------------------------- #

p5_layout <- (plot_spacer() / p5) + plot_layout(heights = c(1, 1.4))

p_combined <- (p1 | p2 | p3 | p4 | p5_layout) +
  plot_layout(widths = c(1.5, 1, 1, 1, 1))

ggsave(p_combined, filename = "figures/Fig1c.pdf", width = 32, height = 10)

writeLines(capture.output(sessionInfo()), "../envs/Fig1c.sessionInfo.txt")
