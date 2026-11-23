############################################################
## Figure S16: BGC recovery within shared species
##
## Inputs:
##   - sources/FigS16_BGC.tsv
##   - sources/FigS16_summary.tsv
##
## Outputs:
##   - figures/FigS16.pdf
##   #- figures/FigS16.png (600 dpi)
##   - ../envs/FigS16.sessionInfo.txt
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)
})

# ----------------------------- #
# settings
# ----------------------------- #

group_colors <- c("MiDAS-HQ" = "#ed665dFF", "GWASMC" = "#59a14f")
group_levels <- c("GWASMC", "MiDAS-HQ")
base_text_size <- 9
value_text_size <- 2.7
median_y_offset <- 0.50
point_size <- 0.35
point_alpha <- 0.17
#png_dpi <- 600

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("../envs", showWarnings = FALSE, recursive = TRUE)

theme_base <- function(show_y_text = TRUE) {
  theme_classic(base_size = base_text_size, base_family = "sans") +
    theme(
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
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(6, 9, 5, 5),
      plot.tag = element_text(size = 11, face = "plain"),
      plot.tag.position = "topleft"
    )
}

figure_width <- 8.6
figure_height <- 3.8

# ----------------------------- #
# Load prepared figure inputs
# ----------------------------- #

data <- read.delim("sources/FigS16_BGC.tsv", stringsAsFactors = FALSE)
summary <- read.delim("sources/FigS16_summary.tsv", stringsAsFactors = FALSE)

data$Catalog <- factor(data$Catalog, levels = group_levels)
summary$Catalog <- factor(summary$Catalog, levels = group_levels)
data_all <- data %>% filter(Cohort == "HQ_shared")
data_quality <- data %>% filter(Cohort == "HQ_quality")
lengths_all <- data_all %>% filter(Length_eligible == 1)
lengths_quality <- data_quality %>% filter(Length_eligible == 1)
summary_all <- summary %>% filter(Cohort == "HQ_shared")
summary_quality <- summary %>% filter(Cohort == "HQ_quality")


# ----------------------------- #
# Panels a-b: all matched pairs
# ----------------------------- #

# Panel a
p_a <- ggplot(data_all, aes(x = Complete_regions, y = Catalog, fill = Catalog)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
    width = 0.62, linewidth = 0.25, colour = "black", alpha = 0.88) +
  geom_point(position = position_jitter(width = 0, height = 0.055, seed = 20261113),
    colour = "black", size = point_size, alpha = point_alpha) +
  geom_boxplot(orientation = "y", width = 0.16, outlier.shape = NA,
    fill = "white", linewidth = 0.3) +
  geom_text(data = summary_all %>% filter(Metric == "Complete_regions"),
    aes(x = Median, y = Catalog, label = format(round(Median), trim = TRUE, scientific = FALSE)),
    inherit.aes = FALSE, nudge_y = median_y_offset, size = value_text_size) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.82), name = NULL) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40),
    expand = expansion(mult = c(0.02, 0.06))) +
  coord_cartesian(xlim = c(-0.8, 41), clip = "off") +
  labs(x = "Complete predicted BGC regions per MAG") +
  theme_base(show_y_text = TRUE)

# Panel b
p_b <- ggplot(lengths_all, aes(x = Median_region_length_kb, y = Catalog, fill = Catalog)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
    width = 0.62, linewidth = 0.25, colour = "black", alpha = 0.88) +
  geom_point(position = position_jitter(width = 0, height = 0.055, seed = 20261113),
    colour = "black", size = point_size, alpha = point_alpha) +
  geom_boxplot(orientation = "y", width = 0.16, outlier.shape = NA,
    fill = "white", linewidth = 0.3) +
  geom_text(data = summary_all %>% filter(Metric == "Median_region_length_kb"),
    aes(x = Median, y = Catalog, label = sprintf("%.2f", Median)),
    inherit.aes = FALSE, nudge_y = median_y_offset, size = value_text_size) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.82), name = NULL) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80),
    expand = expansion(mult = c(0.02, 0.06))) +
  coord_cartesian(xlim = c(0, 82), clip = "off") +
  labs(x = "Median region length per MAG (kbp)") +
  theme_base(show_y_text = FALSE)

# ----------------------------- #
# Panels c-d: quality-comparable pairs
# ----------------------------- #

# Panel c
p_c <- ggplot(data_quality, aes(x = Complete_regions, y = Catalog, fill = Catalog)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
    width = 0.62, linewidth = 0.25, colour = "black", alpha = 0.88) +
  geom_point(position = position_jitter(width = 0, height = 0.055, seed = 20261113),
    colour = "black", size = point_size, alpha = point_alpha) +
  geom_boxplot(orientation = "y", width = 0.16, outlier.shape = NA,
    fill = "white", linewidth = 0.3) +
  geom_text(data = summary_quality %>% filter(Metric == "Complete_regions"),
    aes(x = Median, y = Catalog, label = format(round(Median), trim = TRUE, scientific = FALSE)),
    inherit.aes = FALSE, nudge_y = median_y_offset, size = value_text_size) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.82), name = NULL) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40),
    expand = expansion(mult = c(0.02, 0.06))) +
  coord_cartesian(xlim = c(-0.8, 41), clip = "off") +
  labs(x = "Complete predicted BGC regions per MAG") +
  theme_base(show_y_text = TRUE)

# Panel d
p_d <- ggplot(lengths_quality, aes(x = Median_region_length_kb, y = Catalog, fill = Catalog)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
    width = 0.62, linewidth = 0.25, colour = "black", alpha = 0.88) +
  geom_point(position = position_jitter(width = 0, height = 0.055, seed = 20261113),
    colour = "black", size = point_size, alpha = point_alpha) +
  geom_boxplot(orientation = "y", width = 0.16, outlier.shape = NA,
    fill = "white", linewidth = 0.3) +
  geom_text(data = summary_quality %>% filter(Metric == "Median_region_length_kb"),
    aes(x = Median, y = Catalog, label = sprintf("%.2f", Median)),
    inherit.aes = FALSE, nudge_y = median_y_offset, size = value_text_size) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.82), name = NULL) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80),
    expand = expansion(mult = c(0.02, 0.06))) +
  coord_cartesian(xlim = c(0, 82), clip = "off") +
  labs(x = "Median region length per MAG (kbp)") +
  theme_base(show_y_text = FALSE)

p_combined <- ((p_a | p_b) / (p_c | p_d)) +
  plot_annotation(tag_levels = "a")

# ----------------------------- #
# Assemble and export
# ----------------------------- #

ggsave("figures/FigS16.pdf", p_combined, device = cairo_pdf,
  width = figure_width, height = figure_height, units = "in", bg = "white")
#ggsave("figures/FigS16.png", p_combined,
#  width = figure_width, height = figure_height, units = "in",
#  dpi = png_dpi, type = "cairo", bg = "white")
if (interactive()) print(p_combined)
writeLines(capture.output(sessionInfo()), "../envs/FigS16.sessionInfo.txt")
