############################################################
## Figure S14: Genome contiguity within shared species
##
## Inputs:
##   - sources/FigS14_contiguity.tsv
##   - sources/FigS14_summary.tsv
##   - sources/FigS14_N50_ratios.tsv
##
## Outputs:
##   - figures/FigS14.pdf
##   ##- figures/FigS14.png (600 dpi)
##   - ../envs/FigS14.sessionInfo.txt
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)
})

# ----------------------------- #
# Shared settings
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
# Load the prepared figure inputs
# ----------------------------- #

data <- read.delim("sources/FigS14_contiguity.tsv", stringsAsFactors = FALSE)
summary <- read.delim("sources/FigS14_summary.tsv", stringsAsFactors = FALSE)
ratios <- read.delim("sources/FigS14_N50_ratios.tsv", stringsAsFactors = FALSE)

data$Catalog <- factor(data$Catalog, levels = group_levels)
summary$Catalog <- factor(summary$Catalog, levels = c(group_levels, "MiDAS/GWASMC"))
data_all <- data %>% filter(Cohort == "HQ_shared")
data_quality <- data %>% filter(Cohort == "HQ_quality")
summary_all <- summary %>% filter(Cohort == "HQ_shared")
summary_quality <- summary %>% filter(Cohort == "HQ_quality")
ratios_all <- ratios %>% filter(Cohort == "HQ_shared") %>% mutate(Catalog = "MiDAS/GWASMC")
ratios_quality <- ratios %>% filter(Cohort == "HQ_quality") %>% mutate(Catalog = "MiDAS/GWASMC")


# ----------------------------- #
# Panels a-c: all matched pairs
# ----------------------------- #

# Panel a
p_a <- ggplot(data_all, aes(x = N50_Mbp, y = Catalog, fill = Catalog)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
    width = 0.62, linewidth = 0.25, colour = "black", alpha = 0.88) +
  geom_point(position = position_jitter(width = 0, height = 0.055, seed = 20261113),
    colour = "black", size = point_size, alpha = point_alpha) +
  geom_boxplot(orientation = "y", width = 0.16, outlier.shape = NA,
    fill = "white", linewidth = 0.3) +
  geom_text(data = summary_all %>% filter(Metric == "N50_Mbp"),
    aes(x = Median, y = Catalog, label = formatC(Median, format = "fg", digits = 3, flag = "#")),
    inherit.aes = FALSE, nudge_y = median_y_offset, size = value_text_size) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.82), name = NULL) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10),
    labels = c("0.001", "0.01", "0.1", "1", "10"),
    expand = expansion(mult = c(0.03, 0.06))) +
  coord_cartesian(xlim = c(0.0005, 15), clip = "off") +
  labs(x = "N50 (Mbp)") +
  theme_base(show_y_text = TRUE)

# Panel b
p_b <- ggplot(data_all, aes(x = Contig_count, y = Catalog, fill = Catalog)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
    width = 0.62, linewidth = 0.25, colour = "black", alpha = 0.88) +
  geom_point(position = position_jitter(width = 0, height = 0.055, seed = 20261113),
    colour = "black", size = point_size, alpha = point_alpha) +
  geom_boxplot(orientation = "y", width = 0.16, outlier.shape = NA,
    fill = "white", linewidth = 0.3) +
  geom_text(data = summary_all %>% filter(Metric == "Contig_count"),
    aes(x = Median, y = Catalog, label = format(round(Median), trim = TRUE, big.mark = ",", scientific = FALSE)),
    inherit.aes = FALSE, nudge_y = median_y_offset, size = value_text_size) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.82), name = NULL) +
  scale_x_log10(breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1,000"),
    expand = expansion(mult = c(0.03, 0.06))) +
  coord_cartesian(xlim = c(0.8, 6000), clip = "off") +
  labs(x = "Contigs per MAG") +
  theme_base(show_y_text = FALSE)

# Panel c
p_c <- ggplot(ratios_all, aes(x = N50_ratio, y = Catalog, fill = Catalog)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
    width = 0.62, linewidth = 0.25, colour = "black", alpha = 0.88) +
  geom_point(position = position_jitter(width = 0, height = 0.055, seed = 20261113),
    colour = "black", size = point_size, alpha = point_alpha) +
  geom_boxplot(orientation = "y", width = 0.16, outlier.shape = NA,
    fill = "white", linewidth = 0.3) +
  geom_text(data = summary_all %>% filter(Metric == "N50_ratio"),
    aes(x = Median, y = Catalog, label = format(round(Median), trim = TRUE, big.mark = ",", scientific = FALSE)),
    inherit.aes = FALSE, nudge_y = median_y_offset, size = value_text_size) +
  scale_fill_manual(values = c("MiDAS/GWASMC" = "#ed665dFF"), guide = "none") +
  scale_y_discrete(limits = "MiDAS/GWASMC", expand = expansion(add = 0.82), name = NULL) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000),
    labels = c("0.1", "1", "10", "100", "1,000", "10,000"),
    expand = expansion(mult = c(0.03, 0.06))) +
  coord_cartesian(clip = "off") +
  labs(x = "MiDAS / GWASMC") +
  theme_base(show_y_text = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.3)

# ----------------------------- #
# Panels d-f: quality-comparable pairs
# ----------------------------- #

# Panel d
p_d <- ggplot(data_quality, aes(x = N50_Mbp, y = Catalog, fill = Catalog)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
    width = 0.62, linewidth = 0.25, colour = "black", alpha = 0.88) +
  geom_point(position = position_jitter(width = 0, height = 0.055, seed = 20261113),
    colour = "black", size = point_size, alpha = point_alpha) +
  geom_boxplot(orientation = "y", width = 0.16, outlier.shape = NA,
    fill = "white", linewidth = 0.3) +
  geom_text(data = summary_quality %>% filter(Metric == "N50_Mbp"),
    aes(x = Median, y = Catalog, label = formatC(Median, format = "fg", digits = 3, flag = "#")),
    inherit.aes = FALSE, nudge_y = median_y_offset, size = value_text_size) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.82), name = NULL) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10),
    labels = c("0.001", "0.01", "0.1", "1", "10"),
    expand = expansion(mult = c(0.03, 0.06))) +
  coord_cartesian(xlim = c(0.0005, 15), clip = "off") +
  labs(x = "N50 (Mbp)") +
  theme_base(show_y_text = TRUE)

# Panel e
p_e <- ggplot(data_quality, aes(x = Contig_count, y = Catalog, fill = Catalog)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
    width = 0.62, linewidth = 0.25, colour = "black", alpha = 0.88) +
  geom_point(position = position_jitter(width = 0, height = 0.055, seed = 20261113),
    colour = "black", size = point_size, alpha = point_alpha) +
  geom_boxplot(orientation = "y", width = 0.16, outlier.shape = NA,
    fill = "white", linewidth = 0.3) +
  geom_text(data = summary_quality %>% filter(Metric == "Contig_count"),
    aes(x = Median, y = Catalog, label = format(round(Median), trim = TRUE, big.mark = ",", scientific = FALSE)),
    inherit.aes = FALSE, nudge_y = median_y_offset, size = value_text_size) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.82), name = NULL) +
  scale_x_log10(breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1,000"),
    expand = expansion(mult = c(0.03, 0.06))) +
  coord_cartesian(xlim = c(0.8, 6000), clip = "off") +
  labs(x = "Contigs per MAG") +
  theme_base(show_y_text = FALSE)

# Panel f
p_f <- ggplot(ratios_quality, aes(x = N50_ratio, y = Catalog, fill = Catalog)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
    width = 0.62, linewidth = 0.25, colour = "black", alpha = 0.88) +
  geom_point(position = position_jitter(width = 0, height = 0.055, seed = 20261113),
    colour = "black", size = point_size, alpha = point_alpha) +
  geom_boxplot(orientation = "y", width = 0.16, outlier.shape = NA,
    fill = "white", linewidth = 0.3) +
  geom_text(data = summary_quality %>% filter(Metric == "N50_ratio"),
    aes(x = Median, y = Catalog, label = format(round(Median), trim = TRUE, big.mark = ",", scientific = FALSE)),
    inherit.aes = FALSE, nudge_y = median_y_offset, size = value_text_size) +
  scale_fill_manual(values = c("MiDAS/GWASMC" = "#ed665dFF"), guide = "none") +
  scale_y_discrete(limits = "MiDAS/GWASMC", expand = expansion(add = 0.82), name = NULL) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000),
    labels = c("0.1", "1", "10", "100", "1,000", "10,000"),
    expand = expansion(mult = c(0.03, 0.06))) +
  coord_cartesian(clip = "off") +
  labs(x = "MiDAS / GWASMC") +
  theme_base(show_y_text = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.3)

p_combined <- ((p_a | p_b | p_c) / (p_d | p_e | p_f)) +
  plot_layout(widths = c(1.1, 1, 1.05)) + plot_annotation(tag_levels = "a")

# ----------------------------- #
# Assemble and export
# ----------------------------- #

ggsave("figures/FigS14.pdf", p_combined, device = cairo_pdf,
  width = figure_width, height = figure_height, units = "in", bg = "white")
#ggsave("figures/FigS14.png", p_combined,
#  width = figure_width, height = figure_height, units = "in",
#  dpi = png_dpi, type = "cairo", bg = "white")
if (interactive()) print(p_combined)
writeLines(capture.output(sessionInfo()), "../envs/FigS14.sessionInfo.txt")
