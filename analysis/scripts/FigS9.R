############################################################
## Figure S9: sensitivity of genus-pathway classifications
##
## Inputs:
##   - sources/FigS9_category_counts.tsv
##   - sources/FigS9_stability_counts.tsv
##   - sources/FigS9_transition_counts.tsv
##
## Outputs:
##   - figures/FigS9.pdf
##   #- figures/FigS9.png (600 dpi)
##   - ../envs/FigS9.sessionInfo.txt
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

# ----------------------------- #
# Shared settings
# ----------------------------- #

cutoff_levels <- c("70/30", "80/20", "90/10")
cutoff_labels <- c("70/30" = "70/30", "80/20" = "80/20\n(reference)",
  "90/10" = "90/10")
category_levels <- c("Conserved present", "Mixed", "Conserved absent")
category_colors <- c("Conserved present" = "#4C78A8", "Mixed" = "#F2A541",
  "Conserved absent" = "#D9DDE2")
stability_colors <- c("Unchanged from 80/20" = "#5A93B2", "Changed" = "#ED665D")
transition_colors <- c("Same category" = "#5A93B2", "Different category" = "#ED665D",
  "Zero" = "#F4F5F6")
base_text_size <- 12
value_text_size <- 3.4
matrix_text_size <- 3.5
figure_width <- 9.1
figure_height <- 8.2
#png_dpi <- 600

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("../envs", showWarnings = FALSE, recursive = TRUE)

theme_base <- function() {
  theme_classic(base_size = base_text_size, base_family = "sans") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.3, "cm"),
      axis.text = element_text(colour = "black", size = 10),
      axis.title.x = element_text(margin = margin(t = 7)),
      axis.ticks.x = element_blank(),
      plot.tag = element_text(size = 14, face = "plain"),
      plot.tag.position = "topleft",
      plot.margin = margin(7, 9, 5, 5)
    )
}

format_count <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)

# ----------------------------- #
# Load prepared figure inputs
# ----------------------------- #

categories <- read.delim("sources/FigS9_category_counts.tsv", stringsAsFactors = FALSE)
stability <- read.delim("sources/FigS9_stability_counts.tsv", stringsAsFactors = FALSE)
transitions <- read.delim("sources/FigS9_transition_counts.tsv", stringsAsFactors = FALSE)

# Category labels: conserved absent at the bottom, conserved present at the top.
categories$Cutoff <- factor(categories$Cutoff, levels = cutoff_levels)
categories$Status <- factor(categories$Status, levels = rev(category_levels))
categories <- categories[order(categories$Cutoff, categories$Status), ]
categories$Label_y <- ave(categories$Percentage, categories$Cutoff,
  FUN = function(x) cumsum(x) - x / 2)
categories$Label <- sprintf("%.1f%%", categories$Percentage)
categories$Label_colour <- ifelse(categories$Status == "Conserved present", "white", "#202020")

# Aggregate labels: percentages of all 6,094 genus-pathway combinations.
stability$Cutoff <- factor(stability$Cutoff, levels = cutoff_levels)
stability <- stability[order(stability$Cutoff), ]
stability$Unchanged_label <- paste0(format_count(stability$Unchanged), "\n(",
  sprintf("%.1f", stability$Unchanged_pct), "%)")
stability$Changed_label <- paste0(format_count(stability$Changed), " (",
  sprintf("%.1f", stability$Changed_pct), "%)")
stability_long <- rbind(
  data.frame(Cutoff = stability$Cutoff, Stability = "Unchanged from 80/20", Count = stability$Unchanged),
  data.frame(Cutoff = stability$Cutoff, Stability = "Changed", Count = stability$Changed)
)
stability_long$Stability <- factor(stability_long$Stability, levels = names(stability_colors))

# Matrix labels: percentages within each original 80/20 category (row).
transitions$Baseline_status <- factor(transitions$Baseline_status, levels = rev(category_levels))
transitions$Alternative_status <- factor(transitions$Alternative_status, levels = category_levels)
transitions$Tile <- ifelse(transitions$Count == 0, "Zero", transitions$Change)
transitions$Label <- ifelse(transitions$Count == 0, "0",
  paste0(format_count(transitions$Count), "\n(",
    sprintf("%.1f", transitions$Percentage_of_baseline_category), "%)"))
transitions$Label_colour <- ifelse(transitions$Count == 0, "#777777",
  ifelse(transitions$Change == "Same category", "white", "#202020"))

# ----------------------------- #
# Panel a: category proportions
# ----------------------------- #

p_a <- ggplot(categories, aes(x = Cutoff, y = Percentage, fill = Status)) +
  geom_col(width = 0.58, position = position_stack(reverse = TRUE),
    colour = "white", linewidth = 0.3) +
  geom_text(aes(y = Label_y, label = Label, colour = Label_colour),
    size = value_text_size, show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_manual(values = category_colors, breaks = category_levels) +
  scale_x_discrete(labels = cutoff_labels, expand = expansion(add = 0.65)) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 20),
    expand = expansion(mult = 0)) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Conservation cutoffs (present/absent)",
    y = "Genus-pathway combinations (%)", tag = "a") +
  theme_base()

# ----------------------------- #
# Panel b: aggregate classification changes
# ----------------------------- #

p_b <- ggplot(stability_long, aes(x = Cutoff, y = Count, fill = Stability)) +
  geom_col(width = 0.58, position = position_stack(reverse = TRUE),
    colour = "white", linewidth = 0.3) +
  geom_text(data = stability, aes(x = Cutoff, y = Unchanged / 2, label = Unchanged_label),
    inherit.aes = FALSE, colour = "white", size = 3.1, lineheight = 1.2) +
  geom_text(data = stability, aes(x = Cutoff, y = Total + 185, label = Changed_label),
    inherit.aes = FALSE, colour = "#B9423B", size = value_text_size, vjust = 0) +
  scale_fill_manual(values = stability_colors, breaks = names(stability_colors)) +
  scale_x_discrete(labels = cutoff_labels, expand = expansion(add = 0.65)) +
  scale_y_continuous(limits = c(0, 6990), breaks = c(0, 2000, 4000, 6000),
    labels = format_count, expand = expansion(mult = 0)) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Conservation cutoffs (present/absent)",
    y = "Genus-pathway combinations (n)", tag = "b") +
  theme_base()

# ----------------------------- #
# Panels c-d: category transition matrices
# ----------------------------- #

matrix_panel <- function(cutoff, tag) {
  axis_labels <- c("Conserved present" = "Conserved\npresent", "Mixed" = "Mixed",
    "Conserved absent" = "Conserved\nabsent")
  ggplot(transitions[transitions$Cutoff == cutoff, ],
    aes(x = Alternative_status, y = Baseline_status, fill = Tile)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = Label, colour = Label_colour), size = matrix_text_size,
      lineheight = 1.2, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_manual(values = transition_colors,
      breaks = c("Same category", "Different category"), drop = FALSE) +
    scale_x_discrete(labels = axis_labels, limits = category_levels,
      expand = expansion(mult = 0)) +
    scale_y_discrete(labels = axis_labels, limits = rev(category_levels),
      expand = expansion(mult = 0)) +
    coord_fixed() +
    guides(fill = guide_legend(nrow = 1)) +
    labs(x = paste("Classification at", cutoff),
      y = "Original classification (80/20)", tag = tag) +
    theme_base() +
    theme(axis.line = element_blank(), axis.ticks = element_blank())
}

# Panel c: 80/20 to 70/30.
p_c <- matrix_panel("70/30", "c")

# Panel d: 80/20 to 90/10.
p_d <- matrix_panel("90/10", "d")

# ----------------------------- #
# Assemble and export
# ----------------------------- #

p_upper <- (p_a | p_b) + plot_layout(guides = "keep")
p_lower <- ((p_c | p_d) + plot_layout(guides = "collect")) &
  theme(legend.position = "bottom")
p_combined <- (p_upper / p_lower) + plot_layout(heights = c(1, 1))

ggsave("figures/FigS9.pdf", p_combined,
  width = figure_width, height = figure_height, units = "in", bg = "white")
#ggsave("figures/FigS9.png", p_combined,
#  width = figure_width, height = figure_height, units = "in",
#  dpi = png_dpi, type = "cairo", bg = "white")
if (interactive()) print(p_combined)
writeLines(capture.output(sessionInfo()), "../envs/FigS9.sessionInfo.txt")
