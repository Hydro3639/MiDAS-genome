############################################################
## Figure S15: rRNA recovery within shared species
##
## Input:
##   - sources/FigS15_rRNA.tsv
##     Required columns:
##       Cohort, Catalog, Feature, Positive, Total, Percentage, Length_rule
##       length creteria, same with Fig1d
##
## Output:
##   - figures/FigS15.pdf
##   - ../envs/FigS15.sessionInfo.txt
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
      plot.title = element_text(hjust = 0.5, size = 10, margin = margin(b = 6)),
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

data <- read.delim("sources/FigS15_rRNA.tsv", stringsAsFactors = FALSE)

cohort_sizes <- c(HQ_shared = 2289L, HQ_quality = 738L)
features <- c("5S rRNA", "16S rRNA", "23S rRNA", "Intact operon")
expected_columns <- c("Cohort", "Catalog", "Feature", "Positive", "Total", "Percentage", "Length_rule")


data <- data %>% mutate(
  Catalog = factor(Catalog, levels = group_levels),
  # Show the extra digit when a non-universal result would otherwise read 100.0%.
  Label = ifelse(Positive < Total & round(Percentage, 1) == 100,
                 sprintf("%.2f%%", Percentage), sprintf("%.1f%%", Percentage)),
  Label_x = ifelse(Percentage >= 85, Percentage - 2, Percentage + 2),
  Label_hjust = ifelse(Percentage >= 85, 1, 0)
)

# ----------------------------- #
# Panels a-d: all matched pairs
# ----------------------------- #

# Panel a: 5S rRNA
data_a <- data %>% filter(Cohort == "HQ_shared", Feature == "5S rRNA")
p_a <- ggplot(data_a, aes(x = Percentage, y = Catalog, fill = Catalog)) +
  geom_col(width = 0.55, orientation = "y") +
  geom_text(aes(x = Label_x, label = Label, hjust = Label_hjust),
    size = value_text_size, colour = "black") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 112), breaks = c(0, 50, 100),
    expand = expansion(mult = 0), name = "MAGs (%)") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.65), name = NULL) +
  labs(title = "5S rRNA") +
  theme_base(show_y_text = TRUE)

# Panel b: 16S rRNA
data_b <- data %>% filter(Cohort == "HQ_shared", Feature == "16S rRNA")
p_b <- ggplot(data_b, aes(x = Percentage, y = Catalog, fill = Catalog)) +
  geom_col(width = 0.55, orientation = "y") +
  geom_text(aes(x = Label_x, label = Label, hjust = Label_hjust),
    size = value_text_size, colour = "black") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 112), breaks = c(0, 50, 100),
    expand = expansion(mult = 0), name = "MAGs (%)") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.65), name = NULL) +
  labs(title = "16S rRNA") +
  theme_base(show_y_text = FALSE)

# Panel c: 23S rRNA
data_c <- data %>% filter(Cohort == "HQ_shared", Feature == "23S rRNA")
p_c <- ggplot(data_c, aes(x = Percentage, y = Catalog, fill = Catalog)) +
  geom_col(width = 0.55, orientation = "y") +
  geom_text(aes(x = Label_x, label = Label, hjust = Label_hjust),
    size = value_text_size, colour = "black") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 112), breaks = c(0, 50, 100),
    expand = expansion(mult = 0), name = "MAGs (%)") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.65), name = NULL) +
  labs(title = "23S rRNA") +
  theme_base(show_y_text = FALSE)

# Panel d: Intact operon
data_d <- data %>% filter(Cohort == "HQ_shared", Feature == "Intact operon")
p_d <- ggplot(data_d, aes(x = Percentage, y = Catalog, fill = Catalog)) +
  geom_col(width = 0.55, orientation = "y") +
  geom_text(aes(x = Label_x, label = Label, hjust = Label_hjust),
    size = value_text_size, colour = "black") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 112), breaks = c(0, 50, 100),
    expand = expansion(mult = 0), name = "MAGs (%)") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.65), name = NULL) +
  labs(title = "Intact rRNA operon") +
  theme_base(show_y_text = FALSE)

# ----------------------------- #
# Panels e-h: quality-comparable pairs
# ----------------------------- #

# Panel e: 5S rRNA
data_e <- data %>% filter(Cohort == "HQ_quality", Feature == "5S rRNA")
p_e <- ggplot(data_e, aes(x = Percentage, y = Catalog, fill = Catalog)) +
  geom_col(width = 0.55, orientation = "y") +
  geom_text(aes(x = Label_x, label = Label, hjust = Label_hjust),
    size = value_text_size, colour = "black") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 112), breaks = c(0, 50, 100),
    expand = expansion(mult = 0), name = "MAGs (%)") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.65), name = NULL) +
  labs(title = "5S rRNA") +
  theme_base(show_y_text = TRUE)

# Panel f: 16S rRNA
data_f <- data %>% filter(Cohort == "HQ_quality", Feature == "16S rRNA")
p_f <- ggplot(data_f, aes(x = Percentage, y = Catalog, fill = Catalog)) +
  geom_col(width = 0.55, orientation = "y") +
  geom_text(aes(x = Label_x, label = Label, hjust = Label_hjust),
    size = value_text_size, colour = "black") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 112), breaks = c(0, 50, 100),
    expand = expansion(mult = 0), name = "MAGs (%)") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.65), name = NULL) +
  labs(title = "16S rRNA") +
  theme_base(show_y_text = FALSE)

# Panel g: 23S rRNA
data_g <- data %>% filter(Cohort == "HQ_quality", Feature == "23S rRNA")
p_g <- ggplot(data_g, aes(x = Percentage, y = Catalog, fill = Catalog)) +
  geom_col(width = 0.55, orientation = "y") +
  geom_text(aes(x = Label_x, label = Label, hjust = Label_hjust),
    size = value_text_size, colour = "black") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 112), breaks = c(0, 50, 100),
    expand = expansion(mult = 0), name = "MAGs (%)") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.65), name = NULL) +
  labs(title = "23S rRNA") +
  theme_base(show_y_text = FALSE)

# Panel h: Intact operon
data_h <- data %>% filter(Cohort == "HQ_quality", Feature == "Intact operon")
p_h <- ggplot(data_h, aes(x = Percentage, y = Catalog, fill = Catalog)) +
  geom_col(width = 0.55, orientation = "y") +
  geom_text(aes(x = Label_x, label = Label, hjust = Label_hjust),
    size = value_text_size, colour = "black") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 112), breaks = c(0, 50, 100),
    expand = expansion(mult = 0), name = "MAGs (%)") +
  scale_y_discrete(limits = group_levels, expand = expansion(add = 0.65), name = NULL) +
  labs(title = "Intact rRNA operon") +
  theme_base(show_y_text = FALSE)

p_combined <- ((p_a | p_b | p_c | p_d) / (p_e | p_f | p_g | p_h)) +
  plot_annotation(tag_levels = "a")

# ----------------------------- #
# Assemble and export
# ----------------------------- #

ggsave("figures/FigS15.pdf", p_combined, device = cairo_pdf,
  width = figure_width, height = figure_height, units = "in", bg = "white")
if (interactive()) print(p_combined)
writeLines(capture.output(sessionInfo()), "../envs/FigS15.sessionInfo.txt")
