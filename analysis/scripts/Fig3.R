############################################################
## Figure 3 — MAG recovery, community recovery, and core genera
##
## Panels:
##   a) Continental distribution of HQ and MQ MAGs
##   b) Microbial community recovery estimates based on
##      datasets from this study and external datasets
##   c) Number of HQ MAGs recovered for MiDAS core genera
##
## Inputs:
##   - sources/Fig3b.xlsx
##       Sheet: GTDB.r226
##   - sources/Fig3b.xlsx
##       Sheet: GWMC.DK
##   - sources/Fig3c.xlsx
##       Sheet: Genus
##
## Output:
##   - figures/Fig.3.pdf
##
## Notes:
##   - The script should be run from the analysis directory.
##   - Panel labels and the shared legend are constructed
##     explicitly to preserve the manuscript layout.
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(scales)
  library(patchwork)
  library(grid)
  library(cowplot)
})

# ----------------------------- #
# Input and output paths
# ----------------------------- #

input_recovery <- "sources/Fig3b.xlsx"
input_external <- "sources/Fig3b.xlsx"
input_core     <- "sources/Fig3c.xlsx"

output_figure <- "figures/Fig.3.pdf"

if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}

required_files <- c(input_recovery, input_external, input_core)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input file(s): ",
    paste(missing_files, collapse = ", ")
  )
}

# ============================================================
# Panel a: Continental distribution of HQ and MQ MAGs
# ============================================================

continent_mag_data <- data.frame(
  Category = c(
    "Asia",
    "Europe",
    "North America",
    "Oceania",
    "South America"
  ),
  hqmag_num = c(5110, 7451, 5726, 2083, 1907),
  mqmag_num = c(6498, 11595, 7872, 2967, 2292)
)

continent_mag_data$Category <- factor(
  continent_mag_data$Category,
  levels = c(
    "Asia",
    "Europe",
    "North America",
    "Oceania",
    "South America"
  )
)

continent_mag_long <- tidyr::gather(
  continent_mag_data,
  key = "Variable",
  value = "Value",
  -Category
)

continent_mag_long$Variable <- factor(
  continent_mag_long$Variable,
  levels = rev(c("mqmag_num", "hqmag_num"))
)

continent_labels <- c(
  "Asia" = "Asia\n(n=18)",
  "Europe" = "Europe\n(n=30)",
  "North America" = "North America\n(n=21)",
  "Oceania" = "Oceania\n(n=7)",
  "South America" = "South America\n(n=7)"
)

p_continents <- ggplot(
  continent_mag_long,
  aes(x = Category, y = Value, fill = Variable)
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    width = 0.6,
    color = NA
  ) +
  scale_fill_manual(
    values = c("#f4b3ad", "#ed665dFF"),
    labels = c("MQ", "HQ")
  ) +
  geom_text(
    aes(label = scales::comma(Value), y = Value),
    color = "white",
    position = position_stack(vjust = 0.5),
    size = 6
  ) +
  scale_x_discrete(labels = continent_labels) +
  labs(
    x = "",
    y = "Number of MAGs",
    fill = "Quality"
  ) +
  scale_y_continuous(labels = scales::comma) +
  theme_bw(base_size = 26) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    text = element_text(size = 26),
    axis.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    legend.position = c(0.8, 0.75),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(1.0, "cm"),
    legend.key.spacing.y = unit(0.6, "cm"),
    legend.background = element_rect(
      fill = "white",
      color = NA
    ),
    legend.key = element_blank()
  )

# ============================================================
# Shared dataset labels and colors for Panel b
# ============================================================

recovery_levels <- c(
  "GTDB-HQ species-rep",
  "GTDB-HQ+MQ species-rep",
  "MiDAS-HQ species-rep",
  "MiDAS-HQ+MQ species-rep",
  "MiDAS-HQ+MQ MAGs"
)

recovery_colors <- c(
  "GTDB-HQ species-rep"     = "#1E90FF",
  "GTDB-HQ+MQ species-rep"  = "#00468B",
  "MiDAS-HQ species-rep"    = "#C2A5CF",
  "MiDAS-HQ+MQ species-rep" = "#7B3294",
  "MiDAS-HQ+MQ MAGs"        = "#FFD700"
)

# ============================================================
# Panel b1: Recovery estimates for datasets from this study
# ============================================================

recovery_raw <- readxl::read_excel(
  input_recovery,
  sheet = "GTDB.r226",
  skip = 1
)

required_recovery_columns <- c(
  "SA.ID",
  "GTDB.hq_ratio",
  "GTDB_all_ratio",
  "MiDAS.hq.mOTUs_ratio",
  "MiDAS.all.mOTUs_ratio",
  "MiDAS.all.MAGs_ratio",
  "Continent"
)

missing_recovery_columns <- setdiff(
  required_recovery_columns,
  colnames(recovery_raw)
)

if (length(missing_recovery_columns) > 0) {
  stop(
    "Missing columns in ",
    input_recovery,
    ": ",
    paste(missing_recovery_columns, collapse = ", ")
  )
}

recovery_selected <- recovery_raw[
  1:83,
  required_recovery_columns
]

recovery_long <- tidyr::gather(
  recovery_selected,
  key = "dataset",
  value = "Ratio",
  -SA.ID,
  -Continent
)

recovery_long$dataset <- factor(
  recovery_long$dataset,
  levels = c(
    "GTDB.hq_ratio",
    "GTDB_all_ratio",
    "MiDAS.hq.mOTUs_ratio",
    "MiDAS.all.mOTUs_ratio",
    "MiDAS.all.MAGs_ratio"
  ),
  labels = recovery_levels
)

p_minimap2.mapping <- ggplot(
  recovery_long,
  aes(x = Continent, y = Ratio, fill = dataset)
) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.9,
    outlier.shape = NA
  ) +
  geom_point(
    aes(fill = dataset),
    alpha = 0.7,
    shape = 21,
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.25,
      dodge.width = 0.9
    )
  ) +
  scale_fill_manual(
    values = recovery_colors,
    breaks = recovery_levels,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, by = 15)
  ) +
  scale_x_discrete(labels = continent_labels) +
  labs(
    x = "",
    y = "Microbial community\nrecovery estimation (%)",
    fill = ""
  ) +
  theme_bw(base_size = 26) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 26),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(1.0, "cm"),
    legend.spacing.x = unit(1.0, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.title.align = 0,
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# ============================================================
# Panel b2: Recovery estimates for external datasets
# ============================================================

external_data <- readxl::read_excel(
  input_external,
  sheet = "GWMC.DK",
  skip = 0
)

required_external_columns <- c(
  "Continent",
  "Group",
  "ratio"
)

missing_external_columns <- setdiff(
  required_external_columns,
  colnames(external_data)
)

if (length(missing_external_columns) > 0) {
  stop(
    "Missing columns in ",
    input_external,
    ": ",
    paste(missing_external_columns, collapse = ", ")
  )
}

external_data <- external_data %>%
  mutate(
    Continent = factor(
      Continent,
      levels = sort(unique(Continent))
    )
  )

# Each external sample is represented by five dataset groups.
external_continent_counts <- external_data %>%
  count(Continent) %>%
  arrange(match(Continent, levels(external_data$Continent))) %>%
  mutate(
    n = n / 5,
    label = paste0(
      as.character(Continent),
      "\n(n=",
      n,
      ")"
    )
  )

external_data <- external_data %>%
  left_join(
    external_continent_counts %>%
      select(Continent, label),
    by = "Continent"
  ) %>%
  mutate(
    continent_label = factor(
      label,
      levels = external_continent_counts$label
    )
  )

external_data$Group <- factor(
  external_data$Group,
  levels = c(
    "GTDB-HQ_spRep",
    "GTDB-HQ-MQ_spRep",
    "MiDAS-HQ_spRep",
    "MiDAS-HQ-MQ_spRep",
    "MiDAS-HQ-MQ_MAGs"
  ),
  labels = recovery_levels
)

p_sylph.mapping <- ggplot(
  external_data,
  aes(
    x = continent_label,
    y = ratio,
    fill = Group,
    group = interaction(continent_label, Group)
  )
) +
  geom_boxplot(
    width = 0.7,
    alpha = 0.7,
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    alpha = 0.7,
    shape = 21,
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.25,
      dodge.width = 0.8
    )
  ) +
  scale_fill_manual(
    values = recovery_colors,
    breaks = recovery_levels,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Microbial community\nrecovery estimation (%)",
    fill = ""
  ) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, 15)
  ) +
  theme_bw(base_size = 26) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(1.0, "cm"),
    legend.spacing.x = unit(1.0, "cm"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

# ============================================================
# Panel c: HQ MAG recovery across MiDAS core genera
# ============================================================

core_raw <- readxl::read_excel(
  input_core,
  sheet = "Genus"
)

required_core_columns <- c(
  "Genus",
  "No.of_HQ_MAGs",
  "Category"
)

missing_core_columns <- setdiff(
  required_core_columns,
  colnames(core_raw)
)

if (length(missing_core_columns) > 0) {
  stop(
    "Missing columns in ",
    input_core,
    ": ",
    paste(missing_core_columns, collapse = ", ")
  )
}

core_data <- core_raw[
  1:250,
  required_core_columns
]

core_data$Category <- factor(
  core_data$Category,
  levels = c(
    "Strict core",
    "General core",
    "Loose core"
  )
)

core_medians <- aggregate(
  No.of_HQ_MAGs ~ Category,
  data = core_data,
  FUN = median
)

p_ge <- ggplot(
  core_data,
  aes(
    x = Category,
    y = No.of_HQ_MAGs,
    fill = Category
  )
) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.7,
    outlier.shape = NA
  ) +
  geom_point(
    aes(fill = Category),
    alpha = 0.5,
    shape = 21,
    size = 4,
    position = position_jitter(width = 0.25)
  ) +
  scale_fill_manual(
    values = c(
      "Strict core" = "#ED0000FF",
      "General core" = "#00468BFF",
      "Loose core" = "#42B540FF"
    )
  ) +
  scale_y_log10() +
  coord_flip() +
  geom_text(
    data = core_medians,
    aes(
      x = Category,
      y = No.of_HQ_MAGs,
      label = round(No.of_HQ_MAGs, 3)
    ),
    vjust = -3.5,
    size = 6
  ) +
  labs(
    x = "",
    y = "HQ MAGs/core genus",
    fill = ""
  ) +
  theme_bw(base_size = 26) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    text = element_text(size = 26),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "none"
  )

# ============================================================
# Panel-label helper
# ============================================================

make_panel_label <- function(
    letter,
    size = 28,
    family = "sans",
    face = "plain"
) {
  wrap_elements(
    full = textGrob(
      letter,
      x = 0.02,
      y = 0.98,
      hjust = 0,
      vjust = 1,
      gp = gpar(
        fontsize = size,
        fontfamily = family,
        fontface = face
      )
    )
  )
}

label_a <- make_panel_label("a")
label_b <- make_panel_label("b")
label_c <- make_panel_label("c")

# ============================================================
# Construct the shared legend for Panel b
# ============================================================

extract_real_legend <- function(plot_object) {
  guide_grobs <- cowplot::get_plot_component(
    plot_object,
    "guide-box",
    return_all = TRUE
  )
  
  real_guide <- NULL
  
  for (guide in guide_grobs) {
    if (!inherits(guide, "zeroGrob")) {
      real_guide <- guide
      break
    }
  }
  
  if (is.null(real_guide)) {
    stop("No non-empty legend was found.")
  }
  
  real_guide
}

legend_plot <- ggplot(
  data.frame(
    Dataset = factor(
      recovery_levels,
      levels = recovery_levels
    )
  ),
  aes(x = Dataset, fill = Dataset)
) +
  geom_bar() +
  scale_fill_manual(
    values = recovery_colors,
    breaks = recovery_levels
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(1.2, "cm"),
      keyheight = unit(1.0, "cm")
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(1.0, "cm"),
    legend.key.spacing.x = unit(1.2, "cm")
  )

shared_legend <- extract_real_legend(legend_plot)

legend_block <- cowplot::ggdraw() +
  cowplot::draw_grob(
    shared_legend,
    x = 0,
    y = 0,
    width = 1,
    height = 1
  )

# ============================================================
# Prepare Panel b subplots
# ============================================================

p_b1_noleg <- p_minimap2.mapping +
  theme(legend.position = "none") +
  labs(title = "Datasets from this study") +
  theme(
    plot.title = element_text(
      size = 24,
      hjust = 0.5
    )
  )

p_b2_noleg <- p_sylph.mapping +
  theme(legend.position = "none") +
  labs(title = "External datasets") +
  theme(
    plot.title = element_text(
      size = 24,
      hjust = 0.5
    )
  )

p_bottom_plots <- (p_b1_noleg | p_b2_noleg) +
  plot_layout(widths = c(4.5, 5.5))

p_bottom <- (p_bottom_plots / legend_block) +
  plot_layout(heights = c(1, 0.1))

# ============================================================
# Assemble Figure 3
# ============================================================

top_left_block <- label_a / p_continents +
  plot_layout(heights = c(0.04, 1))

top_right_block <- label_c / p_ge +
  plot_layout(heights = c(0.04, 1))

p_top <- top_left_block | top_right_block

bottom_block <- label_b / p_bottom +
  plot_layout(heights = c(0.04, 1))

p_combined <- (p_top / bottom_block) +
  plot_layout(heights = c(1.5, 2.7))

p_combined

# ============================================================
# Export
# ============================================================

ggsave(
  filename = output_figure,
  plot = p_combined,
  width = 24,
  height = 16
)

# Record the R environment used to generate this figure.
if (!dir.exists("../envs")) {
  dir.create("../envs", recursive = TRUE)
}

writeLines(
  capture.output(sessionInfo()),
  "../envs/Fig.3.sessionInfo.txt"
)
