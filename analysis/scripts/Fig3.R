############################################################
## Input:
##   - sources/Fig3b.xlsx
##     Sheet: "GTDB.r226" (skip = 1)
##     Required columns:
##       SA.ID, Continent, GTDB.hq.ratio, GTDB.all.ratio,
##       midas.hq.mOTUs.ratio, midas.all.mOTUs.ratio,
##       midas.allMAGs.ratio
##     Sheet: "GWMC.DK"
##     Required columns:
##       Accession, ratio, Group, Continent
##   - sources/Fig3b_GWASMC_recovery.tsv
##     Required columns:
##       Dataset, SampleID, GWASMC_HQ, GWASMC_HQMQ
##   - sources/Fig3c.xlsx
##     Sheet: "Genus"
##     Required columns:
##       Genus, No.of_HQ_MAGs, Category
##
## Output:
##   - figures/Fig3.pdf
##   - ../envs/Fig3.sessionInfo.txt
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(patchwork)
  library(tidyr)
  library(scales)
  library(grid)
})

input_reference <- "sources/Fig3b.xlsx"
input_gwasmc <- "sources/Fig3b_GWASMC_recovery.tsv"
input_core <- "sources/Fig3c.xlsx"
output_figure <- "figures/Fig3.pdf"

required_files <- c(input_reference, input_gwasmc, input_core)
if (any(!file.exists(required_files))) {
  stop("Missing input file(s): ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}
for (directory in c("figures", "../envs")) {
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
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
  levels = c("mqmag_num", "hqmag_num")
)

## should auto number generation like what i did in the other R script,
## but is fine here right now
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
    values = c("mqmag_num" = "#f4b3ad", "hqmag_num" = "#ed665dFF"),
    breaks = c("mqmag_num", "hqmag_num"),
    labels = c("mqmag_num" = "MQ", "hqmag_num" = "HQ")
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
# Panel b: Community recovery
# include gwasmc sp (HQ, HQ+MQ) for the ms revision
# GWASMC-HQ was followed the original GWASMC paper,
# no tRNA/rRNA genes were required.
# GTDB and MiDAS follow the stringent MIMAG standard.
# ============================================================

require_columns <- function(data, columns, source) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop("Missing columns in ", source, ": ", paste(missing, collapse = ", "))
}
require_unique <- function(values, source) {
  if (anyNA(values) || any(!nzchar(values)) || anyDuplicated(values)) {
    stop("Missing or duplicate identifiers in ", source)
  }
}
write_tsv <- function(data, path) {
  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

recovery_levels <- c(
  "GTDB-HQ species-rep",
  "GTDB-HQ+MQ species-rep",
  "GWASMC-HQ species-rep",
  "GWASMC-HQ+MQ species-rep",
  "MiDAS-HQ species-rep",
  "MiDAS-HQ+MQ species-rep",
  "MiDAS-HQ+MQ MAGs"
)
recovery_colors <- c(
  "GTDB-HQ species-rep" = "#1E90FF",
  "GTDB-HQ+MQ species-rep" = "#00468B",
  "GWASMC-HQ species-rep" = "#8DC884",
  "GWASMC-HQ+MQ species-rep" = "#2E7D32",
  "MiDAS-HQ species-rep" = "#C2A5CF",
  "MiDAS-HQ+MQ species-rep" = "#7B3294",
  "MiDAS-HQ+MQ MAGs" = "#FFD700"
)
continent_levels <- c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")

# Select the prepared Sylph recovery columns for the 83 metagenomes generated in this study.
study_raw <- as.data.frame(read_excel(input_reference, sheet = "GTDB.r226", skip = 1))
study_columns <- c(
  "SA.ID", "Continent", "GTDB.hq.ratio", "GTDB.all.ratio",
  "midas.hq.mOTUs.ratio", "midas.all.mOTUs.ratio", "midas.allMAGs.ratio"
)
require_columns(study_raw, study_columns, "GTDB.r226")
study <- study_raw[!is.na(study_raw$SA.ID) & grepl("^[0-9]{3}-", study_raw$SA.ID), study_columns]
require_unique(study$SA.ID, "study workbook")
stopifnot(nrow(study) == 83L, !anyNA(study), all(study$Continent %in% continent_levels))

external <- as.data.frame(read_excel(input_reference, sheet = "GWMC.DK"))
require_columns(external, c("Accession", "ratio", "Group", "Continent"), "GWMC.DK")
require_unique(paste(external$Accession, external$Group, sep = "::"), "external sample/group pairs")
external_groups <- c(
  "GTDB-HQ_spRep", "GTDB-HQ-MQ_spRep", "MiDAS-HQ_spRep",
  "MiDAS-HQ-MQ_spRep", "MiDAS-HQ-MQ_MAGs"
)
stopifnot(setequal(unique(external$Group), external_groups))
external_meta <- unique(external[c("Accession", "Continent")])
require_unique(external_meta$Accession, "external metadata")
stopifnot(nrow(external_meta) == 295L, all(external_meta$Continent %in% continent_levels))
for (group in external_groups) {
  stopifnot(setequal(external$Accession[external$Group == group], external_meta$Accession))
}

# GWASMC percentages have already been calculated (sylph-based) before running R.
gwasmc <- read.delim(input_gwasmc, check.names = FALSE, stringsAsFactors = FALSE)
require_columns(gwasmc, c("Dataset", "SampleID", "GWASMC_HQ", "GWASMC_HQMQ"), input_gwasmc)

sample_metadata <- rbind(
  data.frame(Dataset = "Study", SampleID = study$SA.ID, Continent = study$Continent),
  data.frame(Dataset = "External", SampleID = external_meta$Accession, Continent = external_meta$Continent)
)
sample_key <- function(data) paste(data$Dataset, data$SampleID, sep = "::")
require_unique(sample_key(gwasmc), input_gwasmc)
gwasmc$Continent <- sample_metadata$Continent[match(sample_key(gwasmc), sample_key(sample_metadata))]

study_map <- c(
  "GTDB.hq.ratio" = "GTDB-HQ species-rep",
  "GTDB.all.ratio" = "GTDB-HQ+MQ species-rep",
  "midas.hq.mOTUs.ratio" = "MiDAS-HQ species-rep",
  "midas.all.mOTUs.ratio" = "MiDAS-HQ+MQ species-rep",
  "midas.allMAGs.ratio" = "MiDAS-HQ+MQ MAGs"
)
study_long <- do.call(rbind, lapply(names(study_map), function(column) {
  data.frame(Dataset = "Study", SampleID = study$SA.ID, Continent = study$Continent,
             Group = unname(study_map[column]), ratio = as.numeric(study[[column]]),
             stringsAsFactors = FALSE)
}))
external_map <- setNames(
  c("GTDB-HQ species-rep", "GTDB-HQ+MQ species-rep", "MiDAS-HQ species-rep",
    "MiDAS-HQ+MQ species-rep", "MiDAS-HQ+MQ MAGs"), external_groups
)
external_long <- data.frame(
  Dataset = "External", SampleID = external$Accession, Continent = external$Continent,
  Group = unname(external_map[external$Group]), ratio = as.numeric(external$ratio),
  stringsAsFactors = FALSE
)
gw_long <- rbind(
  data.frame(gwasmc[c("Dataset", "SampleID", "Continent")], Group = "GWASMC-HQ species-rep", ratio = gwasmc$GWASMC_HQ),
  data.frame(gwasmc[c("Dataset", "SampleID", "Continent")], Group = "GWASMC-HQ+MQ species-rep", ratio = gwasmc$GWASMC_HQMQ)
)
plot_data <- rbind(study_long, external_long, gw_long)

require_unique(paste(plot_data$Dataset, plot_data$SampleID, plot_data$Group, sep = "::"), "plot data")

plot_data$Group <- factor(plot_data$Group, levels = recovery_levels)
plot_data$Continent <- factor(plot_data$Continent, levels = continent_levels)
plot_data <- plot_data[order(match(plot_data$Dataset, c("Study", "External")),
                             plot_data$Continent, plot_data$Group, plot_data$SampleID), ]

summarize_recovery <- function(data, by) {
  groups <- interaction(data[by], drop = TRUE, lex.order = TRUE)
  do.call(rbind, lapply(split(data, groups), function(part) {
    cbind(part[1L, by, drop = FALSE],
          N = nrow(part), Mean_pct = mean(part$ratio), Median_pct = median(part$ratio),
          Q25_pct = unname(quantile(part$ratio, 0.25, type = 7)),
          Q75_pct = unname(quantile(part$ratio, 0.75, type = 7)),
          Min_pct = min(part$ratio), Max_pct = max(part$ratio))
  }))
}
summary_statistics <- summarize_recovery(plot_data, c("Dataset", "Group"))
continent_statistics <- summarize_recovery(plot_data, c("Dataset", "Continent", "Group"))
#write_tsv(plot_data, "results/Fig3b_plot_data.tsv")
#write_tsv(summary_statistics, "results/Fig3b_summary_statistics.tsv")
#write_tsv(continent_statistics, "results/Fig3b_continent_statistics.tsv")

# Grouped boxplots and sample points as the preprint did.
make_recovery_plot <- function(data, title) {
  samples <- unique(data[c("SampleID", "Continent")])
  counts <- table(factor(samples$Continent, levels = continent_levels))
  labels <- setNames(paste0(names(counts), "\n(n=", as.integer(counts), ")"), names(counts))
  ggplot(data, aes(x = Continent, y = ratio, fill = Group,
                   group = interaction(Continent, Group))) +
    geom_boxplot(width = 0.78, alpha = 0.7, outlier.shape = NA,
                 linewidth = 0.45, position = position_dodge(width = 0.86)) +
    geom_point(shape = 21, size = 1.65, stroke = 0.3, alpha = 0.6,
               position = position_jitterdodge(jitter.width = 0.22, dodge.width = 0.86, seed = 20261111)) +
    scale_fill_manual(values = recovery_colors, breaks = recovery_levels, drop = FALSE) +
    scale_x_discrete(labels = labels, drop = TRUE) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20),
                       expand = expansion(mult = c(0.01, 0.025))) +
    labs(x = NULL, y = "Microbial community\nrecovery estimation (%)", fill = NULL, title = title) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_bw(base_size = 22) +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.3),
      axis.text = element_text(color = "black"), axis.text.x = element_text(size = 19),
      plot.title = element_text(size = 24, hjust = 0.5, margin = margin(b = 14)),
      legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
      legend.title = element_blank(), legend.text = element_text(size = 16),
      legend.key.width = grid::unit(0.5, "cm"), legend.key.height = grid::unit(0.5, "cm"),
      legend.spacing.x = grid::unit(0.15, "cm"), legend.box.margin = margin(t = 12),
      legend.background = element_blank(), legend.key = element_blank(),
      plot.margin = margin(12, 12, 6, 12)
    )
}
p_study <- make_recovery_plot(plot_data[plot_data$Dataset == "Study", ], "Datasets from this study")
p_external <- make_recovery_plot(plot_data[plot_data$Dataset == "External", ], "External datasets")
p_external_shared <- p_external + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())

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
    vjust = -3,
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
# Assemble Figure 3
# ============================================================

# Panel b keeps its own shared legend within the full figure.
p_recovery <- (p_study | p_external_shared) +
  plot_layout(widths = c(5, 6), guides = "collect") &
  theme(legend.position = "bottom")

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

top_left_block <- label_a / p_continents +
  plot_layout(heights = c(0.04, 1))
top_right_block <- label_c / p_ge +
  plot_layout(heights = c(0.04, 1))
p_top <- top_left_block | top_right_block
bottom_block <- label_b / p_recovery +
  plot_layout(heights = c(0.04, 1))
p_combined <- (p_top / bottom_block) +
  plot_layout(heights = c(1.5, 2.7))

# ============================================================
# Export
# ============================================================

ggsave(
  filename = output_figure,
  plot = p_combined,
  width = 24,
  height = 16,
  device = grDevices::cairo_pdf,
  bg = "white"
)
writeLines(capture.output(sessionInfo()), "../envs/Fig3.sessionInfo.txt")
