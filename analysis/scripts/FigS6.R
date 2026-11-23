############################################################
## Figure S6 — Community recovery from individual sample catalogs
##
## Description:
##   This script compares microbial community recovery estimates
##   obtained using individual-sample HQ and HQ+MQ species
##   representatives across globally collected WWTP samples.
##
## Input:
##   - sources/FigS3.xlsx
##       Sheet: GTDB.r226
##
## Required columns:
##   - SA.ID
##   - indi_hq_ratio
##   - indi.all_ratio
##   - Continent
##
## Output:
##   - figures/FigS6.pdf
##
## Notes:
##   - The first 83 data rows are used, as in the original analysis.
##   - The script should be run from the analysis directory.
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(readxl)
  library(grid)
})

# ----------------------------- #
# Input and output paths
# ----------------------------- #
# the input is Fig3b.xlsx

input_file  <- "sources/Fig3b.xlsx"
sheet_name  <- "GTDB.r226"
output_file <- "figures/FigS6.pdf"

if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# ----------------------------- #
# Read and validate input data
# ----------------------------- #

available_sheets <- readxl::excel_sheets(input_file)

if (!sheet_name %in% available_sheets) {
  stop(
    "Sheet '",
    sheet_name,
    "' was not found in ",
    input_file,
    ". Available sheets: ",
    paste(available_sheets, collapse = ", ")
  )
}

raw_data <- readxl::read_excel(
  input_file,
  sheet = sheet_name,
  skip = 1
)

required_columns <- c(
  "SA.ID",
  "indi_hq_ratio",
  "indi.all_ratio",
  "Continent"
)

missing_columns <- setdiff(
  required_columns,
  colnames(raw_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns in ",
    input_file,
    ": ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(raw_data) < 83) {
  stop(
    input_file,
    " contains only ",
    nrow(raw_data),
    " data rows; at least 83 rows are required."
  )
}

# Retain the first 83 samples and the columns used in Figure S6.
selected_data <- raw_data[
  1:83,
  required_columns
]

# ----------------------------- #
# Reshape data for plotting
# ----------------------------- #

plot_data <- selected_data %>%
  tidyr::pivot_longer(
    cols = c(
      "indi_hq_ratio",
      "indi.all_ratio"
    ),
    names_to = "dataset",
    values_to = "Ratio"
  )

plot_data$dataset <- factor(
  plot_data$dataset,
  levels = c(
    "indi_hq_ratio",
    "indi.all_ratio"
  ),
  labels = c(
    "Individual-HQ",
    "Individual-HQ+MQ"
  )
)

# Preserve the original dataset order and color mapping.
dataset_colors <- c(
  "Individual-HQ"    = "#C2A5CF",
  "Individual-HQ+MQ" = "#7B3294"
)

# Sample numbers are fixed to the 83 globally collected WWTP
# can also be extracted automatically, but not necessary
# samples represented in the input table.
continent_labels <- c(
  "Asia"          = "Asia\n(n=18)",
  "Europe"        = "Europe\n(n=30)",
  "North America" = "North America\n(n=21)",
  "Oceania"       = "Oceania\n(n=7)",
  "South America" = "South America\n(n=7)"
)

# ----------------------------- #
# Generate Figure S6
# ----------------------------- #

p_figs6 <- ggplot(
  plot_data,
  aes(
    x = Continent,
    y = Ratio,
    fill = dataset
  )
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
      jitter.width = 0.4,
      dodge.width = 0.9
    )
  ) +
  scale_fill_manual(
    values = dataset_colors,
    breaks = names(dataset_colors)
  ) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, by = 15)
  ) +
  scale_x_discrete(
    labels = continent_labels
  ) +
  labs(
    x = "",
    y = "Microbial community\nrecovery estimation (%)",
    fill = "Dataset"
  ) +
  theme_bw(base_size = 26) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 26),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
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

p_figs6

# ----------------------------- #
# Export
# ----------------------------- #

ggsave(
  filename = output_file,
  plot = p_figs6,
  width = 12,
  height = 8
)

# ----------------------------- #
# Record the R environment
# ----------------------------- #

if (!dir.exists("../envs")) {
  dir.create("../envs", recursive = TRUE)
}

writeLines(
  capture.output(sessionInfo()),
  "../envs/FigS6.sessionInfo.txt"
)
