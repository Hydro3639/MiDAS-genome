############################################################
## Figure S10 — WWTP-level distribution of nitrite oxidizers
##
## Description:
##   This script generates a discrete heatmap showing the
##   relative abundance of seven target nitrite-oxidizing
##   genera across individual wastewater treatment plants.
##
##   Relative abundance is summed at the genus level within
##   each WWTP. Missing genus–WWTP combinations are assigned
##   a relative abundance of zero.
##
## Inputs:
##   - sources/Fig5b1.tsv
##       Species-level abundance table containing:
##       gtdb.r226, abundance, and Sample
##
##   - sources/Fig5b2.tsv
##       Tab-delimited WWTP metadata without a header:
##       key, continent, country, and city
##
## Output:
##   - figures/FigS10.pdf
##
## Notes:
##   - Relative abundance values are assumed to be expressed
##     as percentages (e.g. 2.3 represents 2.3%).
##   - WWTPs are ordered by continent, country, and sample key.
##   - Genera are ordered by their mean relative abundance
##     across all WWTPs.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggtext)
})

# ----------------------------- #
# Input and output paths
# ----------------------------- #

abundance_file <- "sources/Fig5b1.tsv"
metadata_file  <- "sources/Fig5b2.tsv"
output_file    <- "figures/FigS10.pdf"

if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}

required_files <- c(
  abundance_file,
  metadata_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input file(s): ",
    paste(missing_files, collapse = ", ")
  )
}

# ----------------------------- #
# Target genera
# ----------------------------- #

target_genera <- c(
  "g__BQWY01",
  "g__JAHBWI01",
  "g__Nitrospira_A",
  "g__Nitrospira_D",
  "g__Nitrospira_G",
  "g__Nitrospira_C",
  "g__Nitrotoga"
)

# ----------------------------- #
# Read and validate abundance data
# ----------------------------- #

abundance_data <- read.table(
  abundance_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_abundance_columns <- c(
  "gtdb.r226",
  "abundance",
  "Sample"
)

missing_abundance_columns <- setdiff(
  required_abundance_columns,
  colnames(abundance_data)
)

if (length(missing_abundance_columns) > 0) {
  stop(
    "Missing required columns in ",
    abundance_file,
    ": ",
    paste(missing_abundance_columns, collapse = ", ")
  )
}

abundance_data <- abundance_data %>%
  mutate(
    Sample = trimws(Sample),
    gtdb.r226 = trimws(gtdb.r226),
    abundance = as.numeric(abundance),
    
    # Extract genus-level GTDB taxonomy while retaining
    # suffixes such as Nitrospira_A, Nitrospira_D, and
    # Nitrospira_G.
    genus = stringr::str_extract(
      gtdb.r226,
      "g__[^;\\s]+"
    ),
    genus = trimws(genus)
  )

if (anyNA(abundance_data$abundance)) {
  warning(
    "Some abundance values could not be converted to numeric ",
    "and will be treated as missing values."
  )
}

# ----------------------------- #
# Read and validate WWTP metadata
# ----------------------------- #

geo_metadata <- read.table(
  metadata_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (ncol(geo_metadata) != 4) {
  stop(
    metadata_file,
    "must contain exactly four columns: ",
    "key, continent, country, and city."
  )
}

colnames(geo_metadata) <- c(
  "key",
  "continent",
  "country",
  "city"
)

geo_metadata <- geo_metadata %>%
  mutate(
    across(everything(), trimws)
  )

# ----------------------------- #
# Construct metadata for all WWTPs
# ----------------------------- #

# Derive the complete WWTP list from the full abundance table.
# This retains WWTPs where none of the target genera were found.
wwtp_metadata <- abundance_data %>%
  distinct(WWTP = Sample) %>%
  mutate(
    # Convert identifiers such as 076-IN-18 to IN-18.
    key = sub("^[0-9]+-", "", WWTP),
    key = trimws(key)
  ) %>%
  left_join(
    geo_metadata,
    by = "key"
  ) %>%
  mutate(
    continent = ifelse(
      is.na(continent) | continent == "",
      "Unknown",
      continent
    ),
    country = ifelse(
      is.na(country) | country == "",
      "Unknown",
      country
    )
  ) %>%
  distinct(
    WWTP,
    key,
    continent,
    country,
    city
  )

unmatched_wwtps <- wwtp_metadata %>%
  filter(
    continent == "Unknown" |
      country == "Unknown"
  ) %>%
  pull(WWTP)

if (length(unmatched_wwtps) > 0) {
  warning(
    "Geographic metadata could not be resolved for: ",
    paste(unmatched_wwtps, collapse = ", ")
  )
}

# ----------------------------- #
# Sum relative abundance by genus and WWTP
# ----------------------------- #

genus_sum <- abundance_data %>%
  filter(
    !is.na(genus),
    genus != "",
    genus %in% target_genera
  ) %>%
  group_by(
    genus,
    WWTP = Sample
  ) %>%
  summarise(
    RA = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------- #
# Define WWTP order
# ----------------------------- #

continent_levels <- c(
  "Asia",
  "Europe",
  "North America",
  "South America",
  "Oceania",
  "Africa",
  "Unknown"
)

wwtp_metadata <- wwtp_metadata %>%
  mutate(
    continent = factor(
      continent,
      levels = continent_levels
    )
  ) %>%
  arrange(
    continent,
    country,
    key
  )

wwtp_order <- wwtp_metadata$WWTP

wwtp_metadata <- wwtp_metadata %>%
  mutate(
    WWTP = factor(
      WWTP,
      levels = wwtp_order
    )
  )

# ----------------------------- #
# Complete the genus × WWTP grid
# ----------------------------- #

# Every target genus is represented for every WWTP. Missing
# combinations are assigned a relative abundance of zero.
plot_data <- tidyr::expand_grid(
  genus = target_genera,
  WWTP = wwtp_order
) %>%
  left_join(
    genus_sum,
    by = c(
      "genus",
      "WWTP"
    )
  ) %>%
  left_join(
    wwtp_metadata %>%
      mutate(WWTP = as.character(WWTP)),
    by = "WWTP"
  ) %>%
  mutate(
    RA = replace_na(RA, 0),
    WWTP = factor(
      WWTP,
      levels = wwtp_order
    )
  )

# ----------------------------- #
# Order genera by mean abundance
# ----------------------------- #

genus_order <- plot_data %>%
  group_by(genus) %>%
  summarise(
    mean_RA = mean(RA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_RA)) %>%
  pull(genus)

# Reverse the factor levels so that the most abundant genus
# appears at the top of the heatmap.
plot_data <- plot_data %>%
  mutate(
    genus = factor(
      genus,
      levels = rev(genus_order)
    )
  )

# ----------------------------- #
# Define abundance bins and labels
# ----------------------------- #

plot_data <- plot_data %>%
  mutate(
    RA_bin = case_when(
      RA > 1    ~ ">1%",
      RA > 0.1  ~ "0.1%-1%",
      RA > 0.01 ~ "0.01%-0.1%",
      RA > 0    ~ "<0.01%",
      TRUE      ~ NA_character_
    ),
    RA_bin = factor(
      RA_bin,
      levels = c(
        ">1%",
        "0.1%-1%",
        "0.01%-0.1%",
        "<0.01%"
      )
    ),
    RA_label = case_when(
      RA == 0 ~ "",
      TRUE ~ sprintf(
        "%.1f",
        round(RA, 1)
      )
    )
  )

# ----------------------------- #
# Fill colors
# ----------------------------- #

fill_values <- c(
  ">1%"        = "#d7301f",
  "0.1%-1%"    = "#fc8d59",
  "0.01%-0.1%" = "#fdd49e",
  "<0.01%"     = "#fef0d9"
)

# ----------------------------- #
# Construct continent-based axis-label colors
# ----------------------------- #

# Assign alternating red and black colors to consecutive
# continent blocks.
continent_color_table <- wwtp_metadata %>%
  arrange(
    continent,
    country,
    key
  ) %>%
  distinct(continent) %>%
  mutate(
    label_color = rep(
      c("#b2182b", "black"),
      length.out = n()
    )
  )

wwtp_label_data <- wwtp_metadata %>%
  left_join(
    continent_color_table,
    by = "continent"
  ) %>%
  arrange(WWTP)

wwtp_color_vector <- wwtp_label_data$label_color
names(wwtp_color_vector) <- as.character(
  wwtp_label_data$WWTP
)

axis_label_colors <- wwtp_color_vector[wwtp_order]
axis_label_colors[is.na(axis_label_colors)] <- "black"

# Use HTML labels so that each WWTP label can be colored
# independently and rendered consistently by ggtext.
wwtp_axis_labels <- setNames(
  paste0(
    "<span style='color:",
    axis_label_colors,
    ";'>",
    wwtp_order,
    "</span>"
  ),
  wwtp_order
)

# ----------------------------- #
# Generate Figure S10
# ----------------------------- #

p_figS10 <- ggplot(
  plot_data,
  aes(
    x = WWTP,
    y = genus
  )
) +
  geom_tile(
    aes(fill = RA_bin),
    color = "white",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = RA_label),
    size = 2
  ) +
  scale_fill_manual(
    values = fill_values,
    name = "Relative abundance",
    drop = FALSE,
    na.translate = FALSE,
    na.value = "white"
  ) +
  scale_x_discrete(
    labels = wwtp_axis_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub("^g__", "", x)
    },
    drop = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    axis.text.x = ggtext::element_markdown(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7
    ),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

p_figS10

# ----------------------------- #
# Export
# ----------------------------- #

ggsave(
  filename = output_file,
  plot = p_figS10,
  width = 16,
  height = 2.9,
  useDingbats = FALSE
)

# ----------------------------- #
# Record the R environment
# ----------------------------- #

if (!dir.exists("../envs")) {
  dir.create(
    "../envs",
    recursive = TRUE
  )
}

writeLines(
  capture.output(sessionInfo()),
  "../envs/FigS10.sessionInfo.txt"
)
