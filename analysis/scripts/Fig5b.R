############################################################
## Figure 5b - Country/Region-level distribution of nitrite oxidizers
##
## Description:
##   This script generates a heatmap showing the mean relative
##   abundance of seven target nitrite-oxidizing genera (NOB)
##   across countries/regions.
##
##   Relative abundance is first summed at the genus level for
##   each WWTP and then averaged across all WWTPs within each
##   country. A genus not detected in a WWTP is assigned a
##   relative abundance of zero before country-level averaging.
##
## Inputs:
##   - sources/Fig5b1.tsv
##       Species-level abundance table with the columns:
##       species, gtdb.r226, abundance, Sample
##
##   - sources/Fig5b2.tsv
##       Tab-delimited WWTP metadata without a header:
##       key, continent, country, city
##
## Output:
##   - figures/Fig5b.pdf
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
output_file    <- "figures/Fig5b.pdf"

if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}

required_files <- c(abundance_file, metadata_file)
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
  "species",
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
    "Missing columns in ",
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
    
    # Extract the genus-level GTDB taxonomy while retaining
    # suffixes such as Nitrospira_A, Nitrospira_C, and Nitrospira_D.
    genus = stringr::str_extract(
      gtdb.r226,
      "g__[^;\\s]+"
    ),
    genus = trimws(genus)
  )

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
    " must contain exactly four columns: ",
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

# The complete WWTP list is derived from the full abundance table,
# rather than only from WWTPs where a target genus was detected.
wwtp_geo <- abundance_data %>%
  distinct(WWTP = Sample) %>%
  mutate(
    # Convert identifiers such as 076-IN-18 to IN-18.
    key = sub("^[0-9]+-", "", WWTP),
    key = trimws(key)
  ) %>%
  left_join(geo_metadata, by = "key") %>%
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

# Warn when abundance samples cannot be matched to metadata.
unmatched_wwtps <- wwtp_geo %>%
  filter(continent == "Unknown" | country == "Unknown") %>%
  pull(WWTP)

if (length(unmatched_wwtps) > 0) {
  warning(
    "WWTP metadata could not be resolved for: ",
    paste(unmatched_wwtps, collapse = ", ")
  )
}

# ----------------------------- #
# Sum abundance by genus and WWTP
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
# Complete the WWTP × genus grid
# ----------------------------- #

# Missing genus–WWTP combinations are assigned RA = 0.
wwtp_genus_full <- tidyr::expand_grid(
  WWTP = unique(wwtp_geo$WWTP),
  genus = target_genera
) %>%
  left_join(
    wwtp_geo %>%
      select(WWTP, continent, country),
    by = "WWTP"
  ) %>%
  left_join(
    genus_sum,
    by = c("WWTP", "genus")
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
    ),
    RA = replace_na(RA, 0)
  )

# ----------------------------- #
# Calculate country-level means
# ----------------------------- #

# Mean relative abundance is calculated across all WWTPs in
# each country, including WWTPs in which the genus is absent.
country_mean <- wwtp_genus_full %>%
  group_by(
    genus,
    continent,
    country
  ) %>%
  summarise(
    RA = mean(RA, na.rm = TRUE),
    n_wwtp = n_distinct(WWTP),
    .groups = "drop"
  )

# ----------------------------- #
# Define plotting order
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

country_mean <- country_mean %>%
  mutate(
    continent = factor(
      continent,
      levels = continent_levels
    )
  )

# Countries are ordered first by continent and then alphabetically.
country_order <- country_mean %>%
  distinct(continent, country) %>%
  arrange(continent, country) %>%
  pull(country)

# Genera are initially ordered by mean abundance across countries.
genus_order <- country_mean %>%
  group_by(genus) %>%
  summarise(
    mean_RA = mean(RA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_RA)) %>%
  pull(genus)

# ----------------------------- #
# Complete the country × genus grid
# ----------------------------- #

country_to_continent <- country_mean %>%
  distinct(country, continent)

plot_data <- country_mean %>%
  select(
    genus,
    country,
    RA,
    n_wwtp
  ) %>%
  tidyr::complete(
    genus,
    country,
    fill = list(
      RA = 0,
      n_wwtp = 0
    )
  ) %>%
  left_join(
    country_to_continent,
    by = "country"
  ) %>%
  mutate(
    continent = ifelse(
      is.na(continent),
      "Unknown",
      as.character(continent)
    ),
    continent = factor(
      continent,
      levels = continent_levels
    ),
    country = factor(
      country,
      levels = country_order
    ),
    genus = factor(
      genus,
      levels = genus_order
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
    RA_round = round(RA, 1),
    RA_label = case_when(
      # True absence is shown as a blank cell.
      RA == 0 ~ "",
      
      # Values greater than zero that round to 0.0 are shown as 0.
      RA_round == 0 ~ "0",
      
      TRUE ~ sprintf("%.1f", RA_round)
    )
  )

# Recalculate genus order after completing the plotting matrix.
genus_order <- plot_data %>%
  group_by(genus) %>%
  summarise(
    mean_RA = mean(RA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_RA)) %>%
  pull(genus)

plot_data <- plot_data %>%
  mutate(
    genus = factor(
      genus,
      levels = rev(genus_order)
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
# Construct country-axis labels
# ----------------------------- #

# Count the total number of WWTPs represented in each country.
country_wwtp_counts <- wwtp_geo %>%
  distinct(country, WWTP) %>%
  count(
    country,
    name = "n_total"
  )

country_levels <- levels(plot_data$country)

country_n <- country_wwtp_counts$n_total[
  match(
    country_levels,
    country_wwtp_counts$country
  )
]

country_n[is.na(country_n)] <- 0

# Display the country on the first line and the WWTP count
# on the second line.
country_labels_plain <- setNames(
  paste0(
    country_levels,
    "<br/><span style='font-size:9pt;'>(n=",
    country_n,
    ")</span>"
  ),
  country_levels
)

# ----------------------------- #
# Assign continent-based label colors
# ----------------------------- #

# Alternate red and black across consecutive continent blocks.
continent_color_map <- tibble(
  continent = factor(
    continent_levels,
    levels = continent_levels
  ),
  label_color = rep(
    c("#b2182b", "black"),
    length.out = length(continent_levels)
  )
)

country_colors <- plot_data %>%
  distinct(country, continent) %>%
  mutate(
    country_character = as.character(country)
  ) %>%
  left_join(
    continent_color_map,
    by = "continent"
  )

country_color_vector <- country_colors$label_color
names(country_color_vector) <- country_colors$country_character

axis_label_colors <- country_color_vector[country_levels]
axis_label_colors[is.na(axis_label_colors)] <- "black"

country_labels <- setNames(
  paste0(
    "<span style='color:",
    axis_label_colors,
    ";'>",
    country_labels_plain[country_levels],
    "</span>"
  ),
  country_levels
)

# ----------------------------- #
# Generate Figure 5b
# ----------------------------- #

p_fig5b <- ggplot(
  plot_data,
  aes(
    x = country,
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
    size = 3
  ) +
  scale_fill_manual(
    values = fill_values,
    name = "Relative abundance",
    drop = FALSE,
    na.translate = FALSE,
    na.value = "white"
  ) +
  scale_x_discrete(
    labels = country_labels
  ) +
  scale_y_discrete(
    labels = function(x) sub("^g__", "", x)
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(
      size = 10,
      color = "black",
      face = "plain"
    ),
    axis.text.x = ggtext::element_markdown(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10,
      face = "plain"
    ),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

p_fig5b

# ----------------------------- #
# Export
# ----------------------------- #

ggsave(
  filename = output_file,
  plot = p_fig5b,
  width = 10,
  height = 4.1,
  useDingbats = FALSE
)

# ----------------------------- #
# Record the R environment
# ----------------------------- #

if (!dir.exists("../envs")) {
  dir.create("../envs", recursive = TRUE)
}

writeLines(
  capture.output(sessionInfo()),
  "../envs/Fig5b.sessionInfo.txt"
)
