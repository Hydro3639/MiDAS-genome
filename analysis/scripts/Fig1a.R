############################################################
## Fig. 1a — Global distribution of MiDAS WWTPs and process types
##
## Input:
##   - sources/Fig1a.xlsx (sheet: "83samples")
##     Required columns:
##       WWTP_id, Continent, Country, City, Latitude, Longitude,
##       Process_type, Plant_type
##
## Output:
##   - figures/Fig1a.pdf
##
## Notes:
##   - Continent summary labels are computed from the input table.
##   - Legend labels (e.g., "C (17)") are computed from Process_type2.
##   - Label positions are fixed to match the manuscript layout.
############################################################

library(ggplot2)
library(tidyverse)
library(readxl)

# ----------------------------- #
# Load metadata
# ----------------------------- #

Geo_info <- read_excel(
  "sources/Fig1a.xlsx",
  sheet = "83samples"
)

# ----------------------------- #
# Column checks (fail early with a clear message)
# ----------------------------- #

required_cols <- c(
  "WWTP_id", "Continent", "Country", "City",
  "Latitude", "Longitude",
  "Process_type", "Plant_type"
)

missing_cols <- setdiff(required_cols, colnames(Geo_info))
if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Missing required columns in sources/Fig1a.xlsx: ",
      paste(missing_cols, collapse = ", "),
      "\nPlease check column names in the Excel sheet or update `required_cols`."
    )
  )
}

# ----------------------------- #
# Preprocess
# ----------------------------- #

df <- Geo_info %>%
  mutate(
    # Collapse rare process types into "Other"
    Process_type2 = fct_other(
      Process_type,
      keep = c("C", "C,N", "C,N,DN", "C,N,DN,P")
    ),
    Process_type2 = fct_explicit_na(Process_type2, na_level = "Other"),
    Plant_type = fct_explicit_na(as.factor(Plant_type), na_level = "Other")
  )

# ----------------------------- #
# Continent summary labels (computed)
# ----------------------------- #

continent_summary <- df %>%
  group_by(Continent) %>%
  summarise(
    n_country = n_distinct(Country),
    n_city    = n_distinct(City),
    n_sample  = n(),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      Continent, "\n",
      n_country, " countries\n",
      n_city, " cities\n",
      n_sample, " samples"
    )
  )

# Fixed label positions to preserve the manuscript layout
continent_pos <- tibble(
  Continent = c("Asia", "Europe", "North America", "Oceania", "South America"),
  x = c(65, -10, -170, 90, -105),
  y = c(50, 80, 28, -28, -10)
)

continent_labels <- continent_summary %>%
  inner_join(continent_pos, by = "Continent")

# ----------------------------- #
# Legend labels (computed from Process_type2)
# ----------------------------- #

process_levels <- c("C", "C,N", "C,N,DN", "C,N,DN,P", "Other")

process_counts <- df %>%
  mutate(Process_type2 = as.character(Process_type2)) %>%
  count(Process_type2, name = "n")

# Ensure all expected levels exist (fill missing with 0)
process_counts_complete <- tibble(Process_type2 = process_levels) %>%
  left_join(process_counts, by = "Process_type2") %>%
  mutate(n = replace_na(n, 0))

process_labels <- paste0(process_counts_complete$Process_type2, " (", process_counts_complete$n, ")")

# Sanity check: legend counts should sum to number of rows
stopifnot(sum(process_counts_complete$n) == nrow(df))

# ----------------------------- #
# Plot: world map + WWTP locations
# ----------------------------- #

MiDAS_geo_dis <-
  ggplot() +
  borders("world", fill = "#DAD7C5", colour = "grey60", alpha = 0.9) +
  ylim(-90, 90) +
  geom_point(
    data = df,
    aes(x = Longitude, y = Latitude, color = Process_type2),
    size = 4,
    stroke = 0.3
  ) +
  scale_color_manual(
    values = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#6A6599FF"),
    breaks = process_levels,
    labels = process_labels
  ) +
  labs(x = "Longitude", y = "Latitude", color = "Process type") +
  geom_label(
    data = continent_labels,
    aes(x = x, y = y, label = label),
    hjust = 0,
    size = 5,
    fill = "#f5f5dc",
    color = "black"
  ) +
  theme(
    legend.justification = c(1, 1),
    legend.margin = margin(t = 0.2, r = 0.2, b = 0.2, l = 0.2, unit = "cm"),
    legend.key = element_blank(),
    legend.key.size = unit(1.5, "lines"),
    legend.background = element_rect(fill = "#f5f5dc", colour = "black", linewidth = 0.4),
    text = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    legend.position = c(0.18, 0.35)
  )

MiDAS_geo_dis

# ----------------------------- #
# Export
# ----------------------------- #

ggsave(MiDAS_geo_dis, filename = "figures/Fig1a.pdf", width = 16, height = 10)

writeLines(
  capture.output(sessionInfo()),
  "../envs/Fig1a.sessionInfo.txt"
)