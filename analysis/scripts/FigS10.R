# =========================
# Discrete (binned) heatmap with numbers in each cell
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})
setwd("~/Library/CloudStorage/Dropbox/LL-dropbox/AAU/002-After-Arrival-AAU/001-MiDAS-Project/005-Paper-prep/000-Final-version-83samples/010-Main-SI-figures/FigS9-Novel_NOB.abun.dist.indi.wwtp")
abun_file <- "7.genera.abun.tsv"
geo_file  <- "geo.info.tsv"
out_pdf   <- "genus_heatmap.binned.indi.wwtp.pdf"

# ---------- load abundance ----------
df <- read.table(abun_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    Sample = trimws(Sample),
    gtdb.r226 = trimws(gtdb.r226),
    abundance = as.numeric(abundance)
  )

# ---------- extract genus (keeps Nitrospira_A/D/G etc.) ----------
df <- df %>%
  mutate(
    genus = str_extract(gtdb.r226, "g__[^;\\s]+"),
    genus = trimws(genus)
  )

# ---------- sum RA within genus per WWTP ----------
genus_sum <- df %>%
  group_by(genus, WWTP = Sample) %>%
  summarise(RA = sum(abundance, na.rm = TRUE), .groups = "drop")

# ---------- optional: keep only target genera ----------
target_genera <- c(
  "g__BQWY01",
  "g__JAHBWI01",
  "g__Nitrospira_A",
  "g__Nitrospira_D",
  "g__Nitrospira_G",
  "g__Nitrospira_C",
  "g__Nitrotoga"
)

genus_sum <- genus_sum %>%
  filter(genus %in% target_genera)

# ---------- load geo info ----------
geo <- read.table(geo_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(geo) <- c("key", "continent", "country", "city")
geo <- geo %>% mutate(across(everything(), trimws))

# ---------- attach geo metadata ----------
genus_sum_geo <- genus_sum %>%
  mutate(
    key = sub("^[0-9]+-", "", WWTP),   # 076-IN-18 -> IN-18
    key = trimws(key)
  ) %>%
  left_join(geo, by = "key") %>%
  mutate(
    continent = ifelse(is.na(continent) | continent == "", "Unknown", continent),
    country   = ifelse(is.na(country)   | country   == "", "Unknown", country)
  )

# ---------- order continents ----------
continent_levels <- c("Asia","Europe","North America","South America","Oceania","Africa","Unknown")
genus_sum_geo$continent <- factor(genus_sum_geo$continent, levels = continent_levels)

# ---------- order WWTPs by continent -> country -> key ----------
wwtp_order <- genus_sum_geo %>%
  distinct(WWTP, continent, country, key) %>%
  arrange(continent, country, key) %>%
  pull(WWTP)

genus_sum_geo$WWTP <- factor(genus_sum_geo$WWTP, levels = wwtp_order)

# ---------- order genera (optional) ----------
genus_levels <- c("g__JAHBWI01","g__BQWY01","g__Nitrospira_A","g__Nitrospira_D","g__Nitrotoga","g__Nitrospira_G", "g__Nitrospira_C")
genus_sum_geo$genus <- factor(genus_sum_geo$genus, levels = genus_levels)

# ---------- binning to match your legend ----------
# Note: assumes RA is already in percent units (e.g., 2.3 means 2.3%)
genus_sum_geo <- genus_sum_geo %>%
  mutate(
    RA_bin = case_when(
      RA > 1 ~ ">1%",
      RA > 0.1 ~ "0.1%-1%",
      RA > 0.01 ~ "0.01%-0.1%",
      RA > 0 ~ "<0.01%",
      TRUE ~ NA_character_
    ),
    RA_bin = factor(RA_bin, levels = c(">1%","0.1%-1%","0.01%-0.1%","<0.01%")),
    RA_label = case_when(
      RA == 0 ~ "",
      TRUE ~ sprintf("%.1f", round(RA,1))
    )
  )

# ---------- make a complete tile grid (so zeros show as empty tiles) ----------
# This ensures every genus × WWTP combination exists.
plot_df <- genus_sum_geo %>%
  select(genus, WWTP, RA, RA_bin, RA_label, continent, country, key) %>%
  tidyr::complete(genus, WWTP, fill = list(RA = 0, RA_bin = NA, RA_label = ""))


# ---------- compute average RA per genus within each country ----------
genus_levels <- plot_df %>%
  group_by(genus) %>%
  summarise(avg_RA = mean(RA, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(avg_RA)) %>%
  pull(genus)

plot_df$genus <- factor(plot_df$genus, levels = rev(genus_levels))

## ordered by average abun
plot_df$genus <- factor(plot_df$genus, levels = rev(genus_levels))

# ---------- discrete colors (close to your legend style) ----------
# You can tweak these if you want a specific palette.
fill_vals <- c(
  ">1%" = "#d7301f",
  "0.1%-1%" = "#fc8d59",
  "0.01%-0.1%" = "#fdd49e",
  "<0.01%" = "#fef0d9"
)

# ---------- continent separators (optional but nice) ----------
wwtp_meta <- plot_df %>%
  distinct(WWTP) %>%
  mutate(idx = as.integer(WWTP)) %>%
  left_join(
    genus_sum_geo %>% distinct(WWTP, continent) %>% mutate(WWTP = factor(WWTP, levels = wwtp_order)),
    by = "WWTP"
  ) %>%
  arrange(idx)

cont_breaks <- wwtp_meta %>%
  group_by(continent) %>%
  summarise(last_idx = max(idx), .groups = "drop") %>%
  arrange(last_idx)

# ---------- plot ----------

# add continent separators
#if (nrow(cont_breaks) > 1) {
#  p <- p + geom_vline(
#    xintercept = cont_breaks$last_idx[-nrow(cont_breaks)] + 0.5,
#    color = "black",
#    linewidth = 0.25
#  )
#}

# Create continent-level alternating colors
continent_df <- genus_sum_geo %>%
  distinct(WWTP, continent) %>%
  arrange(WWTP) %>%
  distinct(continent) %>%
  mutate(cont_color = rep(c("#b2182b","black"), length.out = n()))

# Join back to get color per WWTP
wwtp_colors <- genus_sum_geo %>%
  distinct(WWTP, continent) %>%
  left_join(continent_df, by = "continent") %>%
  arrange(WWTP)

p <- ggplot(plot_df, aes(x = WWTP, y = genus)) +
  geom_tile(aes(fill = RA_bin), color = "white", linewidth = 0.2) +
  geom_text(aes(label = RA_label), size = 2) +
  scale_fill_manual(
    values = fill_vals,
    name = "Relative abundance",
    drop = FALSE,
    na.translate = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7,
      color = wwtp_colors$cont_color   # <-- here is the trick
    )
  )

p <- p +
  scale_y_discrete(labels = function(x) sub("^g__", "", x))

ggsave(out_pdf, p, width = 16, height = 2.9, useDingbats = FALSE)
print(p)




