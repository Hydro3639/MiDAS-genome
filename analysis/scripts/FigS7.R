############################################################
## Fig S7 — Top 10 phyla by recovered MAGs (top) and spReps (bottom, mirrored)
##
## Input:
##   - sources/FigS7.tsv (tab-delimited; header required)
##     Expected columns:
##       phylum, NumberofMAGs, NumberofspRep
##
## Output:
##   - figures/FigS7.pdf
##
## Notes:
##   - the dataset was extracted from Supplementary files
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(scales)
  library(dplyr)
  library(patchwork)
})

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

# ----------------------------- #
# 1) Read input file
# ----------------------------- #
data_raw <- read.table(
  "sources/FigS7.tsv",
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("phylum", "NumberofMAGs", "NumberofspRep")
missing_cols <- setdiff(required_cols, colnames(data_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sources/FigS7.tsv: ", paste(missing_cols, collapse = ", "))
}

# ----------------------------- #
# 2) Clean and standardize columns
# ----------------------------- #
data_clean <- data_raw %>%
  rename(
    phylum_raw = phylum,
    Genome_number = NumberofMAGs,
    species_number = NumberofspRep
  ) %>%
  mutate(
    Phylum = sub("^p__", "", phylum_raw),
    Phylum = trimws(Phylum)
  )

# ----------------------------- #
# 3) Select top 10 by recovered MAGs
# ----------------------------- #
top10 <- data_clean %>%
  arrange(desc(Genome_number)) %>%
  slice(1:10)

# Keep x-order identical to the selected table order
phylum_levels <- top10$Phylum
top10 <- top10 %>%
  mutate(Phylum = factor(Phylum, levels = phylum_levels))

# Dynamic y-limits
ymax_genome <- ceiling(max(top10$Genome_number, na.rm = TRUE) * 1.05)
ymax_species <- ceiling(max(top10$species_number, na.rm = TRUE) * 1.25)

# ----------------------------- #
# 4) Top panel: number of MAGs
# ----------------------------- #
p_phylum.genome <- ggplot(top10, aes(x = Phylum, y = Genome_number, fill = Phylum)) +
  geom_col() +
  geom_text(
    aes(label = scales::comma(Genome_number)),
    vjust = -0.5,
    color = "black",
    size = 6
  ) +
  scale_fill_aaas() +
  coord_cartesian(ylim = c(0, ymax_genome)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "", y = "Number of MAGs", title = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    text = element_text(size = 22),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )

# ----------------------------- #
# 5) Bottom panel: number of species (mirrored)
# ----------------------------- #
top10 <- top10 %>%
  mutate(species_number_neg = -species_number)

p_phylum.species <- ggplot(top10, aes(x = Phylum, y = species_number_neg, fill = Phylum)) +
  geom_col() +
  geom_text(
    aes(label = scales::comma(species_number)),
    vjust = 1.5,
    color = "black",
    size = 6
  ) +
  scale_fill_aaas() +
  coord_cartesian(ylim = c(-ymax_species, 0)) +
  scale_y_continuous(labels = function(x) scales::comma(abs(x))) +
  labs(x = "Phylum", y = "Number of species") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    text = element_text(size = 22),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )

# Hide x-axis elements on the bottom panel
p_phylum.species <- p_phylum.species +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# ----------------------------- #
# 6) Combine and export
# ----------------------------- #
p_combined <- p_phylum.genome / p_phylum.species +
  plot_layout(ncol = 1, heights = c(2, 1))

print(p_combined)

ggsave("figures/FigS7.pdf", p_combined, width = 14, height = 10)


writeLines(capture.output(sessionInfo()), "../envs/FigS7.sessionInfo.txt")
