############################################################
## Fig S4 — Species improvement across bacterial lineages
##
## Input:
##   - sources/FigS4.txt (tab-delimited; header required)
##     Required columns:
##       Bacteria, Known_in_MiDAS, Num_of_first_HQ_rep,
##       impr_20per, impr_50per, impr_100per
##
## Output:
##   - figures/FigS4.pdf
##
## Notes:
##   - the dataset was extracted from Supplementary files
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggsci)
  library(ggpubr)
})

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

df <- read.table("sources/FigS4.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

required_cols <- c(
  "Bacteria", "Known_in_MiDAS", "Num_of_first_HQ_rep",
  "impr_20per", "impr_50per", "impr_100per"
)
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sources/FigS4.txt: ", paste(missing_cols, collapse = ", "))
}

# ordering
bacteria_levels <- rev(unique(df$Bacteria))

# ----------------------------- #
# Panel a: improvement proportions (faceted stacked bars)
# ----------------------------- #
df_long <- df %>%
  pivot_longer(
    cols = starts_with("impr_"),
    names_to = "improvement",
    values_to = "count"
  ) %>%
  mutate(
    improved  = count / Known_in_MiDAS,
    remaining = 1 - improved
  ) %>%
  pivot_longer(
    cols = c(improved, remaining),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(
    # For the "remaining" part, show the remaining count (Known_in_MiDAS - improved_count)
    count = ifelse(type == "remaining", Known_in_MiDAS - count, count),
    improvement = factor(
      improvement,
      levels = c("impr_20per", "impr_50per", "impr_100per"),
      labels = c("20%", "50%", "100%")
    ),
    type = factor(type, levels = c("remaining", "improved")),
    Bacteria = factor(Bacteria, levels = bacteria_levels)
  )

plot1 <- ggplot(df_long, aes(x = value, y = Bacteria, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = ifelse(value > 0.01, scales::comma(count), "")),
    position = position_stack(vjust = 0.5),
    size = 8,
    color = "black"
  ) +
  facet_wrap(~ improvement, nrow = 1, scales = "free_x") +
  scale_x_continuous(
    name = "Proportion",
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_y_discrete(name = NULL) +
  scale_fill_npg(
    name = "",
    labels = c("Unimproved lineages", "Lineages with improved species diversity")
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(hjust = 0.5, color = "black"),
    axis.text.y = element_text(hjust = 0.5, color = "black"),
    text = element_text(size = 28)
  )

# ----------------------------- #
# Panel b: First HQ representatives (log10 x-axis)
# ----------------------------- #
plot2 <- ggplot(df, aes(x = Num_of_first_HQ_rep, y = Bacteria)) +
  geom_bar(stat = "identity", fill = "#f4b3ad") +
  geom_text(
    aes(label = scales::comma(Num_of_first_HQ_rep)),
    hjust = 1.5,
    size = 8,
    color = "black"
  ) +
  scale_x_continuous(
    name = "First HQ representatives",
    trans = "log10"
  ) +
  scale_y_discrete(name = NULL, limits = bacteria_levels) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(hjust = 0.5, color = "black"),
    text = element_text(size = 28)
  ) +
  labs(title = "")

# ----------------------------- #
# Combine and export
# ----------------------------- #
blank_space <- ggplot() + theme_void()

right_col <- ggpubr::ggarrange(
  blank_space,
  plot2,
  ncol = 1,
  heights = c(0.2, 2.0),
  labels = c("", "b"),
  font.label = list(size = 24, face = "bold")
)

p_combined <- ggpubr::ggarrange(
  plot1,
  right_col,
  ncol = 2,
  widths = c(6, 2),
  labels = c("a"),
  label.x = 0.01,
  label.y = 0.9,
  font.label = list(size = 24, color = "black")
)

p_combined

ggsave(p_combined, filename = "figures/FigS4.pdf", width = 18, height = 9)

writeLines(capture.output(sessionInfo()), "../envs/FigS4.sessionInfo.txt")
