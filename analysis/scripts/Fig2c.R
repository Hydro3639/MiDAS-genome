############################################################
## Fig. 2c — Classified vs unclassified HQ species across taxonomic ranks
##
## Description:
##   This figure summarizes the number of classified and unclassified
##   high-quality (HQ) species representatives across taxonomic ranks
##   (phylum to species).
##
## Data source:
##   Numbers were calculated based on the GTDB r226 taxonomy file
##   for HQ species representatives derived from 83 activated sludge samples.
##
## Notes:
##   - "Unclassified" refers to taxa lacking an assigned name at the
##     corresponding taxonomic rank in GTDB r226.
##   - Percentages are calculated relative to the total number of HQ
##     species representatives (n = 12,047).
##
## Output:
##   - Fig2c.pdf
############################################################

library(ggplot2)
library(tidyverse)
library(scales)

# ----------------------------- #
# Input data (precomputed counts)
# ----------------------------- #

df <- data.frame(
  Category = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  Classified = c(12047, 12044, 12015, 11844, 10337, 2135),
  Unclassified = c(0, 3, 32, 203, 1710, 9912)
)

# Ensure correct taxonomic order
df$Category <- factor(
  df$Category,
  levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
)

# Long format
df_long <- df %>%
  pivot_longer(
    cols = c(Classified, Unclassified),
    names_to = "Variable",
    values_to = "Value"
  )

# ----------------------------- #
# Plot: percentage of unclassified taxa
# ----------------------------- #

df_unclassified <- df_long %>%
  filter(Variable == "Unclassified")

p_unclassified <- ggplot(
  df_unclassified,
  aes(x = Category, y = Value / 12047 * 100)
) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    fill = "#45818e",
    color = "black"
  ) +
  geom_text(
    aes(
      label = scales::comma(Value),
      y = Value / 12047 * 100 + 4
    ),
    color = "black",
    size = 6
  ) +
  labs(
    x = NULL,
    y = "Percentage of novel species (%)"
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    text = element_text(size = 24),
    legend.position = "none"
  )

p_unclassified

ggsave(
  p_unclassified,
  filename = "figures/Fig2c.pdf",
  width = 10,
  height = 6
)

writeLines(capture.output(sessionInfo()), "../envs/Fig2c.sessionInfo.txt")
