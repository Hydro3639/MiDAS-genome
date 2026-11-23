############################################################
## Fig S4 — Bacteria PD/PG summary (three aligned bar panels)
##
## Input:
##   - sources/FigS4.txt (tab-delimited; header required)
##     Required columns:
##       Clade, Percent_PD, gtdb_PD, midas_PD, midas_PG
##
## Output:
##   - figures/FigS4.pre.pdf
##
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(ggpubr)
})

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

data <- read.table("sources/FigS4.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

required_cols <- c("Clade", "Percent_PD", "gtdb_PD", "midas_PD", "midas_PG")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sources/FigS4.txt: ", paste(missing_cols, collapse = ", "))
}

clade_levels <- rev(unique(data$Clade))

# ----------------------------- #
# Panel 1: Total PD (Percent_PD)
# ----------------------------- #
data_selected_TPD <- data[, c("Clade", "Percent_PD")]
data_selected_TPD$Clade <- factor(data_selected_TPD$Clade, levels = clade_levels)

data_long1 <- tidyr::pivot_longer(
  data_selected_TPD,
  cols = "Percent_PD",
  names_to = "variable",
  values_to = "value"
)

p1 <- ggplot(data_long1, aes(x = value, y = Clade)) +
  geom_bar(
    aes(y = Clade, x = 100, fill = variable),
    stat = "identity",
    position = position_dodge2(width = 1, preserve = "single"),
    color = "#595959",
    fill = "white",
    alpha = 1
  ) +
  geom_bar(
    stat = "identity",
    position = "identity",
    color = "#595959",
    fill = "#A6A6A6",
    alpha = 1
  ) +
  geom_text(
    aes(x = 20, label = sprintf("%.1f%%", value)),
    position = position_dodge2(width = 1, preserve = "single"),
    hjust = 0.2,
    vjust = 0.5,
    color = "black"
  ) +
  labs(x = "", y = "") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

# ----------------------------- #
# Panel 2: Relative taxon PD (gtdb_PD vs midas_PD)
# ----------------------------- #
data_selected_RPD <- data[, c("Clade", "gtdb_PD", "midas_PD")]
data_selected_RPD$Clade <- factor(data_selected_RPD$Clade, levels = clade_levels)

data_long2 <- tidyr::pivot_longer(
  data_selected_RPD,
  cols = c("gtdb_PD", "midas_PD"),
  names_to = "variable",
  values_to = "value"
)

p2 <- ggplot(data_long2, aes(y = Clade, x = value, fill = variable)) +
  geom_bar(
    aes(y = Clade, x = 100, fill = variable),
    stat = "identity",
    position = position_dodge2(width = 1, preserve = "single"),
    color = "#595959",
    fill = "white",
    alpha = 1
  ) +
  geom_bar(
    stat = "identity",
    position = position_dodge2(width = 1, preserve = "single")
  ) +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("gtdb_PD" = "#E1872799", "midas_PD" = "#7876B1FF")) +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none"
  )

# ----------------------------- #
# Panel 3: MiDAS PG (midas_PG)
# ----------------------------- #
data_selected_PG <- data[, c("Clade", "midas_PG")]
data_selected_PG$Clade <- factor(data_selected_PG$Clade, levels = clade_levels)

data_long3 <- tidyr::pivot_longer(
  data_selected_PG,
  cols = "midas_PG",
  names_to = "variable",
  values_to = "value"
)

p3 <- ggplot(data_long3, aes(x = value, y = Clade, fill = variable)) +
  geom_bar(
    aes(y = Clade, x = 100, fill = variable),
    stat = "identity",
    position = position_dodge2(width = 1, preserve = "single"),
    color = "#595959",
    fill = "white",
    alpha = 1
  ) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    aes(x = 20, label = sprintf("%.1f%%", value)),
    position = position_dodge2(width = 1, preserve = "single"),
    hjust = 0.2,
    vjust = 0.5,
    color = "black"
  ) +
  scale_fill_manual(values = c("midas_PG" = "#7876B1FF")) +
  labs(x = "", y = "") +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none"
  )

# ----------------------------- #
# Combine and export
# ----------------------------- #
p_combined <- ggpubr::ggarrange(
  p1, p2, p3,
  ncol = 3,
  common.legend = TRUE,
  legend = "bottom",
  widths = c(5, 2.5, 1.5)
)

p_combined

ggsave("figures/FigS4.pre.pdf", p_combined, width = 9, height = 12)

writeLines(capture.output(sessionInfo()), "../envs/FigS4.sessionInfo.txt")
