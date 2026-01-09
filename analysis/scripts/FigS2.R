############################################################
## Fig S2 — Distributions of genome size, coding density, GC content
##
## Input:
##   - sources/FigS2.tsv (tab-delimited; header required)
##     Required columns:
##       Genome_Size, Coding_Density, GC_Content
##
## Output:
##   - figures/FigS2.pdf
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(ggpubr)
})

# ----------------------------- #
# I/O
# ----------------------------- #
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

metadata <- read.table("sources/FigS2.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

required_cols <- c("Genome_Size", "Coding_Density", "GC_Content")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sources/FigS2.tsv: ", paste(missing_cols, collapse = ", "))
}

# ----------------------------- #
# colors
# ----------------------------- #
fill_color <- "#8491B4FF"
vline_color <- "orange"

## median value
v_gsize_mb <- median(metadata$Genome_Size, na.rm = TRUE) / 1e6
v_coding_pct <- median(metadata$Coding_Density, na.rm = TRUE) * 100
v_gc_pct <- median(metadata$GC_Content, na.rm = TRUE) * 100

# ----------------------------- #
# Panel 1: Genome size (Mb)
# ----------------------------- #
p1_gSize <- ggplot(metadata, aes(x = Genome_Size / 1e6)) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = fill_color) +
  theme_bw() +
  labs(x = "Genome size (Mb)", y = "HQ MAGs (22,277)") +
  geom_vline(xintercept = v_gsize_mb, linetype = "dashed", color = vline_color, linewidth = 1) +
  scale_x_continuous(
    limits = c(0, 26),
    breaks = seq(1, 26, 3),
    labels = seq(1, 26, 3)
  ) +
  scale_y_continuous(
    limits = c(0, 0.35),
    labels = scales::percent
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    text = element_text(size = 20)
  )

# ----------------------------- #
# Panel 2: Coding density (%)
# ----------------------------- #
p2_codingDensity <- ggplot(metadata, aes(x = Coding_Density * 100)) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = fill_color, bins = 20) +
  theme_bw() +
  labs(x = "Coding density (%)", y = "") +
  geom_vline(xintercept = v_coding_pct, linetype = "dashed", color = vline_color, linewidth = 1) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    labels = seq(0, 100, 10)
  ) +
  scale_y_continuous(
    limits = c(0, 0.20),
    labels = scales::percent
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    text = element_text(size = 20)
  )

# ----------------------------- #
# Panel 3: GC content (%)
# ----------------------------- #
p3_gc <- ggplot(metadata, aes(x = GC_Content * 100)) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = fill_color, bins = 20) +
  theme_bw() +
  geom_vline(xintercept = v_gc_pct, linetype = "dashed", color = vline_color, linewidth = 1) +
  labs(x = "GC content (%)", y = "") +
  scale_x_continuous(
    breaks = seq(20, 80, 12),
    labels = seq(20, 80, 12)
  ) +
  scale_y_continuous(
    limits = c(0, 0.05),
    labels = scales::percent
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    text = element_text(size = 20)
  )

# ----------------------------- #
# Combine and export
# ----------------------------- #
p <- ggpubr::ggarrange(p1_gSize, p2_codingDensity, p3_gc, ncol = 3, nrow = 1)

p

ggsave(p, filename = "figures/FigS2.pdf", width = 12, height = 4)

writeLines(capture.output(sessionInfo()), "../envs/FigS2.sessionInfo.txt")


