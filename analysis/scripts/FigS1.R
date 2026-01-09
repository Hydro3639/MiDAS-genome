############################################################
## Fig S1 — Nanopore sequencing metrics across sequencing phases
##
## Input:
##   - sources/FigS1.txt (tab-delimited; header required)
##     Required columns:
##       Pore_Num, Yield_Gb, SeqPhase, SeqOrder, MedianQ, MedianLength
##
## Output:
##   - figures/FigS1.pdf
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(dplyr)
  library(ggpubr)
  library(grid)       # for grobs
})

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

data <- read.table("sources/FigS1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

required_cols <- c("Pore_Num", "Yield_Gb", "SeqPhase", "SeqOrder", "MedianQ", "MedianLength")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sources/FigS1.txt: ", paste(missing_cols, collapse = ", "))
}

# ----------------------------- #
# Parameters
# ----------------------------- #
pore_threshold  <- 8138
yield_threshold <- 138

acc_phase1_y  <- 17.7
acc_phase2_y  <- 19.1
acc_phase1_x1 <- 1
acc_phase1_x2 <- 46.5
acc_phase2_x1 <- 47.5
acc_phase2_x2 <- 83

# ----------------------------- #
# Panel a: Pore number vs yield
# ----------------------------- #
p_yield <- ggplot(data, aes(x = Pore_Num, y = Yield_Gb, fill = SeqPhase)) +
  geom_point(alpha = 0.7, shape = 21, size = 8) +
  ylim(20, 200) +
  geom_text(aes(label = SeqOrder), vjust = 0.5, size = 2, alpha = 0.6, color = "white") +
  geom_vline(xintercept = pore_threshold, linetype = "dashed", color = "orange", linewidth = 1.5) +
  geom_hline(yintercept = yield_threshold, linetype = "dashed", color = "orange", linewidth = 1.5) +
  annotate("label", x = pore_threshold, y = 80, label = as.character(pore_threshold),
           size = 6, fill = "white", color = "black") +
  annotate("label", x = 6000, y = yield_threshold, label = paste0(yield_threshold, " Gbp"),
           size = 6, fill = "white", color = "black") +
  theme_bw() +
  scale_fill_aaas() +
  scale_color_aaas() +
  theme(
    text = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  ) +
  labs(y = "Yield (Gbp; Len>1K and >Q10)", x = "Flowcell pore number")

# ----------------------------- #
# Panel b: Accuracy distribution
#   - Build a base version WITH legend (for extraction)
#   - Then build the display version WITHOUT legend
# ----------------------------- #
p_acc_base <- ggplot(data, aes(x = SeqOrder, y = MedianQ, fill = SeqPhase)) +
  geom_point(alpha = 0.7, shape = 21, size = 8, show.legend = TRUE) +
  geom_segment(x = acc_phase1_x1, xend = acc_phase1_x2, y = acc_phase1_y, yend = acc_phase1_y,
               linetype = "dashed", color = "orange", linewidth = 1.5) +
  geom_segment(x = acc_phase2_x1, xend = acc_phase2_x2, y = acc_phase2_y, yend = acc_phase2_y,
               linetype = "dashed", color = "orange", linewidth = 1.5) +
  annotate("label", x = 65, y = acc_phase1_y - 0.02, label = "17.7 (98.3%)",
           size = 6, fill = "white", color = "orange") +
  annotate("label", x = 30, y = acc_phase2_y, label = "19.1 (98.8%)",
           size = 6, fill = "white", color = "orange") +
  theme_bw() +
  labs(x = "Sequencing order", y = "Median accuracy (Qvalue)", fill = "Sequencing phase") +
  scale_fill_aaas() +
  theme(
    text = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    axis.text = element_text(color = "black")
  )

p_acc <- p_acc_base + theme(legend.position = "none")

# ----------------------------- #
# Panel c: Read length distribution
# ----------------------------- #
p_len <- ggplot(data, aes(x = "", y = MedianLength)) +
  geom_violin(width = 0.5) +
  geom_boxplot(width = 0.35) +
  geom_jitter(aes(color = SeqPhase), width = 0.2, height = 0, alpha = 0.7, size = 4) +
  theme_bw() +
  labs(x = "MiDAS samples", y = "Sequencing length (Kbp)") +
  scale_color_aaas() +
  theme(
    text = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  ) +
  coord_flip()

# ----------------------------- #
# Legend extraction: pull "guide-box" from ggplotGrob
# ----------------------------- #
p_acc_for_legend <- p_acc_base +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "Sequencing Phase"))

g <- ggplotGrob(p_acc_for_legend)

# Find guide-box grob
guide_idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")

if (length(guide_idx) == 0) {
  stop("No 'guide-box' found in ggplotGrob. Legend was not generated in this environment.")
}

legend_grob <- g$grobs[[guide_idx[1]]]

# Turn grob into a ggplot-like object that ggarrange can place
legend_plot <- ggpubr::as_ggplot(grid::grobTree(legend_grob))

# ----------------------------- #
# Margins
# ----------------------------- #
p_yield <- p_yield + theme(plot.margin = margin(5, 20, 5, 20))
p_acc   <- p_acc   + theme(plot.margin = margin(5, 20, 5, 20))
p_len   <- p_len   + theme(plot.margin = margin(5, 20, 5, 20))

# ----------------------------- #
# Combine panels and export
# ----------------------------- #
p_combined <- ggpubr::ggarrange(
  p_yield,
  p_acc,
  ggpubr::ggarrange(
    p_len,
    legend_plot,
    ncol = 1,
    heights = c(2, 1.2)
  ),
  ncol = 3,
  widths = c(2, 2, 2),
  labels = c("a", "b", "c"),
  label.x = 0.01,
  label.y = 0.99,
  font.label = list(size = 16)
)

p_combined

ggsave(p_combined, filename = "figures/FigS1.pdf", width = 14, height = 6)

writeLines(capture.output(sessionInfo()), "../envs/FigS1.sessionInfo.txt")

