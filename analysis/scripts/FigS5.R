############################################################
## Fig S5 — RED-based novel lineages
##
## Input:
##   - The numbers was extracted from the RED-based evaluation
##
## Output:
##   - figures/FigS5.pdf
##
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

df <- data.frame(
  Category = c("Class", "Order", "Family", "Genus"),
  Novel.number = c(2, 21, 104, 1040),
  stringsAsFactors = FALSE
)

df$Category <- factor(df$Category, levels = c("Class", "Order", "Family", "Genus"))

p_novel <- ggplot(df, aes(x = Category, y = Novel.number)) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    color = "#45818e",
    fill = "#45818e",
    alpha = 0.8
  ) +
  geom_text(
    aes(label = scales::comma(Novel.number)),
    vjust = -0.5,
    color = "black",
    size = 8
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "", y = "RED-based novel lineages") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    text = element_text(size = 36)
  )

p_novel

ggsave(p_novel, filename = "figures/FigS5.pdf", width = 10, height = 9)

writeLines(capture.output(sessionInfo()), "../envs/FigS5.sessionInfo.txt")
