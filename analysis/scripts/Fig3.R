############################################################
## Figure 3 - Global MAG recovery and core taxa summaries
## representation
##
## Outputs:
##   - figures/Fig3.pdf
##
## Inputs:
##   - sources/Fig3b.xlsx (sheet: "GTDB.r226")
##   - sources/Fig3c.tsv
##   - sources/Fig3d.xlsx (sheet: "Genus")
##
## Panel tags order (required):
##   a (top-left), b (top-right), d (bottom-left), c (bottom-right)
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(scales)
  library(patchwork)
})

# ----------------------------- #
# Helper: ensure output directory exists
# ----------------------------- #
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

# ----------------------------- #
# Panel a: Continental distribution of HQ/MQ MAGs
# (counts precomputed from Supplementary files)
# ----------------------------- #
panel_a_data <- tibble::tibble(
  Category  = c("Asia", "Europe", "North America", "Oceania", "South America"),
  hqmag_num = c(5110, 7451, 5726, 2083, 1907),
  mqmag_num = c(6498, 11595, 7872, 2967, 2292)
) %>%
  mutate(Category = factor(Category, levels = c("Asia", "Europe", "North America", "Oceania", "South America"))) %>%
  pivot_longer(cols = c(hqmag_num, mqmag_num), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = factor(Variable, levels = rev(c("mqmag_num", "hqmag_num"))))

p_continents <- ggplot(panel_a_data, aes(x = Category, y = Value, fill = Variable)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6, color = NA) +
  scale_fill_manual(values = c("#f4b3ad", "#ed665dFF"), labels = c("MQ", "HQ")) +
  geom_text(
    aes(label = scales::comma(Value), y = Value),
    color = "white",
    position = position_stack(vjust = 0.5),
    size = 6
  ) +
  labs(x = NULL, y = "Number of MAGs", fill = "Quality") +
  scale_y_continuous(labels = scales::comma) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    text = element_text(size = 26),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_blank()
  )

# ----------------------------- #
# Panel b: Recovery estimation across continents (Fig3b.xlsx)
# ----------------------------- #
excel_file_b <- "sources/Fig3b.xlsx"
sheet_b <- "GTDB.r226"

raw_b <- readxl::read_excel(excel_file_b, sheet = sheet_b, skip = 1)

selected_rows <- 1:83
selected_columns <- c(
  "SA.ID",
  "GTDB.hq_ratio",
  "GTDB_all_ratio",
  "indi_hq_ratio",
  "indi.all_ratio",
  "MiDAS.hq.mOTUs_ratio",
  "MiDAS.all.mOTUs_ratio",
  "MiDAS.all.MAGs_ratio",
  "Continent"
)
b <- raw_b[selected_rows, selected_columns]

b_long <- b %>%
  pivot_longer(
    cols = -c(`SA.ID`, Continent),
    names_to = "dataset",
    values_to = "Ratio"
  ) %>%
  mutate(
    dataset = factor(
      dataset,
      levels = c(
        "GTDB.hq_ratio", "GTDB_all_ratio",
        "indi_hq_ratio", "indi.all_ratio",
        "MiDAS.hq.mOTUs_ratio", "MiDAS.all.mOTUs_ratio", "MiDAS.all.MAGs_ratio"
      ),
      labels = c(
        "GTDB-HQ SR", "GTDB-HQ+MQ SR",
        "Individual-HQ SR", "Individual-HQ+MQ SR",
        "MiDAS-HQ SR", "MiDAS-HQ+MQ SR",
        "MiDAS-HQ+MQ MAGs"
      )
    )
  )

# Auto compute (n=) labels from the same table
continent_n <- b %>%
  count(Continent, name = "n") %>%
  mutate(label = paste0(Continent, "\n(n=", n, ")"))

# Keep the continent order as it appears in data
continent_levels <- unique(b$Continent)
continent_label_map <- continent_n %>%
  mutate(Continent = factor(Continent, levels = continent_levels)) %>%
  arrange(Continent)

x_labels_b <- setNames(continent_label_map$label, as.character(continent_label_map$Continent))

p_minimap2.mapping <- ggplot(b_long, aes(x = Continent, y = Ratio, fill = dataset)) +
  geom_boxplot(alpha = 0.7, width = 0.9) +
  geom_point(alpha = 0.7, shape = 21, size = 2,
             position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = setNames(
      c("#1E90FF", "#00468B", "#C2A5CF", "#7B3294", "#66C08F", "#0F4C30", "#FFD700"),
      levels(b_long$dataset)
    )
  ) +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, by = 15)) +
  scale_x_discrete(labels = x_labels_b) +
  theme_bw() +
  labs(x = NULL, y = "Microbial community\nrecovery estimation (%)", fill = NULL) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 26),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(1.0, "cm"),
    legend.spacing.x = unit(1.0, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.title.align = 0,
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# ----------------------------- #
# Panel c: Evaluation based on third-party datasets (Fig3c.tsv)
# ----------------------------- #
samples_per_site <- 5

c_dt <- read.delim("sources/Fig3c.tsv", sep = "\t", header = TRUE)

continent_levels_c <- unique(c_dt$Continent)
c_dt <- c_dt %>%
  mutate(Continent = factor(Continent, levels = continent_levels_c))

continent_counts <- c_dt %>%
  count(Continent, name = "n_samples") %>%
  mutate(
    n = n_samples / samples_per_site,
    label = paste0(as.character(Continent), "\n(n=", n, ")")
  )

c_dt <- c_dt %>%
  left_join(continent_counts %>% select(Continent, label), by = "Continent") %>%
  mutate(continent_label = factor(label, levels = continent_counts$label))

c_dt$Group <- factor(
  c_dt$Group,
  levels = c("GTDB-HQ_spRep", "GTDB-HQ-MQ_spRep",
             "MiDAS-HQ_spRep", "MiDAS-HQ-MQ_spRep",
             "MiDAS-HQ-MQ_MAGs"),
  labels = c("GTDB-HQ SR", "GTDB-HQ+MQ SR",
             "MiDAS-HQ SR", "MiDAS-HQ+MQ SR",
             "MiDAS-HQ+MQ MAGs")
)

p_sylph.mapping <- ggplot(
  c_dt,
  aes(
    x = continent_label,
    y = ratio,
    fill = Group,
    group = interaction(continent_label, Group)
  )
) +
  geom_boxplot(
    width = 0.7,
    alpha = 0.7,
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    alpha = 0.7,
    shape = 21,
    size = 2,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)
  ) +
  scale_fill_manual(
    values = setNames(
      c("#1E90FF", "#00468B", "#66C08F", "#0F4C30", "#FFD700"),
      levels(c_dt$Group)
    )
  ) +
  labs(x = NULL, y = "Microbial community\nrecovery estimation (%)", fill = "Dataset") +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 15)) +
  theme_bw(base_size = 26) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "none",         # IMPORTANT: must be "none"
    panel.background = element_rect(fill = "white"),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_blank()
  )

# ----------------------------- #
# Panel d: HQ MAGs per core genus (Fig3d.xlsx; sheet: "Genus")
# ----------------------------- #
excel_file_d <- "sources/Fig3d.xlsx"
sheet_d <- "Genus"

raw_d <- readxl::read_excel(excel_file_d, sheet = sheet_d)

d <- raw_d[1:250, c("Genus", "No.of_HQ_MAGs", "Category")]

# Keep ordering with "Strict/General/Loose"
d$Category <- factor(d$Category, levels = c("Strict core", "General core", "Loose core"))

medians <- aggregate(`No.of_HQ_MAGs` ~ Category, data = d, FUN = median)

p_ge <- ggplot(d, aes(x = Category, y = `No.of_HQ_MAGs`, fill = Category)) +
  geom_boxplot(alpha = 0.7, width = 0.7) +
  geom_point(alpha = 0.5, shape = 21, size = 4,
             position = position_jitter(width = 0.25)) +
  scale_fill_manual(values = c(
    "Strict core"  = "#ED0000FF",
    "General core" = "#00468BFF",
    "Loose core"   = "#42B540FF"
  )) +
  scale_y_log10() +
  coord_flip() +
  geom_text(
    data = medians,
    aes(x = Category, y = `No.of_HQ_MAGs`, label = round(`No.of_HQ_MAGs`, 3)),
    vjust = -4.5,
    size = 6
  ) +
  labs(x = NULL, y = "HQ MAGs/core genus", fill = NULL) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    text = element_text(size = 26),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "none"
  )

## add tags
p_a <- p_continents + labs(tag = "a")
p_b <- p_minimap2.mapping + labs(tag = "b")
p_d <- p_ge + labs(tag = "d")
p_c <- p_sylph.mapping + labs(tag = "c")

p_combined <- (
  (p_a / p_d) |
    (p_b / p_c)
) +
  plot_layout(widths = c(2.8, 3.5),
              heights = c(1.5, 2.5)) &
  theme(plot.tag = element_text(size = 26, face = "bold"))

p_combined

ggsave(p_combined,
       filename = "figures/Fig3.pdf",
       width = 24, height = 16)

writeLines(capture.output(sessionInfo()), "../envs/Fig3.sessionInfo.txt")

