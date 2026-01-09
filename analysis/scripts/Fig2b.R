############################################################
## Fig. 2b — Rarefaction curves (HQ vs All MAGs)
##
## Inputs:
##   - sources/Fig2b-I.txt   (HQ; two columns: mOTU, Count; no header)
##   - sources/Fig2b-II.txt  (All; two columns: mOTU, Count; no header)
##
## Output:
##   - figures/Fig2b.pdf
############################################################

library(data.table)   # robust fast I/O
library(ggplot2)
library(scales)
library(parallel)

n_replicates <- 1000
step_size <- 100

# ----------------------------- #
# Rarefaction for ONE replicate
# ----------------------------- #
calculate_rarefaction <- function(sequence, max_subsample, step = 100) {
  subsample_sizes <- seq(step, max_subsample, by = step)
  n_points <- length(subsample_sizes)
  
  n_motus <- integer(n_points)
  for (i in seq_len(n_points)) {
    size <- subsample_sizes[i]
    sampled <- sample(sequence, size, replace = FALSE)
    n_motus[i] <- length(unique(sampled))
  }
  data.frame(Reads = subsample_sizes, mOTUs = n_motus)
}

# ----------------------------- #
# Run rarefaction for an input file
# ----------------------------- #
run_rarefaction <- function(path, n_reps, step) {
  input_dt <- data.table::fread(
    path, header = FALSE, sep = "\t",
    col.names = c("mOTU", "Count")
  )
  
  # Map each mOTU to integer IDs and expand to MAG sequence
  motu_id <- seq_len(nrow(input_dt))
  mag_sequence <- rep.int(motu_id, times = input_dt$Count)
  total_mags <- length(mag_sequence)
  
  # Use all cores - 1
  n_cores <- max(1, parallel::detectCores() - 1)
  
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # Export objects/functions to workers
  parallel::clusterExport(
    cl,
    varlist = c("calculate_rarefaction", "mag_sequence", "total_mags", "step"),
    envir = environment()
  )
  
  # Ensure reproducible RNG across workers
  parallel::clusterSetRNGStream(cl, iseed = 1)
  
  # Run replicates in parallel
  res_list <- parallel::parLapply(cl, seq_len(n_reps), function(i) {
    # Keep your original “seed per replicate” behavior
    set.seed(i)
    out <- calculate_rarefaction(mag_sequence, max_subsample = total_mags, step = step)
    out$Replicate <- i
    out
  })
  
  rarefaction_df <- do.call(rbind, res_list)
  
  # Base-R summary by Reads (no data.table syntax)
  reads_vals <- sort(unique(rarefaction_df$Reads))
  
  mean_vals <- tapply(rarefaction_df$mOTUs, rarefaction_df$Reads, mean)
  low_vals  <- tapply(rarefaction_df$mOTUs, rarefaction_df$Reads, quantile, probs = 0.025)
  up_vals   <- tapply(rarefaction_df$mOTUs, rarefaction_df$Reads, quantile, probs = 0.975)
  
  stats <- data.frame(
    Reads = reads_vals,
    Mean  = as.numeric(mean_vals[as.character(reads_vals)]),
    Lower = as.numeric(low_vals[as.character(reads_vals)]),
    Upper = as.numeric(up_vals[as.character(reads_vals)])
  )
  
  list(stats = stats, total_mags = total_mags)
}

# ----------------------------- #
# HQ (Fig2b-I)
# ----------------------------- #
hq <- run_rarefaction("sources/Fig2b-I.txt", n_replicates, step_size)
stats_hq <- hq$stats
total_mags_hq <- hq$total_mags

p1 <- ggplot(stats_hq, aes(x = Reads, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#9bca3c", alpha = 0.9) +
  geom_line(color = "#2C3E50", linewidth = 1) +
  labs(x = "Recovered HQ MAGs", y = "Number of species observed") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme_bw() +
  theme(
    text = element_text(size = 22),
    axis.title.y = element_text(size = 24),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )
p1

# ----------------------------- #
# All (Fig2b-II)
# ----------------------------- #
all <- run_rarefaction("sources/Fig2b-II.txt", n_replicates, step_size)
stats_all <- all$stats
total_mags_all <- all$total_mags

stats_hq$Type <- "HQ"
stats_all$Type <- "All"
combined_stats <- rbind(stats_hq, stats_all)

# ----------------------------- #
# Combined plot (Fig2b)
# ----------------------------- #
p2 <- ggplot(combined_stats, aes(x = Reads, y = Mean, color = Type, fill = Type)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.5, color = NA) +
  geom_line(aes(y = Lower), linewidth = 0) +
  geom_line(aes(y = Upper), linewidth = 0) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = c("HQ" = "black", "All" = "black")) +
  scale_fill_manual(values = c("HQ" = "#FFA500", "All" = "#800080")) +
  labs(x = "Number of recovered MAGs", y = "Number of species observed") +
  scale_x_continuous(labels = scales::comma, limits = c(0, total_mags_all)) +
  scale_y_continuous(labels = scales::comma) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    text = element_text(size = 20),
    axis.title.y = element_text(size = 24),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  ) +
  geom_text(aes(x = 15000, y = 5000, label = "HQ"), size = 8, color = "#FFA500") +
  geom_text(aes(x = 23000, y = 20000, label = "All (HQ+MQ)"), size = 8, color = "#800080")

p2

ggsave(p2, filename = "figures/Fig2b.pdf", width = 10, height = 6)

writeLines(capture.output(sessionInfo()), "../envs/Fig2b.sessionInfo.txt")

