# TCR gene-segment usage diversity / entropy analysis.
#
# For each group and each V/J segment type, summarises how diverse the gene
# usage is using several standard metrics, then draws a bar plot (mean across
# groups) with the individual group values overlaid as points, faceted by
# segment (chain) type and with a chosen column on the x axis.
#
# Dependencies: dplyr, ggplot2

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Diversity / entropy metrics from a vector of gene counts.
#   richness            number of distinct genes observed
#   shannon_entropy     Shannon entropy H = -sum(p * log(p))  (natural log)
#   normalized_entropy  Pielou's evenness H / log(richness), in [0, 1]
#   inv_simpson         inverse Simpson index 1 / sum(p^2) (effective # genes)
#   gini                Gini coefficient of usage inequality, in [0, 1]
#                       (0 = every gene used equally, ->1 = a few genes dominate)
diversity_metrics <- function(counts) {
  p <- counts / sum(counts)
  H <- -sum(p * log(p))
  richness <- length(p)
  # Gini coefficient: sort counts, then the standard mean-difference formula.
  x <- sort(counts)
  n <- length(x)
  gini <- if (n > 1) (2 * sum(seq_len(n) * x)) / (n * sum(x)) - (n + 1) / n else NA_real_
  data.frame(
    richness           = richness,
    shannon_entropy    = H,
    normalized_entropy = if (richness > 1) H / log(richness) else NA_real_,
    inv_simpson        = 1 / sum(p^2),
    gini               = gini
  )
}

# Calculate gene-segment usage diversity across groups.
#
#   data       data.frame with one row per clonotype / cell
#   x_col      column to place on the x axis (e.g. "antigen.epitope")
#   group_cols NULL, or extra column(s) that split each x value into replicate
#              units (e.g. c("meta.subject.id")); each (x, group) combination
#              becomes one point, and the bar is their mean
#   gene_cols  gene-segment columns to summarise (defaults to V/J of both chains)
#   unique_by  NULL (default) to count every row, or a vector of clone-defining
#              columns (e.g. c("v.alpha", "cdr3.alpha")) to collapse duplicate
#              clones so each unique clone counts once per unit
#   metric     which metric to plot on y (one of the columns of diversity_metrics);
#              the returned table always contains all of them
#
# Returns a list with:
#   diversity  data.frame: one row per (x, group, segment) with all metrics
#   plot       ggplot object; x = x_col, y = chosen metric, bar = mean across
#              groups + points for individual groups, faceted by segment type
tcr_gene_diversity <- function(data,
                               x_col,
                               group_cols = NULL,
                               gene_cols = c("v.alpha", "j.alpha", "v.beta", "j.beta"),
                               unique_by = NULL,
                               metric = "shannon_entropy") {
  gene_cols <- intersect(gene_cols, names(data))
  cell_cols <- c(x_col, group_cols)  # columns identifying one diversity unit
  stopifnot(length(gene_cols) > 0, all(cell_cols %in% names(data)))
  # Optionally keep only unique clones: one row per distinct clone (as defined
  # by unique_by) within each unit. NULL means count every row as-is.
  if (!is.null(unique_by)) {
    unique_by <- intersect(unique_by, names(data))
    data <- distinct(data, across(all_of(c(cell_cols, unique_by))), .keep_all = TRUE)
  }
  # Gene counts per (unit, segment, gene), then collapse to diversity metrics
  # per (unit, segment).
  diversity <- lapply(gene_cols, function(col) {
    data %>%
      filter(!is.na(.data[[col]]), .data[[col]] != "") %>%
      count(across(all_of(cell_cols)), gene = .data[[col]], name = "n") %>%
      group_by(across(all_of(cell_cols))) %>%
      group_modify(~ diversity_metrics(.x$n)) %>%
      ungroup() %>%
      mutate(segment = col)
  }) %>%
    bind_rows()
  metric_cols <- c("richness", "shannon_entropy", "normalized_entropy",
                   "inv_simpson", "gini")
  stopifnot(metric %in% metric_cols)
  diversity$segment <- factor(diversity$segment, levels = gene_cols)
  plot <- ggplot(diversity, aes(x = .data[[x_col]], y = .data[[metric]])) +
    stat_summary(fun = mean, geom = "col", fill = "grey80", width = 0.7) +
    geom_point(position = position_jitter(width = 0.15, height = 0),
               size = 1, alpha = 0.6) +
    facet_wrap(~ segment) +
    labs(x = x_col, y = metric) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  list(diversity = diversity, plot = plot)
}
