# TCR chain (gene-segment) usage analysis.
#
# Computes how often each V/J gene segment is used, separately within each
# group (e.g. donor, sample, or epitope), and visualises the per-group
# frequencies as a boxplot + points, one facet per segment type.
#
# Dependencies: dplyr, ggplot2

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Order gene-segment names by gene position rather than alphabetically, so
# e.g. TRBV2 comes before TRBV10 (and TRBV6-1 before TRBV6-2).
order_genes <- function(genes) {
  genes  <- unique(genes)
  base   <- sub("\\*.*$", "", genes)        # drop the "*01" allele suffix
  prefix <- sub("[0-9].*$", "", base)       # leading letters, e.g. "TRBV"
  nums   <- regmatches(base, gregexpr("[0-9]+", base))  # numeric parts
  family <- vapply(nums, function(x) if (length(x) >= 1) as.integer(x[1]) else 0L, 1L)
  member <- vapply(nums, function(x) if (length(x) >= 2) as.integer(x[2]) else 0L, 1L)
  genes[order(prefix, family, member)]
}

# Calculate gene-segment usage frequencies across groups.
#
#   data       data.frame with one row per clonotype / cell
#   group_cols character vector of column names defining a group (the unit
#              across which frequencies are compared, e.g. c("donor"))
#   gene_cols  gene-segment columns to summarise (defaults to V/J of both chains)
#   unique_by  NULL (default) to count every row, or a vector of clone-defining
#              columns (e.g. c("v.alpha", "cdr3.alpha")) to collapse duplicate
#              clones so each unique clone counts once per group
#
# Returns a list with:
#   usage  long data.frame: one row per (group, segment, gene) with its
#          frequency within that group
#   plot   ggplot object; x = gene segment (ordered by gene position),
#          y = frequency, boxplot + points, faceted by segment type
tcr_chain_usage <- function(data,
                            group_cols,
                            gene_cols = c("v.alpha", "j.alpha", "v.beta", "j.beta"),
                            unique_by = NULL) {
  gene_cols <- intersect(gene_cols, names(data))
  stopifnot(length(gene_cols) > 0, all(group_cols %in% names(data)))
  # Optionally keep only unique clones: one row per distinct clone (as defined
  # by unique_by) within each group. NULL means count every row as-is.
  if (!is.null(unique_by)) {
    unique_by <- intersect(unique_by, names(data))
    data <- distinct(data, across(all_of(c(group_cols, unique_by))), .keep_all = TRUE)
  }
  # For each segment column, count gene occurrences per group and divide by the
  # group's total to get a frequency; then stack all segments into one table.
  usage <- lapply(gene_cols, function(col) {
    data %>%
      filter(!is.na(.data[[col]]), .data[[col]] != "") %>%
      group_by(across(all_of(group_cols))) %>%
      mutate(.group_total = n()) %>%
      group_by(across(all_of(c(group_cols, col)))) %>%
      summarise(
        segment = col,
        gene    = first(.data[[col]]),
        count   = n(),
        freq    = n() / first(.group_total),
        .groups = "drop"
      ) %>%
      select(all_of(group_cols), segment, gene, count, freq)
  }) %>%
    bind_rows() %>%
    mutate(gene = factor(gene, levels = order_genes(gene)))
  plot <- ggplot(usage, aes(x = gene, y = freq)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.15, height = 0),
               size = 0.8, alpha = 0.6) +
    facet_wrap(~ segment, scales = "free_x") +
    labs(x = "Gene segment", y = "Frequency") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  list(usage = usage, plot = plot)
}
