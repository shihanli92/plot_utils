# CDR3 sub-sequence composition comparison.
#
# Pulls a window of amino acids from each CDR3 (e.g. positions 3:5) and compares
# the amino-acid composition at those positions across groups. Returns the
# per-group, per-position frequency table plus two plots: a position x amino-acid
# heatmap (faceted by group) and a sequence logo per group.
#
# Positions are 1-based from the start of the CDR3 (the conserved C is
# position 1). Sequences shorter than max(positions) are dropped.
#
# Dependencies: dplyr, ggplot2, ggseqlogo

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Compare CDR3 sub-sequence amino-acid composition across groups.
#
#   data       data.frame with one row per clonotype / cell
#   group_cols character vector of column names defining a group to compare
#   positions  amino-acid positions to pull from each CDR3 (e.g. 3:5)
#   cdr3_col   the CDR3 column to read (default "cdr3.beta")
#   unique_by  NULL (default) to count every row, or a vector of clone-defining
#              columns to collapse duplicate clones so each counts once per group
#
# Returns a list with:
#   composition  data.frame: one row per (group, position, amino acid) with freq
#   heatmap      ggplot tile plot; position x amino acid, fill = frequency,
#                faceted by group
#   logo         ggseqlogo plot; one sequence logo per group (NULL if ggseqlogo
#                is not installed)
tcr_cdr3_composition <- function(data,
                                 group_cols,
                                 positions = 3:5,
                                 cdr3_col = "cdr3.beta",
                                 unique_by = NULL) {
  stopifnot(all(group_cols %in% names(data)), cdr3_col %in% names(data))
  # Optionally collapse duplicate clones (one row per clone per group).
  if (!is.null(unique_by)) {
    unique_by <- intersect(unique_by, names(data))
    data <- distinct(data, across(all_of(c(group_cols, unique_by))), .keep_all = TRUE)
  }
  # Keep sequences long enough to cover the requested window.
  seqs <- data[[cdr3_col]]
  keep <- !is.na(seqs) & nchar(seqs) >= max(positions)
  data <- data[keep, , drop = FALSE]
  seqs <- seqs[keep]
  # One group label per row (joins multiple group_cols with " | ").
  group <- do.call(paste, c(data[group_cols], sep = " | "))
  # Pull the amino acid at each position: rows = sequences, cols = positions.
  aa_mat <- vapply(positions, function(p) substr(seqs, p, p),
                   character(length(seqs)))
  # Long table of (group, position, amino acid), then frequencies per
  # (group, position).
  long <- data.frame(
    group    = rep(group, times = length(positions)),
    position = rep(positions, each = length(seqs)),
    aa       = as.vector(aa_mat),
    stringsAsFactors = FALSE
  )
  composition <- long %>%
    count(group, position, aa, name = "n") %>%
    group_by(group, position) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()
  # --- Heatmap: position x amino acid, fill = frequency, faceted by group ----
  comp_plot <- composition
  comp_plot$position <- factor(comp_plot$position, levels = positions)
  comp_plot$aa <- factor(comp_plot$aa, levels = rev(sort(unique(comp_plot$aa))))
  heatmap <- ggplot(comp_plot, aes(x = position, y = aa, fill = freq)) +
    geom_tile() +
    facet_wrap(~ group) +
    scale_fill_viridis_c(name = "Frequency", limits = c(0, NA)) +
    labs(x = "CDR3 position", y = "Amino acid") +
    theme_bw()
  # --- Sequence logo: one logo per group ------------------------------------
  logo <- NULL
  if (requireNamespace("ggseqlogo", quietly = TRUE)) {
    # The pulled window as one string per sequence, split into per-group vectors.
    windows <- do.call(paste0, lapply(seq_along(positions), function(i) aa_mat[, i]))
    by_group <- split(windows, group)
    logo <- ggseqlogo::ggseqlogo(by_group, method = "prob") +
      scale_x_continuous(breaks = seq_along(positions), labels = positions) +
      labs(x = "CDR3 position", y = "Frequency")
  } else {
    warning("Package 'ggseqlogo' not installed; returning logo = NULL.")
  }
  list(composition = composition, heatmap = heatmap, logo = logo)
}
