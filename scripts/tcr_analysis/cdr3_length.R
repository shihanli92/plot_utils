# CDR3 length distribution comparison.
#
# Computes the CDR3 length (number of amino acids) of each TCR and compares the
# length distribution across groups. Returns the per-group length-frequency
# table, a per-group summary (mean/median/sd/range), and a distribution plot
# (frequency polygon, one line per group).
#
# Dependencies: dplyr, ggplot2

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Compare CDR3 length distributions across groups.
#
#   data       data.frame with one row per clonotype / cell
#   group_cols character vector of column names defining a group to compare
#   cdr3_col   the CDR3 column to read (default "cdr3.beta")
#   unique_by  NULL (default) to count every row, or a vector of clone-defining
#              columns to collapse duplicate clones so each counts once per group
#
# Returns a list with:
#   lengths  data.frame: one row per (group, length) with count and frequency
#            within that group
#   summary  data.frame: one row per group with n, mean, median, sd, min, max
#   plot     ggplot object; x = CDR3 length, y = frequency, one line + points
#            per group (frequency polygon)
tcr_cdr3_length <- function(data,
                            group_cols,
                            cdr3_col = "cdr3.beta",
                            unique_by = NULL) {
  stopifnot(all(group_cols %in% names(data)), cdr3_col %in% names(data))
  # Optionally collapse duplicate clones (one row per clone per group).
  if (!is.null(unique_by)) {
    unique_by <- intersect(unique_by, names(data))
    data <- distinct(data, across(all_of(c(group_cols, unique_by))), .keep_all = TRUE)
  }
  # Keep populated sequences and measure each one's length.
  seqs <- data[[cdr3_col]]
  keep <- !is.na(seqs) & seqs != ""
  data <- data[keep, , drop = FALSE]
  seqs <- seqs[keep]
  # One group label per row (joins multiple group_cols with " | ").
  group <- do.call(paste, c(data[group_cols], sep = " | "))
  long <- data.frame(
    group  = group,
    length = nchar(seqs),
    stringsAsFactors = FALSE
  )
  # Length frequencies per group (each group's frequencies sum to 1).
  lengths <- long %>%
    count(group, length, name = "count") %>%
    group_by(group) %>%
    mutate(freq = count / sum(count)) %>%
    ungroup()
  # Per-group summary statistics of the raw lengths.
  summary <- long %>%
    group_by(group) %>%
    summarise(
      n      = n(),
      mean   = mean(length),
      median = median(length),
      sd     = sd(length),
      min    = min(length),
      max    = max(length),
      .groups = "drop"
    )
  plot <- ggplot(lengths, aes(x = length, y = freq, colour = group)) +
    geom_line() +
    geom_point(size = 1, alpha = 0.8) +
    scale_x_continuous(breaks = function(lim) seq(floor(lim[1]), ceiling(lim[2]), by = 2)) +
    labs(x = "CDR3 length (aa)", y = "Frequency", colour = NULL) +
    theme_bw()
  list(lengths = lengths, summary = summary, plot = plot)
}
