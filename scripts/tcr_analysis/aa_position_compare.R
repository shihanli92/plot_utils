# Amino-acid composition at a CDR3 position, compared across groups.
#
# Pulls a k-mer starting at one specified CDR3 position (k = 1 single amino acid,
# 2 doublet, 3 triplet) and compares its frequency across groups. Two modes:
#   * all groups (groups = NULL)      -> heatmap: group x motif, fill = freq
#   * two groups (groups = c(g1, g2)) -> scatter: one point per motif,
#                                        x = freq in g1, y = freq in g2
#
# Position is 1-based from the start of the CDR3 (the conserved C is position 1).
# Sequences too short to cover position .. position+k-1 are dropped.
#
# Dependencies: dplyr, ggplot2

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})


# Compare amino-acid (k-mer) composition at a CDR3 position across groups.
#
#   data       data.frame with one row per clonotype / cell
#   group_cols character vector of column names defining a group
#   position   single CDR3 position to start reading from (e.g. 5)
#   k          motif length: 1 = single amino acid, 2 = doublet, 3 = triplet
#   cdr3_col   the CDR3 column to read (default "cdr3.beta")
#   groups     NULL for a heatmap of all groups, or a length-2 vector of group
#              labels for a g1-vs-g2 scatter
#   unique_by  NULL (default) to count every row, or a vector of clone-defining
#              columns to collapse duplicate clones so each counts once per group
#   top_n      NULL to plot every motif, or keep only the top_n most frequent
#              motifs (by total count across groups) in the plot. Useful for
#              doublets/triplets, where the number of motifs explodes. The
#              returned composition table always contains every motif.
#
# Returns a list with:
#   composition  data.frame: one row per (group, motif) with count and freq
#   plot         ggplot object; heatmap (all groups) or scatter (two groups)
tcr_aa_position_compare <- function(data,
                                    group_cols,
                                    position,
                                    k = 1,
                                    cdr3_col = "cdr3.beta",
                                    groups = NULL,
                                    unique_by = NULL,
                                    top_n = NULL) {

  stopifnot(all(group_cols %in% names(data)), cdr3_col %in% names(data),
            length(position) == 1, k >= 1)

  # Optionally collapse duplicate clones (one row per clone per group).
  if (!is.null(unique_by)) {
    unique_by <- intersect(unique_by, names(data))
    data <- distinct(data, across(all_of(c(group_cols, unique_by))), .keep_all = TRUE)
  }

  # Keep sequences long enough to cover the k-mer, then pull positions
  # position .. position + k - 1 as a single motif.
  last <- position + k - 1
  seqs <- data[[cdr3_col]]
  keep <- !is.na(seqs) & nchar(seqs) >= last
  data <- data[keep, , drop = FALSE]
  seqs <- seqs[keep]

  group <- do.call(paste, c(data[group_cols], sep = " | "))

  # Motif frequency per group at this position.
  composition <- data.frame(group = group, motif = substr(seqs, position, last),
                            stringsAsFactors = FALSE) %>%
    count(group, motif, name = "n") %>%
    group_by(group) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()

  # Labels that adapt to the motif length.
  unit_lab  <- if (k == 1) "Amino acid" else paste0(k, "-mer")
  title_lab <- if (k == 1) paste0("CDR3 position ", position)
               else paste0("CDR3 positions ", position, "-", last)

  # Optionally restrict the plot to the most frequent motifs (table keeps all).
  comp_plot <- composition
  if (!is.null(top_n)) {
    keep_motifs <- composition %>%
      group_by(motif) %>%
      summarise(total = sum(n), .groups = "drop") %>%
      slice_max(total, n = top_n, with_ties = FALSE) %>%
      pull(motif)
    comp_plot <- composition[composition$motif %in% keep_motifs, ]
  }

  if (is.null(groups)) {
    # --- Heatmap across all groups ------------------------------------------
    comp_plot$motif <- factor(comp_plot$motif, levels = rev(sort(unique(comp_plot$motif))))
    plot <- ggplot(comp_plot, aes(x = group, y = motif, fill = freq)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Frequency", limits = c(0, NA)) +
      labs(x = NULL, y = unit_lab, title = title_lab) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  } else {
    # --- Scatter of two groups ----------------------------------------------
    stopifnot(length(groups) == 2, all(groups %in% composition$group))
    # Build a wide table with one row per motif (missing -> frequency 0).
    motifs <- sort(unique(comp_plot$motif))
    g1 <- comp_plot[comp_plot$group == groups[1], ]
    g2 <- comp_plot[comp_plot$group == groups[2], ]
    wide <- data.frame(
      motif = motifs,
      x     = g1$freq[match(motifs, g1$motif)],
      y     = g2$freq[match(motifs, g2$motif)]
    )
    wide$x[is.na(wide$x)] <- 0
    wide$y[is.na(wide$y)] <- 0

    plot <- ggplot(wide, aes(x = x, y = y, label = motif)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
      geom_point() +
      geom_text(vjust = -0.6, size = 3, check_overlap = TRUE) +
      labs(x = paste0("Frequency (", groups[1], ")"),
           y = paste0("Frequency (", groups[2], ")"),
           title = title_lab) +
      coord_equal() +
      theme_bw()
  }

  list(composition = composition, plot = plot)
}
