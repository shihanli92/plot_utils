# CDR3 self-reactivity from the position-6/7 doublet hydrophobicity score.
#
# Each CDR3 is scored by the amino-acid doublet at positions 6 and 7 using a
# fixed lookup matrix: 2 = self-reactive, 1 = neutral, 0 = hydrophilic. The
# proportion of CDR3s in each group falling into each class is returned, along
# with a stacked-bar plot (proportion self-reactive is the labelled segment).
#
# Positions are 1-based from the start of the CDR3 (the conserved C is
# position 1). Sequences too short to cover both positions are dropped, as are
# doublets containing a residue not present in the scoring matrix.
#
# The scoring matrix is defined for the position-6/7 doublet; `positions` is
# exposed only so the same matrix can be applied to a different residue pair.
#
# Dependencies: dplyr, ggplot2

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Doublet self-reactivity scores. Rows = position-6 amino acid, columns =
# position-7 amino acid; entry = 2 (self-reactive) / 1 (neutral) / 0 (hydrophilic).
p6_order <- c("W", "C", "F", "Y", "L", "I", "H", "V", "M", "A",
              "P", "S", "N", "R", "D", "G", "T", "Q", "E", "K")

selfreactivity_scores <- as.matrix(data.frame(
  C = c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
  W = c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1),
  Y = c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1),
  F = c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1),
  L = c(2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1),
  P = c(2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1),
  V = c(2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0),
  G = c(2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0),
  I = c(2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0),
  H = c(2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0),
  N = c(2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0),
  M = c(2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0,0,0),
  A = c(2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,0,0,0),
  R = c(2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0),
  S = c(2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0),
  T = c(2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0),
  D = c(2,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0),
  Q = c(2,2,2,2,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0),
  K = c(2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0),
  E = c(1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  row.names = p6_order, check.names = FALSE
))

# Score CDR3 self-reactivity from the position-6/7 doublet and summarise per group.
#
#   data       data.frame with one row per clonotype / cell
#   group_cols character vector of column names defining a group to compare
#   positions  the two CDR3 positions forming the doublet (default c(6, 7))
#   cdr3_col   the CDR3 column to read (default "cdr3.beta")
#   unique_by  NULL (default) to count every row, or a vector of clone-defining
#              columns to collapse duplicate clones so each counts once per group
#   scores     scoring matrix (rows = first-position aa, columns = second-position
#              aa); defaults to the built-in position-6/7 doublet matrix
#
# Returns a list with:
#   scored       data.frame: one row per scored CDR3 with group, the two residues,
#                its numeric score, and its class
#   proportions  data.frame: one row per (group, class) with n, group total, and
#                proportion within the group
#   plot         ggplot stacked bar; x = group, y = proportion, fill = class
tcr_cdr3_selfreactivity <- function(data,
                                    group_cols,
                                    positions = c(6, 7),
                                    cdr3_col = "cdr3.beta",
                                    unique_by = NULL,
                                    scores = selfreactivity_scores) {
  stopifnot(all(group_cols %in% names(data)), cdr3_col %in% names(data),
            length(positions) == 2)
  # Optionally collapse duplicate clones (one row per clone per group).
  if (!is.null(unique_by)) {
    unique_by <- intersect(unique_by, names(data))
    data <- distinct(data, across(all_of(c(group_cols, unique_by))), .keep_all = TRUE)
  }
  # Keep sequences long enough to cover both doublet positions.
  seqs <- data[[cdr3_col]]
  keep <- !is.na(seqs) & nchar(seqs) >= max(positions)
  data <- data[keep, , drop = FALSE]
  seqs <- seqs[keep]
  # One group label per row (joins multiple group_cols with " | ").
  group <- do.call(paste, c(data[group_cols], sep = " | "))
  aa1 <- substr(seqs, positions[1], positions[1])
  aa2 <- substr(seqs, positions[2], positions[2])
  # Drop doublets whose residues are not in the scoring matrix (e.g. "*", "X").
  scorable <- aa1 %in% rownames(scores) & aa2 %in% colnames(scores)
  if (any(!scorable)) {
    warning(sum(!scorable), " CDR3(s) dropped: residue not in scoring matrix.")
  }
  group <- group[scorable]; aa1 <- aa1[scorable]; aa2 <- aa2[scorable]
  # Look up the doublet score, then label each class.
  score <- scores[cbind(aa1, aa2)]
  class <- factor(c("hydrophilic", "neutral", "self-reactive")[score + 1],
                  levels = c("self-reactive", "neutral", "hydrophilic"))
  scored <- data.frame(group = group, p1 = aa1, p2 = aa2,
                       score = score, class = class, stringsAsFactors = FALSE)
  # Proportion of each class within each group (all three classes always present).
  proportions <- scored %>%
    count(group, class, name = "n", .drop = FALSE) %>%
    group_by(group) %>%
    mutate(total = sum(n), prop = n / total) %>%
    ungroup()
  plot <- ggplot(proportions, aes(x = group, y = prop, fill = class)) +
    geom_col() +
    scale_fill_manual(
      values = c("self-reactive" = "#b2182b",
                 "neutral"       = "grey75",
                 "hydrophilic"   = "#2166ac"),
      name = NULL
    ) +
    labs(x = NULL, y = "Proportion of CDR3s") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  list(scored = scored, proportions = proportions, plot = plot)
}
