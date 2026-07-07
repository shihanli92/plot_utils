# conga-style neighbourhood scoring of TCR features.
#
# Given a table of per-TCR feature scores (e.g. from tcr_score() in
# tcr_scoring.R) and a k-nearest-neighbour graph over the same TCRs, test -- for
# every TCR and every feature -- whether that TCR's neighbourhood is enriched
# for high or low feature values relative to the rest of the repertoire.
#
# This re-derives conga's graph-vs-feature step (Bradley lab,
# conga/correlations.py: gex_nbrhood_rank_tcr_scores / _get_split_mean_var):
#   * scores are used raw (no standardisation)
#   * foreground = a TCR's neighbours (+ itself); background = all others
#   * per-feature fg/bg mean & variance via the second moment
#     (var = mean(x^2) - mean(x)^2), background derived from the global moments
#   * significance = Welch's two-sample t-test (unequal variance):
#         t = (mean_fg - mean_bg) / sqrt(var_fg/n_fg + var_bg/n_bg)
#   * crude Bonferroni: p * (n_nonempty_nbrhoods * n_features)
#
# The neighbour matrix is supplied by the caller, so no distance library is
# required here. It pairs naturally with tcrdistR::tcrdist_knn(), whose
# $knn_indices is exactly the expected format (see the example below).
#
# Dependencies: base R only.

# Neighbourhood enrichment z-scores / p-values for TCR features.
#
#   scores       numeric matrix or data.frame, N rows (TCRs) x F feature columns
#                (non-numeric columns are dropped with a warning)
#   nbrs         integer matrix, N x K, of 1-based row indices into `scores`
#                giving each TCR's K nearest neighbours (NA entries allowed for
#                ragged neighbourhoods; a TCR must not list itself)
#   include_self TRUE (default) to include each TCR in its own foreground, as
#                conga does; FALSE to score neighbours only
#   var_floor    variance floor for the t denominator (conga uses 1e-12)
#
# Returns a list with (all N x F, same dimnames as the feature columns):
#   z          Welch t-statistic per (TCR, feature): >0 = neighbourhood enriched
#              for high values, <0 = enriched for low values
#   pval       two-sided Welch p-value
#   pval_adj   Bonferroni-adjusted p (x n_nonempty_nbrhoods x n_features, capped 1)
#   n_fg       foreground size per TCR (neighbours + self)
tcr_nbrhood_zscore <- function(scores, nbrs, include_self = TRUE,
                               var_floor = 1e-12) {
  # Coerce scores to a numeric matrix, dropping non-numeric columns.
  if (is.data.frame(scores)) {
    num <- vapply(scores, is.numeric, logical(1))
    if (any(!num)) warning("dropping non-numeric column(s): ",
                           paste(names(scores)[!num], collapse = ", "))
    scores <- as.matrix(scores[, num, drop = FALSE])
  }
  stopifnot(is.matrix(scores), is.numeric(scores))
  nbrs <- as.matrix(nbrs)
  N <- nrow(scores); F <- ncol(scores)
  stopifnot(nrow(nbrs) == N, N >= 3L)
  if (any(nbrs < 1L | nbrs > N, na.rm = TRUE)) {
    stop("nbrs contains indices outside 1..", N)
  }

  # Foreground size per TCR: valid neighbours (+ self). Constant when nbrs is a
  # full K-column matrix, but computed per-row to allow ragged neighbourhoods.
  k_valid <- rowSums(!is.na(nbrs))
  n_fg <- k_valid + as.integer(include_self)
  n_bg <- N - n_fg
  if (any(n_fg < 1L) || any(n_bg < 1L)) {
    stop("every TCR needs >=1 foreground and >=1 background member; ",
         "check neighbourhood size vs N")
  }
  wt_fg <- n_fg / N; wt_bg <- 1 - wt_fg

  z <- pval <- matrix(NA_real_, N, F, dimnames = list(NULL, colnames(scores)))

  for (j in seq_len(F)) {
    x <- scores[, j]
    xm <- matrix(x[nbrs], nrow = N)            # neighbour values, N x K (NA-safe)
    # Foreground sums (optionally include self), then means.
    fg_sum   <- rowSums(xm, na.rm = TRUE)   + if (include_self) x   else 0
    fg_sqsum <- rowSums(xm^2, na.rm = TRUE) + if (include_self) x^2 else 0
    fg_mean   <- fg_sum   / n_fg
    fg_meansq <- fg_sqsum / n_fg
    # Background derived from the global moments (conga's _get_split_mean_var).
    M   <- mean(x); MSQ <- mean(x^2)
    bg_mean   <- (M   - wt_fg * fg_mean)   / wt_bg
    bg_meansq <- (MSQ - wt_fg * fg_meansq) / wt_bg
    var_fg <- pmax(fg_meansq - fg_mean^2, 0)
    var_bg <- pmax(bg_meansq - bg_mean^2, 0)
    # Welch's t and its Satterthwaite dof (variance floored like conga).
    vf <- pmax(var_fg, var_floor) / n_fg
    vb <- pmax(var_bg, var_floor) / n_bg
    se <- sqrt(vf + vb)
    tval <- (fg_mean - bg_mean) / se
    df <- (vf + vb)^2 / (vf^2 / (n_fg - 1) + vb^2 / (n_bg - 1))
    p <- 2 * pt(-abs(tval), df)
    # Degenerate neighbourhoods (no variance, no signal) -> t=0, p=1.
    bad <- !is.finite(tval); tval[bad] <- 0; p[bad] <- 1
    z[, j] <- tval; pval[, j] <- p
  }

  # Crude Bonferroni over all non-empty neighbourhoods x features.
  # (Assign into a matrix so dim/dimnames are preserved -- pmin() would drop them.)
  n_nonempty <- sum(n_fg > 0L)
  pval_adj <- pval * n_nonempty * F
  pval_adj[pval_adj > 1] <- 1

  list(z = z, pval = pval, pval_adj = pval_adj, n_fg = n_fg)
}

# Convenience: pull the significant (TCR, feature) hits out of a result.
#
#   res      output of tcr_nbrhood_zscore()
#   alpha    adjusted-p threshold (default 0.05)
#   Returns a data.frame (tcr, feature, z, pval, pval_adj), sorted by pval_adj.
tcr_nbrhood_hits <- function(res, alpha = 0.05) {
  idx <- which(res$pval_adj < alpha, arr.ind = TRUE)
  if (!nrow(idx)) return(data.frame(tcr = integer(), feature = character(),
                                    z = numeric(), pval = numeric(),
                                    pval_adj = numeric()))
  feats <- colnames(res$z)
  out <- data.frame(
    tcr      = idx[, 1],
    feature  = feats[idx[, 2]],
    z        = res$z[idx],
    pval     = res$pval[idx],
    pval_adj = res$pval_adj[idx],
    stringsAsFactors = FALSE
  )
  out[order(out$pval_adj, -abs(out$z)), ]
}

# ---------------------------------------------------------------------------
# Example (requires the tcrdistR package for the neighbour graph):
#
#   source("scripts/tcr_analysis/tcr_scoring.R")
#   library(tcrdistR)
#   # 1. feature scores
#   st  <- tcr_score(tcrs, features = c("hydropathy","charge","cdr3len","old_imhc"))
#   # 2. kNN graph over the same TCRs (columns: va, cdr3a, vb, cdr3b)
#   knn <- tcrdist_knn(tcrs[, c("va","cdr3a","vb","cdr3b")], "human", K = 50L)
#   # 3. neighbourhood enrichment
#   res <- tcr_nbrhood_zscore(st, knn$knn_indices)
#   head(tcr_nbrhood_hits(res))
# ---------------------------------------------------------------------------
