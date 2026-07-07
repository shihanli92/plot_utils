# conga-style TCR feature scoring.
#
# Computes a table of per-TCR "feature" scores, re-derived in R from the conga
# package (Bradley lab, github.com/phbradley/conga, conga/tcr_scoring.py). A set
# of private per-feature helpers each compute one feature; the master function
# tcr_score() dispatches a requested list of features into a single score table
# (one row per TCR, one column per feature), mirroring conga's
# make_tcr_score_table().
#
# Only the self-contained features are implemented -- the trained-model scores
# (cd8, imhc, oldcd8) and neighbour-graph scores (nndists_tcr, N_ins) are not.
#
# Embedded reference data, taken verbatim from conga/data (cite: phbradley/conga):
#   * aa_props        -- conga/data/aa_props.tsv (20 amino acids x 28 properties)
#   * trav_order/traj_order -- conga/data/imgt_tra_locus_order.txt (locus order)
#
# CDR3 positions are 1-based; the conserved C is position 1. The "FG loop" is
# conga's cdr3[4:-4] (drop 4 residues from each end). Rows with an NA/empty CDR3
# are dropped before scoring.
#
# Dependencies: base R only.

# ---------------------------------------------------------------------------
# Embedded reference data
# ---------------------------------------------------------------------------

# Amino-acid property table (rows = amino acid, columns = property).
aa_props <- structure(c(-0.591, -1.343, 1.05, 1.357, -1.006, -0.384, 0.336,
-1.239, 1.831, -1.019, -0.663, 0.945, 0.189, 0.931, 1.538, -0.228,
-0.032, -1.337, -0.595, 0.26, -1.302, 0.465, 0.302, -1.453, -0.59,
1.652, -0.417, -0.547, -0.561, -0.987, -1.524, 0.828, 2.081,
-0.179, -0.055, 1.399, 0.326, -0.279, 0.009, 0.83, -0.733, -0.862,
-3.656, 1.477, 1.891, 1.33, -1.673, 2.131, 0.533, -1.505, 2.219,
1.299, -1.628, -3.005, 1.502, -4.76, 2.213, -0.544, 0.672, 3.097,
1.57, -1.02, -0.259, 0.113, -0.397, 1.045, -1.474, 0.393, -0.277,
1.266, -1.005, -0.169, 0.421, -0.503, 0.44, 0.67, 0.908, 1.242,
-2.128, -0.838, -0.146, -0.255, -3.242, -0.837, 0.412, 2.064,
-0.078, 0.816, 1.648, -0.912, 1.212, 0.933, -1.392, -1.853, 2.897,
-2.647, 1.313, -1.262, -0.184, 1.512, 1.29, 1.11, 1.04, 1.44,
1.07, 0.56, 1.22, 0.97, 1.23, 1.3, 1.47, 0.9, 0.52, 1.27, 0.96,
0.82, 0.82, 0.91, 0.99, 0.72, 0.9, 0.74, 0.72, 0.75, 1.32, 0.92,
1.08, 1.45, 0.77, 1.02, 0.97, 0.76, 0.64, 0.8, 0.99, 0.95, 1.21,
1.49, 1.14, 1.25, 0, 0, -1, -1, 0, 0, 0.5, 0, 1, 0, 0, 0, 0,
0, 1, 0, 0, 0, 0, 0, 0.049, 0.02, 0.051, 0.051, 0.051, 0.06,
0.034, 0.047, 0.05, 0.078, 0.027, 0.058, 0.051, 0.051, 0.066,
0.057, 0.064, 0.049, 0.022, 0.07, 0, -1, 1, 1, -1, 1, -1, -1,
1, -1, 1, 1, 1, 1, 1, 1, 0, -1, -1, -1, 1.8, 2.5, -3.5, -3.5,
2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5,
-0.8, -0.7, 4.2, -0.9, -1.3, -1.56, 0.12, 0.58, -1.45, -0.21,
1.46, -0.41, -0.73, -0.34, -1.04, -1.4, 1.14, 2.06, -0.47, 0.22,
0.81, 0.26, -0.74, 0.3, 1.38, -0.48, 1.1, 0.7, -0.12, -0.44,
0.46, 1.63, -1.78, 0.6, 0.93, 0.27, -1.73, -0.28, -2.33, 0.93,
-0.23, 0.19, 0.65, -0.6, 0.53, -1.67, -0.89, -0.22, 0.19, 0.98,
-1.96, 0.52, -0.16, 0.82, 0, 0.18, -0.07, -0.33, 0.24, 1.27,
-1.08, -0.7, -0.71, 2.1, 1.48, -0.97, 0.45, -1.58, -1.61, -0.36,
-0.23, -0.28, 1.79, -0.23, -0.24, -0.42, -0.12, -1.15, 0.07,
1.37, 0.16, 1.21, 2.04, -0.72, 0.8, -0.27, -1.05, 0.81, 1.17,
-1.43, -0.16, 0.28, -0.77, 1.7, -1.1, -0.73, 0.81, -0.75, 1.1,
1.87, 0.42, 0.63, -0.4, -1.57, -0.56, -0.93, -0.71, -0.92, -1.31,
0.22, 0.1, 1.61, -0.54, 1.54, -0.55, 2, 0.18, 0.88, 1.1, -1.7,
-0.21, -0.1, 0.5, -1.16, 0, -0.78, 2.41, 0.15, 0.4, -0.81, -0.11,
1.01, 0.03, -1.62, -2.05, 1.52, 0.37, -0.45, 0.59, 0.46, -0.43,
0.21, -0.81, 0.57, -0.68, -0.2, 1.52, -1.52, 0.04, 0.67, 1.32,
-1.85, -0.83, 1.15, 0.96, 0.26, -0.09, 0.3, 0.84, 0.92, -1.89,
0.24, -1.07, -0.48, -0.31, -0.08, -0.69, 0.47, 0.38, 1.1, 2.36,
0.47, 0.51, -0.08, -0.76, 0.11, 1.23, -2.3, -0.71, -0.39, -1.15,
-1.15, 0.06, -0.4, 1.03, 0.21, 1.13, 0.76, -0.35, 1.71, -1.66,
1.13, 0.66, -0.48, 0.45, -1.27, 1.1, 0.74, -0.03, 0.23, -0.97,
-0.56, -0.46, -2.3, -0.05, -2.8455, -3.782, -2.116, -2.141, -5.017,
-2.499, -2.927, -4.641, -1.789, -5.023, -4.1915, -2.349, -2.443,
-2.2505, -2.402, -2.308, -2.6145, -4.093, -4.1375, -3.7505, 0,
0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0.047,
0.015, 0.071, 0.094, 0.021, 0.071, 0.022, 0.032, 0.105, 0.052,
0.017, 0.062, 0.052, 0.053, 0.068, 0.072, 0.064, 0.048, 0.007,
0.032, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1,
1, 1, 0.065, 0.015, 0.074, 0.089, 0.029, 0.07, 0.025, 0.035,
0.08, 0.063, 0.016, 0.053, 0.054, 0.051, 0.059, 0.071, 0.065,
0.048, 0.012, 0.033, 0.78, 0.8, 1.41, 1, 0.58, 1.64, 0.69, 0.51,
0.96, 0.59, 0.39, 1.28, 1.91, 0.97, 0.88, 1.33, 1.03, 0.47, 0.75,
1.05, 67, 86, 91, 109, 135, 48, 118, 124, 135, 124, 124, 96,
90, 114, 148, 73, 93, 105, 163, 141), dim = c(20L, 28L), dimnames = list(
    c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
    "N", "P", "Q", "R", "S", "T", "V", "W", "Y"), c("af1", "af2",
    "af3", "af4", "af5", "alpha", "beta", "charge", "core", "disorder",
    "hydropathy", "kf1", "kf10", "kf2", "kf3", "kf4", "kf5",
    "kf6", "kf7", "kf8", "kf9", "mjenergy", "polarity", "rim",
    "strength", "surface", "turn", "volume")))

amino_acids <- rownames(aa_props)                     # A..Y
property_names <- colnames(aa_props)                  # 28 property columns

# TRAV / TRAJ gene names in IMGT locus order (conga/data/imgt_tra_locus_order.txt).
trav_order <- c("TRAV1-1", "TRAV1-2", "TRAV2", "TRAV3", "TRAV4", "TRAV5", "TRAV6",
  "TRAV7", "TRAVA", "TRAV8-1", "TRAV9-1", "TRAV10", "TRAV11", "TRAV12-1", "TRAV8-2",
  "TRAV8-3", "TRAVB", "TRAV13-1", "TRAV14-1", "TRAV11-1", "TRAV12-2", "TRAV8-4",
  "TRAV8-5", "TRAV13-2", "TRAV14/DV4", "TRAV9-2", "TRAV15", "TRAV12-3", "TRAV8-6",
  "TRAV16", "TRAV17", "TRAV18", "TRAV19", "TRAVC", "TRAV20", "TRAV21", "TRAV8-6-1",
  "TRAV22", "TRAV23/DV6", "TRAV24", "TRAV25", "TRAV26-1", "TRAV8-7", "TRAV27",
  "TRAV28", "TRAV29/DV5", "TRAV30", "TRAV31", "TRAV32", "TRAV33", "TRAV26-2",
  "TRAV34", "TRAV35", "TRAV36/DV7", "TRAV37", "TRAV38-1", "TRAV38-2/DV8", "TRAV39",
  "TRAV40", "TRAV41", "TRAV46")
traj_order <- c("TRAJ61", "TRAJ60", "TRAJ59", "TRAJ58", "TRAJ57", "TRAJ56", "TRAJ55",
  "TRAJ54", "TRAJ53", "TRAJ52", "TRAJ51", "TRAJ50", "TRAJ49", "TRAJ48", "TRAJ47",
  "TRAJ46", "TRAJ45", "TRAJ44", "TRAJ43", "TRAJ42", "TRAJ41", "TRAJ40", "TRAJ39",
  "TRAJ38", "TRAJ37", "TRAJ36", "TRAJ35", "TRAJ34", "TRAJ33", "TRAJ32", "TRAJ31",
  "TRAJ30", "TRAJ29", "TRAJ28", "TRAJ27", "TRAJ26", "TRAJ25", "TRAJ24", "TRAJ23",
  "TRAJ22", "TRAJ21", "TRAJ20", "TRAJ19", "TRAJ18", "TRAJ17", "TRAJ16", "TRAJ15",
  "TRAJ14", "TRAJ13", "TRAJ12", "TRAJ11", "TRAJ10", "TRAJ9", "TRAJ8", "TRAJ7",
  "TRAJ6", "TRAJ5", "TRAJ4", "TRAJ3", "TRAJ2", "TRAJ1")

# conga constants
.fg_trim <- 4L         # FG loop = cdr3[fg_trim : -fg_trim]
.center_len <- 5L      # center window length
.cdr3_score_modes <- c("fg", "cen")
.default_score_mode <- "fg"

# ---------------------------------------------------------------------------
# Small shared helpers
# ---------------------------------------------------------------------------

# Strip an IMGT allele suffix: "TRAV1-2*01" -> "TRAV1-2".
.strip_allele <- function(g) sub("\\*.*$", "", g)

# FG-loop substring cdr3[4:-4] (empty string if length <= 8).
.fg_loop <- function(cdr3) {
  n <- nchar(cdr3)
  ifelse(n > 2L * .fg_trim, substr(cdr3, .fg_trim + 1L, n - .fg_trim), "")
}

# Center window: drop the leading C, then the central .center_len residues.
# Returns "" when the (post-C) sequence is shorter than the window.
.center_window <- function(cdr3) {
  s <- substr(cdr3, 2L, nchar(cdr3))          # drop conserved leading C
  n <- nchar(s)
  ntrim <- (n - .center_len) %/% 2L
  ifelse(n < .center_len, "", substr(s, ntrim + 1L, ntrim + .center_len))
}

# ---------------------------------------------------------------------------
# Private per-feature scorers (each vectorised over a vector of TCR fields)
# ---------------------------------------------------------------------------

# cdr3len: alpha length + 2 * beta length (beta double-weighted).
.score_cdr3len <- function(cdr3a, cdr3b) nchar(cdr3a) + 2L * nchar(cdr3b)

# old_imhc: hand-tuned FG-loop score, beta double-weighted. CDR3 <= 8 aa -> 0.
.old_imhc_cdr3 <- function(cdr3) {
  fg <- .fg_loop(cdr3)
  cnt <- function(s, aa) lengths(regmatches(s, gregexpr(aa, s, fixed = TRUE)))
  ifelse(nchar(cdr3) <= 8L, 0,
         nchar(fg) + 3 * cnt(fg, "C") + 2 * cnt(fg, "W") + cnt(fg, "R") +
           cnt(fg, "K") + 0.5 * cnt(fg, "H") - cnt(fg, "D") - cnt(fg, "E"))
}
.score_old_imhc <- function(cdr3a, cdr3b) {
  .old_imhc_cdr3(cdr3a) + 2 * .old_imhc_cdr3(cdr3b)
}

# mait: invariant MAIT alpha chain (human/rhesus rules; TRAV1 for mouse).
.score_mait <- function(va, ja, cdr3a, organism = "human") {
  va <- .strip_allele(va); ja <- .strip_allele(ja)
  if (grepl("mouse", organism)) {
    hit <- va == "TRAV1" & grepl("^TRAJ33", ja) & nchar(cdr3a) == 12L
  } else {  # human / rhesus
    hit <- grepl("^TRAV1-2", va) &
      (grepl("^TRAJ33", ja) | grepl("^TRAJ20", ja) | grepl("^TRAJ12", ja)) &
      nchar(cdr3a) == 12L
  }
  as.numeric(hit)
}

# inkt: invariant NKT TCR (human/rhesus needs beta TRBV25; mouse uses alpha only).
.score_inkt <- function(va, ja, cdr3a, vb, organism = "human") {
  va <- .strip_allele(va); ja <- .strip_allele(ja); vb <- .strip_allele(vb)
  if (grepl("mouse", organism)) {
    hit <- grepl("^TRAV11", va) & grepl("^TRAJ18", ja) & nchar(cdr3a) == 15L
  } else {  # human / rhesus
    hit <- grepl("^TRAV10", va) & grepl("^TRAJ18", ja) &
      nchar(cdr3a) %in% c(14L, 15L, 16L) & grepl("^TRBV25", vb)
  }
  as.numeric(hit)
}

# alphadist: combined TRAV + TRAJ locus position (missing gene -> list midpoint).
.score_alphadist <- function(va, ja) {
  va <- .strip_allele(va); ja <- .strip_allele(ja)
  iv <- match(va, trav_order)                 # 1-based, NA if absent
  ij <- match(ja, traj_order)
  va_dist <- ifelse(is.na(iv), 0.5 * (length(trav_order) - 1L),
                    length(trav_order) - iv)   # conga: len-1 - index0
  ja_dist <- ifelse(is.na(ij), 0.5 * (length(traj_order) - 1L), ij - 1L)  # index0
  va_dist + ja_dist
}

# gene-presence: 1 if the chain's V/J (allele-stripped) equals `name`, else 0.
# The chain/region is inferred from the gene-name prefix (TRAV/TRAJ/TRBV/TRBJ).
.score_gene <- function(name, va, ja, vb, jb) {
  col <- switch(substr(name, 1L, 4L),
                TRAV = va, TRAJ = ja, TRBV = vb, TRBJ = jb,
                stop("gene '", name, "' must start with TRAV/TRAJ/TRBV/TRBJ"))
  as.numeric(.strip_allele(col) == name)
}

# Mean of a per-residue lookup over the FG-loop (mode "fg") or center (mode "cen")
# window of one CDR3. Empty window -> mean of the lookup over all 20 residues.
.window_mean <- function(cdr3, lookup, mode) {
  win <- if (mode == "cen") .center_window(cdr3) else .fg_loop(cdr3)
  vapply(win, function(w) {
    if (!nzchar(w)) return(mean(lookup))
    mean(lookup[strsplit(w, "", fixed = TRUE)[[1]]], na.rm = TRUE)
  }, numeric(1), USE.NAMES = FALSE)
}

# property_score for one CDR3: an aa_props column, an "<AA>_frac" indicator, or
# "arofrac" (fraction of aromatic F/Y/W/H; an extension beyond conga's list).
.property_cdr3 <- function(cdr3, score_name, mode) {
  if (nchar(score_name) == 6L && substr(score_name, 1L, 1L) %in% amino_acids &&
      substr(score_name, 2L, 6L) == "_frac") {
    aa <- substr(score_name, 1L, 1L)
    lookup <- setNames(as.numeric(amino_acids == aa), amino_acids)
  } else if (score_name == "arofrac") {
    lookup <- setNames(as.numeric(amino_acids %in% c("F", "Y", "W", "H")), amino_acids)
  } else {
    lookup <- setNames(aa_props[, score_name], amino_acids)
  }
  .window_mean(cdr3, lookup, mode)
}

# property_score for a TCR: weighted sum of the alpha and beta CDR3 scores.
.score_property <- function(cdr3a, cdr3b, score_name, mode,
                            alpha_weight = 1, beta_weight = 1) {
  alpha_weight * .property_cdr3(cdr3a, score_name, mode) +
    beta_weight * .property_cdr3(cdr3b, score_name, mode)
}

# ---------------------------------------------------------------------------
# Feature-name handling
# ---------------------------------------------------------------------------

# Fixed (non-property) feature keys implemented here.
.fixed_features <- c("cdr3len", "old_imhc", "mait", "inkt", "alphadist")

# Every property-style feature name (bare, plus explicit _fg / _cen variants).
.property_features <- local({
  base <- c(property_names, paste0(amino_acids, "_frac"), "arofrac")
  c(base, paste0(base, "_fg"), paste0(base, "_cen"))
})

# Canonical list of all directly-nameable features (gene-presence names, i.e.
# any TRAV/TRAJ/TRBV/TRBJ gene, are additionally accepted by tcr_score()).
all_tcr_features <- c(.fixed_features, .property_features)

# Is `name` a V/J gene-presence request?
.is_gene_feature <- function(name) grepl("^TR[AB][VJ]", name)

# Split a property feature name into (score_name, mode), matching conga's parser:
# a trailing _fg/_cen sets the mode, otherwise the default mode applies.
.parse_property <- function(name, default_mode) {
  parts <- strsplit(name, "_", fixed = TRUE)[[1]]
  tail <- parts[length(parts)]
  if (tail %in% .cdr3_score_modes) {
    list(score_name = paste(parts[-length(parts)], collapse = "_"), mode = tail)
  } else {
    list(score_name = name, mode = default_mode)
  }
}

# ---------------------------------------------------------------------------
# Master function
# ---------------------------------------------------------------------------

# Compute a conga-style TCR feature score table.
#
#   data          data.frame with one row per clonotype / cell
#   features      character vector of feature keys to compute (see
#                 all_tcr_features, plus any TRAV/TRAJ/TRBV/TRBJ gene name for a
#                 gene-presence indicator). AA-property features may carry a
#                 trailing "_fg" (FG loop, default) or "_cen" (central window).
#   cols          named vector mapping the six TCR fields to column names:
#                 c(cdr3a=, va=, ja=, cdr3b=, vb=, jb=). Defaults to VDJdb naming.
#   organism      "human" (default), "rhesus", or "mouse" (affects mait / inkt)
#   mode          default window for AA-property features ("fg" or "cen")
#   alpha_weight,
#   beta_weight   chain weights for AA-property features (conga defaults 1, 1)
#   keep_cols     optional column names from `data` to carry through unchanged
#                 (e.g. c("antigen.epitope","meta.tissue"))
#
# Returns a data.frame with one row per retained TCR (rows with an NA/empty
# alpha or beta CDR3 are dropped) and one column per requested feature, plus any
# keep_cols. Feature columns appear in the order requested.
tcr_score <- function(data,
                      features = c("cdr3len", "old_imhc", "mait", "inkt",
                                   "alphadist", "hydropathy", "charge", "volume"),
                      cols = c(cdr3a = "cdr3.alpha", va = "v.alpha", ja = "j.alpha",
                               cdr3b = "cdr3.beta", vb = "v.beta", jb = "j.beta"),
                      organism = "human",
                      mode = "fg",
                      alpha_weight = 1, beta_weight = 1,
                      keep_cols = NULL) {
  need <- c("cdr3a", "va", "ja", "cdr3b", "vb", "jb")
  stopifnot(all(need %in% names(cols)), all(cols %in% names(data)),
            mode %in% .cdr3_score_modes,
            is.null(keep_cols) || all(keep_cols %in% names(data)))

  # Validate every requested feature up front.
  bad <- features[!(features %in% all_tcr_features | vapply(features, .is_gene_feature, logical(1)))]
  if (length(bad)) {
    stop("unknown feature(s): ", paste(bad, collapse = ", "),
         "\nSee all_tcr_features, or pass a TRAV/TRAJ/TRBV/TRBJ gene name.")
  }

  # Pull the TCR fields, drop rows with an NA/empty CDR3 on either chain.
  cdr3a <- as.character(data[[cols[["cdr3a"]]]]); cdr3b <- as.character(data[[cols[["cdr3b"]]]])
  va <- as.character(data[[cols[["va"]]]]); ja <- as.character(data[[cols[["ja"]]]])
  vb <- as.character(data[[cols[["vb"]]]]); jb <- as.character(data[[cols[["jb"]]]])
  ok <- !is.na(cdr3a) & cdr3a != "" & !is.na(cdr3b) & cdr3b != ""
  data <- data[ok, , drop = FALSE]
  cdr3a <- cdr3a[ok]; cdr3b <- cdr3b[ok]
  va <- va[ok]; ja <- ja[ok]; vb <- vb[ok]; jb <- jb[ok]

  # Dispatch each requested feature to its scorer.
  score_one <- function(name) {
    if (name == "cdr3len")        return(.score_cdr3len(cdr3a, cdr3b))
    if (name == "old_imhc")       return(.score_old_imhc(cdr3a, cdr3b))
    if (name == "mait")           return(.score_mait(va, ja, cdr3a, organism))
    if (name == "inkt")           return(.score_inkt(va, ja, cdr3a, vb, organism))
    if (name == "alphadist")      return(.score_alphadist(va, ja))
    if (.is_gene_feature(name))   return(.score_gene(name, va, ja, vb, jb))
    p <- .parse_property(name, mode)   # AA-property / _frac / arofrac
    .score_property(cdr3a, cdr3b, p$score_name, p$mode, alpha_weight, beta_weight)
  }
  out <- as.data.frame(lapply(features, score_one), stringsAsFactors = FALSE)
  names(out) <- features

  if (!is.null(keep_cols)) out <- cbind(data[, keep_cols, drop = FALSE], out)
  rownames(out) <- NULL
  out
}
