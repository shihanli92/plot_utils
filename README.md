# plot_utils

R and Python utilities for analysing and plotting T-cell receptor (TCR)
repertoire data. Standalone and self-contained — no external package required.

## Contents

| Script | Language | Purpose |
| --- | --- | --- |
| `scripts/extract_vdjdb_paired.py` | Python (pandas, requests) | Download the latest [VDJdb](https://github.com/antigenomics/vdjdb-db) release, extract fully-paired human TCR records (both `cdr3.alpha` and `cdr3.beta` populated), and write `vdjdb_paired.csv`. Falls back to the slim table and reconstructs pairs by `complex.id` when no `*_full.txt` is present. |
| `scripts/tcr_analysis/chain_usage.R` | R (dplyr, ggplot2) | V/J gene-segment usage frequencies per group; boxplot + points, faceted by segment. |
| `scripts/tcr_analysis/gene_diversity.R` | R (dplyr, ggplot2) | Gene-usage diversity per group: richness, Shannon entropy, evenness, inverse Simpson, Gini; bar (mean) + points, faceted by segment. |
| `scripts/tcr_analysis/cdr3_composition.R` | R (dplyr, ggplot2, ggseqlogo) | Compare CDR3 sub-sequence amino-acid composition across groups; position x AA heatmap + per-group sequence logo. |
| `scripts/tcr_analysis/cdr3_length.R` | R (dplyr, ggplot2) | Compare CDR3 length (aa) distributions across groups; per-group length-frequency table + summary stats, frequency-polygon plot. |
| `scripts/tcr_analysis/aa_position_compare.R` | R (dplyr, ggplot2) | Compare the amino-acid (k-mer) at one CDR3 position across groups; heatmap of all groups, or a two-group g1-vs-g2 frequency scatter. |

All three R scripts share the same options: a `group_cols` vector defining the
comparison unit, a `gene_cols` / `cdr3_col` selector, and a `unique_by` argument
(NULL to count every row, or a vector of clone-defining columns to collapse
duplicate clones so each unique clone counts once per group).

## Usage

### Extract paired TCRs from VDJdb

```bash
pip install pandas requests
python scripts/extract_vdjdb_paired.py
```

### Analyses

```r
data <- read.csv("vdjdb_paired.csv", check.names = FALSE)

source("scripts/tcr_analysis/chain_usage.R")
u <- tcr_chain_usage(data, group_cols = "antigen.epitope")
u$usage; u$plot

source("scripts/tcr_analysis/gene_diversity.R")
d <- tcr_gene_diversity(data, x_col = "antigen.epitope",
                        group_cols = "reference.id", metric = "gini")
d$diversity; d$plot

source("scripts/tcr_analysis/cdr3_composition.R")
c3 <- tcr_cdr3_composition(data, group_cols = "antigen.epitope", positions = 3:5)
c3$composition; c3$heatmap; c3$logo

source("scripts/tcr_analysis/cdr3_length.R")
l <- tcr_cdr3_length(data, group_cols = "antigen.epitope", cdr3_col = "cdr3.beta")
l$lengths; l$summary; l$plot

source("scripts/tcr_analysis/aa_position_compare.R")
# heatmap of the amino acid at position 5 across all epitopes
p <- tcr_aa_position_compare(data, group_cols = "antigen.epitope", position = 5)
p$composition; p$plot
# two-group scatter of doublets (k = 2) at position 5
p2 <- tcr_aa_position_compare(data, group_cols = "antigen.epitope", position = 5,
                              k = 2, groups = c("epitopeA", "epitopeB"), top_n = 30)
p2$plot
```
