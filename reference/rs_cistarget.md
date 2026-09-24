# Run CisTarget motif enrichment analysis

**\[experimental\]** Core Rust function for motif enrichment analysis
using recovery curves. Gene sets are processed in parallel; only motifs
with an NES at or above `nes_threshold` are returned.

## Usage

``` r
rs_cistarget(
  rankings,
  gs_list,
  auc_threshold,
  nes_threshold,
  max_rank,
  method,
  n_mean,
  verbose
)
```

## Arguments

- rankings:

  Integer matrix with motif rankings for genes (genes in rows, motifs in
  columns). Lower ranks indicate higher regulatory potential.

- gs_list:

  List of integer vectors. Each element contains 1-based indices of
  genes in the gene set (matching row indices in rankings).

- auc_threshold:

  Integer. Absolute number of top-ranked genes to use for the AUC
  calculation (e.g., for 5% of 10000 genes, use 500).

- nes_threshold:

  Numeric. Normalised Enrichment Score threshold for filtering
  significant motifs.

- max_rank:

  Integer. Maximum rank to consider for the recovery curves (at most
  `nrow(rankings)`).

- method:

  String. Recovery curve calculation method, one of
  `c("approx", "icistarget")`. Anything else falls back to `"approx"`.

- n_mean:

  Integer. Window size for the smoothing in the approximate method.

- verbose:

  Boolean. Report progress per decile of gene sets.

## Value

List of lists, one per gene set, each containing

- motif_idx - 1-based column index of the motif in `rankings`.

- nes - Normalised enrichment score.

- auc - Area under the recovery curve.

- rank_at_max - Rank at which the leading edge is reached.

- n_enriched - Number of genes in the leading edge.

- leading_edge - List of 1-based row indices of the leading edge genes,
  one element per motif.
