# Run CisTarget motif enrichment analysis

Core Rust function for motif enrichment analysis using recovery curves.

## Usage

``` r
rs_cistarget(
  rankings,
  gs_list,
  auc_threshold,
  nes_threshold,
  max_rank,
  method,
  n_mean
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

  Absolute number of top-ranked genes to use for AUC calculation (e.g.,
  for 5% of 10000 genes, use 500).

- nes_threshold:

  Normalised Enrichment Score threshold for filtering significant motifs

- max_rank:

  Maximum rank to consider (typically nrow(rankings)).

- method:

  Recovery curve calculation method: "approx" or "icistarget".

- n_mean:

  Number of points for averaging in approximate method.

## Value

List of lists, one per gene set, each containing motif_idx, nes, auc,
rank_at_max, n_enriched, and leading_edge.
