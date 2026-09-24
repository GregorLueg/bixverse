# Calculate VISION pathway scores in Rust with auto-correlation

**\[experimental\]** The function will take in a list of gene sets that
contains lists of `"pos"` and `"neg"` gene indices (0-indexed). You
don't have to provide the `"neg"`, but it can be useful to classify the
delta of two stats (EMT, Th1; Th2) etc. Additionally, it will take a
random gene list and calculate an auto-correlation score based on
Geary's C to identify pathways that show significant patterns on the kNN
graph generated on the provided embedding.

## Usage

``` r
rs_vision_with_autocorrelation(
  f_path,
  embd,
  knn_data,
  gs_list,
  random_gs_list,
  vision_params,
  cells_to_keep,
  cluster_membership,
  streaming,
  verbose,
  seed
)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- embd:

  Numerical matrix. The embedding matrix to use to generate the kNN
  graph. Rows must align with `cells_to_keep`.

- knn_data:

  Optional list. This contains pre-computed kNN data (`indices`
  (0-indexed), `dist`, `dist_metric` and `k`). The user has to ensure
  consistency! If provided, this will be used rather than a graph built
  from the parameter list.

- gs_list:

  Nested list. Each sublist contains the (0-indexed!) `pos` and `neg`
  gene indices of that specific gene set.

- random_gs_list:

  Double-nested list. The outer list represents the gene set clusters
  and the inner list the permuted gene sets (same structure as
  `gs_list`) of that cluster.

- vision_params:

  List. Contains various parameters to use in terms of the kNN
  generation.

- cells_to_keep:

  Integer. Vector of indices (0-indexed) of the cells to keep.

- cluster_membership:

  Integer. Vector that indicates to which of the permuted gene set
  clusters (1-indexed) the given gene set belongs.

- streaming:

  Boolean. Shall the data be streamed.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

- seed:

  Integer. Random seed for reproducibility.

## Value

A list with the following items:

- autocor_res - List with `auto_cor` (1 - Geary's C), `p_val`
  (empirical, against the permuted gene sets of the same cluster) and
  `fdr`.

- vision_mat - A matrix of cells x vision scores per gene set.
