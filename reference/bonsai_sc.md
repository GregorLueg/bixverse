# Build a Bonsai tree over the cells

Bonsai reconstructs a tree over the cells in which every cell is a leaf
and the internal nodes are inferred ancestral states, with branch
lengths that carry the amount of change between them. Unlike a kNN graph
it uses each measurement's error bar, which is what Sanity provides: the
raw counts go through Sanity first, for posterior log fold changes with
error bars, and Bonsai builds the tree on those. The tree is then laid
out in 2D. For details, please refer to de Groot, et al. and Breda, et
al.

No HVG selection needed. By default every gene goes in, streamed through
Sanity in chunks, and only the genes with enough signal over their own
noise (`min_signal_to_noise` in
[`params_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/params_sc_bonsai.md))
are kept for the tree. That is the gene selection of the Bonsai paper,
and it keeps memory at one chunk plus the survivors. Pass `hvg` to
restrict the candidates.

Runtime grows a little faster than linearly with the number of cells. In
bonsai-rs's own benchmarks the search took about 70 seconds at 10,000
cells and 210 seconds at 25,000 on ten cores. Sanity comes on top,
linear in the number of genes it has to fit; on the CPU that is the
larger share once all genes go in.

The object itself is not touched. Store the leaf coordinates with
[`set_bonsai_embedding()`](https://gregorlueg.github.io/bixverse/reference/set_bonsai_embedding.md)
if you want them next to the other embeddings.

On `MetaCells` every metacell is a leaf. Their raw counts are sums over
their cells, and a sum of Poisson counts is Poisson, so Sanity treats a
metacell exactly as it treats a cell: its total counts are the library
size, and a bigger metacell gets tighter error bars, which Bonsai
weights by. That is the way to large data: 100k cells in metacells of 50
is a 2,000-leaf tree. What you give up is any structure inside a
metacell.

## Usage

``` r
bonsai_sc(
  object,
  hvg = NULL,
  bonsai_params = params_sc_bonsai(),
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `MetaCells` class.

- hvg:

  Optional integer. Restrict the candidate genes to these, e.g. the
  output of
  [`get_hvg()`](https://gregorlueg.github.io/bixverse/reference/get_hvg.md)
  plus one. Please provide 1-indexed genes here! If `NULL`, every gene
  in the object is a candidate.

- bonsai_params:

  List. See
  [`params_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/params_sc_bonsai.md).

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A `BonsaiTree` S3 object with:

- nodes - data.table with `node`, `parent` (`NA` for the root),
  `branch`, `is_leaf`, `cell_id` (`NA` for the inferred ancestors),
  `n_cells` (the cells behind a leaf: 1 for single cells, the
  originating cells for a metacell; `NA` for the ancestors), `x` and
  `y`. The leaves come first, in the order of the object's cells or
  metacells.

- loglik - The loglikelihood of the final tree.

- steps - data.table with the loglikelihood after each search step and
  the step's wall time in seconds.

- timings - data.table with the wall time in seconds of each stage:
  `sanity`, `ingest`, `bonsai` (the whole search), `layout`, and
  `total`, the whole call as R saw it.

- genes_used - The genes the tree was built on.

- genes_dropped - The candidate genes left out: no counts in the cells,
  ill-conditioned Sanity posteriors, or a signal-to-noise ratio below
  `min_signal_to_noise`.

- cell_idx - The cells or metacells the tree was built over (0-indexed).

- layout, hyperbolic - The current layout.

- params - The parameters of the run.

## References

de Groot, et al., Nat. Biotechnol., 2026; Breda, et al., Nat.
Biotechnol., 2021.

## Examples

``` r
# a tree over the demo cells, genes selected by their signal-to-noise
sc <- demo_single_cells(prepped = FALSE)
tree <- bonsai_sc(sc, .verbose = FALSE)
tree
#> BonsaiTree: 500 leaves, 496 inferred ancestors
#>   Genes: 50 used, 0 dropped
#>   Loglikelihood: -5619.8
#>   Layout: equal_angle
#>   Seconds: sanity 0.2 | ingest 0.0 | bonsai 0.9 | layout 0.0 | total 1.2

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
