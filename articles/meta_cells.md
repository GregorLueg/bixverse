# Meta cells with bixverse

## Intro

Meta cells are small groups of transcriptionally similar cells
aggregated into a single representative profile. The motivation is
twofold: firstly, single cell counts are sparse and noisy, and many
downstream methods (correlation networks, GRN inference, archetype-based
modelling) behave poorly on raw single cell data. Aggregating into meta
cells reduces sparsity while preserving the heterogeneity that bulk
pseudo-replicates would average out. Secondly, they reduce the
computational burden substantially.

`bixverse` provides three meta cell algorithms with different underlying
ideas:

- **hdWGCNA-style meta cells** ([Morabito, et al.,
  2023](https://doi.org/10.1016/j.crmeth.2023.100498)) iterate over a
  kNN graph, picking seed cells and aggregating their neighbours, with a
  constraint on how many cells two meta cells are allowed to share.
  Simple, fast and reasonably effective. Works directly on a kNN graph.

- **SEACells** ([Persad, et al.,
  2023](https://doi.org/10.1038/s41587-023-01716-9)) uses kernel
  archetypal analysis. The algorithm finds a set of archetypes (the
  SEACells) such that every cell can be expressed as a convex
  combination of nearby archetypes. Tends to produce purer aggregations
  but is the most computationally expensive of the three.

- **SuperCells** ([Bilous, et al.,
  2022](https://doi.org/10.1186/s12859-022-04861-1)) applies the
  walktrap community detection algorithm to the kNN graph and treats
  each community as a meta cell. Granularity is controlled indirectly
  through a graining factor (the average number of cells per meta cell).

One implementation detail worth flagging upfront. In `bixverse`, the
`SingleCells` class keeps counts on disk in DuckDB and binary Rust-based
files and streams them in as needed. The `MetaCells` class does not.
After aggregation the count matrix is small enough to hold entirely in
memory as a sparse matrix, and all downstream operations (HVG, PCA,
neighbours, clustering, embeddings) run against this in-memory
representation (most of the code still makes usage of Rust’s incredible
performance nonetheless). Practically that means the meta cell pipeline
is substantially faster than the equivalent operations on the parent
`SingleCells` object for two reasons: data is already in memory and more
importantly, it is simply less data.

``` r

library(bixverse)
library(bixverse.plots)
library(ggplot2)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
library(patchwork)
library(magrittr)
```

## Preparing the single cell data

We use the same CD34 cells from the [SEACells
vignette](https://github.com/dpeerlab/SEACells/blob/main/notebooks/SEACell_computation.ipynb)
to set up the parent `SingleCells` object: load (QC not needed - it’s
already filtered data), HVG selection, PCA, kNN graph. We can use the
provided cell type labels to check purity.

``` r

cd34_path <- download_cd34_data()

tempdir_cd34 <- tempdir()

sc_object <- SingleCells(dir_data = tempdir_cd34)

sc_object <- load_h5ad(object = sc_object, h5_path = cd34_path)
#>  Using light streaming for the CSR to CSC conversion.
#> Loading observations data from h5ad into the DuckDB.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Loading variables data from h5ad into the DuckDB.
#> 
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
```

Let’s run HVG detection, PCA and kNN (+ sNN) generation.

``` r

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = 2000L,
  .verbose = FALSE
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L,
  .verbose = FALSE
)

sc_object <- find_neighbours_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(ann_dist = "euclidean", knn_method = "kmknn")
  ),
  .verbose = FALSE
)

# we need the kNN object for the diffusion maps
knn_object <- get_knn_obj(x = sc_object)
```

## Generating meta cells

All three generators take a `SingleCells` object, read counts from the
binary files as needed, and return a `MetaCells` object. They share a
common output structure: an obs table with one row per meta cell
(recording which original cells went into it), a var table mirrored from
the parent object, and an in-memory sparse count matrix.

### hdWGCNA-style (bootstrapped metacells)

The hdWGCNA-style algorithm reuses the kNN graph already on the object.
The two parameters that matter most are `target_no_metacells` and
`max_shared`, which controls how many cells two meta cells may have in
common. The original paper uses `1000L` meta cells, but we want to also
compare the different metrics downstream, so, we will set it to `250L`.
We also recalculate the kNN graph with larger k, to include more cells
per given meta cell.

``` r

hdwgcna <- generate_bt_meta_cells_sc(
  object = sc_object,
  sc_meta_cell_params = params_sc_bt_metacells(
    target_no_metacells = 250L,
    knn = list(k = 25L)
  ),
  regenerate_knn = TRUE, # regenerate kNN graph with 25 neighbours
  .verbose = TRUE
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

hdwgcna
#> Single cell experiment (Meta Cells).
#>   Meta cell method: meta_cells_hdwgcna
#>   Merged: FALSE
#>   No meta cells: 250
#>   No genes: 12464
#>   No cells aggregated: 3973
#>   No obs rows in source: 6881
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
#>   Stale artefacts: none
```

Importantly, not all cells HAVE to be assigned to a meta cell and cells
can occur in several meta cells at once. We can appreciate here that a
number of cells remains unassigned and each MetaCell does contain self +
neighbours, i.e., 26 (`self + k = 25L`) original cells. Let’s check how
often the cells occur… ?

``` r

no_duplicated_cells <- table(unlist(hdwgcna[[]]$original_cell_idx))

hist(
  no_duplicated_cells,
  xlab = "No of times a cell occurs",
  main = "No of times a cell is part of a meta cell",
  breaks = 10L
)
```

![](meta_cells_files/figure-html/hdwgcna%20-%20cells%20per%20meta%20cell-1.png)

As we can appreciate, some of the cells are indeed shared across meta
cells. Let’s add the diffusion coordinates to this one for later
analysis

``` r

hdwgcna <- calc_diffusion_coordinates(object = hdwgcna, knn_data = knn_object)
```

### SEACells

SEACells operates on the PCA embedding (or any batch corrected
embedding) and the kNN graph and runs an iterative kernel archetypal
analysis. `n_sea_cells` sets the number of archetypes;
`min_iter`/`max_iter` bound the optimisation; `convergence_epsilon` sets
the early stopping threshold relative to the initial RSS. The algorithm
was aggressively optimised to be faster. With `pruning = TRUE` (default)
small values are set to zero, removing basically dust and accelerating
the calculations substantially. The rule of thumb by the authors is one
SEACell per 75 cells, i.e., 90 SEACells for this experiment. We will set
this higher here to `250L` to make this a bit more comparable with the
other methods…

``` r

seacells <- generate_seacells_sc(
  object = sc_object,
  seacell_params = params_sc_seacells(
    n_sea_cells = 250L,
    min_iter = 10L,
    convergence_epsilon = 0.001
  ),
  .verbose = TRUE
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

seacells
#> Single cell experiment (Meta Cells).
#>   Meta cell method: seacell
#>   Merged: FALSE
#>   No meta cells: 250
#>   No genes: 12464
#>   No cells aggregated: 6881
#>   No obs rows in source: 6881
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
#>   Stale artefacts: none
```

We can appreciate that each SEACell can contain a varying degree of
original cells contributing to this specific metacell.

``` r

no_cells_seacells <- purrr::map_dbl(seacells[[]]$original_cell_idx, length)

hist(
  no_cells_seacells,
  xlab = "No cells per meta cell",
  main = "Cells per SEACells meta cell",
  breaks = 25L
)
```

![](meta_cells_files/figure-html/seacells%20-%20cells%20per%20meta%20cell-1.png)

We will also add the diffusion coordinates to this one

``` r

seacells <- calc_diffusion_coordinates(object = seacells, knn_data = knn_object)
```

#### SEACells on large data

SEACells is infamous for being very slow and difficult to run on larger
data sets. Most of the hot path has been rewritten in `bixverse`, which
is what makes the algorithm feasible beyond a few tens of thousands of
cells:

1.  The initialisation is the part of the original algorithm that scales
    the worst. The greedy column subset selection needs a full `K^2`
    column for every candidate cell, which is `O(N^2)` and simply will
    not finish on millions of cells. `bixverse` therefore switches
    strategy by size: above `greedy_threshold` (20,000 cells by default)
    it drops the greedy top-up and samples the archetypes at random,
    which is effectively free. When a better-than-random start is wanted
    but the greedy pass is still out of reach, setting `n_landmarks`
    enables a Nyström route: a small set of density-weighted landmarks
    is chosen (5 to 10 times `n_sea_cells` is a sensible range), the
    diffusion operator is built and eigendecomposed on those landmarks
    alone, and the multiscale embedding is carried to all cells via a
    Nyström extension before the usual max-min waypoint sampling. The
    eigendecomposition stays at landmark scale `L x L` instead of
    `N x N`.
2.  The optimisation loop itself runs in bounded memory. `K^2` is never
    formed; every `K^2 X` term is evaluated as `K (K X)`, so memory
    stays bounded by the number of non-zeros in `K` rather than `N^2`.
    Neither Frank-Wolfe gradient is ever materialised. The A update runs
    cell-major: a column of `A` is a convex combination of at most
    `max_fw_iters` one-hot atoms, so the solver carries
    `(index, weight)` pairs and patches the gradient with rank-1
    corrections instead of rebuilding it from the sparse `A` on every
    iteration. The argmin over the two gradient terms is fused and runs
    on SIMD (NEON on Apple Silicon, SSE2, AVX-2 or AVX-512 picked at
    runtime on x86). The B update walks one gradient column per
    archetype in the same spirit. Everything is chunked across rayon
    threads.
3.  `K^2 B` is the one product the inner loop needs on every single
    iteration, so it is cached rather than recomputed. A Frank-Wolfe
    step changes `B` by a rank-k update, and the matching change to
    `K^2 B` is the same scaling plus one weighted column of `K^2` per
    archetype; the pruning corrections fold into the same delta. The
    update is exact, and a full recompute still runs every eight
    iterations to bound floating point drift and the sparsity pattern.
    The B loop stops early once the relative Frank-Wolfe duality gap
    drops below `1e-3`.
4.  `pruning` is on by default and should stay on. The Frank-Wolfe
    updates only ever add atoms, so without pruning the non-zeros in `A`
    and `B` climb monotonically across the whole fit and the time per
    iteration keeps rising instead of settling. The accuracy cost sits
    in `pruning_threshold`, not in the flag. Keep it below the smallest
    weight the schedule can produce, `2 / (T (T + 1))` for
    `T = max_fw_iters`; above that you start deleting live mass and
    shift the solution. The default `1e-7` removes numerical dust only.
5.  The RSS for the convergence check never materialises the `N x N`
    reconstruction `K B A`. The squared residual is expanded with the
    trace identity and evaluated through cyclic reordering, so every
    intermediate is at most `N x k` or `k x k`. Same Frobenius norm,
    just without the dense reconstruction.

Most of that is automatic. The knobs you actually reach for on large
data are `greedy_threshold` and `n_landmarks` for the initialisation,
and `pruning`/`pruning_threshold` for the loop. Together they are what
made it possible to run a million cells locally, which the original
implementation cannot do.

### SuperCells

SuperCells runs walktrap on the kNN graph. The number of meta cells is
set indirectly through `graining_factor`: with a factor of 30 over ~6800
cells you should get roughly ~230 meta cells. The option to run the
kernel-based version of SuperCells 2.0 from [Hérault et
al.](https://www.biorxiv.org/content/10.64898/2026.02.19.706848v1) is
enabled by default. (The multi-modal versions are still to come. Watch
the space…)

``` r

supercells <- generate_supercells_sc(
  object = sc_object,
  sc_supercell_params = params_sc_supercell(
    graining_factor = 30
  ),
  .verbose = TRUE
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

supercells
#> Single cell experiment (Meta Cells).
#>   Meta cell method: supercells
#>   Merged: FALSE
#>   No meta cells: 230
#>   No genes: 12464
#>   No cells aggregated: 6881
#>   No obs rows in source: 6881
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
#>   Stale artefacts: none
```

Let’s plot the number of cells per SuperCell.

``` r

no_cells_supercells <- purrr::map_dbl(supercells[[]]$original_cell_idx, length)

hist(
  no_cells_supercells,
  xlab = "No cells per meta cell",
  main = "Cells per Supercell",
  breaks = 25L
)
```

![](meta_cells_files/figure-html/supercells%20-%20cells%20per%20meta%20cell-1.png)

Similar to SEACells, we get a gradient here.

``` r

supercells <- calc_diffusion_coordinates(
  object = supercells,
  knn_data = knn_object
)
```

#### SuperCells on large data

SuperCell is already far cheaper than SEACells on large data: there is
no archetypal analysis, just Walktrap community detection on the kNN
graph. The one part that does not scale naively is the random-walk
representation. `compute_walk_probabilities()` within the Rust code
gives every cell a sparse vector of landing probabilities for a
`walk_length`-step walk; on a well-connected graph these vectors spread
to touch a large fraction of cells after only a few steps, so the store
drifts towards dense and memory grows as `O(n × support)`.

`max_support` caps this. With `max_support = k` each initial walk vector
keeps only its `k` largest entries (the dropped tail holds negligible
walk mass), bounding the store at roughly `k × n` regardless of how far
the walks spread. This makes the run an approximation: the Ward merge
criterion is driven by distances between walk vectors, and truncating
them can shift the merge order slightly, so you may get marginally
different communities than the exact run. With `max_support = NULL` (the
default) the walks are kept exact. Beyond saving memory, a smaller `k`
also speeds up every distance and merge operation, since both are linear
in the support.

``` r

supercells_large <- generate_supercells_sc(
  object = sc_object,
  sc_supercell_params = params_sc_supercell(
    graining_factor = 30,
    max_support = 256L # would bound the walk vectors to 256
  ),
  .verbose = TRUE
)
```

## Metrics

### Purity

A simple sanity check is to ask, for each meta cell, what fraction of
its constituent cells share the same underlying label. We use the Leiden
clusters from the parent `SingleCells` object as a proxy for cell type
identity. The labels themselves are imperfect, so absolute numbers
should be read with that caveat, but relative differences between
methods are reasonably indicative.

``` r

# memberships are positions in the full obs table, so the labels have to come
# from the unfiltered obs rather than from `sc_object[[...]]`
cell_labels <- as.character(get_sc_obs(sc_object)$celltype)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

hdwgcna <- calc_meta_cell_purity(hdwgcna, original_cell_type = cell_labels)
seacells <- calc_meta_cell_purity(seacells, original_cell_type = cell_labels)
supercells <- calc_meta_cell_purity(
  supercells,
  original_cell_type = cell_labels
)

purity_dt <- rbind(
  data.table(method = "hdWGCNA", purity = hdwgcna[[]]$mc_purity),
  data.table(method = "SEACells", purity = seacells[[]]$mc_purity),
  data.table(method = "SuperCells", purity = supercells[[]]$mc_purity)
)

purity_dt[, .(mean = mean(purity), median = median(purity)), by = method]
#>        method      mean    median
#>        <char>     <num>     <num>
#> 1:    hdWGCNA 0.8838462 0.9615385
#> 2:   SEACells 0.8954497 0.9861111
#> 3: SuperCells 0.8944299 0.9784323
```

``` r

ggplot(purity_dt, aes(x = method, y = purity)) +
  geom_violin(aes(fill = method), alpha = 0.6) +
  geom_boxplot(width = 0.05, outlier.size = 0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Method", y = "Meta cell purity (Leiden)") +
  ylim(0, 1)
```

![](meta_cells_files/figure-html/purity%20plot-1.png)

In terms of cell type purity, the methods are basically the same.
Potentially more interesting is their behaviour in the manifold, see
below…

### Manifold regions

We called the
[`calc_diffusion_coordinates()`](https://gregorlueg.github.io/bixverse/reference/calc_diffusion_coordinates.md)
with the kNN graph on the data. This will tell us which parts of the
manifold are being sampled here, i.e., for example high, medium or low
density regions (based on the diffusion map and distance to the 150k-th
neighbour within that). This gives an idea how many rare cell states the
meta cells capture.

``` r

hdwgcna_regions <- as.data.table(
  table(hdwgcna[[]]$density_region) / nrow(hdwgcna[[]])
)[, method := "hdWGCNA"]
seacells_regions <- as.data.table(
  table(seacells[[]]$density_region) / nrow(seacells[[]])
)[, method := "SEACells"]
supercells_regions <- as.data.table(
  table(supercells[[]]$density_region) / nrow(supercells[[]])
)[, method := "SuperCells"]

region_dt <- rbind(
  hdwgcna_regions,
  seacells_regions,
  supercells_regions
)[, V1 := factor(V1, levels = c("high", "mid", "low"))]
setnames(region_dt, old = c("V1", "N"), new = c("region", "proportion"))
```

``` r

ggplot(region_dt, aes(x = method, y = proportion, fill = region)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Method", y = "Proportion of meta cells", fill = "Density region") +
  scale_fill_manual(
    values = setNames(
      c("#2c2d54", "#969bc7", "#6f9954"),
      c("high", "mid", "low")
    )
  )
```

![](meta_cells_files/figure-html/regions%20plot-1.png)

This is one of the most important differences between the methods.
hdWGCNA will disproportionally sample regions of the manifolds with
medium density and basically proportion-based. SEACells on the other
hand generates meta cells that capture sparse (more heterogenous)
regions of the manifold. SuperCells sits between the two methods.

### Compactness and separation

We can also assess other metrics on the manifold representation.
Compactness (how close are the cells within a given metacell to the
centroid) and separation (how far is the closest centroid). The first
gives an (unbiased) indication of purity; the latter of diversity
captured in the manifold. We can do this trivially via:

``` r

hdwgcna <- calc_manifold_metrics(hdwgcna)
seacells <- calc_manifold_metrics(seacells)
supercells <- calc_manifold_metrics(supercells)
```

Let’s generate a plotting data.table

``` r

metrics_dt <- rbind(
  data.table(
    method = "hdWGCNA",
    separation = hdwgcna[[]]$separation,
    compactness = hdwgcna[[]]$compactness,
    region = hdwgcna[[]]$density_region
  ),
  data.table(
    method = "SEACells",
    separation = seacells[[]]$separation,
    compactness = seacells[[]]$compactness,
    region = seacells[[]]$density_region
  ),
  data.table(
    method = "SuperCells",
    separation = supercells[[]]$separation,
    compactness = supercells[[]]$compactness,
    region = supercells[[]]$density_region
  )
)[, region := factor(region, levels = c("high", "mid", "low"))]

metrics_dt[,
  .(
    mean_separation = mean(separation),
    median_separation = median(separation),
    mean_compactness = mean(compactness),
    median_compactness = median(compactness)
  ),
  .(method)
]
#>        method mean_separation median_separation mean_compactness
#>        <char>           <num>             <num>            <num>
#> 1:    hdWGCNA       0.2230430         0.1587665       0.03730329
#> 2:   SEACells       0.3199538         0.2855716       0.04906456
#> 3: SuperCells       0.3138431         0.2789584       0.03011156
#>    median_compactness
#>                 <num>
#> 1:         0.02216349
#> 2:         0.02163796
#> 3:         0.01738179
```

Some differences are visible here. hdWGCNA has the worst separation and
best compactness (not unexpected given the bootstrapping method only
sampling direct neighbours). SEACells has the worst compactness across
the three methods (not unsurprising as the meta cells here capture more
sparse regions of the manifold) and middling separation; SuperCells
(version 2.0, the default setting here) is in the middle.

``` r

per_region_stats <- metrics_dt[,
  .(
    mean_separation = mean(separation),
    median_separation = median(separation),
    mean_compactness = mean(compactness),
    median_compactness = median(compactness)
  ),
  .(method, region)
]
setorder(per_region_stats, method, region)

per_region_stats[]
#>        method region mean_separation median_separation mean_compactness
#>        <char> <fctr>           <num>             <num>            <num>
#> 1:   SEACells   high      0.10425718        0.05595133       0.04914284
#> 2:   SEACells    mid      0.24700755        0.24671055       0.01891915
#> 3:   SEACells    low      0.44994459        0.40122062       0.07254122
#> 4: SuperCells   high      0.08777602        0.08696582       0.01144669
#> 5: SuperCells    mid      0.25622191        0.24006854       0.01325881
#> 6: SuperCells    low      0.45096919        0.43743290       0.05507259
#> 7:    hdWGCNA   high      0.06090325        0.05360688       0.02073870
#> 8:    hdWGCNA    mid      0.18660221        0.17095762       0.02184195
#> 9:    hdWGCNA    low      0.49851019        0.41214153       0.09340455
#>    median_compactness
#>                 <num>
#> 1:         0.01400404
#> 2:         0.01440683
#> 3:         0.05276826
#> 4:         0.01103607
#> 5:         0.01098484
#> 6:         0.05249530
#> 7:         0.02045949
#> 8:         0.01897707
#> 9:         0.07422031
```

We can also look at this per density region in the manifold where we can
observe that compactness is smallest in the high density regions,
whereas separation is the highest in the low density regions. And below
as plots:

#### Compactness plots

``` r

ggplot(metrics_dt, aes(x = method, y = compactness)) +
  geom_violin(aes(fill = method), alpha = 0.6) +
  geom_boxplot(width = 0.05, outlier.size = 0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Method", y = "Compactness")
```

![](meta_cells_files/figure-html/compactness%20plot-1.png)

``` r

ggplot(metrics_dt, aes(x = method, y = compactness)) +
  geom_boxplot(aes(fill = region), outlier.size = 0.5) +
  theme_bw() +
  labs(x = "Method", y = "Compactness") +
  scale_fill_manual(
    values = setNames(
      c("#2c2d54", "#969bc7", "#6f9954"),
      c("high", "mid", "low")
    )
  )
```

![](meta_cells_files/figure-html/compactness%20per%20region-1.png)

#### Separation plots

``` r

ggplot(metrics_dt, aes(x = method, y = separation)) +
  geom_violin(aes(fill = method), alpha = 0.6) +
  geom_boxplot(width = 0.05, outlier.size = 0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Method", y = "Separation")
```

![](meta_cells_files/figure-html/separation%20plot-1.png)

``` r

ggplot(metrics_dt, aes(x = method, y = separation)) +
  geom_boxplot(aes(fill = region), outlier.size = 0.5) +
  theme_bw() +
  labs(x = "Method", y = "Compactness") +
  scale_fill_manual(
    values = setNames(
      c("#2c2d54", "#969bc7", "#6f9954"),
      c("high", "mid", "low")
    )
  )
```

![](meta_cells_files/figure-html/separation%20per%20region-1.png)

Based on this, we can conclude the following:

- Purity is VERY similar across the methods (for this data set).
- At matched K (~250 metacells), SuperCells achieves the lowest median
  compactness, slightly ahead of hdWGCNA and SEACells. Separation is
  comparable between SuperCells and SEACells, with hdWGCNA noticeably
  lower.
- The methods differ markedly in where they place metacells. SEACells
  assigns 50% of metacells to low-density regions versus 24% for hdWGCNA
  and 41% for SuperCells, consistent with its kernel archetypal
  formulation seeking out underrepresented states.
- Method choice should reflect the downstream question. For tight
  aggregation in dense regions (e.g., bulk-like pseudo-replicates of
  common cell types), SuperCells is efficient and effective. For
  preserving rare populations or transitioning states (e.g.,
  differentiation trajectories), SEACells’ rare-state bias is desirable,
  with the trade-off of more heterogeneous metacells in those regions.
  hdWGCNA is fast and simple but assignments overlap more in embedding
  space.

## Working with the MetaCells class

This is where the in-memory representation pays off. HVG selection, PCA,
neighbour graphs and embeddings on a few hundred meta cells take seconds
and do not touch the Rust binary files on disk. We will just take
forward SEACells here, but the methods below work across all of
`MetaCells`.

### HVGs and PCA

The same HVG and PCA method dispatches you know from single cells work
here…

``` r

seacells <- find_hvg_sc(
  object = seacells,
  hvg_no = 2000L,
  .verbose = FALSE
)

seacells <- calculate_pca_sc(
  object = seacells,
  no_pcs = 30L
)
```

### Neighbours, clustering and UMAP

And also the same interfaces for neighbours, (Leiden) clustering and
UMAP can be used here.

``` r

seacells <- find_neighbours_sc(
  object = seacells,
  neighbours_params = params_sc_neighbours(
    knn = list(k = 10L, knn_method = "exhaustive")
  )
)
#> 
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.

seacells <- find_clusters_sc(seacells, res = 1.0, name = "leiden_clusters")

seacells <- umap_sc(seacells, k = 10L, knn_method = "exhaustive")
#> Running UMAP.
#> Using n_epochs = 500 (dataset <10k samples or adam_parallel optimiser)
#> Using provided kNN graph.
```

``` r

embedding_plot_sc(
  object = seacells,
  embedding = "umap",
  colour_by = "leiden_clusters",
  discrete = TRUE
)
```

![](meta_cells_files/figure-html/mc%20umap%20plot-1.png)

The plot is much sparser than the equivalent UMAP on the original cells,
which is the point: each dot is a denoised aggregate of ~15-20 cells,
and the structure that survives aggregation is the structure that’s
robust to sparsity.

### Co-expression module detection on meta cells

Meta cells compress information and reduce sparsity. That makes them
ideal for co-expression module detection methods and we can use a SCENIC
implementation within bixverse [akin to the single cell
version](https://gregorlueg.github.io/bixverse/articles/analysis_single_cell.html#scenic)
also on MetaCells.

Let’s run this quickly…

``` r

tf_dt <- data.table::fread(
  "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
  header = FALSE,
  col.names = "tf"
)

scenic_res <- scenic_grn_sc(
  object = seacells,
  tf_ids = tf_dt$tf,
  scenic_params = params_scenic(
    learner_type = "randomforest",
    gene_batch_size = 64L,
    # due to the small data set size, we should set 'min_samples_leaf' lower
    # than for a massive single cell data set
    learner_params = list(min_samples_leaf = 10L)
  ),
  .verbose = TRUE
)
#> No target genes supplied, running gene filter...
#> SCENIC gene filter: 12464 / 12464 genes pass.
#> Warning in `method(scenic_grn_sc, bixverse::MetaCells)`(object = <object>, :
#> 610 TF identifier(s) not found in the object and dropped.
#> SCENIC: 12464 target genes, 1282 TFs, 250 cells
```

``` r

scenic_res <- identify_tf_to_genes(
  scenic_res,
  n_sd = 2,
  .verbose = TRUE
)
#> Extracting TF to gene associations via per-gene threshold (mean + 2.0 * SD).

scenic_res <- tf_to_genes_correlations(
  x = scenic_res,
  object = seacells,
  cor_filter = 0.01,
  .verbose = TRUE
)
#> Calculating the pairwise correlations between the TFs and genes
#> Removing TF <> gene pairs with cors <= 0.010
#> Removing self loops (TF controlling its own expression

tf_gene_dt <- get_tf_to_gene(scenic_res)

# we will not do any filtering
tf_to_gene_ls <- split(tf_gene_dt$gene, tf_gene_dt$tf)
```

In a proper situation we would filter down the TF to gene associations
via motifs, please refer to the details in this
[vignette](https://gregorlueg.github.io/bixverse/articles/analysis_single_cell.html#scenic).
We will skip this step for now and just run AUCell to show how it works.
Which statistic you get is controlled by
[`params_sc_aucell()`](https://gregorlueg.github.io/bixverse/reference/params_sc_aucell.md);
the default is the Mann-Whitney AUC, which we stick with here.
`"recovery"` (the original AUCell recovery curve) and `"ap"` (average
precision) are the other two options and weight the top of the ranking
more heavily.

``` r

auc_res <- aucell_sc(
  object = seacells,
  gs_list = tf_to_gene_ls,
  aucell_params = params_sc_aucell(auc_type = "wilcox")
)

umap_dt <- as.data.table(
  get_embedding(seacells, "umap"),
  keep.rownames = "meta_cell_id"
)

umap_dt[, `:=`(IRF1 = auc_res[, "IRF1"], TCF4 = auc_res[, "TCF4"])]

p1 <- ggplot(umap_dt, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(fill = IRF1), size = 2.5, shape = 21) +
  theme_bw() +
  labs(fill = "IRF1 AUC") +
  scale_fill_viridis_c()

p2 <- ggplot(umap_dt, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(fill = TCF4), size = 2.5, shape = 21) +
  theme_bw() +
  labs(fill = "TCF4 AUC") +
  scale_fill_viridis_c()

p1 + p2
```

![](meta_cells_files/figure-html/run%20AUCell%20on%20meta%20cells-1.png)

## Pseudo-bulking

Pseudo-bulks are conceptually adjacent but solve a different problem.
Where meta cells are dense aggregates designed to feed into
co-expression and archetype methods, pseudo-bulks are coarse summaries
(typically one bulk per sample-by-cluster combination) used as input to
bulk-style DGE tools such as limma-voom, edgeR or DESeq2 to get around
pseudo-replication problems and p-value inflation in cell-based DGEs,
see [Zimmerman et
al.](https://www.nature.com/articles/s41467-021-21038-1)

`get_pseudobulked_sc` takes a named list of cell IDs and returns either
a dense matrix or a sparse `dgRMatrix`. With `assay = "raw"` it sums raw
counts (what bulk DGE tools expect); with `assay = "norm"` it averages
normalised counts.

``` r

cl_dt <- sc_object[[c("leiden", "cell_id")]]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpygfqOu/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
cell_list <- split(cl_dt$cell_id, cl_dt$leiden)

pb_counts <- get_pseudobulked_sc(
  object = sc_object,
  cell_list = cell_list,
  return_format = "sparse",
  assay = "raw",
  .verbose = FALSE
)

dim(pb_counts)
#> [1]    14 12464
```

The returned matrix has one row per group and one column per gene. From
here the standard bulk DGE pipeline applies; `bixverse` does not wrap
the DGE call itself.

## Clean up

``` r

unlink(tempdir_cd34, recursive = TRUE, force = TRUE)
```
