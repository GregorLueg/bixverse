# KNN searches for single cell

## Intro

kNN searches are the bedrock of a lot of methods in the single cell
space. This vignette is basically a primer to understand the large
number of (approximate) nearest neighbour search methods available in
`bixverse`, when to use what, when to be careful with them and hopefully
gives you a general idea on what they can (and cannot do) well.

``` r

library(bixverse)
library(ggplot2)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
```

## Loading the data and basic workflow

We start by downloading the PBMC3k data set bundled with the package and
loading it via Cell Ranger-style MTX I/O, run a quick and dirty QC, HVG
and PCA to then move onto the kNN searches…

``` r

# load in the data
pbmc3k_path <- bixverse:::download_pbmc3k()

tempdir_pbmc <- tempdir()

sc_object <- SingleCells(dir_data = tempdir_pbmc)

mtx_io_params <- get_cell_ranger_params(pbmc3k_path)

sc_object <- load_mtx(
  object = sc_object,
  sc_mtx_io_param = mtx_io_params,
  streaming = FALSE,
  .verbose = TRUE
)
#> Loading observations data from flat file into the DuckDB.
#> Loading variable data from flat file into the DuckDB.

# ensembl to gene symbols
setnames_sc(
  object = sc_object,
  table = "var",
  old = "column1",
  new = "gene_symbol"
)

var <- get_sc_var(sc_object)

ensembl_to_symbol <- setNames(var$gene_symbol, var$gene_id)
symbol_to_ensembl <- setNames(var$gene_id, var$gene_symbol)

gs_of_interest <- list(
  MT = var[grepl("^MT-", gene_symbol), gene_id],
  Ribo = var[grepl("^RPS|^RPL", gene_symbol), gene_id]
)

# remove high mitochondrial genes and generally low quality cells
sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  streaming = FALSE,
  .verbose = TRUE
)

qc_df <- sc_object[[c("cell_id", "lib_size", "nnz", "MT")]]

metrics <- list(
  log10_lib_size = log10(qc_df$lib_size),
  log10_nnz = log10(qc_df$nnz),
  MT = qc_df$MT
)

directions <- c(
  log10_lib_size = "twosided",
  log10_nnz = "twosided",
  MT = "above"
)

qc <- run_cell_qc(metrics, directions, threshold = 3)

sc_object[["outlier"]] <- qc$combined

cells_to_keep <- qc_df[!qc$combined, cell_id]

sc_object <- set_cells_to_keep(sc_object, cells_to_keep)

# HVG
sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = 2000L,
  .verbose = TRUE
)

# PCA
sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L,
  sparse_svd = TRUE
)
#> Using sparse SVD solving on scaled data on 2000 HVG.
```

## kNN searches

Okay, with all of this out of the way, let’s look at the different
methods available in bixverse to generate kNN graphs for single cell.
You will realise that there is a LARGE number of methods that use
something like `knn = list()` with the function
[`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md).
There is a lot of different parameters in there. Let’s first check out
what we have there…

``` r

str(params_knn_defaults())
#> List of 13
#>  $ k              : int 15
#>  $ knn_method     : chr "kmknn"
#>  $ ann_dist       : chr "euclidean"
#>  $ n_trees        : int 50
#>  $ search_budget  : NULL
#>  $ delta          : num 0.001
#>  $ diversify_prob : num 0
#>  $ ef_budget      : NULL
#>  $ m              : int 16
#>  $ ef_construction: int 200
#>  $ ef_search      : int 100
#>  $ n_list         : NULL
#>  $ n_probe        : NULL
```

### The methods

There is a lot of different parameters in here and it is worthwhile
actually understanding them and where they are being called used. Let’s
first give an overview of the different kNN searches implemented/exposed
in this package (all happens via
[ann-search-rs](https://crates.io/crates/ann-search-rs)).

| Method | Features | Class | Use case | Ideal k range |
|----|----|----|----|----|
| Exhaustive | What the name says | Exact | Works well on smaller data sets where the CPU can just crunch through the numbers | Does not matter |
| KmKnn | An exact method that takes a shortcut via k-means cluster | (Fast) exact | On well structured (low dimensional) data, it gives you an exact kNN search that is substantially faster than Exhaustive on large data sets | Does not matter |
| Annoy | The classic approximate nearest neighbour search you know from Seurat | Tree-based | A memory-based version of Annoy that avoids the round trip to the disk-based index of the original. Good default for data sets with up to half a million cells with fast index building (but slower querying) | Ideally ≤ 100 as otherwise the query phase can become quite long… |
| HNSW | A powerful graph-based index for large data sets and vector search index powering a large number of data bases. | Graph-based | Need to tackle large data sets … ? HNSW does that for you. An important note is that due to the benign race condition during graph construction this index is not deterministic. On large data sets this makes basically no difference, but something to keep in mind for smaller ones! | ≤50 or you need to bump the `ef_search` substantially |
| NNDescent | Another graph-based index that is ideal for tackling large data sets | Graph-based | Another graph-based index that is very good at dealing with large data sets. Uses tricks from PyNNDescent and is designed to be very memory-efficient | ≤30. The initial kNN is designed with a node degree of 30, so you might see degradation in the performance when you ask for higher k. If you try this, set `ef_budget` high! |
| IVF | A cluster-based index that partitions the data via k-means into Voronoi cells | Cluster-based | If you read a bit about vector searches, this is a very common algorithm which works well with quantisations in particular (and GPU-acceleration…). For the usual ephemeral kNN graph generation in single cell maybe not the right call, but if you need A LOT of neighbours on large data sets, this is a very good choice! | Does not matter (outside of absurd high numbers) |

There are additional GPU-accelerated kNN searches in the sister package,
see
[here](https://gregorlueg.github.io/bixverse.gpu/articles/gpu_single_cell.html).

### The parameters

Now let’s see which parameter belongs to which method and does what.

| Parameter | Default | Meaning |
|----|----|----|
| k | `15L` | Kinda obvious… Number of neighbours |
| knn_method | `"kmknn"` | The approximate nearest neighbour method. One of `c("kmknn", "exhaustive", "annoy", "ivf", "hnsw", "nndescent")` |
| ann_dist | `"euclidean"` | The distance metric to use. Two options are available… Euclidean distance or Cosine distance. Very important, the underlying library returns the **SQUARED** Euclidean distance for speed purposes. If you use the Euclidean distance for something, keep this in mind! |
| n_trees | `50L` | **Annoy:** The number of trees to use for Annoy. Generally speaking the generation is very fast in terms of indexing and in different (emperical) experiments after 50 trees one gets diminishing returns |
| search_budget | `NULL` | **Annoy:** The search budget during querying. If `NULL` defaults to `k * n_trees * 20`. Higher values allow for better Recall. |
| delta | `0.001` | **NNDescent:** The termination criterium for the NNDescent iterations. When less than 0.001 (0.1%) of the edges are updated, the iterations stop and the algorithm is considered converged. |
| diversify_prob | `0.0` | **NNDescent:** If you want to do additional diversification after the NNDescent algorithm converged. This is based on the original paper, but in reality usually degrades graph quality while not providing too many benefits. |
| ef_budget | `NULL` | **NNDescent:** The beam width for the query phase for NNDescent. If `NULL`, it will use the library default, i.e., `(k * 2).clamp(50, 200)).max(k)`. Generally not really needed to touch, but limits the number of k neighbours you can confidently return. |
| m | `16L` | **HNSW:** The node degree during the index creation. Creates 16 edges per given node in the graph-based index. |
| ef_construction | `200L` | **HNSW:** The budget to generate good connections during construction of the index. Higher values yields a better initial graph-index. |
| ef_search | `100L` | **HNSW:** The beam width during querying of the index. If a high quality graph was generated, usually does not need to be too high. If you set `k` high, you should adopt this parameter. |
| n_list | `NULL` | **IVF and KmKnn:** This parameter is shared between the IVF and KmKnn index. Defines the number of Voronoi cells/clusters to generate. If not provided defaults to `sqrt(N)`. Be careful on smaller data sets and possible set it to more clusters here. |
| n_probe | `NULL` | **IVF:** The number of clusters to probe per query. Defaults to `sqrt(n_list)`. On smaller data sets you need to set this higher. |

### Exploration of the kNN search methods

Okay, still do not get it. When do I use what and how … ? The standard
approach would be to do something like:

``` r

sc_object <- find_neighbours_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "exhaustive")
  )
)
#> 
#> Generating sNN graph (full: FALSE).
#> Transforming sNN data to igraph.
```

You can however also generate the kNN indices like this:

``` r

exhaustive_knn <- generate_knn_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "exhaustive")
  )
)

exhaustive_knn
#> SingleCellNearestNeighbour: 2163 cells, k = 15
#>   Distance metric: euclidean
#>   Index range: [0, 2162]
#>   Distance range: [20.8050, 24447.3223]
```

This is the class that is actually sored in the object after running
[`find_neighbours_sc()`](https://gregorlueg.github.io/bixverse/reference/find_neighbours_sc.md).

``` r

exhaustive_knn_obj <- get_knn_obj(sc_object)

tinytest::expect_equal(
  current = exhaustive_knn,
  target = exhaustive_knn_obj,
  info = "same class"
)
#> ----- PASSED      : <-->
#>  call| tinytest::expect_equal(current = exhaustive_knn, target = exhaustive_knn_obj, 
#>  call| info = "same class")
#>  info| same class
```

Potentially more interesting is what you can do with the approximate
nearest neighbour searches.

``` r

annoy_knn <- generate_knn_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "annoy")
  )
)

annoy_knn
#> SingleCellNearestNeighbour: 2163 cells, k = 15
#>   Distance metric: euclidean
#>   Index range: [0, 2162]
#>   Distance range: [20.8050, 24447.3223]
```

The print message showed this line here:

    Recall of approximate nearest neighbours search in random subset: 1.00

The approximate versions have the option to test the index against a
subset of the data (1000 randomly selected cells) and compare against
the exhaustive search for this 1000 cells.

``` r

annoy_vs_exhaustive <- calc_knn_metrics(exhaustive_knn, annoy_knn)

cat(
  "Annoy vs exhaustive.\n",
  sprintf(
    " Recall@15: %.3f.\n",
    annoy_vs_exhaustive$final_recall
  ),
  sprintf(
    " Distance ratio: %.3f.\n",
    annoy_vs_exhaustive$final_ratio
  )
)
#> Annoy vs exhaustive.
#>   Recall@15: 1.000.
#>   Distance ratio: 1.000.
```

This allows us to actually explore some parameters more in depth:

``` r

annoy_knn_less_trees <- generate_knn_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "annoy", n_trees = 1L)
  )
)

annoy_knn_less_trees
#> SingleCellNearestNeighbour: 2163 cells, k = 15
#>   Distance metric: euclidean
#>   Index range: [0, 2162]
#>   Distance range: [20.8050, 24860.6562]
```

With one tree, we can appreciate that the Recall is dropping:

    Recall of approximate nearest neighbours search in random subset: 0.75

Let’s get some more detailed info here:

``` r

annoy_lt_vs_exhaustive <- calc_knn_metrics(exhaustive_knn, annoy_knn_less_trees)

cat(
  "Annoy vs exhaustive.\n",
  sprintf(
    " Recall@15: %.3f.\n",
    annoy_lt_vs_exhaustive$final_recall
  ),
  sprintf(
    " Distance ratio: %.3f.\n",
    annoy_lt_vs_exhaustive$final_ratio
  )
)
#> Annoy vs exhaustive.
#>   Recall@15: 0.771.
#>   Distance ratio: 1.033.
```

We can appreciate that with one tree we do not get a Recall of 1 and the
distance ratio is 4.2% worse than the exhaustive one. Let’s test another
index

``` r

nndescent_knn <- generate_knn_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "nndescent")
  )
)

nndescent_knn
#> SingleCellNearestNeighbour: 2163 cells, k = 15
#>   Distance metric: euclidean
#>   Index range: [0, 2162]
#>   Distance range: [20.8050, 24452.8691]
```

And benchmark it:

``` r

nndescent_vs_exhaustive <- calc_knn_metrics(exhaustive_knn, nndescent_knn)

cat(
  "Annoy vs exhaustive.\n",
  sprintf(
    " Recall@15: %.3f.\n",
    nndescent_vs_exhaustive$final_recall
  ),
  sprintf(
    " Distance ratio: %.3f.\n",
    nndescent_vs_exhaustive$final_ratio
  )
)
#> Annoy vs exhaustive.
#>   Recall@15: 1.000.
#>   Distance ratio: 1.000.
```

As we can see, this approximate index also reaches perfect recall and
distance ratio.

## Conclusions

We could iterate through the others, but it does not have much of a
point. The goal of this vignette was to check the type of
hyperparameters and type of indices provided and give an idea when to
use which and how to work with them and validate parameters if neededed.

## Clean up

``` r

unlink(tempdir_pbmc, recursive = TRUE, force = TRUE)
```
