# defaults blocks --------------------------------------------------------------

# These are merged into the single cell constructors via p_merge(), so their
# rules are also the rules the merged elements are checked against.

spec_knn_graph_defaults <- param_defaults(
  name = "knn_graph_defaults",
  title = "Helper function to generate the kNN graph defaults",
  description = paste(
    "Controls how a kNN is turned into the undirected, weighted graph that",
    "the clustering runs on. `\"snn\"` is the shared nearest neighbour graph",
    "of bluster/Seurat, `\"umap\"` the fuzzy simplicial set ScanPy builds",
    "(its `obsp/connectivities`) and `\"knn\"` the plain kNN union. Each",
    "method only reads its own parameters."
  ),
  references = paste(
    "Xu & Su, Bioinformatics, 2015; McInnes, et al., arXiv, 2018."
  ),
  checker = "KnnGraph",
  label = "kNN graph params",
  hint = paste(
    "pruning must be in [0, 1]; local_connectivity must be >= 0; mix_ratio",
    "must be in [0, 1]; rescue and binary must be single booleans."
  ),
  fields = list(
    knn_sym = p_choice(
      "snn",
      c("snn", "umap", "knn"),
      doc = "How to symmetrise the kNN into a graph."
    ),
    snn_similarity = p_choice(
      "jaccard",
      c("jaccard", "rank"),
      doc = paste(
        "`\"snn\"` only. `\"jaccard\"` weighs an edge by the Jaccard index",
        "of the two (self-inclusive) neighbour sets. `\"rank\"` takes the",
        "shared neighbour with the lowest combined rank, so closer shared",
        "neighbours give heavier edges. Both are scaled to `[0, 1]`."
      )
    ),
    pruning = p_dbl(
      1 / 16,
      range = "[0, 1]",
      doc = paste(
        "`\"snn\"` only. Edges below this weight are dropped. For the",
        "Jaccard weight a threshold `p` keeps pairs that share at least a",
        "fraction `2p / (1 + p)` of their (self-inclusive) neighbourhoods,",
        "rounded up to whole cells, so the same value means the same thing",
        "at any `k`. The default asks for about 12%, i.e. 2 of 16 shared",
        "cells at `k = 15`, in line with Seurat's `1/15`. Avoid a threshold",
        "that equals one of the discrete weights exactly, or float rounding",
        "decides whether the edge survives."
      )
    ),
    rescue = p_lgl(
      TRUE,
      doc = paste(
        "`\"snn\"` only. Reconnect the small components that pruning cuts",
        "off via their strongest sNN edge, so they do not turn up as tiny",
        "clusters."
      )
    ),
    local_connectivity = p_dbl(
      1,
      range = "[0,)",
      doc = paste(
        "`\"umap\"` only. Number of nearest neighbours assumed to be fully",
        "connected, as in umap-learn."
      )
    ),
    mix_ratio = p_dbl(
      1,
      range = "[0, 1]",
      doc = paste(
        "`\"umap\"` only. `1` is the fuzzy union of the directed kNN",
        "memberships, `0` the fuzzy intersection (umap-learn's",
        "`set_op_mix_ratio`)."
      )
    ),
    binary = p_lgl(
      FALSE,
      doc = paste(
        "`\"knn\"` only. If `TRUE`, every edge weighs 1. Otherwise mutual",
        "neighbours weigh 2 and one-sided ones 1."
      )
    )
  )
)

spec_knn_defaults <- param_defaults(
  name = "knn_defaults",
  title = "Helper function to generate kNN defaults",
  description = paste(
    "This function generates various sensible default parameters",
    "for all of the different approximate nearest neighbours that",
    "are available within this package."
  ),
  checker = "Knn",
  label = "kNN params",
  hint = paste(
    "k must be >= 0; n_trees, m, ef_construction and ef_search must be",
    ">= 1; search_budget, ef_budget, n_list and n_probe must be NULL or",
    ">= 1; delta and diversify_prob must be in [0, 1]; extract_knn must be",
    "a single boolean."
  ),
  fields = list(
    k = p_int(15L, range = "[0,)", doc = "Number of neighbours."),
    knn_method = p_choice(
      "kmknn",
      c("kmknn", "hnsw", "annoy", "nndescent", "ivf", "exhaustive"),
      doc = paste(
        "Which method to use for the approximate nearest neighbour",
        "search."
      )
    ),
    ann_dist = p_choice(
      "euclidean",
      c("euclidean", "cosine"),
      doc = paste(
        "Which distance metric to use for the approximate nearest",
        "neighbour search."
      )
    ),
    n_trees = p_int(
      50L,
      range = "[1,)",
      doc = "Annoy param: number of trees to generate for Annoy."
    ),
    search_budget = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Annoy param: optional search budget per tree for Annoy. If",
        "not provided, it will default to `n_tree * k * 20L`."
      )
    ),
    delta = p_dbl(
      0.001,
      range = "[0,1]",
      doc = "NNDescent param: early termination criterium for NNDescent."
    ),
    diversify_prob = p_dbl(
      0.0,
      range = "[0,1]",
      doc = paste(
        "NNDescent param: diversification probability for the",
        "NNDescent index. This will diversify the index at the end and",
        "identify potentially better edges."
      )
    ),
    ef_budget = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "NNDescent param: optional query budget parameter. Can",
        "accelerate querying, but at the cost of Recall."
      )
    ),
    extract_knn = p_lgl(
      FALSE,
      doc = paste(
        "NNDescent param: hand back the graph the descent already built",
        "instead of beam searching it. Skips the query pass entirely, so",
        "it is much faster, at the cost of some recall. Rows the descent",
        "never filled come back padded with duplicate edges. Ignored by",
        "every other method."
      )
    ),
    m = p_int(
      16L,
      range = "[1,)",
      doc = "HNSW param: number of connections between layers for HNSW."
    ),
    ef_construction = p_int(
      200L,
      range = "[1,)",
      doc = paste(
        "HNSW param: size of dynamic candidate list during",
        "construction."
      )
    ),
    ef_search = p_int(
      100L,
      range = "[1,)",
      doc = paste(
        "HNSW param: size of candidate list (higher = better recall,",
        "slower)."
      )
    ),
    n_list = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "IVF param: number of clusters/centroids to generate. `NULL`",
        "generates `sqrt(n)` lists."
      )
    ),
    n_probe = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "IVF param: number of clusters/centroids to query. `NULL`",
        "queries `sqrt(n_lists)` clusters."
      )
    )
  )
)

spec_hvg_defaults <- param_defaults(
  name = "hvg_defaults",
  title = "Helper function to generate HVG defaults",
  checker = NULL,
  label = "HVG params",
  fields = list(
    min_gene_var_pctl = p_dbl(
      0.7,
      range = "[0,1]",
      doc = "Which percentile of the highly variable genes to include."
    ),
    hvg_method = p_choice(
      "vst",
      c("vst", "meanvarbin", "dispersion", "scran"),
      doc = paste(
        "Which method to use to identify HVG. `\"scran\"` runs with the",
        "default trend parameters, see [params_hvg_scran_defaults()]."
      )
    ),
    loess_span = p_dbl(
      0.3,
      range = "(0,)",
      doc = "In case of `\"vst\"` the span of the loess function."
    ),
    clip_max = p_dbl(
      NULL,
      range = "(0,)",
      null_ok = TRUE,
      doc = "The maximum clipping value (optional)."
    ),
    n_bins = p_int(
      20L,
      range = "[1,)",
      doc = paste(
        "The number of bins to use for the `\"meanvarbin\"` and",
        "`\"dispersion\"` HVG detection."
      )
    ),
    binning_strategy = p_choice(
      "equal_width",
      c("equal_width", "equal_frequency"),
      doc = paste(
        "Which binning strategy to use for `\"meanvarbin\"` and",
        "`\"dispersion\"`."
      )
    )
  )
)

spec_hvg_scran_defaults <- param_defaults(
  name = "hvg_scran_defaults",
  title = "Helper function to generate the scran HVG trend defaults",
  description = paste(
    "Trend parameters for `method = \"scran\"` in [params_sc_hvg()]. They",
    "mirror scrapper's `fitVarianceTrend()` defaults before scrapper 1.8,",
    "which switched `use.min.width` to `TRUE`. The lowess span comes from",
    "`loess_span` in [params_sc_hvg()]."
  ),
  checker = NULL,
  label = "scran HVG trend params",
  fields = list(
    mean_filter = p_lgl(
      TRUE,
      doc = "Shall genes below `min_mean` be left out of the trend fit."
    ),
    min_mean = p_dbl(
      0.1,
      range = "[0,)",
      doc = "Minimum mean log-expression for a gene to enter the trend fit."
    ),
    transform = p_lgl(
      TRUE,
      doc = "Shall the variances be fourth-root transformed before the fit."
    ),
    use_min_width = p_lgl(
      FALSE,
      doc = paste(
        "Shall the lowess window be defined by `min_width` and",
        "`min_window_count` instead of the span."
      )
    ),
    min_width = p_dbl(
      1,
      range = "(0,)",
      doc = "Minimum window width, only used with `use_min_width = TRUE`."
    ),
    min_window_count = p_int(
      200L,
      range = "[1,)",
      doc = paste(
        "Minimum number of genes per window, only used with",
        "`use_min_width = TRUE`."
      )
    )
  )
)

spec_norm_doublets_defaults <- param_defaults(
  name = "norm_doublets_defaults",
  title = paste(
    "Helper function to generate normalisation defaults for doublet",
    "detection."
  ),
  checker = NULL,
  label = "doublet normalisation params",
  fields = list(
    log_transform = p_lgl(TRUE, doc = "Shall the counts be log-normalised."),
    mean_center = p_lgl(FALSE, doc = "Shall mean centring be applied."),
    normalise_variance = p_lgl(
      FALSE,
      doc = "Shall the variance be normalised."
    ),
    target_size = p_dbl(1e6, range = "[0,)", doc = "Target library size.")
  )
)

spec_pca_defaults <- param_defaults(
  name = "pca_defaults",
  title = "Helper function to generate default parameters for PCA",
  checker = NULL,
  label = "PCA params",
  fields = list(
    no_pcs = p_int(30L, range = "[1,)", doc = "Number of PCs to consider."),
    svd_solver = p_choice(
      "randomised",
      c("randomised", "covariance", "exact"),
      doc = paste(
        "Which solver to use. `\"randomised\"` (default) is a randomised",
        "SVD, approximate in the trailing components. `\"covariance\"`",
        "builds the gene x gene cross-product and eigendecomposes it.",
        "`\"exact\"` is Lanczos on the sparse path and a full SVD on the",
        "dense one. See [params_sc_pca()] for the trade-offs."
      )
    )
  )
)

spec_fast_cluster_default <- param_defaults(
  name = "fast_cluster_default",
  title = paste(
    "Helper function to generate default parameters for the fast",
    "clustering for the doublet detection methods"
  ),
  checker = "FastClusterDefault",
  label = "FastCluster params",
  hint = paste(
    "batch_size and k_means_iter must be integers; n_centroids must be an",
    "integer or NULL."
  ),
  fields = list(
    km_type = p_choice(
      "minibatch",
      c("minibatch", "standard"),
      doc = "The type of k-means clustering."
    ),
    n_centroids = p_int(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "The number of centroids to use. `NULL` uses",
        "`sqrt(N_cells) * 4` centroids."
      )
    ),
    k_means_iter = p_int(
      100L,
      doc = "Number of maximum k-means iterations."
    ),
    batch_size = p_int(
      4098L,
      doc = paste(
        "Max batch size, capped at `N_cells / 2` depending on the data",
        "set."
      )
    )
  )
)

spec_kmeans_defaults <- param_defaults(
  name = "kmeans_defaults",
  title = "K-mean parameter defaults.",
  description = paste(
    "Helper function to generate defaults for the k-mean clustering were",
    "more control is needed."
  ),
  checker = "KMeans",
  label = "k-means params",
  hint = paste(
    "k_means_iter must be an integer; gemm and hamerly must be booleans",
    "or NULL."
  ),
  fields = list(
    k_means_iter = p_int(
      30L,
      doc = "The number of iterations to use for the clustering."
    ),
    k_means_init = p_choice(
      "parallel",
      c("parallel", "random"),
      doc = "The initialisation."
    ),
    gemm = p_lgl(
      FALSE,
      null_ok = TRUE,
      doc = paste(
        "Controls which CPU implementation is used by the method. GEMM",
        "is faster with large dimensionality."
      )
    ),
    hamerly = p_lgl(
      TRUE,
      null_ok = TRUE,
      doc = paste(
        "Shall a faster exact method be used leveraging the triangle",
        "inequality. Faster on large data sets with large numbers of",
        "centroids."
      )
    )
  )
)
