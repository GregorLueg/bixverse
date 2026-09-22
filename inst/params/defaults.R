# defaults blocks --------------------------------------------------------------

# These are merged into the single cell constructors via p_merge(), so their
# rules are also the rules the merged elements are checked against.

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
      c("vst", "mvb", "dispersion"),
      doc = "Which method to use to identify HVG."
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
      doc = "The number of bins to use for the `\"mvb\"` HVG detection."
    ),
    binning_strategy = p_choice(
      "equal_width",
      c("equal_width", "equal_frequency"),
      doc = "Which binning strategy to use for `\"mvb\"`."
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
    random_svd = p_lgl(TRUE, doc = "Shall randomised SVD be used."),
    sparse = p_lgl(
      FALSE,
      doc = paste(
        "Shall sparse solvers be used that do not do scaling. If set to",
        "yes, in the case of `random_svd = FALSE`, Lanczos iterations",
        "are used to solve the sparse SVD. With `random_svd = TRUE`, the",
        "sparse initial matrix is multiplied with the random matrix,",
        "yielding a much smaller dense matrix that does not increase the",
        "memory pressure massively."
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
    "batch_size and kmeans_iters must be integers; n_centroids must be an",
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
    kmeans_iters = p_int(
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
