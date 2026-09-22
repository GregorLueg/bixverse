spec_scrublet <- param_spec(
  name = "scrublet",
  title = "Wrapper function for Scrublet doublet detection parameters",
  description = paste(
    "Constructor for the various Scrublet parameters. In this",
    "case, the default for the kNN graph generation was set to",
    "`\"hnsw\"` as this algorithm showed the best performance in",
    "different empirical benchmarks."
  ),
  return_order = c(
    "normalisation",
    "hvg",
    "pca",
    "knn",
    "sim_doublet_ratio",
    "expected_doublet_rate",
    "stdev_doublet_rate",
    "n_bins_histogram",
    "manual_threshold"
  ),
  checker = "ScScrublet",
  label = "Scrublet params",
  hint = paste(
    "no_pcs must be >= 1; n_bins_histogram must be >= 10; n_bins",
    "must be >= 1; min_gene_var_pctl, expected_doublet_rate and",
    "stdev_doublet_rate must be in [0, 1]; loess_span and",
    "sim_doublet_ratio must be > 0; target_size must be >= 0;",
    "log_transform, mean_center, normalise_variance and",
    "random_svd must be booleans; clip_max and manual_threshold",
    "must be NULL or positive numerics."
  ),
  fields = list(
    sim_doublet_ratio = p_dbl(
      1.5,
      range = "(0,)",
      doc = paste(
        "Number of doublets to simulate relative to the number of",
        "observed cells. For example, 2.0 simulates twice as many",
        "doublets as there are cells."
      )
    ),
    expected_doublet_rate = p_dbl(
      0.1,
      range = "[0,1]",
      doc = paste(
        "Expected doublet rate for the experiment, typically",
        "0.05-0.10 depending on cell loading. Must be between 0 and",
        "1."
      )
    ),
    stdev_doublet_rate = p_dbl(
      0.02,
      range = "[0,1]",
      doc = "Uncertainty in the expected doublet rate."
    ),
    n_bins_histogram = p_int(
      100L,
      range = "[10,)",
      doc = paste(
        "Number of bins for histogram-based automatic threshold",
        "detection. Typically 50-100."
      )
    ),
    manual_threshold = p_dbl(
      NULL,
      range = "[0,)",
      null_ok = TRUE,
      doc = paste(
        "Manual doublet score threshold. If `NULL` (default),",
        "threshold is automatically detected from simulated doublet",
        "score distribution."
      )
    ),
    normalisation = p_merge(
      "norm_doublets_defaults",
      doc = paste(
        "Optional overrides for normalisation parameters. See",
        "[bixverse::params_norm_doublets_defaults()] for available",
        "parameters: `log_transform`, `mean_center`,",
        "`normalise_variance`, `target_size`."
      )
    ),
    hvg = p_merge(
      "hvg_defaults",
      doc = paste(
        "Optional overrides for highly variable gene selection",
        "parameters. See [bixverse::params_hvg_defaults()] for",
        "available parameters: `min_gene_var_pctl`, `hvg_method`,",
        "`loess_span`, `clip_max`."
      )
    ),
    pca = p_merge(
      "pca_defaults",
      doc = paste(
        "Optional overrides for PCA parameters. See",
        "[bixverse::params_pca_defaults()] for available parameters:",
        "`no_pcs`, `random_svd`, `sparse` and `skip_first_pc`."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      default = list(k = 0L),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`.",
        "Note: this function defaults to `k = 0L` (automatic",
        "neighbour detection)."
      )
    )
  )
)

spec_boost <- param_spec(
  name = "boost",
  title = "Wrapper function for Boost parameters",
  return_order = c(
    "normalisation",
    "hvg",
    "pca",
    "knn",
    "fast_cluster_params",
    "boost_rate",
    "replace",
    "resolution",
    "fast_cluster",
    "n_iters",
    "p_thresh",
    "voter_thresh"
  ),
  checker = "ScBoost",
  label = "Boost params",
  hint = paste(
    "no_pcs, n_bins and n_iters must be >= 1; min_gene_var_pctl,",
    "boost_rate and voter_thresh must be in [0, 1]; loess_span,",
    "resolution and p_thresh must be > 0; target_size must be >",
    "0; log_transform, mean_center, normalise_variance, replace,",
    "random_svd and fast_cluster must be booleans; clip_max must",
    "be NULL or a positive numeric."
  ),
  fields = list(
    boost_rate = p_dbl(
      0.25,
      range = "[0,1]",
      doc = "Boosting rate for the algorithm. Must be between 0 and 1."
    ),
    replace = p_lgl(FALSE, doc = "Whether to use replacement during boosting."),
    resolution = p_dbl(
      1,
      range = "(0,)",
      doc = paste(
        "Resolution parameter for graph-based clustering. Higher",
        "values lead to more clusters."
      )
    ),
    n_iters = p_int(
      20L,
      range = "[1,)",
      doc = "Number of iterations to run the algorithm."
    ),
    p_thresh = p_dbl(
      1e-07,
      range = "(0,)",
      doc = "P-value threshold for significance testing."
    ),
    voter_thresh = p_dbl(
      0.9,
      range = "[0,1]",
      doc = paste(
        "Voter threshold across iterations. Proportion of iterations",
        "a cell must be assigned to a cluster to be considered a",
        "member. Must be between 0 and 1."
      )
    ),
    fast_cluster = p_lgl(
      FALSE,
      doc = paste(
        "Shall fast Louvain clustering be applied, i.e., k-means",
        "clustering and use the centroids for kNN graph generation",
        "and Louvain clustering with then backpropagating the",
        "membership based on centroid proximity."
      )
    ),
    normalisation = p_merge(
      "norm_doublets_defaults",
      doc = paste(
        "Optional overrides for normalisation parameters. See",
        "[bixverse::params_norm_doublets_defaults()] for available",
        "parameters: `log_transform`, `mean_center`,",
        "`normalise_variance`, `target_size`."
      )
    ),
    hvg = p_merge(
      "hvg_defaults",
      doc = paste(
        "Optional overrides for highly variable gene selection",
        "parameters. See [bixverse::params_hvg_defaults()] for",
        "available parameters: `min_gene_var_pctl`, `hvg_method`,",
        "`loess_span`, `clip_max`."
      )
    ),
    pca = p_merge(
      "pca_defaults",
      doc = paste(
        "Optional overrides for PCA parameters. See",
        "[bixverse::params_pca_defaults()] for available parameters:",
        "`no_pcs`, `random_svd`."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      default = list(k = 0L),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`.",
        "Note: this function defaults to `k = 0L` (automatic",
        "neighbour detection)."
      )
    ),
    fast_cluster_params = p_merge(
      "fast_cluster_default",
      doc = paste(
        "Optional overrides for the fast clustering parameters. Only",
        "relevant if `fast_cluster = TRUE`. See",
        "[params_fast_cluster_default()] for available parameters:",
        "`km_type`, `n_centroids`, `kmeans_iters` and `batch_size`."
      )
    )
  )
)

spec_scdblfinder <- param_spec(
  name = "scdblfinder",
  title = paste(
    "Wrapper function for scDblFinder doublet detection",
    "parameters"
  ),
  description = paste(
    "Constructor for the scDblFinder parameters. This method",
    "combines cluster-aware doublet simulation with a",
    "gradient-boosted classifier trained on engineered features."
  ),
  return_order = c(
    "normalisation",
    "pca",
    "knn",
    "fast_cluster_params",
    "n_genes",
    "doublet_ratio",
    "heterotypic_bias",
    "cluster_resolution",
    "cluster_iters",
    "fast_cluster",
    "n_iterations",
    "gbm_n_trees",
    "max_depth",
    "learning_rate",
    "min_samples_leaf",
    "subsample_rate",
    "cv_folds",
    "cv_early_stop",
    "se_fraction",
    "include_pcs",
    "expected_doublet_rate",
    "manual_threshold",
    "cxds_genes"
  ),
  checker = "ScDblFinder",
  label = "scDblFinder params",
  fields = list(
    n_genes = p_int(
      1352L,
      range = "[1,)",
      doc = "Number of top-expressed genes to use as features."
    ),
    doublet_ratio = p_dbl(
      1,
      range = "(0,)",
      doc = "Ratio of simulated doublets to observed cells."
    ),
    heterotypic_bias = p_dbl(
      1,
      range = "[0,1]",
      doc = paste(
        "Fraction of simulated pairs forced to come from different",
        "clusters (0-1)."
      )
    ),
    cluster_resolution = p_dbl(
      1,
      range = "(0,)",
      doc = "Resolution for the initial Louvain clustering."
    ),
    cluster_iters = p_int(
      10L,
      range = "[1,)",
      doc = "Number of Louvain iterations per clustering step."
    ),
    fast_cluster = p_lgl(
      FALSE,
      doc = paste(
        "Shall fast Louvain clustering be applied, i.e., k-means",
        "clustering and use the centroids for kNN graph generation",
        "and Louvain clustering with then backpropagating the",
        "membership based on centroid proximity."
      )
    ),
    n_iterations = p_int(
      3L,
      range = "[1,)",
      doc = "Number of refinement iterations. Typically 2-3."
    ),
    gbm_n_trees = p_int(
      200L,
      range = "[1,)",
      doc = "Maximum number of boosting rounds for the GBM classifier."
    ),
    max_depth = p_int(
      4L,
      range = "[1,)",
      doc = "Maximum tree depth. Shallow trees (3-5) work best."
    ),
    learning_rate = p_dbl(
      0.3,
      range = "(0,)",
      doc = "Shrinkage applied to each tree."
    ),
    min_samples_leaf = p_int(
      20L,
      range = "[1,)",
      doc = "Minimum training samples per leaf."
    ),
    subsample_rate = p_dbl(
      0.75,
      range = "(0,1]",
      doc = "Fraction of samples used per tree."
    ),
    cv_folds = p_int(
      5L,
      range = "[2,)",
      doc = paste(
        "Number of cross-validation folds for boosting round",
        "selection."
      )
    ),
    cv_early_stop = p_int(
      2L,
      range = "[1,)",
      doc = "Early stopping patience per CV fold."
    ),
    se_fraction = p_dbl(
      1,
      range = "[0,)",
      doc = paste(
        "Multiplier on the standard error for the SE rule used in",
        "round selection."
      )
    ),
    include_pcs = p_free(
      19L,
      doc = paste(
        "Number of leading principal components to include as",
        "classifier features."
      )
    ),
    expected_doublet_rate = p_dbl(
      NULL,
      range = "(0,1]",
      null_ok = TRUE,
      doc = paste(
        "Expected doublet rate as a percentage. If not provided, will",
        "be calculated internally."
      )
    ),
    cxds_genes = p_int(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "Number of CXDS genes to consider. If not provided, defaults",
        "to `500L`."
      )
    ),
    manual_threshold = p_dbl(
      NULL,
      range = "[0,)",
      null_ok = TRUE,
      doc = paste(
        "Manual score threshold. If `NULL` (default), expected-rate",
        "thresholding is used."
      )
    ),
    normalisation = p_merge(
      "norm_doublets_defaults",
      default = list(mean_center = TRUE),
      doc = paste(
        "Optional overrides for normalisation parameters. See",
        "[bixverse::params_norm_doublets_defaults()]."
      )
    ),
    pca = p_merge(
      "pca_defaults",
      doc = paste(
        "Optional overrides for PCA parameters. See",
        "[bixverse::params_pca_defaults()]."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      default = list(k = 0L),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()]. NNDescent works better",
        "for the larger k-values often used here."
      )
    ),
    fast_cluster_params = p_merge(
      "fast_cluster_default",
      doc = paste(
        "Optional overrides for the fast clustering parameters. Only",
        "relevant if `fast_cluster = TRUE`. See",
        "[params_fast_cluster_default()] for available parameters:",
        "`km_type`, `n_centroids`, `kmeans_iters` and `batch_size`."
      )
    )
  )
)

spec_sc_neighbours <- param_spec(
  name = "sc_neighbours",
  title = paste(
    "Wrapper function for parameters for neighbour identification",
    "in single cell"
  ),
  checker = "ScNeighbours",
  label = "neighbour params",
  hint = "full_snn must be a boolean; pruning must be in [0, 1].",
  fields = list(
    full_snn = p_lgl(
      TRUE,
      doc = paste(
        "Shall the full shared nearest neighbour graph be generated",
        "that generates edges between all cells instead of between",
        "only neighbours."
      )
    ),
    pruning = p_dbl(
      1 / 12,
      range = "[0, 1]",
      doc = paste(
        "Weights below this threshold will be set to 0 in the",
        "generation of the sNN graph. Seurat uses for example `1/15`",
        "with `k = 20`. As the default k is set to 15, we set it to",
        "`1/12`. Track this against `k` rather than leaving it: the",
        "threshold is a share of the neighbourhood, so the same value",
        "prunes far harder at a larger `k`. Over-pruning fails",
        "quietly, in that you still get a clustering, but cells left",
        "with too few shared neighbours drop out as singleton",
        "communities, which then show up downstream as one-cell",
        "clusters with inflated [bixverse::run_paga_sc()]",
        "connectivities."
      )
    ),
    snn_similarity = p_choice(
      "jaccard",
      c("rank", "jaccard"),
      doc = paste(
        "The Jaccard similarity calculates the Jaccard between the",
        "neighbours, whereas the rank method calculates edge weights",
        "based on the ranking of shared neighbours. For the rank",
        "method, the weight is determined by finding the shared",
        "neighbour with the lowest combined rank across both cells,",
        "where lower-ranked (closer) shared neighbours result in",
        "higher edge weights Both methods produce weights normalised",
        "to the range `[0, 1]`."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`."
      )
    )
  )
)

spec_sc_fast_cluster <- param_spec(
  name = "sc_fast_cluster",
  title = "Fast single cell clustering parameters",
  return_order = c(
    "kmeans_iters",
    "batch_size",
    "drift_threshold",
    "lr_alpha",
    "louvain_iters",
    "full_snn",
    "pruning",
    "snn_similarity",
    "knn"
  ),
  checker = "ScFastCluster",
  label = "fast clustering params",
  hint = paste(
    "kmeans_iters, batch_size and louvain_iters must be integers",
    ">= 1; drift_threshold and lr_alpha must be single numerics;",
    "full_snn must be a boolean; pruning must be NULL or a single",
    "numeric."
  ),
  fields = list(
    kmeans_iters = p_int(
      100L,
      range = "[1,)",
      doc = "Number of iterations for k-means clustering."
    ),
    batch_size = p_int(
      4096L,
      range = "[1,)",
      doc = "Batch size for mini batch k-means clustering."
    ),
    drift_threshold = p_dbl(
      1e-04,
      doc = paste(
        "The drift for the mini batch k-means clustering. If the",
        "centroid drift is below this, the mini batch k-means",
        "terminates."
      )
    ),
    lr_alpha = p_dbl(
      1,
      doc = "Learning rate alpha parameter for mini batch k-means."
    ),
    full_snn = p_lgl(
      FALSE,
      doc = paste(
        "Shall the full shared nearest neighbour graph be generated",
        "that generates edges between all cells instead of between",
        "only neighbours."
      )
    ),
    pruning = p_dbl(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "Weights below this threshold will be set to 0 in the",
        "generation of the sNN graph. If not provided, defaults to `1",
        "/ ceil(k * 0.8)`."
      )
    ),
    snn_similarity = p_choice(
      "jaccard",
      c("jaccard", "rank"),
      doc = paste(
        "The Jaccard similarity calculates the Jaccard between the",
        "neighbours, whereas the rank method calculates edge weights",
        "based on the ranking of shared neighbours. For the rank",
        "method, the weight is determined by finding the shared",
        "neighbour with the lowest combined rank across both cells,",
        "where lower-ranked (closer) shared neighbours result in",
        "higher edge weights Both methods produce weights normalised",
        "to the range `[0, 1]`."
      )
    ),
    louvain_iters = p_int(
      10L,
      range = "[1,)",
      doc = "Number of iterations for Louvain clustering."
    ),
    knn = p_merge(
      "knn_defaults",
      default = list(k = 5L),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`. Sets",
        "the default `k = 5L`."
      )
    )
  )
)

spec_sc_vision <- param_spec(
  name = "sc_vision",
  title = paste(
    "Wrapper function for parameters for VISION with",
    "auto-correlation"
  ),
  checker = "ScVision",
  label = "VISION params",
  hint = "n_perm and n_cluster must be integers >= 1.",
  fields = list(
    n_perm = p_int(
      500L,
      range = "[1,)",
      doc = "Number of random gene sets to generate per cluster."
    ),
    n_cluster = p_int(
      5L,
      range = "[1,)",
      doc = paste(
        "Number of clusters for the random gene set clustering",
        "generation."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      overrides = list(k = 15L, nn_max_iter = 15L),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`."
      )
    )
  )
)

spec_sc_aucell <- param_spec(
  name = "sc_aucell",
  title = "Wrapper function for parameters for AUCell",
  description = paste(
    "The three statistics consume the same within-cell ranking,",
    "but weight it very differently. `\"wilcox\"` is a pure",
    "function of the gene set's rank sum, so a gene at rank 2 and",
    "one at rank 200 count for almost the same thing.",
    "`\"recovery\"` and `\"ap\"` are top-heavy. For SCENIC use",
    "`\"recovery\"`, the default. `\"wilcox\"` is a bixverse",
    "addition and its flatter score does not separate into on/off",
    "populations, so it binarises badly. Ranking ties are",
    "averaged (midranks) rather than broken at random the way",
    "AUCell does it, which is why there is no need for SCENIC's",
    "trick of setting the cutoff to the 1st percentile of genes",
    "detected per cell. Undetected genes all collapse onto one",
    "rank well outside any sensible `max_rank`. Note the recovery",
    "AUC is normalised by `max_rank * length(gene_set)`, matching",
    "pySCENIC and AUCell's `normAUC = FALSE`. Modern AUCell",
    "divides by the attainable maximum instead, so absolute",
    "values differ by a per-gene-set constant. Cell ordering",
    "within a regulon is unaffected."
  ),
  references = "Aibar, et al., Nat Methods, 2017",
  checker = "ScAucell",
  label = "AUCell params",
  hint = paste(
    "max_rank must be NULL or a single numeric >= 1; standardise must",
    "be a single logical."
  ),
  # Rust parses max_rank with as_real(), so an integer would silently be
  # dropped and fall back to the automatic cutoff.
  extra_ctor = quote(
    if (!is.null(max_rank)) {
      max_rank <- as.double(max_rank)
    }
  ),
  fields = list(
    auc_type = p_choice(
      "recovery",
      c("recovery", "wilcox", "ap"),
      doc = paste(
        "Which statistic to calculate. `\"wilcox\"` is the AUC",
        "derived from the Mann-Whitney U statistic over the full",
        "ranking, with the null at 0.5 for any gene set size.",
        "`\"recovery\"` is the recovery-curve AUC under a rank",
        "cutoff, i.e. the actual AUCell statistic of Aibar, et al.",
        "`\"ap\"` is average precision, the most top-heavy of the",
        "three, but its null tracks the gene set prevalence so raw",
        "values are not comparable across gene sets of different size",
        "unless `standardise` is on."
      )
    ),
    max_rank = p_dbl(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Rank cutoff for `\"recovery\"`, counted from the top of the",
        "within-cell ranking. If `NULL`, resolves to the top 5% of",
        "the gene universe, following Aibar, et al. Ignored by the",
        "other two statistics."
      )
    ),
    standardise = p_lgl(
      FALSE,
      doc = paste(
        "Shall each gene set's scores be z-scored across the cells.",
        "This is what makes `\"ap\"` comparable across gene sets of",
        "different size."
      )
    )
  )
)

spec_scenic_binarise <- param_spec(
  name = "scenic_binarise",
  title = "Wrapper function for parameters for the SCENIC binarisation",
  description = paste(
    "Each regulon gets its own threshold. A two-component",
    "Gaussian mixture is fitted and compared against a single",
    "Gaussian by BIC; if the mixture wins, the threshold is the",
    "kernel density minimum between the two component means,",
    "otherwise it falls back to `mean + 2 * sd`. This follows",
    "pySCENIC. AUCell fits six candidates and then lets the",
    "density trough override all of them whenever one exists, so",
    "the two land in much the same place. Turn `bw_adjust` up if",
    "shallow wobbles in the density are being picked up as",
    "troughs. AUCell effectively runs at `2`."
  ),
  references = "Aibar, et al., Nat Methods, 2017",
  checker = "ScenicBinarise",
  label = "SCENIC binarisation params",
  hint = paste(
    "bw_adjust must be a positive numeric; n_grid must be numeric >=",
    "3; n_bins must be numeric >= 2."
  ),
  # Rust parses these with as_real(), so integers need to go over as doubles.
  extra_ctor = quote({
    bw_adjust <- as.double(bw_adjust)
    n_grid <- as.double(n_grid)
    n_bins <- as.double(n_bins)
  }),
  fields = list(
    bw_adjust = p_dbl(
      1,
      range = "(0,)",
      doc = paste(
        "Multiplier on the Silverman bandwidth of the kernel density",
        "estimate. Higher values smooth more."
      )
    ),
    n_grid = p_int(
      512L,
      range = "[3,)",
      check_as = "N1[3,)",
      doc = paste(
        "Number of points at which the density is evaluated between",
        "the two component means."
      )
    ),
    n_bins = p_int(
      512L,
      range = "[2,)",
      check_as = "N1[2,)",
      doc = "Number of histogram bins used to approximate the density."
    )
  )
)

spec_sc_hotspot <- param_spec(
  name = "sc_hotspot",
  title = "Wrapper function for parameters for HotSpot",
  description = paste(
    "`weighted_graph` controls how the kNN distances become edge",
    "weights. The default of `FALSE` follows the reference",
    "implementation: the distances only decide who is a neighbour",
    "and every retained edge weighs one. Set it to `TRUE` for the",
    "Gaussian kernel, whose width is the `ceil(k /",
    "neighborhood_factor)`-th neighbour distance."
  ),
  references = "DeTomaso and Yosef, Cell Systems, 2021",
  checker = "ScHotspot",
  label = "HotSpot params",
  hint = paste(
    "normalise and weighted_graph must be booleans,",
    "neighborhood_factor a positive number."
  ),
  fields = list(
    model = p_choice(
      "danb",
      c("danb", "normal", "bernoulli"),
      doc = "Model to use for modelling the GEX."
    ),
    normalise = p_lgl(TRUE, doc = "Shall the data be normalised."),
    weighted_graph = p_lgl(
      FALSE,
      doc = paste(
        "Shall the Gaussian kernel be applied to the neighbour",
        "distances."
      )
    ),
    neighborhood_factor = p_dbl(
      3,
      range = "(0,)",
      strict = TRUE,
      doc = paste(
        "Kernel width is the `ceil(k / neighborhood_factor)`-th",
        "neighbour distance. Only read when `weighted_graph = TRUE`."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`."
      )
    )
  )
)

spec_dialogue_pmd <- param_spec(
  name = "dialogue_pmd",
  title = "Wrapper function for the DIALOGUE decomposition parameters",
  description = paste(
    "Stage one of DIALOGUE: the penalised matrix decomposition",
    "that turns the per-cell-type features into multicellular",
    "programmes, and the provisional gene signatures that come",
    "off it."
  ),
  details = "The defaults follow upstream's `DLG.get.param`. Two knobs are worth thinking\nabout before anything else. `k` is how many programmes you are asking for,\nand there is no sweep to help you pick it. `n_permutations` sets the\nresolution of the empirical p-value: with the default of `100` the smallest\np you can observe is `0.01`, so lower it for a quick look and leave it alone\nfor anything you intend to believe.\n\n`averaging` is exposed and honoured here. Upstream takes the same argument\nand then ignores it, hard-coding column medians, so `\"median\"` is what every\npublished DIALOGUE run actually used.",
  references = "Jerby-Arnon & Regev, Nature Biotechnology, 2022",
  checker = "DialoguePmd",
  label = "DIALOGUE decomposition params",
  hint = paste(
    "k >= 1, n_permutations >= 2, cap in [0, 0.5), the p-value",
    "cutoffs in (0, 1] and min_ci in [0, 1]."
  ),
  fields = list(
    k = p_int(
      2L,
      range = "[1,)",
      doc = paste(
        "Number of multicellular programmes to extract. Must be at",
        "least 1."
      )
    ),
    n_permutations = p_int(
      100L,
      range = "[2,)",
      doc = paste(
        "Permutations backing the empirical p-value per programme.",
        "Must be at least 2."
      )
    ),
    extra_sparse = p_lgl(
      FALSE,
      doc = paste(
        "Tune the L1 bound by permutation instead of fixing it at",
        "`sqrt(p_1) / 2`. Costs ten more fits per permutation."
      )
    ),
    abn_c = p_int(
      15L,
      range = "[0,)",
      doc = paste(
        "Minimum cells a sample must contribute, within a cell type,",
        "before it counts towards the feature-level ANOVA."
      )
    ),
    p_anova = p_dbl(
      0.05,
      range = "(0,1]",
      doc = paste(
        "BH-adjusted ANOVA cutoff for keeping a feature. Must be in",
        "`(0, 1]`."
      )
    ),
    centre = p_lgl(
      TRUE,
      doc = paste(
        "Centre and scale the sample-level feature matrix, then",
        "winsorise it."
      )
    ),
    cap = p_dbl(
      0.01,
      range = "[0,0.5)",
      doc = paste(
        "Winsorising tail fraction applied to each column. Must be in",
        "`[0, 0.5)`."
      )
    ),
    spatial = p_lgl(
      FALSE,
      doc = paste(
        "Spatial data: skip the ANOVA feature filter entirely. Niches",
        "are small, so a feature need not vary across them to be",
        "real."
      )
    ),
    n_genes = p_int(
      200L,
      range = "[1,)",
      doc = paste(
        "Genes taken per programme per direction when building a",
        "signature."
      )
    ),
    min_ci = p_dbl(
      0.05,
      range = "[0,1]",
      doc = paste(
        "Minimum absolute correlation for a gene to enter a",
        "signature. Must be in `[0, 1]`."
      )
    ),
    averaging = p_choice(
      "median",
      c("median", "mean"),
      doc = "How cell-level features are collapsed per sample."
    ),
    mcp_assignment_p = p_dbl(
      0.1,
      range = "(0,1]",
      doc = paste(
        "Empirical p below which a cell type pair counts as connected",
        "when deciding which cell types a programme spans. Must be in",
        "`(0, 1]`."
      )
    ),
    seed = p_int(1234L, doc = "Seed for the permutation null.")
  )
)

spec_dialogue_hlm <- param_spec(
  name = "dialogue_hlm",
  title = "Wrapper function for the DIALOGUE mixed model parameters",
  description = paste(
    "Stage two of DIALOGUE: for every ordered pair of cell types",
    "and every candidate gene, a random-intercept mixed model",
    "over samples asking whether a cell's own programme score",
    "tracks the partner cell type's expression of that gene in",
    "the same sample."
  ),
  details = "This stage dominates the runtime, and `satterthwaite` is the knob that\ndecides how badly. Turning it off falls back to the residual count for the\ndenominator degrees of freedom, which is far cheaper and barely differs once\na cell type has thousands of cells.\n\n`use_cell_quality` conditions on the cell's own quality covariate, which\nstage one has already regressed out of the scores by ordinary least squares.\nThe default conditions on it twice, because upstream does.",
  references = "Jerby-Arnon & Regev, Nature Biotechnology, 2022",
  checker = "DialogueHlm",
  label = "DIALOGUE mixed model params",
  hint = paste(
    "min_cells_per_sample must be a non-negative integer, the",
    "rest are booleans."
  ),
  fields = list(
    min_cells_per_sample = p_int(
      2L,
      range = "[0,)",
      doc = paste(
        "Minimum cells a sample must contribute, in *both* cell types",
        "of a pair, before it takes part in that pair's models."
      )
    ),
    use_tme_qc = p_lgl(
      TRUE,
      doc = paste(
        "Include the partner cell type's mean quality in that sample",
        "as a fixed effect. Upstream's `tme.qc`."
      )
    ),
    use_cell_quality = p_lgl(
      TRUE,
      doc = paste(
        "Include the responding cell's own quality as a fixed effect.",
        "Upstream's `cellQ`."
      )
    ),
    satterthwaite = p_lgl(
      TRUE,
      doc = paste(
        "Compute Satterthwaite denominator degrees of freedom, as",
        "`lmerTest` does."
      )
    )
  )
)

spec_dialogue_refine <- param_spec(
  name = "dialogue_refine",
  title = "Wrapper function for the DIALOGUE refinement parameters",
  description = paste(
    "Stage three of DIALOGUE: the cross-partner meta-analysis",
    "that decides which genes survive, and the non-negative refit",
    "of the programme scores onto them."
  ),
  details = "Two gene lists come out. The permissive one asks only for a Fisher-combined\np below `permissive_p`; the strict one is looser on the p-value but also\ndemands that *every* partner supports the gene. They are not nested by\nthreshold, they are nested by evidence, and the strict list is the one to\nquote.",
  references = "Jerby-Arnon & Regev, Nature Biotechnology, 2022",
  checker = "DialogueRefine",
  label = "DIALOGUE refinement params",
  hint = paste(
    "the p-value cutoffs and early_stop_cor must be in (0, 1],",
    "min_support_fraction in [0, 1] and min_stratum a",
    "non-negative integer."
  ),
  fields = list(
    support_p = p_dbl(
      0.1,
      range = "(0,1]",
      doc = paste(
        "Adjusted p below which one partner counts as supporting a",
        "gene. Must be in `(0, 1]`."
      )
    ),
    min_support_fraction = p_dbl(
      1 / 3,
      range = "[0,1]",
      doc = paste(
        "Minimum supporting fraction for a stratum to enter the",
        "staged fit. Must be in `[0, 1]`."
      )
    ),
    min_stratum = p_int(
      5L,
      range = "[0,)",
      doc = "Minimum genes in a stratum before it is worth fitting."
    ),
    early_stop_cor = p_dbl(
      0.95,
      range = "(0,1]",
      doc = paste(
        "Correlation between the original score and the running fit",
        "at which the staged fit stops early. Must be in `(0, 1]`."
      )
    ),
    permissive_p = p_dbl(
      0.001,
      range = "(0,1]",
      doc = paste(
        "Fisher-combined p for the permissive gene list, where a gene",
        "is carried by partner support rather than by a positive",
        "coefficient. Must be in `(0, 1]`."
      )
    ),
    strict_p = p_dbl(
      0.05,
      range = "(0,1]",
      doc = paste(
        "Fisher-combined p for the strict gene list, which also",
        "demands that every partner supports the gene. Must be in",
        "`(0, 1]`."
      )
    )
  )
)


spec_sc_bt_metacells <- param_spec(
  name = "sc_bt_metacells",
  title = paste(
    "Wrapper function for parameters for bootstrapped meta cell",
    "generation"
  ),
  description = paste(
    "This function generates parameters for the bootstrapped meta",
    "cell generation based on hdWGCNA, see Morabito, et al., Cell",
    "Rep. Methods, 2023."
  ),
  checker = "ScBootstrappedMetacells",
  label = "bootstrapped metacell params",
  hint = paste(
    "max_shared, target_no_metacells and max_iter must be",
    "integers >= 1."
  ),
  fields = list(
    max_shared = p_int(
      15L,
      range = "[1,)",
      doc = paste(
        "Maximum number of allowed shared neighbours for the meta",
        "cell to be considered."
      )
    ),
    target_no_metacells = p_int(
      1000L,
      range = "[1,)",
      doc = "Target number of meta-cells to generate."
    ),
    max_iter = p_int(
      5000L,
      range = "[1,)",
      doc = "Maximum number of iterations for the algorithm."
    ),
    knn = p_merge(
      "knn_defaults",
      overrides = list(k = 25L, ann_dist = "cosine"),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`."
      )
    )
  )
)

spec_sc_seacells <- param_spec(
  name = "sc_seacells",
  title = "Wrapper function for the SEACells parameters",
  checker = "ScSeacells",
  label = "SEACells params",
  hint = paste(
    "n_sea_cells, max_fw_iters, max_iter, min_iter and",
    "greedy_threshold must be integers >= 1; convergence_epsilon",
    "and pruning_threshold must be numeric; pruning must be a",
    "boolean; graph_building must be a string; n_landmarks must",
    "be an integer or NULL."
  ),
  fields = list(
    n_sea_cells = p_int(range = "[1,)", doc = "Number of SEA cells to detect."),
    max_fw_iters = p_int(
      50L,
      range = "[1,)",
      doc = "Maximum iterations for the Franke-Wolfe algorithm."
    ),
    convergence_epsilon = p_dbl(
      0.001,
      doc = paste(
        "Convergence threshold. Algorithm stops when RSS change <",
        "epsilon * RSS(0)."
      )
    ),
    max_iter = p_int(
      100L,
      range = "[1,)",
      doc = "Maximum iterations to run SEACells for."
    ),
    min_iter = p_int(
      10L,
      range = "[1,)",
      doc = "Minimum iterations to run SEACells for."
    ),
    greedy_threshold = p_int(
      20000L,
      range = "[1,)",
      doc = paste(
        "Maximum number of cells before defaulting to rapid random",
        "selection of archetypes."
      )
    ),
    graph_building = p_chr("union", doc = "Graph building method."),
    pruning = p_lgl(
      TRUE,
      doc = "Shall tiny values be pruned during Franke-Wolfe updates."
    ),
    pruning_threshold = p_dbl(
      1e-07,
      doc = paste(
        "If `pruning = TRUE` values below which threshold shall be",
        "pruned."
      )
    ),
    n_landmarks = p_int(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "If provided, it will use the Nystroem extension during the",
        "archetype finding. Useful for larger data sets."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      overrides = list(k = 25L),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`."
      )
    )
  )
)

spec_sc_supercell <- param_spec(
  name = "sc_supercell",
  title = "Wrapper function for parameters for SuperCell generation",
  checker = "ScSupercell",
  label = "SuperCell params",
  hint = paste(
    "walk_length must be an integer >= 1; k_ith and max_support",
    "must be an integer or NULL; graining_factor must be numeric;",
    "use_kernel must be a boolean."
  ),
  fields = list(
    walk_length = p_int(
      3L,
      range = "[1,)",
      doc = "Walk length for the Walktrap algorithm."
    ),
    graining_factor = p_dbl(
      20,
      doc = paste(
        "Graining level of data (proportion of number of single cells",
        "in the initial dataset to the number of metacells in the",
        "final dataset). (One meta cell per 20 cells.)"
      )
    ),
    use_kernel = p_lgl(
      TRUE,
      doc = paste(
        "Shall a kernel function akin to MAGIC be applied akin to the",
        "approach in SuperCell2, see Hérault, et al., bioRxiv, 2026",
        "and van Dijk, et al., Cell, 2018."
      )
    ),
    k_ith = p_int(
      NULL,
      null_ok = TRUE,
      doc = "The k-ith neighbour to use for the kernel."
    ),
    max_support = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Caps each cell's walk-probability vector to its top entries",
        "by mass, bounding memory at ~`max_support * n_cells` on",
        "large data. Makes the result an approximation. `NULL`",
        "(default) keeps the walks exact."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      overrides = list(k = 5L, ann_dist = "cosine"),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`."
      )
    )
  )
)

spec_sc_bbknn <- param_spec(
  name = "sc_bbknn",
  title = "Wrapper function for the BBKNN parameters",
  checker = "ScBbknn",
  label = "BBKNN params",
  hint = paste(
    "neighbours_within_batch must be an integer >= 1; trim must",
    "be NULL or an integer >= 1; set_op_mix_ratio must be a",
    "numeric in [0, 1]; local_connectivity must be numeric."
  ),
  fields = list(
    neighbours_within_batch = p_int(
      3L,
      range = "[1,)",
      doc = "Number of neighbours to consider per batch."
    ),
    set_op_mix_ratio = p_dbl(
      1,
      range = "[0,1]",
      doc = "Mixing ratio between union (1.0) and intersection (0.0)."
    ),
    local_connectivity = p_dbl(
      1,
      doc = paste(
        "UMAP connectivity computation parameter, how many nearest",
        "neighbours of each cell are assumed to be fully connected."
      )
    ),
    trim = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Trim the neighbours of each cell to these many top",
        "connectivities. May help with population independence and",
        "improve the tidiness of clustering. If `NULL`, it defaults",
        "to `10 * neighbours_within_batch`."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      overrides = quote(list(k = neighbours_within_batch * 2L)),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`."
      )
    )
  )
)

spec_sc_fastmnn <- param_spec(
  name = "sc_fastmnn",
  title = "Wrapper function for the fastMNN parameters",
  checker = "ScFastmnn",
  label = "fastMNN params",
  hint = paste(
    "no_pcs must be an integer >= 1; ndist must be a positive",
    "numeric; size_factor must be numeric; cos_norm, randomised,",
    "sparse_svd, mean_center, normalise_variance and clr must be",
    "booleans."
  ),
  fields = list(
    ndist = p_dbl(
      3,
      range = "(0,)",
      doc = "Number of median distances for the tricube kernel bandwidth."
    ),
    cos_norm = p_lgl(
      TRUE,
      doc = "Apply cosine normalisation before computing distances."
    ),
    no_pcs = p_int(
      30L,
      range = "[1,)",
      doc = "Number of PCs to use for MNN calculations."
    ),
    sparse_svd = p_lgl(TRUE, doc = "Shall the sparse SVD be used."),
    knn = p_merge(
      "knn_defaults",
      default = list(k = 20L),
      overrides = list(k = 20L, ann_dist = "cosine"),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`."
      )
    ),
    pca = p_merge(
      "sc_pca",
      default = quote(params_sc_pca()),
      doc = paste(
        "Parameters to feed through to the optional recalculation of",
        "the PCA, see [params_sc_pca()]."
      )
    )
  )
)

spec_sc_harmony <- param_spec(
  name = "sc_harmony",
  title = "Default parameters for Harmony batch correction",
  checker = "ScHarmony",
  label = "Harmony params",
  hint = paste(
    "max_iter_kmeans, max_iter_harmony and window_size must be",
    "integers >= 1; k must be NULL or an integer; sigma, theta and",
    "lambda must be numeric vectors with non-negative values;",
    "block_size must be in (0, 1]; epsilon_kmeans and epsilon_harmony",
    "must be > 0."
  ),
  class_tag = "params_sc_harmony",
  fields = list(
    k = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Number of clusters for k-means clustering. If not provided,",
        "it will be automatically determined as `min(round(N / 30),",
        "100)`."
      )
    ),
    sigma = p_dbl(
      0.1,
      range = "[0,)",
      len = "+",
      doc = paste(
        "Per-cluster diversity weights. Either a single value",
        "(broadcast to all clusters) or a vector of length k."
      )
    ),
    theta = p_dbl(
      2,
      range = "[0,)",
      len = "+",
      doc = paste(
        "Per-variable diversity penalties. Either a single value",
        "(broadcast to all variables) or a vector of length equal to",
        "the number of batch variables."
      )
    ),
    lambda = p_dbl(
      1,
      range = "[0,)",
      len = "+",
      doc = paste(
        "Ridge regression penalty for the linear model. Typically a",
        "single value that is broadcast to all design matrix columns."
      )
    ),
    block_size = p_dbl(
      0.2,
      range = "(0,1]",
      doc = paste(
        "Fraction of cells to update per block during optimisation",
        "(0.0-1.0). Lower values reduce memory usage but increase",
        "computation time."
      )
    ),
    max_iter_kmeans = p_int(
      20L,
      range = "[1,)",
      doc = "Maximum number of k-means iterations per Harmony round."
    ),
    max_iter_harmony = p_int(
      10L,
      range = "[1,)",
      doc = "Maximum number of Harmony outer iterations."
    ),
    epsilon_kmeans = p_dbl(
      1e-05,
      range = "(0,)",
      doc = paste(
        "Convergence threshold for k-means clustering. Stops when the",
        "relative change in cluster assignments falls below this",
        "value."
      )
    ),
    epsilon_harmony = p_dbl(
      1e-04,
      range = "(0,)",
      doc = paste(
        "Convergence threshold for Harmony. Stops when the relative",
        "change in the objective function falls below this value."
      )
    ),
    window_size = p_int(
      2L,
      range = "[1,)",
      doc = paste(
        "Number of previous iterations to consider when checking",
        "convergence."
      )
    ),
    kmeans = p_merge(
      "kmeans_defaults",
      doc = paste(
        "Optional overrides for the k-means clustering algorithm",
        "Possible parameters are `\"k_means_iter\"`,",
        "`\"k_means_init\"`, `\"gemm\"` and `\"hamerly\"`, see",
        "[params_kmeans_defaults()]."
      )
    )
  )
)

spec_sc_harmony_v2 <- param_spec(
  name = "sc_harmony_v2",
  title = "Default parameters for Harmony v2 batch correction",
  checker = "ScHarmonyV2",
  label = "Harmony v2 params",
  hint = paste(
    "max_iter_kmeans, max_iter_harmony and window_size must be",
    "integers >= 1; k must be NULL or an integer; sigma, theta and",
    "lambda must be numeric vectors with non-negative values;",
    "block_size must be in (0, 1]; epsilon_kmeans, epsilon_harmony",
    "and batch_proportion_cutoff must be > 0; alpha must be in (0,",
    "1); tau must be >= 0; use_dynamic_lambda must be a single",
    "logical."
  ),
  class_tag = "params_sc_harmony_v2",
  fields = list(
    k = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Number of clusters for k-means clustering. If not provided,",
        "it will be automatically determined as `min(round(N / 30),",
        "100)`."
      )
    ),
    sigma = p_dbl(
      0.1,
      range = "[0,)",
      len = "+",
      doc = paste(
        "Per-cluster diversity weights. Either a single value",
        "(broadcast to all clusters) or a vector of length k."
      )
    ),
    theta = p_dbl(
      2,
      range = "[0,)",
      len = "+",
      doc = paste(
        "Per-variable diversity penalties. Either a single value",
        "(broadcast to all variables) or a vector of length equal to",
        "the number of batch variables."
      )
    ),
    lambda = p_dbl(
      1,
      range = "[0,)",
      len = "+",
      doc = paste(
        "Ridge regression penalty for the linear model. Typically a",
        "single value that is broadcast to all design matrix columns.",
        "Ignored when `use_dynamic_lambda = TRUE`."
      )
    ),
    block_size = p_dbl(
      0.2,
      range = "(0,1]",
      doc = paste(
        "Fraction of cells to update per block during optimisation",
        "(0.0-1.0). Lower values reduce memory usage but increase",
        "computation time."
      )
    ),
    max_iter_kmeans = p_int(
      4L,
      range = "[1,)",
      doc = "Maximum number of k-means iterations per Harmony round."
    ),
    max_iter_harmony = p_int(
      10L,
      range = "[1,)",
      doc = "Maximum number of Harmony outer iterations."
    ),
    epsilon_kmeans = p_dbl(
      0.001,
      range = "(0,)",
      doc = paste(
        "Convergence threshold for k-means clustering. Stops when the",
        "relative change in cluster assignments falls below this",
        "value."
      )
    ),
    epsilon_harmony = p_dbl(
      0.01,
      range = "(0,)",
      doc = paste(
        "Convergence threshold for Harmony. Stops when the relative",
        "change in the objective function falls below this value."
      )
    ),
    window_size = p_int(
      3L,
      range = "[1,)",
      doc = paste(
        "Number of previous iterations to consider when checking",
        "convergence."
      )
    ),
    alpha = p_dbl(
      0.2,
      range = "(0,1)",
      doc = paste(
        "Scaling factor for dynamic lambda estimation. Must be in (0,",
        "1). Only relevant when `use_dynamic_lambda = TRUE`."
      )
    ),
    tau = p_dbl(
      0,
      range = "[0,)",
      doc = paste(
        "Scaling factor for theta based on batch size. A value of 0",
        "disables batch-size scaling of theta."
      )
    ),
    batch_proportion_cutoff = p_dbl(
      1e-05,
      range = "(0,)",
      doc = paste(
        "Cutoff for pruning batches with small proportions during",
        "ridge regression."
      )
    ),
    use_dynamic_lambda = p_lgl(
      FALSE,
      doc = paste(
        "If `TRUE`, lambda is estimated dynamically per cluster",
        "instead of using the fixed `lambda` value."
      )
    ),
    kmeans = p_merge(
      "kmeans_defaults",
      doc = paste(
        "Optional overrides for the k-means clustering algorithm",
        "Possible parameters are `\"k_means_iter\"`,",
        "`\"k_means_init\"`, `\"gemm\"` and `\"hamerly\"`, see",
        "[params_kmeans_defaults()]."
      )
    )
  )
)

spec_sc_seurat_cca <- param_spec(
  name = "sc_seurat_cca",
  title = "Wrapper function for the Seurat CCA parameters",
  references = "Stuart, et al., Cell, 2019",
  checker = "ScSeuratCca",
  label = "Seurat CCA params",
  fields = list(
    num_cc = p_int(
      30L,
      range = "[1,)",
      doc = paste(
        "Number of canonical correlation dimensions to compute for",
        "the anchor space. The effective rank used is `max(num_cc,",
        "dims)`."
      )
    ),
    dims = p_int(
      30L,
      range = "[1,)",
      doc = paste(
        "Number of dimensions used for the anchor kNN queries and the",
        "size of the returned embedding."
      )
    ),
    k_anchor = p_int(
      5L,
      range = "[1,)",
      doc = paste(
        "Neighbourhood size for the mutual nearest neighbour anchor",
        "search."
      )
    ),
    k_filter = p_int(
      200L,
      range = "[0,)",
      doc = "Neighbourhood size for the gene-space anchor filter."
    ),
    k_score = p_int(
      30L,
      range = "[1,)",
      doc = "Neighbourhood size for the shared-neighbour anchor scoring."
    ),
    k_weight = p_int(
      100L,
      range = "[1,)",
      doc = paste(
        "Neighbourhood size for the kernel weights applied during the",
        "correction."
      )
    ),
    n_top_features = p_int(
      200L,
      range = "[1,)",
      doc = paste(
        "Number of top-loading genes used for the gene-space anchor",
        "filter."
      )
    ),
    l2_norm = p_lgl(
      TRUE,
      doc = paste(
        "Shall the canonical correlation embedding be L2-normalised",
        "per cell."
      )
    ),
    sd = p_dbl(
      1,
      range = "(0,)",
      doc = paste(
        "Bandwidth divisor of the Gaussian kernel used for the anchor",
        "weights."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      overrides = list(ann_dist = "cosine"),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`. Note",
        "that `k` is unused here, the neighbourhood sizes come from",
        "`k_anchor`, `k_filter`, `k_score` and `k_weight`."
      )
    ),
    pca = p_merge(
      "sc_pca",
      default = quote(params_sc_pca()),
      doc = paste(
        "Parameters to feed through to the optional recalculation of",
        "the PCA, see [params_sc_pca()]."
      )
    )
  )
)

spec_sc_seurat_rpca <- param_spec(
  name = "sc_seurat_rpca",
  title = "Wrapper function for the Seurat rPCA parameters",
  references = "Stuart, et al., Cell, 2019",
  checker = "ScSeuratRpca",
  label = "Seurat rPCA params",
  fields = list(
    dims = p_int(
      30L,
      range = "[1,)",
      doc = paste(
        "Number of dimensions used for the per-batch PCA projections,",
        "the anchor kNN queries and the size of the returned",
        "embedding."
      )
    ),
    k_anchor = p_int(
      5L,
      range = "[1,)",
      doc = paste(
        "Neighbourhood size for the mutual nearest neighbour anchor",
        "search."
      )
    ),
    k_score = p_int(
      30L,
      range = "[1,)",
      doc = "Neighbourhood size for the shared-neighbour anchor scoring."
    ),
    k_weight = p_int(
      100L,
      range = "[1,)",
      doc = paste(
        "Neighbourhood size for the kernel weights applied during the",
        "correction."
      )
    ),
    l2_norm = p_lgl(
      TRUE,
      doc = "Shall the projected embeddings be L2-normalised per cell."
    ),
    sd = p_dbl(
      1,
      range = "(0,)",
      doc = paste(
        "Bandwidth divisor of the Gaussian kernel used for the anchor",
        "weights."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      overrides = list(ann_dist = "cosine"),
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`. Note",
        "that `k` is unused here, the neighbourhood sizes come from",
        "`k_anchor`, `k_score` and `k_weight`."
      )
    ),
    pca = p_merge(
      "sc_pca",
      default = quote(params_sc_pca()),
      doc = paste(
        "Parameters to feed through to the optional recalculation of",
        "the PCA, see [params_sc_pca()]."
      )
    )
  )
)

spec_scenic_random_forest_defaults <- param_defaults(
  name = "scenic_random_forest_defaults",
  title = paste(
    "Default parameters for the SCENIC RandomForest regression",
    "learner"
  ),
  checker = NULL,
  fields = list(
    n_trees = p_free(250L, doc = "TODO"),
    min_samples_leaf = p_free(50L, doc = "TODO"),
    n_features_split = p_free(0L, doc = "TODO"),
    subsample_rate = p_free(0.632, doc = "TODO"),
    bootstrap = p_free(FALSE, doc = "TODO"),
    max_depth = p_free(8L, doc = "TODO"),
    subsample_frac = p_free(NULL, doc = "TODO")
  )
)

spec_scenic_extra_trees_defaults <- param_defaults(
  name = "scenic_extra_trees_defaults",
  title = paste(
    "Default parameters for the SCENIC ExtraTrees regression",
    "learner"
  ),
  checker = NULL,
  fields = list(
    n_trees = p_free(500L, doc = "TODO"),
    min_samples_leaf = p_free(50L, doc = "TODO"),
    n_features_split = p_free(0L, doc = "TODO"),
    n_thresholds = p_free(1L, doc = "TODO"),
    max_depth = p_free(8L, doc = "TODO"),
    subsample_frac = p_free(NULL, doc = "TODO")
  )
)

spec_scenic_gradient_boosting_defaults <- param_defaults(
  name = "scenic_gradient_boosting_defaults",
  title = paste(
    "Default parameters for the SCENIC GradientBoosting",
    "(GRNBoost2) regression learner"
  ),
  checker = NULL,
  fields = list(
    n_trees_max = p_free(1000L, doc = "TODO"),
    learning_rate = p_free(0.01, doc = "TODO"),
    max_depth = p_free(3L, doc = "TODO"),
    min_samples_leaf = p_free(50L, doc = "TODO"),
    early_stop_window = p_free(25L, doc = "TODO"),
    subsample_rate = p_free(0.9, doc = "TODO"),
    n_features_split = p_free(0L, doc = "TODO")
  )
)

spec_scenic <- param_spec(
  name = "scenic",
  title = "Constructor for SCENIC parameters",
  checker = "Scenic",
  label = "SCENIC params",
  hint = paste(
    "min_counts, n_pcs, n_subsample, min_samples_leaf and max_depth",
    "must be integers >= 1; n_features_split must be an integer >= 0;",
    "min_cells must be in (0, 1]."
  ),
  extra_check = quote({
    res <- check_list_shape(
      x,
      c("min_samples_leaf", "n_features_split", "max_depth")
    )
    if (!isTRUE(res)) {
      return(res)
    }
    res <- apply_qtest_rules(
      x,
      list(
        min_samples_leaf = "I1[1,)",
        n_features_split = "I1[0,)",
        max_depth = "I1[1,)"
      ),
      label = "SCENIC params"
    )
    if (!isTRUE(res)) {
      return(res)
    }
    if (x$learner_type == "randomforest") {
      if (is.null(x$n_trees) || !checkmate::qtest(x$n_trees, "I1[1,)")) {
        return("n_trees must be a positive integer for randomforest.")
      }
      if (
        is.null(x$subsample_rate) ||
          !checkmate::qtest(x$subsample_rate, "N1(0,1]")
      ) {
        return("subsample_rate must be a numeric in (0, 1] for randomforest.")
      }
      if (is.null(x$bootstrap) || !checkmate::qtest(x$bootstrap, "B1")) {
        return("bootstrap must be a single logical for randomforest.")
      }
      if (
        !is.null(x$subsample_frac) &&
          !checkmate::qtest(x$subsample_frac, "N1(0,1]")
      ) {
        return("subsample_frac must be a numeric in (0, 1] or NULL.")
      }
    }
    if (x$learner_type == "extratrees") {
      if (is.null(x$n_trees) || !checkmate::qtest(x$n_trees, "I1[1,)")) {
        return("n_trees must be a positive integer for extratrees.")
      }
      if (
        is.null(x$n_thresholds) ||
          !checkmate::qtest(x$n_thresholds, "I1[1,)")
      ) {
        return("n_thresholds must be a positive integer for extratrees.")
      }
      if (
        !is.null(x$subsample_frac) &&
          !checkmate::qtest(x$subsample_frac, "N1(0,1]")
      ) {
        return("subsample_frac must be a numeric in (0, 1] or NULL.")
      }
    }
    if (x$learner_type == "grnboost2") {
      if (
        is.null(x$n_trees_max) ||
          !checkmate::qtest(x$n_trees_max, "I1[1,)")
      ) {
        return("n_trees_max must be a positive integer for grnboost2.")
      }
      if (
        is.null(x$learning_rate) ||
          !checkmate::qtest(x$learning_rate, "N1(0,1]")
      ) {
        return("learning_rate must be a numeric in (0, 1] for grnboost2.")
      }
      if (
        is.null(x$early_stop_window) ||
          !checkmate::qtest(x$early_stop_window, "I1[1,)")
      ) {
        return("early_stop_window must be a positive integer for grnboost2.")
      }
      if (
        is.null(x$subsample_rate) ||
          !checkmate::qtest(x$subsample_rate, "N1(0,1]")
      ) {
        return("subsample_rate must be a numeric in (0, 1] for grnboost2.")
      }
    }
  }),
  fields = list(
    min_counts = p_int(
      50L,
      range = "[1,)",
      doc = paste(
        "Minimum total counts a gene needs to be included in the",
        "analysis."
      )
    ),
    min_cells = p_dbl(
      0.03,
      range = "(0,1]",
      doc = paste(
        "Minimum proportion of cells (between 0 and 1) that must",
        "express a gene for it to be considered."
      )
    ),
    learner_type = p_choice(
      "randomforest",
      c("randomforest", "extratrees", "grnboost2"),
      doc = "Regression learner to use."
    ),
    gene_batch_strategy = p_choice(
      "correlated",
      c("random", "correlated"),
      doc = paste(
        "Strategy for grouping target genes into batches. Only used",
        "for `\"randomforest\"` and `\"extratrees\"` learners;",
        "ignored for `\"grnboost2\"`."
      )
    ),
    gene_batch_size = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Number of genes per batch. If `NULL` (default), the batch",
        "size is determined automatically. Ignored for",
        "`\"grnboost2\"`."
      )
    ),
    n_pcs = p_int(
      50L,
      range = "[1,)",
      doc = "Number of PCs to use for the correlated gene batch strategy."
    ),
    n_subsample = p_int(
      100000L,
      range = "[1,)",
      doc = paste(
        "Cell subsampling threshold for the correlated gene batch",
        "strategy. If the number of cells meets or exceeds this",
        "value, `n_subsample` cells are randomly selected prior to",
        "running randomised SVD."
      )
    ),
    learner_params = p_merge(
      quote(switch(
        learner_type,
        extratrees = params_scenic_extra_trees_defaults(),
        grnboost2 = params_scenic_gradient_boosting_defaults(),
        params_scenic_random_forest_defaults()
      )),
      doc = paste(
        "Optional overrides for the regression learner parameters.",
        "For `\"randomforest\"`, see",
        "[bixverse::params_scenic_random_forest_defaults()]. For",
        "`\"extratrees\"`, see",
        "[bixverse::params_scenic_extra_trees_defaults()]. For",
        "`\"grnboost2\"`, see",
        "[bixverse::params_scenic_gradient_boosting_defaults()]."
      )
    )
  )
)

spec_meld <- param_spec(
  name = "meld",
  title = "Constructor for MELD parameters",
  return_order = c(
    "knn",
    "beta",
    "offset",
    "order",
    "filter",
    "chebyshev_order",
    "lap_type",
    "normalise_indicators"
  ),
  checker = "Meld",
  label = "MELD params",
  hint = paste(
    "beta and order must be positive numerics; offset must be in",
    "[0, 1]; chebyshev_order must be an integer >= 2;",
    "normalise_indicators must be a single logical."
  ),
  fields = list(
    beta = p_dbl(
      60,
      range = "(0,)",
      doc = paste(
        "Smoothing strength; larger values produce smoother",
        "densities. Must be strictly positive."
      )
    ),
    offset = p_dbl(
      0,
      range = "[0,1]",
      doc = paste(
        "Shift of the filter centre in the rescaled spectrum. Must be",
        "in `[0, 1]`."
      )
    ),
    order = p_dbl(
      1,
      range = "(0,)",
      doc = paste(
        "Filter falloff sharpness; larger values approach a square",
        "low-pass. Must be strictly positive."
      )
    ),
    filter = p_choice(
      "heat",
      c("heat", "laplacian"),
      doc = "Filter family to use."
    ),
    chebyshev_order = p_int(
      50L,
      range = "[2,)",
      doc = paste(
        "Number of Chebyshev coefficients (polynomial terms). Must be",
        ">= 2."
      )
    ),
    lap_type = p_choice(
      "combinatorial",
      c("combinatorial", "normalised"),
      doc = "Type of Laplacian to use for spectral filtering."
    ),
    normalise_indicators = p_lgl(
      TRUE,
      doc = paste(
        "If `TRUE`, each column of the indicator matrix is divided by",
        "its column sum before filtering, making cross-condition",
        "densities comparable regardless of cells-per-condition."
      )
    ),
    knn = p_merge(
      "knn_defaults",
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`."
      )
    )
  )
)

spec_sc_wnn <- param_spec(
  name = "sc_wnn",
  title = "Wrapper function for WNN parameters",
  checker = "ScWnn",
  label = "WNN params",
  hint = paste(
    "k_nn, knn_range and s_nn must be integers >= 1; sigma_idx must",
    "be an integer >= 0; sd_scale, kernel_power and sigma_floor must",
    "be > 0; cross_const must be >= 0."
  ),
  fields = list(
    k_nn = p_int(
      20L,
      range = "[1,)",
      doc = "Final number of multimodal neighbours per cell."
    ),
    knn_range = p_int(
      200L,
      range = "[1,)",
      doc = paste(
        "Candidate pool size per modality. Each cell's kNN input must",
        "contain at least this many neighbours."
      )
    ),
    sigma_method = p_choice(
      "snn_farthest",
      c("snn_farthest", "sigma_idx"),
      doc = "Bandwidth method."
    ),
    sigma_idx = p_int(
      19L,
      range = "[0,)",
      doc = paste(
        "`\"sigma_idx\"` only: 0-based kNN index for bandwidth.e.",
        "`k_nn - 1`)."
      )
    ),
    snn_type = p_choice(
      "full_connection",
      c("full_connection", "limited"),
      doc = paste(
        "sNN type. The limited version only considers edges that",
        "exist in the kNN."
      )
    ),
    s_nn = p_int(
      20L,
      range = "[1,)",
      doc = paste(
        "`\"snn_farthest\"` only: kNN size used to build the SNN",
        "graph."
      )
    ),
    sd_scale = p_dbl(1, range = "(0,)", doc = "Multiplier on sigma."),
    kernel_power = p_dbl(1, range = "(0,)", doc = "Kernel exponent power."),
    cross_const = p_dbl(
      1e-04,
      range = "[0,)",
      doc = "Cross-modality kernel stabiliser."
    ),
    sigma_floor = p_dbl(
      1e-08,
      range = "(0,)",
      doc = "Minimum sigma value (avoids division by zero)."
    ),
    knn = p_merge(
      "knn_defaults",
      doc = paste(
        "Optional overrides for kNN parameters. See",
        "[bixverse::params_knn_defaults()] for available parameters:",
        "`k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,",
        "`delta`, `diversify_prob`, `ef_budget`, `extract_knn`, `m`,",
        "`ef_construction`, `ef_search`, `n_list` and `n_probe`."
      )
    )
  )
)

spec_sc_palantir <- param_spec(
  name = "sc_palantir",
  title = "Wrapper function for Palantir parameters",
  description = paste(
    "Parameters controlling the Palantir trajectory inference.",
    "The kNN graph you hand to [bixverse::run_palantir_sc()]",
    "feeds the diffusion kernel. The geodesics are measured over",
    "a second kNN graph that Palantir builds internally on the",
    "multiscale space, and `knn` is what controls that one. It is",
    "a different knob from `k` in the kNN parameter block, which",
    "sizes the backend index. Palantir overrides `k` and",
    "`ann_dist` for its internal search, so only `knn_method` and",
    "the backend tuning parameters have an effect."
  ),
  references = "Setty, et al., Nat. Biotechnol., 2019.",
  checker = "ScPalantir",
  label = "Palantir params",
  hint = paste(
    "n_dcs must be >= 3; n_eigs must be NULL or >= 3; knn must be >=",
    "6; num_waypoints must be >= 1; max_iterations must be >= 2;",
    "branch_prob_threshold must be in [0, 1]; scale_components and",
    "use_early_cell_as_start must be booleans."
  ),
  fields = list(
    n_dcs = p_int(
      10L,
      range = "[3,)",
      doc = paste(
        "Diffusion components to extract before the multiscale",
        "scaling."
      )
    ),
    n_eigs = p_int(
      NULL,
      range = "[3,)",
      null_ok = TRUE,
      doc = paste(
        "Eigenvectors to retain, not components: the trivial leading",
        "eigenvector is counted here and then dropped, so `3L` leaves",
        "two multiscale components. If `NULL`, the count is picked",
        "from the largest eigengap, as the reference does."
      )
    ),
    knn = p_int(
      30L,
      range = "[6,)",
      doc = paste(
        "Neighbours for the geodesic graph over the multiscale space,",
        "in the reference's self-inclusive convention."
      )
    ),
    num_waypoints = p_int(
      1200L,
      range = "[1,)",
      doc = "Target waypoint count for the max-min sampler."
    ),
    scale_components = p_lgl(
      TRUE,
      doc = paste(
        "Min-max scale each multiscale component to `[0, 1]` before",
        "any distance is taken."
      )
    ),
    use_early_cell_as_start = p_lgl(
      TRUE,
      doc = paste(
        "Use the provided early cell directly rather than snapping it",
        "to the nearest diffusion-map boundary cell."
      )
    ),
    max_iterations = p_int(
      25L,
      range = "[2,)",
      doc = "Iteration cap for the pseudotime refinement."
    ),
    branch_prob_threshold = p_dbl(
      0.01,
      range = "[0,1]",
      doc = "Fate probabilities below this are zeroed."
    ),
    lanczos_basis_size = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Krylov basis vectors held at once during the diffusion-map",
        "eigendecomposition. If `NULL`, derived from the requested",
        "component count."
      )
    ),
    lanczos_max_restarts = p_int(
      16L,
      range = "[1,)",
      doc = "Maximum restart cycles for the Lanczos solver."
    ),
    lanczos_tol = p_dbl(
      1e-08,
      range = "(0,)",
      doc = "Relative residual tolerance for the Lanczos solver."
    ),
    knn_params = p_merge(
      "knn_defaults",
      doc = paste(
        "Optional overrides for the kNN parameters of the internal",
        "multiscale search. See [bixverse::params_knn_defaults()] for",
        "available parameters: `k`, `knn_method`, `ann_dist`,",
        "`search_budget`, `n_trees`, `delta`, `diversify_prob`,",
        "`ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`",
        "and `n_probe`."
      )
    )
  )
)

spec_sc_magic <- param_spec(
  name = "sc_magic",
  title = "Wrapper function for MAGIC imputation parameters",
  description = paste(
    "Parameters controlling the MAGIC imputation run by",
    "[bixverse::run_magic_sc()]. The defaults mirror the",
    "reference implementation."
  ),
  references = "van Dijk, et al., Cell, 2018.",
  checker = "ScMagic",
  label = "MAGIC params",
  hint = "layer must be one of 'norm' or 'raw'.",
  fields = list(
    n_steps = p_int(
      3L,
      range = "[0,)",
      doc = paste(
        "Diffusion steps applied to the counts. Zero is legal and",
        "hands back the un-imputed values, which is a cheap way to",
        "compare the two."
      )
    ),
    clip_threshold = p_dbl(
      0.01,
      range = "[0,)",
      doc = "Imputed values below this are zeroed after the last step."
    ),
    gene_batch_size = p_int(
      1000L,
      range = "[1,)",
      doc = paste(
        "Genes streamed off the binary store per block. Bounds the",
        "scratch memory and is clamped to the number of requested",
        "genes."
      )
    ),
    layer = p_choice(
      "norm",
      c("norm", "raw"),
      doc = paste(
        "Which stored layer to impute. The operator preserves",
        "per-cell mass, so imputed values sit on the scale of",
        "whatever went in: imputing raw counts and imputing",
        "log-normalised counts are different operations rather than",
        "the same one rescaled."
      )
    ),
    allow_large = p_lgl(
      FALSE,
      doc = paste(
        "Skip the output size guard. The dense output is capped at",
        "1e9 elements, i.e. 4 GB of `f32`."
      )
    )
  )
)

spec_sc_branch_selection <- param_spec(
  name = "sc_branch_selection",
  title = "Wrapper function for the branch cell selection parameters",
  description = paste(
    "Parameters controlling which cells",
    "[bixverse::run_gene_trends_sc()] assigns to each branch. The",
    "threshold on a fate's probability is an expanding quantile",
    "over the pseudotime-sorted cells, made monotone with a",
    "cumulative maximum, so a fate's bar can only rise as",
    "differentiation proceeds. The defaults are the reference",
    "ones."
  ),
  references = "Setty, et al., Nat. Biotechnol., 2019.",
  checker = "ScBranchSelection",
  label = "branch selection params",
  hint = paste(
    "q must be in [0, 1]; eps must be >= 0; resolution must be >=",
    "1."
  ),
  fields = list(
    q = p_dbl(
      0.01,
      range = "[0,1]",
      doc = paste(
        "Upper-tail quantile of the fate probability used as the",
        "threshold."
      )
    ),
    eps = p_dbl(
      0.01,
      range = "[0,)",
      doc = "Slack subtracted from the threshold before the comparison."
    ),
    resolution = p_int(
      500L,
      range = "[1,)",
      doc = "Number of pseudotime buckets, capped at the cell count."
    )
  )
)

spec_sc_gene_trends <- param_spec(
  name = "sc_gene_trends",
  title = "Wrapper function for gene trend parameters",
  description = paste(
    "Parameters controlling the landmark Gaussian process that",
    "[bixverse::run_gene_trends_sc()] fits per branch. The kernel",
    "is a Matern-5/2 one and the prediction grid doubles as the",
    "landmark set. The defaults come from the reference and are",
    "prior-dominated. Palantir's pseudotime is min-max scaled to",
    "`[0, 1]`, so a `length_scale` of `1.0` spans the entire",
    "domain and a `sigma` of `1.0` sits at roughly the signal",
    "scale of log-normalised expression. The posterior will",
    "flatten genuine transient structure and resolve almost any",
    "gene into a smooth monotone or single-peaked curve. That is",
    "a presentation choice, not inference. Shorten `length_scale`",
    "before believing a bump."
  ),
  references = "Setty, et al., Nat. Biotechnol., 2019.",
  checker = "ScGeneTrend",
  label = "gene trend params",
  hint = "weighting must be one of 'hard_mask' or 'fate_probability'.",
  fields = list(
    resolution = p_int(
      500L,
      range = "[2,)",
      doc = paste(
        "Grid points per branch. Kept at the default even when a",
        "branch holds fewer cells, as the reference does."
      )
    ),
    weighting = p_choice(
      "hard_mask",
      c("hard_mask", "fate_probability"),
      doc = paste(
        "With `\"hard_mask\"` every selected cell enters its branch's",
        "fit with equal weight, which is what the reference does.",
        "With `\"fate_probability\"` every cell enters every fit",
        "weighted by its fate probability, which is more defensible:",
        "a cell at 0.6 is not a member."
      )
    ),
    length_scale = p_dbl(1, range = "(0,)", doc = "Matern-5/2 length scale."),
    sigma = p_dbl(1, range = "(0,)", doc = "Noise standard deviation."),
    jitter = p_dbl(
      1e-06,
      range = "[0,)",
      doc = paste(
        "Added to the landmark covariance diagonal before the",
        "Cholesky."
      )
    ),
    max_jitter_retries = p_int(
      3L,
      range = "[0,)",
      doc = paste(
        "Times the jitter is raised and the Cholesky retried before",
        "giving up."
      )
    ),
    chunk_size = p_int(
      2048L,
      range = "[1,)",
      doc = paste(
        "Training points held at once when accumulating the",
        "cross-covariance."
      )
    )
  )
)

spec_lda <- param_spec(
  name = "lda",
  title = "Wrapper function for the LDA parameters",
  description = paste(
    "Solver options for the variational Bayes latent Dirichlet",
    "allocation, see [bixverse::run_lda()]."
  ),
  details = "The defaults follow pycisTopic, so the knobs mean the same thing on both\nsides. `alpha_by_topic = TRUE` turns `alpha` into the Griffiths and\nSteyvers `50 / k` heuristic that cisTopic defaults to; set it to `FALSE` if\nyou want `alpha` taken literally.\n\n`learning = \"batch\"` sweeps every document once per iteration and is\nmonotone in the bound. `\"online\"` takes decaying steps from shuffled\nmini-batches, which reaches a usable fit in far fewer passes on a large\ncorpus at the cost of that guarantee. `batch_size` and `n_epochs` are only\nread by the online variant.",
  references = paste(
    "Hoffman, Blei and Bach, NIPS, 2010; Bravo Gonzalez-Blas, et",
    "al., Nat Methods, 2019"
  ),
  checker = "Lda",
  label = "LDA params",
  hint = paste(
    "alpha, eta, tol and inner_tol must be positive numerics;",
    "max_iter, inner_max_iter, check_every, batch_size and",
    "n_epochs must be positive integers; alpha_by_topic and",
    "eta_by_topic must be booleans."
  ),
  fields = list(
    alpha = p_dbl(
      50,
      range = "(0,)",
      doc = "Dirichlet prior on the document-topic distributions."
    ),
    alpha_by_topic = p_lgl(
      TRUE,
      doc = "Shall `alpha` be divided by the topic count."
    ),
    eta = p_dbl(
      0.1,
      range = "(0,)",
      doc = "Dirichlet prior on the topic-term distributions."
    ),
    eta_by_topic = p_lgl(
      FALSE,
      doc = "Shall `eta` be divided by the topic count."
    ),
    max_iter = p_int(
      150L,
      range = "[1,)",
      doc = paste(
        "Maximum outer iterations. Ignored by the online variant,",
        "which counts epochs instead."
      )
    ),
    tol = p_dbl(
      0.001,
      range = "(0,)",
      doc = "Relative change in the bound below which the solver stops."
    ),
    inner_max_iter = p_int(
      100L,
      range = "[1,)",
      doc = "Maximum fixed-point iterations of the per-document E-step."
    ),
    inner_tol = p_dbl(
      0.001,
      range = "(0,)",
      doc = paste(
        "Relative L1 change in the variational parameters below which",
        "the per-document E-step stops."
      )
    ),
    check_every = p_int(
      10L,
      range = "[1,)",
      doc = "Iterations between bound evaluations."
    ),
    learning = p_choice(
      "batch",
      c("batch", "online"),
      doc = "Batch or online variational inference."
    ),
    batch_size = p_int(
      1024L,
      range = "[1,)",
      doc = "Documents per mini-batch. Online only."
    ),
    n_epochs = p_int(
      10L,
      range = "[1,)",
      doc = "Passes over the corpus. Online only."
    )
  )
)

spec_edger_ql <- param_spec(
  name = "edger_ql",
  title = paste(
    "Wrapper function for parameters for the edgeR",
    "quasi-likelihood workflow"
  ),
  description = paste(
    "Parameters for the edgeR quasi-likelihood chain, implemented",
    "in Rust via the `edge-rs` crate and gated against edgeR",
    "4.8.2. Defaults are edgeR's own. The `legacy` switch picks",
    "between two genuinely different pipelines. The current route",
    "estimates its own dispersion from the most abundant genes",
    "and skips `estimateDisp()`, which is where most of the",
    "runtime went and is edgeR 4's own recommendation. The legacy",
    "route shrinks the raw residual deviance, needs a dispersion",
    "handed to it, and is the only one where the Poisson bound",
    "bites."
  ),
  references = "Chen, Lun and Smyth, F1000Research, 2016",
  checker = "EdgeRQl",
  label = "edgeR QL params",
  hint = paste(
    "min_mean must be a non-negative number; the rest are",
    "booleans."
  ),
  fields = list(
    norm_method = p_choice(
      "TMM",
      c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
      doc = paste(
        "Library size normalisation. `\"none\"` leaves every factor",
        "at one, which is what Milo's `logMS` amounts to."
      )
    ),
    filter = p_lgl(
      TRUE,
      doc = paste(
        "Run `filterByExpr()` before fitting. Turn this off for",
        "anything that is not gene expression, e.g. Milo",
        "neighbourhood counts, where the heuristic means nothing."
      )
    ),
    min_mean = p_dbl(
      0,
      range = "[0,)",
      doc = paste(
        "Drop features whose mean count across samples is below this.",
        "Applied on top of `filter`."
      )
    ),
    robust = p_lgl(
      FALSE,
      doc = paste(
        "Robust empirical Bayes squeezing, giving outlier features",
        "their own smaller prior degrees of freedom."
      )
    ),
    legacy = p_lgl(
      FALSE,
      doc = "Take edgeR's pre-4.0 quasi-likelihood pipeline."
    )
  )
)

spec_limma_voom <- param_spec(
  name = "limma_voom",
  title = "Wrapper function for parameters for the limma-voom workflow",
  description = paste(
    "Parameters for the limma linear model chain, implemented in",
    "Rust via the `edge-rs` crate and gated against limma 3.66.0.",
    "Defaults are limma's own, except for `filter`: inside",
    "[bixverse::BulkDge()] the genes were already filtered by",
    "[bixverse::qc_bulk_dge()]. `route = \"voom\"` is",
    "`voomLmFit()`: precision weights from the mean-variance",
    "trend, then weighted least squares. `route = \"trend\"` is",
    "limma-trend: log-CPM straight into `lmFit()`, with the trend",
    "absorbed by `eBayes(trend = TRUE)`. The empirical Bayes",
    "trend follows the route."
  ),
  references = paste(
    "Law, et al., Genome Biol, 2014; Smyth, Stat Appl Genet Mol",
    "Biol, 2004"
  ),
  checker = "LimmaVoom",
  label = "limma-voom params",
  hint = paste(
    "min_mean must be non-negative, prior_count positive or NULL,",
    "span in (0, 1], proportion in (0, 1); the rest are booleans."
  ),
  fields = list(
    route = p_choice(
      "voom",
      c("voom", "trend"),
      doc = "Whether to run the voom or the limma-trend route."
    ),
    norm_method = p_choice(
      "TMM",
      c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
      doc = "Library size normalisation."
    ),
    filter = p_lgl(FALSE, doc = "Run `filterByExpr()` before fitting."),
    min_mean = p_dbl(
      0,
      range = "[0,)",
      doc = paste(
        "Drop genes whose mean count across samples is below this.",
        "Applied on top of `filter`."
      )
    ),
    robust = p_lgl(
      FALSE,
      doc = "Robust empirical Bayes, `eBayes(robust = TRUE)`."
    ),
    prior_count = p_dbl(
      NULL,
      range = "(0,)",
      null_ok = TRUE,
      doc = paste(
        "Count added before the log. `NULL` takes the route's own",
        "default, `0.5` for voom and `2` for trend."
      )
    ),
    adaptive_span = p_lgl(
      TRUE,
      doc = paste(
        "Derive the lowess span from the number of genes, as limma",
        "does since 3.56. Only used by voom."
      )
    ),
    span = p_dbl(
      0.5,
      range = "(0,1]",
      doc = paste(
        "Lowess span for the voom trend, only read if `adaptive_span",
        "= FALSE`."
      )
    ),
    proportion = p_dbl(
      0.01,
      range = "(0,1)",
      doc = paste(
        "Assumed proportion of differentially expressed genes, only",
        "used for the B-statistic."
      )
    )
  )
)

spec_nebula <- param_spec(
  name = "nebula",
  title = "Wrapper function for parameters for NEBULA",
  description = paste(
    "Parameters for the NEBULA negative binomial gamma mixed",
    "model, implemented in Rust via the `edge-rs` crate and",
    "ported from the `nebula` package's own C++. Defaults are the",
    "R package's own. NEBULA splits the variance into a",
    "subject-level random effect and a cell-level overdispersion.",
    "Run it on meta cells and the cell-level term becomes the",
    "spread between aggregates within a subject rather than",
    "between cells, so it is smaller and absorbs whatever the",
    "aggregation smoothed away. The subject-level term keeps its",
    "meaning either way."
  ),
  references = "He, et al., Commun Biol, 2021",
  checker = "Nebula",
  label = "NEBULA params",
  hint = paste(
    "The overdispersion bounds and `eps` must be strictly",
    "positive; `gene_batch_size` must be at least 1."
  ),
  extra_ctor = quote({
    if (min_sigma >= max_sigma) {
      stop("`min_sigma` needs to be below `max_sigma`.")
    }
    if (min_phi >= max_phi) {
      stop("`min_phi` needs to be below `max_phi`.")
    }
  }),
  extra_check = quote({
    if (x[["min_sigma"]] >= x[["max_sigma"]]) {
      return("`min_sigma` in NEBULA params is not below `max_sigma`.")
    }
    if (x[["min_phi"]] >= x[["max_phi"]]) {
      return("`min_phi` in NEBULA params is not below `max_phi`.")
    }
  }),
  fields = list(
    nebula_method = p_choice(
      "ln",
      c("ln", "hl"),
      doc = paste(
        "Which variant to run. NEBULA downgrades `\"ln\"` to `\"hl\"`",
        "below 30 cells per subject, as the R package does."
      )
    ),
    min_sigma = p_dbl(
      1e-04,
      range = "(0,)",
      doc = "Lower bound on the subject-level overdispersion."
    ),
    min_phi = p_dbl(
      1e-04,
      range = "(0,)",
      doc = "Lower bound on the cell-level overdispersion."
    ),
    max_sigma = p_dbl(
      10,
      range = "(0,)",
      doc = "Upper bound on the subject-level overdispersion."
    ),
    max_phi = p_dbl(
      1000,
      range = "(0,)",
      doc = "Upper bound on the cell-level overdispersion."
    ),
    cutoff_cell = p_dbl(
      20,
      range = "[0,)",
      doc = paste(
        "Refit both overdispersions when the product of the cells per",
        "subject and the estimated `phi` falls below this."
      )
    ),
    kappa = p_dbl(
      800,
      range = "[0,)",
      doc = paste(
        "Threshold on NEBULA's `kappa_obs` above which the",
        "subject-level overdispersion from stage one is trusted as",
        "is."
      )
    ),
    cpc = p_dbl(
      0.005,
      range = "[0,)",
      doc = "Drop a gene whose mean count per cell is at most this."
    ),
    mincp = p_int(
      5L,
      range = "[0,)",
      doc = "Drop a gene expressed in fewer than this many cells."
    ),
    reml = p_lgl(
      FALSE,
      doc = paste(
        "Estimate the overdispersions by restricted maximum",
        "likelihood. The R package only honours this for `NBLMM`,",
        "which the Rust port does not implement, so this arm has not",
        "been validated against an R reference. Leave it off unless",
        "you know why you want it."
      )
    ),
    eps = p_dbl(
      1e-06,
      range = "(0,)",
      doc = "Absolute stopping tolerance for the optimiser."
    ),
    gene_batch_size = p_int(
      1000L,
      range = "[1,)",
      doc = paste(
        "Genes read and fitted per batch. Bounds how much of the",
        "store is resident at once and changes nothing about the",
        "answer, since NEBULA is gene-independent."
      )
    ),
    shrink_dispersion = p_lgl(
      TRUE,
      doc = paste(
        "Shrink the cell-level overdispersions towards an empirical",
        "Bayes prior once the sweep is done."
      )
    )
  )
)

spec_sc_cellsweep <- param_spec(
  name = "sc_cellsweep",
  title = "Default parameters for CellSweep denoising",
  description = paste(
    "Mirrors the CellSweep reference implementation's defaults.",
    "The pseudocounts (`celltype_lambda`, `ambient_lambda`,",
    "`bulk_lambda`) are given on the scale you see here and",
    "divided by the gene count internally. Two of these are worth",
    "knowing about before you touch anything else.",
    "`freeze_ambient_profile = TRUE` keeps the ambient profile at",
    "its empty-droplet estimate, which is the recommended path",
    "and the only one where `alpha_cap`, the repulsion terms and",
    "cell-type reassignment are live. And `freeze_empties` only",
    "accepts `TRUE`: the reference gives empty droplets a",
    "cell-type component they have no label for, which indexes",
    "past the end of the profile matrix and wraps onto the last",
    "cell type."
  ),
  checker = "ScCellsweep",
  label = "CellSweep params",
  class_tag = "params_sc_cellsweep",
  extra_ctor = quote(
    if (!isTRUE(freeze_empties)) {
      stop(paste(
        "`freeze_empties = FALSE` is not supported. The reference gives",
        "empty droplets a cell-type component they have no label for,",
        "which indexes past the end of the profile matrix."
      ))
    }
  ),
  extra_check = quote(
    if (!isTRUE(x$freeze_empties)) {
      return("`freeze_empties = FALSE` is not supported by CellSweep.")
    }
  ),
  fields = list(
    freeze_empties = p_lgl(
      TRUE,
      doc = paste(
        "Keep the contamination fraction of empty droplets pinned at",
        "1. Only `TRUE` is supported, see the description."
      )
    ),
    freeze_ambient_profile = p_lgl(
      TRUE,
      doc = paste(
        "Keep the ambient profile at its empty-droplet estimate",
        "rather than re-estimating it as a mixture over the cell-type",
        "profiles."
      )
    ),
    init_alpha = p_dbl(
      0.9,
      range = "[0.1,0.9]",
      doc = paste(
        "Starting ambient fraction for every real barcode. With",
        "`freeze_ambient_profile = TRUE` the final result barely",
        "depends on it, so it sits at `alpha_cap`."
      )
    ),
    init_beta = p_dbl(
      0.1,
      range = "[0.1,0.9]",
      doc = paste(
        "Starting bulk contamination fraction. Set below `init_alpha`",
        "on purpose: bulk and ambient are not fully separable, so",
        "this biases unassignable contamination towards ambient."
      )
    ),
    alpha_cap = p_dbl(
      0.9,
      range = "[0,1]",
      doc = paste(
        "Ceiling on the per-cell ambient fraction before the",
        "log-likelihood converges. Barcodes wanting to exceed it are",
        "excluded from the cell-type profile update and allowed to",
        "switch cell type."
      )
    ),
    repulsion_strength = p_dbl(
      1e-04,
      range = "[0,1e-3]",
      doc = paste(
        "Strength of the repulsion pushing cell-type profiles away",
        "from the ambient profile. Scales with cluster mass, so it is",
        "inert on small data and only bites at realistic cell counts."
      )
    ),
    max_frac_gene_repulsion = p_dbl(
      0.2,
      range = "(0,1]",
      doc = paste(
        "Ceiling on the fraction of any single profile entry that",
        "repulsion may remove."
      )
    ),
    celltype_lambda = p_dbl(
      50,
      range = "[0,)",
      doc = paste(
        "Pseudocount smoothing the cell-type profile update. Higher",
        "values give smoother profiles."
      )
    ),
    ambient_lambda = p_dbl(
      50,
      range = "[0,)",
      doc = "Pseudocount smoothing the ambient profile estimate."
    ),
    bulk_lambda = p_dbl(
      10,
      range = "[0,)",
      doc = "Pseudocount smoothing the bulk profile estimate."
    ),
    eps = p_dbl(1e-12, range = "(0,)", doc = "Floor on denominators."),
    log_eps = p_dbl(
      1e-300,
      range = "(0,)",
      doc = "Floor on the argument of `log`."
    ),
    max_iter = p_int(2000L, range = "[2,)", doc = "Hard cap on EM iterations."),
    del0_ll_tol = p_dbl(
      0.001,
      range = "(0,)",
      doc = paste(
        "Log-likelihood change, as a fraction of the first EM step's",
        "change, below which stage one ends and parameter convergence",
        "starts being checked."
      )
    ),
    min_ll_tol = p_dbl(
      1e-06,
      range = "(0,)",
      doc = paste(
        "Floor on the adaptive tolerance, relative to the current",
        "log-likelihood. Stops `del0_ll_tol` chasing floating point",
        "noise."
      )
    ),
    tol_p = p_dbl(
      1e-04,
      range = "(0,)",
      doc = paste(
        "Convergence threshold on the maximum row-wise L1 change in",
        "the cell-type profiles."
      )
    ),
    tol_f = p_dbl(
      1e-04,
      range = "(0,)",
      doc = paste(
        "Convergence threshold on the change in the total",
        "contamination fraction."
      )
    ),
    norm_from_rounded = p_lgl(
      FALSE,
      doc = paste(
        "Derive the normalised layer from the integerised counts",
        "rather than the denoised floats. Consistent across the two",
        "layers at the cost of the sub-integer signal, which is where",
        "CellSweep is most informative."
      )
    ),
    seed = p_int(
      42L,
      range = "[0,)",
      doc = "Seed for the stochastic rounding of the denoised counts."
    )
  )
)

spec_sc_empty_droplets <- param_spec(
  name = "sc_empty_droplets",
  title = "Parameters for identifying empty droplets",
  description = paste(
    "CellSweep trains its ambient profile on the empty droplets,",
    "so it needs them called before it runs. Four ways to do",
    "that, in descending order of how much you should trust them.",
    "`\"supplied\"` takes an existing logical column from the obs",
    "table. This is the common case: if you have CellRanger's",
    "filtered barcode list you already know which barcodes are",
    "empty, and nothing here will beat that. `\"umi_cutoff\"`",
    "calls everything below an absolute library size empty.",
    "`\"expected_cells\"` turns a cell count into that cutoff via",
    "the sorted library sizes. `\"knee\"` finds the cutoff from",
    "the curvature of the rank / log-count curve; it is",
    "experimental in the reference too, and smoothing puts the",
    "curvature minimum a few ranks ahead of the actual cliff, so",
    "it sweeps up the real barcodes nearest the transition."
  ),
  checker = "ScEmptyDroplets",
  label = "empty droplet params",
  class_tag = "params_sc_empty_droplets",
  extra_ctor = quote({
    required <- switch(
      method,
      supplied = "is_empty_column",
      umi_cutoff = "umi_cutoff",
      expected_cells = "expected_cells",
      knee = NULL
    )
    if (!is.null(required) && is.null(get(required))) {
      stop(sprintf("`method = '%s'` needs `%s`.", method, required))
    }
  }),
  extra_check = quote({
    required <- switch(
      x$method,
      supplied = "is_empty_column",
      umi_cutoff = "umi_cutoff",
      expected_cells = "expected_cells",
      knee = NULL
    )
    if (!is.null(required) && is.null(x[[required]])) {
      return(sprintf("`method = '%s'` needs `%s`.", x$method, required))
    }
  }),
  fields = list(
    method = p_choice(
      "supplied",
      c("supplied", "umi_cutoff", "expected_cells", "knee"),
      doc = "How the empty droplets are called, see the description."
    ),
    is_empty_column = p_chr(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "Name of the logical obs column holding the mask. Required",
        "for `method = \"supplied\"`, ignored otherwise."
      )
    ),
    umi_cutoff = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Barcodes with a library size strictly below this are empty.",
        "Required for `method = \"umi_cutoff\"`, ignored otherwise."
      )
    ),
    expected_cells = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Number of real cells expected in the run. Required for",
        "`method = \"expected_cells\"`, ignored otherwise."
      )
    )
  )
)
