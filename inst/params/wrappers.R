spec_ica_general <- param_spec(
  name = "ica_general",
  title = "Wrapper function for standard ICA parameters",
  checker = "Ica",
  label = "ICA params",
  hint = paste(
    "maxit must be an integer; alpha must be in [1, 2]; 0 <",
    "max_tol < 1; verbose must be a boolean."
  ),
  fields = list(
    maxit = p_int(200L, doc = "Maximum number of iterations for ICA."),
    alpha = p_dbl(
      1,
      range = "[1, 2]",
      strict = TRUE,
      doc = paste(
        "The alpha parameter for the logcosh version of ICA. Should",
        "be between 1 to 2."
      )
    ),
    max_tol = p_dbl(
      1e-04,
      range = "(0, 1)",
      strict = TRUE,
      doc = paste(
        "Should be `0 < max_tol < 1`. Maximum tolerance of the",
        "algorithm."
      )
    ),
    verbose = p_lgl(FALSE, doc = "Controls verbosity of the function.")
  )
)

spec_ica_ncomp <- param_spec(
  name = "ica_ncomp",
  title = "Wrapper function for ICA ncomp iterations",
  description = paste(
    "Wrapper function to provide parameters through which ncomps",
    "to iterate through."
  ),
  checker = "IcaNcomps",
  label = "ICA n-components params",
  hint = paste(
    "max_no_comp and steps must be integers; custom_seq must be",
    "NULL or a vector of integers."
  ),
  fields = list(
    max_no_comp = p_int(75L, doc = "Maximum number of ncomp to test."),
    steps = p_int(5L, doc = "In which steps to move from 5 onwards."),
    custom_seq = p_int(
      NULL,
      null_ok = TRUE,
      len = "+",
      doc = paste(
        "If you wish to provide a custom version of no_comp to",
        "iterate through. If NULL, you will iterate through `c(2, 3,",
        "4, 5, 5 + step, ... max_no_comp - step, max_no_comp)`"
      )
    )
  )
)

spec_ica_randomisation <- param_spec(
  name = "ica_randomisation",
  title = "Wrapper function for ICA randomisation",
  checker = "IcaIter",
  label = "ICA randomisation params",
  hint = paste(
    "random_init and folds must be integers; cross_validate must",
    "be a boolean."
  ),
  fields = list(
    cross_validate = p_lgl(
      FALSE,
      doc = paste(
        "Do you want to apply a cross-validation type approach and",
        "split the data into `folds` folds to assess within data",
        "stability of the component."
      )
    ),
    random_init = p_int(50L, doc = "Number of random initialisations to use."),
    folds = p_int(
      10L,
      doc = paste(
        "Number of folds to use if `cross_validate` is set to `TRUE`.",
        "To note, you will be running `random_init * folds` ICA runs."
      )
    )
  )
)

spec_cor_graph <- param_spec(
  name = "cor_graph",
  title = "Wrapper function for graph generation",
  checker = "CorGraph",
  label = "correlation graph params",
  hint = paste(
    "min_cor and fdr_threshold must be in [0, 1]; epsilon must be",
    "a double; verbose must be a boolean."
  ),
  fields = list(
    epsilon = p_dbl(
      2,
      strict = TRUE,
      doc = "Defines the epsilon parameter for the radial basis function."
    ),
    min_cor = p_dbl(
      0.2,
      range = "[0, 1]",
      strict = TRUE,
      doc = paste(
        "Minimum absolute correlation that needs to be observed in",
        "either data set. Only relevant for differential",
        "correlation-based graphs."
      )
    ),
    fdr_threshold = p_dbl(
      0.05,
      range = "[0, 1]",
      strict = TRUE,
      doc = "Maximum FDR for the differential correlation p-value."
    ),
    verbose = p_lgl(
      TRUE,
      doc = "Controls verbosity of the graph generation function."
    )
  )
)

spec_graph_resolution <- param_spec(
  name = "graph_resolution",
  title = paste(
    "Wrapper function to generate resolution parameters for",
    "Leiden or Louvain clustering."
  ),
  checker = "GraphRes",
  label = "resolution params",
  hint = paste(
    "min_res and max_res must be doubles; number_res must be an",
    "integer."
  ),
  fields = list(
    min_res = p_dbl(0.1, strict = TRUE, doc = "Minimum resolution to test."),
    max_res = p_dbl(10, strict = TRUE, doc = "Maximum resolution to test."),
    number_res = p_int(
      15L,
      doc = paste(
        "Number of resolutions to test between the `max_res` and",
        "`min_res.`"
      )
    )
  )
)

spec_community_detection <- param_spec(
  name = "community_detection",
  title = "Wrapper function to generate community detection parameters",
  checker = "Community",
  label = "community params",
  hint = paste(
    "min_nodes, max_nodes and min_seed_nodes must be integers",
    "(with max_nodes >= min_nodes); initial_res must be a double;",
    "network_threshold and pval_threshold must be in (0, 1]."
  ),
  extra_ctor = quote(
    checkmate::qassert(max_nodes, sprintf("I1[%i,)", min_nodes))
  ),
  extra_check = quote(
    if (!checkmate::qtest(x$max_nodes, sprintf("I1[%i,)", x$min_nodes))) {
      return("The element `max_nodes` in community params is below min_nodes.")
    }
  ),
  fields = list(
    max_nodes = p_int(
      300L,
      doc = "Maximum number of nodes in a given community."
    ),
    min_nodes = p_int(
      10L,
      doc = "Minimum number of nodes in a given community."
    ),
    min_seed_nodes = p_int(
      2L,
      doc = "Minimum number of seed nodes within a community."
    ),
    initial_res = p_dbl(
      0.5,
      doc = "Initial resolution parameter to start with."
    ),
    threshold_type = p_choice(
      "prop_based",
      c("prop_based", "pval_based"),
      doc = paste(
        "You can chose to include a certain proportion of the network",
        "with the highest diffusion scores, or use p-values based on",
        "permutations."
      )
    ),
    network_threshold = p_dbl(
      0.5,
      range = "(0, 1]",
      doc = paste(
        "The proportion of the network to include. Used if",
        "`threshold_type = \"prop_based\"`."
      )
    ),
    pval_threshold = p_dbl(
      0.1,
      range = "(0, 1]",
      doc = paste(
        "The maximum p-value for nodes to be included. Used if",
        "`threshold_type = \"pval_based\"`."
      )
    )
  )
)

spec_gsea <- param_spec(
  name = "gsea",
  title = "Wrapper function to generate GSEA parameters",
  checker = "GSEA",
  label = "GSEA params",
  hint = paste(
    "min_size and max_size must be integers (with max_size >",
    "min_size and min_size >= 3); gsea_param must be a double;",
    "sample_size must be an integer; eps must be a float."
  ),
  fields = list(
    min_size = p_int(
      5L,
      range = "[3,)",
      doc = "Minimum number of genes per gene set."
    ),
    max_size = p_int(
      500L,
      range = "[4,)",
      doc = "Maximum number of genes per gene set."
    ),
    gsea_param = p_dbl(1, doc = "GSEA parameter."),
    sample_size = p_int(
      101L,
      doc = paste(
        "Number of samples to iterate through for the multi-level",
        "implementation of fgsea."
      )
    ),
    eps = p_dbl(
      1e-50,
      doc = paste(
        "Boundary for calculating the p-value. Used for the multi-",
        "level implementation of fgsea."
      )
    )
  )
)

spec_blitzgsea <- param_spec(
  name = "blitzgsea",
  title = "Wrapper function to generate blitzGSEA parameters",
  references = "Lachmann, et al., Bioinformatics, 2022",
  checker = "BlitzGsea",
  label = "blitzGSEA params",
  hint = paste(
    "min_size and max_size must be integers (with max_size > min_size and",
    "min_size >= 3); permutations and anchors must be integers >= 2;",
    "symmetric, centre and ks_test must be booleans; seed must be a",
    "non-negative double."
  ),
  # The Rust side takes a 64 bit seed and reads it as a double, since R has no
  # integer type wide enough. An integer would be dropped and the default run.
  extra_ctor = quote({
    if (permutations < BLITZ_MIN_PERMUTATIONS_SPLIT) {
      warning(sprintf(
        paste(
          "%i permutations is below %i, so the positive and negative tails",
          "will be pooled into a single gamma regardless of `symmetric`."
        ),
        permutations,
        BLITZ_MIN_PERMUTATIONS_SPLIT
      ))
    }
    seed <- as.double(seed)
  }),
  fields = list(
    min_size = p_int(
      5L,
      range = "[3,)",
      doc = "Minimum number of genes per gene set."
    ),
    max_size = p_int(
      500L,
      range = "[4,)",
      doc = "Maximum number of genes per gene set."
    ),
    permutations = p_int(
      2000L,
      range = "[2,)",
      doc = paste(
        "Random gene sets drawn per anchor size during calibration.",
        "Below `1000L` the two tails are pooled into a single gamma",
        "regardless of `symmetric`."
      )
    ),
    anchors = p_int(
      40L,
      range = "[2,)",
      doc = paste(
        "Number of log-spaced anchor sizes requested. Sizes that",
        "collide after rounding are collapsed, so the realised grid",
        "is usually a little smaller."
      )
    ),
    symmetric = p_lgl(
      FALSE,
      doc = paste(
        "Pool both tails into one gamma instead of fitting them",
        "separately."
      )
    ),
    centre = p_lgl(
      TRUE,
      doc = paste(
        "Centre the signature on its mean before scoring. The",
        "enrichment score is not invariant to an offset, so the",
        "calibration and the scoring have to agree on this."
      )
    ),
    ks_test = p_lgl(
      TRUE,
      doc = paste(
        "Run the Kolmogorov-Smirnov goodness-of-fit diagnostic at",
        "every anchor. Costs a sort per anchor."
      )
    ),
    seed = p_dbl(42, range = "[0,)", doc = "Random seed for the calibration.")
  )
)

spec_gsva <- param_spec(
  name = "gsva",
  title = "Wrapper function to generate GSVA parameters",
  checker = "GSVA",
  label = "GSVA params",
  hint = paste(
    "min_size and max_size must be integers (with max_size >",
    "min_size and min_size >= 3); tau must be a double; max_diff",
    "and abs_rank must be booleans."
  ),
  fields = list(
    tau = p_dbl(
      1,
      doc = paste(
        "Tau parameter, usual recommendation is to use `1.0` here.",
        "Larger values emphasise the tails more."
      )
    ),
    min_size = p_int(
      5L,
      range = "[3,)",
      doc = "Minimum number of genes per gene set."
    ),
    max_size = p_int(
      500L,
      range = "[4,)",
      doc = "Maximum number of genes per gene set."
    ),
    max_diff = p_lgl(
      TRUE,
      doc = paste(
        "Scoring mode for GSVA, if `TRUE` = difference; if `FALSE` =",
        "larger absolute value."
      )
    ),
    abs_rank = p_lgl(FALSE, doc = "If `TRUE` = pos - neg, `FALSE` = pos + neg.")
  )
)

spec_ssgsea <- param_spec(
  name = "ssgsea",
  title = "Wrapper function to generate ssGSEA parameters",
  checker = "SingleSampleGSEA",
  label = "ssGSEA params",
  hint = paste(
    "min_size and max_size must be integers (with max_size > min_size and",
    "min_size >= 3); alpha must be a double in (0, 1); normalise must be a",
    "boolean."
  ),
  fields = list(
    alpha = p_dbl(
      0.25,
      range = "(0, 1)",
      doc = paste(
        "The exponent defining the weight of the tail in the random",
        "walk performed by ssGSEA."
      )
    ),
    min_size = p_int(
      5L,
      range = "[3,)",
      doc = "Minimum number of genes per gene set."
    ),
    max_size = p_int(
      500L,
      range = "[4,)",
      doc = "Maximum number of genes per gene set."
    ),
    normalise = p_lgl(TRUE, doc = "Shall the scores be normalised.")
  )
)

spec_coremo <- param_spec(
  name = "coremo",
  title = "Wrapper function to generate CoReMo parameters",
  checker = "CoReMo",
  label = "CoReMo params",
  hint = paste(
    "k_min and k_max must be integers; min_size must be an",
    "integer or NULL; junk_module_threshold and epsilon must be",
    "floats."
  ),
  fields = list(
    epsilon = p_dbl(
      2,
      doc = paste(
        "Epsilon parameter for the chosen RBF function, see",
        "`rbf_func`. The higher, the more aggressively low",
        "correlations will be shrunk."
      )
    ),
    k_min = p_int(
      2L,
      doc = paste(
        "Minimum and maximum number of cuts to use for the",
        "hierarchical clustering."
      )
    ),
    k_max = p_int(
      150L,
      doc = paste(
        "Minimum and maximum number of cuts to use for the",
        "hierarchical clustering."
      )
    ),
    min_size = p_int(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "Minimum size of the clusters. Smaller clusters will be",
        "combined together."
      )
    ),
    junk_module_threshold = p_dbl(
      0.05,
      doc = paste(
        "Threshold for the minimum correlation to be observed in a",
        "module."
      )
    ),
    rbf_func = p_choice(
      "gaussian",
      c("gaussian", "inverse_quadratic", "bump"),
      doc = paste(
        "Type of RBF you wish to apply to down-weigh weak",
        "correlations."
      )
    ),
    cor_method = p_choice(
      "spearman",
      c("spearman", "pearson"),
      doc = "The type of correlation to use."
    )
  )
)

spec_dgrdl <- param_spec(
  name = "dgrdl",
  title = "Wrapper function to generate DGRDL parameters",
  checker = "DGRDL",
  label = "DGRDL params",
  hint = paste(
    "sparsity, dict_size, max_iter, k_neighbours and admm_iter",
    "must be integers; alpha, beta and rho must be floats."
  ),
  fields = list(
    sparsity = p_int(
      5L,
      doc = "Sparsity constraint (max non-zero coefficients per signal)"
    ),
    dict_size = p_int(5L, doc = "Dictionary size"),
    alpha = p_dbl(1, doc = "Sample context regularisation weight."),
    beta = p_dbl(1, doc = "Feature effect regularisation weight."),
    max_iter = p_int(
      20L,
      doc = "Maximum number of iterations for the main algorithm."
    ),
    k_neighbours = p_int(5L, doc = "Number of neighbours in the KNN graph."),
    admm_iter = p_int(5L, doc = "ADMM iterations for sparse coding."),
    rho = p_dbl(1, doc = "ADMM step size.")
  )
)

spec_nmf_hals <- param_spec(
  name = "nmf_hals",
  title = "Wrapper function for NMF (HALS) parameters",
  checker = "NmfHals",
  label = "NMF HALS params",
  hint = paste(
    "max_iter and check_every must be positive integers; tol and",
    "eps must be positive numerics."
  ),
  fields = list(
    max_iter = p_int(
      250L,
      range = "[1,)",
      doc = "Maximum number of HALS iterations."
    ),
    tol = p_dbl(
      1e-04,
      range = "(0,)",
      doc = paste(
        "Convergence tolerance on the relative change in",
        "reconstruction loss."
      )
    ),
    eps = p_dbl(
      1e-10,
      range = "(0,)",
      doc = "Numerical floor for non-negativity / division safety."
    ),
    check_every = p_int(
      10L,
      range = "[1,)",
      doc = "Convergence check interval in iterations."
    ),
    nmf_init = p_choice(
      "nndsvd",
      c("nndsvd", "svd", "random"),
      doc = paste(
        "`\"nndsvd\"` and `\"svd\"` both map to deterministic NNDSVD",
        "initialisation; `\"random\"` uses random non-negative draws.",
        "For stabilised (multi-run) NMF this field is ignored and",
        "random init is always used."
      )
    )
  )
)

spec_nmf_consensus <- param_spec(
  name = "nmf_consensus",
  title = "Wrapper function for consensus NMF parameters",
  description = paste(
    "Controls the clustering step of consensus NMF: the",
    "components of every restart are pooled, outliers are dropped",
    "by local density, and the survivors are k-means clustered",
    "into `k` groups whose median becomes the consensus factor."
  ),
  references = "Kotliar et al., eLife, 2019",
  checker = "NmfConsensus",
  label = "NMF consensus params",
  hint = paste(
    "n_neighbours must be a non-negative integer (0 = auto);",
    "density_threshold must be a numeric in [0, 2] (>= 2 disables",
    "the filter); kmeans_iters and kmeans_n_init must be positive",
    "integers."
  ),
  fields = list(
    consensus_target = p_choice(
      "h",
      c("h", "w"),
      doc = paste(
        "`\"h\"` clusters the gene programmes (the spectra), which is",
        "what cNMF does and what you almost always want. `\"w\"`",
        "clusters in sample space instead, which on single cell data",
        "means cell space: it pools a dense `(k * n_runs) x n_cells`",
        "matrix and runs an exhaustive cosine search over it, so it",
        "gets expensive fast."
      )
    ),
    n_neighbours = p_int(
      0L,
      range = "[0,)",
      doc = paste(
        "Neighbours used for the local density estimate. `0L` picks",
        "`ceiling(0.3 * n_runs)` for you."
      )
    ),
    density_threshold = p_dbl(
      0.5,
      range = "[0,2]",
      doc = paste(
        "Components whose mean cosine distance to their neighbours",
        "exceeds this are dropped as unstable. Cosine distance cannot",
        "exceed 2, so any value `>= 2` disables the filter entirely."
      )
    ),
    kmeans_iters = p_int(
      100L,
      range = "[1,)",
      doc = "Maximum k-means iterations."
    ),
    kmeans_n_init = p_int(
      3L,
      range = "[1,)",
      doc = "Number of k-means restarts."
    )
  )
)

spec_snf <- param_spec(
  name = "snf",
  title = "Wrapper function to generate SNF parameters",
  return_order = c("k", "t", "mu", "alpha", "distance_metric", "normalise"),
  checker = "SNF",
  test_fn = TRUE,
  label = "SNF params",
  hint = paste(
    "k and t must be positive integers; mu must be a float in [0,",
    "1]; alpha must be a float; normalise must be a boolean."
  ),
  fields = list(
    k = p_int(20L, doc = "Number of neighbours to consider."),
    t = p_int(20L, doc = "Number of iterations for the SNF algorithm."),
    mu = p_dbl(
      0.5,
      range = "[0,1]",
      doc = "Normalisation factor for the Gaussian kernel width."
    ),
    alpha = p_dbl(
      1,
      doc = "Normalisation parameter controlling the fusion strength."
    ),
    normalise = p_lgl(TRUE, doc = "Shall continuous values be Z-scored."),
    distance_metric = p_choice(
      "euclidean",
      c("euclidean", "manhattan", "canberra", "cosine"),
      doc = paste(
        "Which distance metric to use for the continuous",
        "calculations. In case of pure categorical, Hamming will be",
        "used, for mixed data types Gower distance is used."
      )
    )
  )
)

spec_cistarget <- param_spec(
  name = "cistarget",
  title = "Wrapper function to CisTarget parameters",
  description = paste(
    "`auc_threshold` and `max_rank` are two different cutoffs and",
    "are easy to confuse. The first truncates the recovery curve",
    "used to score motif enrichment, the second sets how deep",
    "into the ranking the background curve (mean + 2 SD across",
    "all motifs) is built, and therefore where the leading edge",
    "cuts. RcisTarget uses `maxRank = 5000` and `nMean = 100`;",
    "the defaults here follow it. `auc_threshold` follows",
    "pySCENIC at 5%, RcisTarget itself uses 3%."
  ),
  checker = "Cistarget",
  label = "CisTarget params",
  hint = paste(
    "auc_threshold must be numeric [0, 1]; nes_threshold must be",
    "numeric; max_rank and n_mean must be positive integers;",
    "rcc_method must be a single string; high_conf_cats and",
    "low_conf_cats must be character vectors."
  ),
  fields = list(
    auc_threshold = p_dbl(
      0.05,
      range = "[0,1]",
      doc = paste(
        "Numeric between 0 and 1. Proportion of genes to use for AUC",
        "threshold calculation. Default is 0.05 (5% of genes)."
      )
    ),
    nes_threshold = p_dbl(
      3,
      doc = paste(
        "Normalised Enrichment Score threshold for significant",
        "motifs. Default is 3.0."
      )
    ),
    max_rank = p_int(
      5000L,
      range = "[1,)",
      doc = paste(
        "Depth of the recovery curves used to derive the background",
        "and the leading edge. Clamped to the number of genes in the",
        "ranking database. Default is 5000, the RcisTarget value."
      )
    ),
    n_mean = p_int(
      100L,
      range = "[1,)",
      doc = paste(
        "Window for the rolling mean smoothing the background",
        "recovery curve. Only read when `rcc_method = \"approx\"`.",
        "Default is 100, the RcisTarget value."
      )
    ),
    rcc_method = p_choice(
      "approx",
      c("approx", "icistarget"),
      doc = paste(
        "Method for recovery curve calculation. Either \"approx\"",
        "(approximate, faster) or \"icistarget\" (exact, slower)."
      )
    ),
    high_conf_cats = p_chr(
      c("directAnnotation", "inferredBy_Orthology"),
      len = "+",
      doc = paste(
        "Annotation categories considered high confidence. Default",
        "includes direct annotations and orthology-based inferences."
      )
    ),
    low_conf_cats = p_chr(
      c("inferredBy_MotifSimilarity", "inferredBy_MotifSimilarity_n_Orthology"),
      len = "+",
      doc = paste(
        "Annotation categories considered lower confidence. Default",
        "includes motif similarity-based inferences."
      )
    )
  )
)

spec_label_propagation <- param_spec(
  name = "label_propagation",
  title = "Wrapper function to generate label propagation parameters",
  checker = "LabelProp",
  label = "label propagation params",
  hint = paste(
    "alpha must be in [0, 1]; iter must be a positive integer; tolerance",
    "must be a double; symmetrise must be a boolean."
  ),
  fields = list(
    alpha = p_dbl(
      0.9,
      range = "[0, 1]",
      strict = TRUE,
      doc = paste(
        "Controls the spreading strength. Higher values anchor",
        "labelled nodes more strongly to their original label."
      )
    ),
    iter = p_int(
      100L,
      range = "[1,]",
      doc = "Maximum number of iterations to run."
    ),
    tolerance = p_dbl(
      1e-06,
      strict = TRUE,
      doc = paste(
        "Convergence threshold. Stops early if the maximum change",
        "across all nodes falls below this value."
      )
    ),
    symmetrise = p_lgl(TRUE, doc = "Whether to symmetrise the graph."),
    symmetry_strategy = p_choice(
      "average",
      c("average", "avg", "min", "max"),
      doc = paste(
        "Strategy to resolve weight conflicts when symmetrising. Only",
        "relevant when `symmetrise = TRUE` and edge weights are",
        "provided."
      )
    ),
    max_hops = p_int(
      NULL,
      range = "[0,]",
      null_ok = TRUE,
      doc = paste(
        "If provided, restricts label spreading to nodes within this",
        "many hops of any labelled node. Nodes beyond this limit",
        "remain all-zeroes."
      )
    )
  )
)

spec_module_membership <- param_spec(
  name = "module_membership",
  title = "Wrapper function to generate module membership parameters",
  description = paste(
    "Controls how a `gene x k` loading matrix from ICA, NMF or",
    "DGRDL is turned into module membership. Genes are kept where",
    "they sit in the tail of a component's loading distribution,",
    "which means membership is not exclusive: a gene loading",
    "strongly on three components belongs to three modules. Genes",
    "in no tail belong to nothing, which is the background",
    "category an argmax assignment cannot give you."
  ),
  references = "Biton, et al., Cell Rep, 2014",
  checker = "ModuleMembership",
  label = "module membership params",
  hint = "cutoff must be a positive float; fdr must sit in (0, 1].",
  fields = list(
    method = p_choice(
      "zscore",
      c("zscore", "fdr"),
      doc = paste(
        "`\"zscore\"` standardises each component and keeps `abs(z) >",
        "cutoff`. `\"fdr\"` converts to two-sided p-values against a",
        "Normal null fitted the same way, Benjamini-Hochberg adjusts,",
        "and keeps `padj < fdr`."
      )
    ),
    cutoff = p_dbl(
      3,
      range = "(0,)",
      doc = "Absolute z threshold for `method = \"zscore\"`."
    ),
    fdr = p_dbl(
      0.05,
      range = "(0,1]",
      doc = "Adjusted p-value threshold for `method = \"fdr\"`."
    ),
    tails = p_choice(
      "auto",
      c("auto", "upper", "both"),
      doc = paste(
        "`\"auto\"` uses an upper-tail-only test when every loading",
        "is non-negative (the NMF case) and a two-sided one",
        "otherwise. `\"upper\"` and `\"both\"` force the choice."
      )
    ),
    scaling = p_choice(
      "robust",
      c("robust", "standard"),
      doc = paste(
        "`\"robust\"` centres and scales each component by its median",
        "and MAD. `\"standard\"` uses the mean and standard deviation",
        "instead, which is stricter and less forgiving of skewed",
        "loadings (e.g. NMF)."
      )
    )
  )
)

spec_synthetic_bulk_rnaseq <- param_spec(
  name = "synthetic_bulk_rnaseq",
  title = paste(
    "Wrapper function to generate synthetic bulk RNAseq",
    "parameters"
  ),
  description = paste(
    "Parameters for [bixverse::synthetic_bulk_cor_matrix()].",
    "Counts come from a negative binomial with a mean-dispersion",
    "trend; co-expression modules are planted by putting each",
    "module's genes on a shared latent factor. The `generator`",
    "picks how loadings and factors are drawn, which is what",
    "makes a given dataset a fair or unfair benchmark for a given",
    "method: \\itemize{ \\item `\"hub_modular\"` - LogNormal",
    "loadings on a Normal factor. Some genes end up far more",
    "connected than others, so this is the WGCNA-style default.",
    "\\item `\"modular\"` - Beta(5, 2) loadings on a Normal",
    "factor. Homogeneous within-module correlation and no hubs.",
    "\\item `\"non_negative_factor\"` - LogNormal loadings on a",
    "Gamma factor. The activity matrix is non-negative by",
    "construction, so NMF has a ground truth it can actually",
    "reach. \\item `\"non_gaussian_factor\"` - LogNormal loadings",
    "on a Laplace factor. Non-Gaussian sources satisfy ICA",
    "identifiability. }"
  ),
  details = "`noise_std` and `factor_std` default to `0.1` and `0.5`, not to the\n`bixverse-rs` values of `0.3` and `0.3`. At the crate defaults the\n`\"modular\"` generator plants modules too weakly to detect at 1000 genes by\n100 samples: the within-module minus cross-module mean absolute Spearman gap\ncomes out around `0.06`, against `0.17` to `0.23` for the other three. The\nvalues here put all four generators in the `0.30` to `0.39` band, so a\ncomparison across generators reflects the method rather than the signal\nstrength it happened to be handed. Pass the crate values explicitly if you\nwant a harder problem.",
  references = "Zhang & Horvath, Stat Appl Genet Mol Biol, 2005",
  checker = "SyntheticBulk",
  extra_ctor = quote(checkmate::assertTRUE(sum(module_sizes) <= num_genes)),
  extra_check = quote(
    if (sum(x[["module_sizes"]]) > x[["num_genes"]]) {
      return("The sum of `module_sizes` must not exceed `num_genes`.")
    }
  ),
  label = "synthetic bulk params",
  hint = paste(
    "num_samples, num_genes and seed must be integers;",
    "module_sizes must be an integer vector (a double vector is",
    "silently ignored downstream); hub_percentile must sit in (0,",
    "1]; the remaining distribution parameters must be positive",
    "floats."
  ),
  fields = list(
    num_samples = p_int(
      100L,
      range = "[1,)",
      doc = "Number of samples (columns) to simulate."
    ),
    num_genes = p_int(
      1000L,
      range = "[1,)",
      doc = "Number of genes (rows) to simulate."
    ),
    module_sizes = p_int(
      c(100L, 100L, 100L),
      len = "*",
      doc = paste(
        "Sizes of the co-expression modules. The sum must be smaller",
        "or equal to `num_genes`. Genes are assigned in contiguous",
        "blocks from the first gene onwards, any remainder is",
        "background. Use `integer(0)` for no modules. Must be an",
        "integer vector, see the note below."
      )
    ),
    generator = p_choice(
      "hub_modular",
      c("hub_modular", "modular", "non_negative_factor", "non_gaussian_factor"),
      doc = paste(
        "Which topology and distribution family to plant., see the",
        "description."
      )
    ),
    seed = p_int(
      123L,
      range = "[0,)",
      doc = "Seed for reproducibility purposes."
    ),
    mean_exp_gamma_shape = p_dbl(
      5,
      range = "(0,)",
      doc = paste(
        "Shape and scale of the Gamma the per-gene mean expression is",
        "drawn from."
      )
    ),
    mean_exp_gamma_scale = p_dbl(
      10,
      range = "(0,)",
      doc = paste(
        "Shape and scale of the Gamma the per-gene mean expression is",
        "drawn from."
      )
    ),
    disp_intercept = p_dbl(
      0.2,
      range = "(0,)",
      doc = paste(
        "Intercept and slope of the negative binomial dispersion",
        "trend `disp = 1 / (a + b * mean)`. This is what gives you",
        "heteroskedasticity: lowly expressed genes show higher",
        "variance."
      )
    ),
    disp_slope = p_dbl(
      0.3,
      range = "(0,)",
      doc = paste(
        "Intercept and slope of the negative binomial dispersion",
        "trend `disp = 1 / (a + b * mean)`. This is what gives you",
        "heteroskedasticity: lowly expressed genes show higher",
        "variance."
      )
    ),
    noise_std = p_dbl(
      0.1,
      range = "[0,)",
      doc = paste(
        "Per-gene per-sample noise standard deviation on the latent",
        "log-signal. Smaller values track the module factor more",
        "tightly and give stronger within-module correlation."
      )
    ),
    factor_std = p_dbl(
      0.5,
      range = "(0,)",
      doc = paste(
        "Standard deviation of the Normal factor. Only used by",
        "`\"hub_modular\"` and `\"modular\"`; the other two",
        "generators draw their factor from",
        "`factor_shape`/`factor_scale` instead."
      )
    ),
    factor_shape = p_dbl(
      2,
      range = "(0,)",
      doc = paste(
        "Shape and scale of the Gamma factor for",
        "`\"non_negative_factor\"`. `factor_scale` doubles as the",
        "Laplace scale for `\"non_gaussian_factor\"`."
      )
    ),
    factor_scale = p_dbl(
      0.3,
      range = "(0,)",
      doc = paste(
        "Shape and scale of the Gamma factor for",
        "`\"non_negative_factor\"`. `factor_scale` doubles as the",
        "Laplace scale for `\"non_gaussian_factor\"`."
      )
    ),
    loading_mu = p_dbl(
      0,
      doc = paste(
        "Location and scale of the LogNormal the loadings are drawn",
        "from. Unused by `\"modular\"`, which draws Beta(5, 2)."
      )
    ),
    loading_sigma = p_dbl(
      0.7,
      range = "(0,)",
      doc = paste(
        "Location and scale of the LogNormal the loadings are drawn",
        "from. Unused by `\"modular\"`, which draws Beta(5, 2)."
      )
    ),
    hub_percentile = p_dbl(
      0.1,
      range = "(0,1]",
      doc = paste(
        "Top fraction of module genes flagged as hubs by loading",
        "rank. Must be in `(0, 1]`."
      )
    )
  )
)

spec_bulk_sparsity <- param_spec(
  name = "bulk_sparsity",
  title = "Wrapper function to generate bulk sparsification parameters",
  description = paste(
    "Parameters for [bixverse::simulate_dropouts()]. Dropout",
    "falls out of the library size rather than an explicit",
    "per-gene dropout curve: a size factor `s_j ~ LogNormal(0,",
    "capture_efficiency_sigma)` is drawn per sample, giving a",
    "target library size of `target_library_size * s_j`, and each",
    "gene is binomially thinned towards that target."
  ),
  references = "Zappia, et al., Genome Biol, 2017",
  checker = "BulkSparsity",
  label = "bulk sparsity params",
  hint = paste(
    "target_library_size and capture_efficiency_sigma must be",
    "positive floats; seed must be a positive integer."
  ),
  fields = list(
    strategy = p_choice(
      "seq_depth",
      "seq_depth",
      doc = paste(
        "Which dropout strategy to apply. Currently only",
        "`\"seq_depth\"`."
      )
    ),
    target_library_size = p_dbl(
      20000,
      range = "(0,)",
      doc = "Reference library size per sample."
    ),
    capture_efficiency_sigma = p_dbl(
      0.5,
      range = "(0,)",
      doc = paste(
        "Standard deviation of the LogNormal size-factor",
        "distribution. Larger values spread the library sizes further",
        "apart."
      )
    ),
    seed = p_int(
      123L,
      range = "[0,)",
      doc = "Seed for reproducibility purposes."
    )
  )
)

spec_sc_synthetic_data <- param_spec(
  name = "sc_synthetic_data",
  title = paste(
    "Default parameters for generation of synthetic single cell",
    "data (RNA)"
  ),
  description = paste(
    "For the generation of synthetic single cell data mostly for",
    "testing or showcasing purposes. The default configurations",
    "generates 1000 cells x 100 genes with genes 1:10 being cell",
    "markers for cell type 1, genes 11:20 for cell type 2 and",
    "genes 21:30 for cell type."
  ),
  return_order = c(
    "n_cells",
    "n_genes",
    "marker_genes",
    "n_batches",
    "batch_effect_strength",
    "n_samples",
    "sample_bias"
  ),
  checker = "ScSyntheticData",
  label = "synthetic data params",
  hint = paste(
    "n_cells, n_genes and n_batches must be integers; n_samples",
    "must be an integer or NULL."
  ),
  extra_ctor = quote({
    checkmate::assertList(marker_genes, types = "list", names = "named")
    checkmate::assert(
      checkmate::testNull(sample_bias),
      checkmate::testChoice(
        sample_bias,
        c("even", "slightly_uneven", "very_uneven")
      )
    )
    if (is.null(n_samples) != is.null(sample_bias)) {
      stop(paste(
        "`n_samples` and `sample_bias` must be provided together.",
        "Supply both to add sample membership, or neither to omit it."
      ))
    }
  }),
  extra_check = quote({
    res <- checkmate::checkList(
      x$marker_genes,
      types = "list",
      names = "named"
    )
    if (!isTRUE(res)) {
      return(res)
    }
    if (!is.null(x[["sample_bias"]])) {
      res <- apply_choice_rules(
        x,
        list(sample_bias = c("even", "slightly_uneven", "very_uneven")),
        label = "synthetic data params"
      )
      if (!isTRUE(res)) {
        return(res)
      }
    }
  }),
  fields = list(
    n_cells = p_int(1000L, doc = "Number of cells."),
    n_genes = p_int(100L, doc = "Number of genes."),
    n_batches = p_int(1L, doc = "Number of batches."),
    marker_genes = p_free(
      list(
        cell_type_1 = list(marker_genes = 0:9L),
        cell_type_2 = list(marker_genes = 10:19L),
        cell_type_3 = list(marker_genes = 20:29L)
      ),
      doc = paste(
        "A nested list that indicates which gene indices are markers",
        "for which cell."
      )
    ),
    batch_effect_strength = p_choice(
      "strong",
      c("strong", "medium", "weak"),
      doc = "The strength of the batch effect to add."
    ),
    n_samples = p_int(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "Shall sample membership be added to the synthetic data. If",
        "you want sample information you need to provide `n_samples`",
        "and `sample_bias`."
      )
    ),
    sample_bias = p_free(
      NULL,
      doc = "One of `c(\"even\", \"slightly_uneven\", \"very_uneven\")`"
    )
  )
)

spec_sc_synthetic_data_adt <- param_spec(
  name = "sc_synthetic_data_adt",
  title = paste(
    "Default parameters for generation of synthetic single cell",
    "data (ADT)"
  ),
  description = paste(
    "For the generation of synthetic single cell data mostly for",
    "testing or showcasing purposes. In this case, ADT counts to",
    "test multi-modal integration. The default configurations",
    "generates 1000 cells x 15 proteins with probes 1:3 being",
    "cell markers for cell type 1, genes 4:6 for cell type 2 and",
    "genes 7:9 for cell type. Columns 13:15 represents isotype",
    "controls."
  ),
  return_order = c(
    "n_cells",
    "n_proteins",
    "marker_genes",
    "n_batches",
    "isotype_controls",
    "batch_effect_strength"
  ),
  checker = "ScSyntheticDataAdt",
  label = "synthetic ADT params",
  hint = paste(
    "n_cells, n_proteins and n_batches must be integers;",
    "isotype_controls must be a non-empty integer vector."
  ),
  extra_ctor = quote(
    checkmate::assertList(marker_genes, types = "list", names = "named")
  ),
  extra_check = quote({
    res <- checkmate::checkList(
      x$marker_genes,
      types = "list",
      names = "named"
    )
    if (!isTRUE(res)) {
      return(res)
    }
  }),
  fields = list(
    n_cells = p_int(1000L, doc = "Number of cells."),
    n_proteins = p_int(15L, doc = "Number of proteins"),
    n_batches = p_int(1L, doc = "Number of batches."),
    marker_genes = p_free(
      list(
        cell_type_1 = list(marker_genes = 0:2L),
        cell_type_2 = list(marker_genes = 3:5L),
        cell_type_3 = list(marker_genes = 6:8L)
      ),
      doc = paste(
        "A nested list that indicates which gene indices are markers",
        "for which cell."
      )
    ),
    isotype_controls = p_int(
      12L:14L,
      len = "+",
      doc = "The columns that defines the isotype controls. (0-indexed!)"
    ),
    batch_effect_strength = p_choice(
      "strong",
      c("strong", "medium", "weak"),
      doc = "The strength of the batch effect to add."
    )
  )
)

spec_sc_synthetic_dialogue <- param_spec(
  name = "sc_synthetic_dialogue",
  title = "Default parameters for generation of synthetic DIALOGUE data",
  description = paste(
    "Shapes the fixture with a planted multicellular programme",
    "that [bixverse::generate_dialogue_test_data()] builds. The",
    "defaults give 14 samples x 25 cells x 3 cell types = 1050",
    "cells over 400 genes."
  ),
  details = "`n_genes` defaults higher here than on the Rust side, which plants the same\nstructure into 90. R runs the counts through the normal ingestion path, and\non a small panel the planted genes are a large enough share of the library\nthat log-normalisation divides the signal back out: at 90 genes the library\nsize tracks the programme at an r of 0.75 and background genes pick up a\nspurious correlation of their own. At 400 the planted block is a few percent\nof the library and the contrast survives. The Rust tests feed the planted\nlayer straight in, so they never meet this.",
  references = "Jerby-Arnon & Regev, Nature Biotechnology, 2022",
  checker = "ScSyntheticDialogue",
  extra_ctor = quote({
    checkmate::assertTRUE(n_sample_features <= n_features)
    checkmate::assertTRUE(n_planted * n_cell_types <= n_genes)
  }),
  extra_check = quote({
    if (x$n_sample_features > x$n_features) {
      return("n_sample_features cannot exceed n_features.")
    }
    if (x$n_planted * x$n_cell_types > x$n_genes) {
      return("The planted gene blocks do not fit into n_genes.")
    }
  }),
  label = "synthetic DIALOGUE params",
  hint = paste(
    "n_cell_types and n_features must be at least 2, the",
    "remaining counts positive integers."
  ),
  fields = list(
    n_samples = p_int(
      14L,
      range = "[1,)",
      doc = "Samples the experiment spans. DIALOGUE needs at least 5."
    ),
    cells_per_sample = p_int(
      25L,
      range = "[1,)",
      doc = "Cells per sample per cell type."
    ),
    n_cell_types = p_int(
      3L,
      range = "[2,)",
      doc = "Number of cell types. Must be at least 2."
    ),
    n_features = p_int(
      8L,
      range = "[2,)",
      doc = "Feature columns per cell type. Must be at least 2."
    ),
    n_sample_features = p_int(
      5L,
      range = "[1,)",
      doc = paste(
        "Feature columns carrying a per-sample component. The first",
        "of those is the shared programme, the rest are",
        "cell-type-specific nuisance; anything past this count is",
        "pure noise and exists so the ANOVA filter has something to",
        "reject."
      )
    ),
    n_genes = p_int(400L, range = "[1,)", doc = "Number of genes."),
    n_planted = p_int(
      8L,
      range = "[0,)",
      doc = paste(
        "Planted genes per cell type. The blocks are contiguous, so",
        "`n_planted * n_cell_types` has to fit into `n_genes`."
      )
    )
  )
)

spec_sc_synthetic_cellsweep <- param_spec(
  name = "sc_synthetic_cellsweep",
  title = paste(
    "Default parameters for generation of synthetic CellSweep",
    "data"
  ),
  description = paste(
    "Shapes the fixture with a planted ambient profile that",
    "[bixverse::generate_cellsweep_test_data()] builds. The",
    "defaults give 600 real barcodes over 3 cell types plus 2000",
    "empty droplets, on 200 genes."
  ),
  details = "The soup is the first cell type plus flat background rather than a mixture\nof every cell type profile. A soup sitting in the span of the cell type\nprofiles makes the contamination fraction unidentifiable, and the fixture\nwould then be testing the repulsion term rather than the model.",
  checker = "ScSyntheticCellsweep",
  extra_ctor = quote(checkmate::assertTRUE(n_markers * n_celltypes <= n_genes)),
  extra_check = quote(
    if (x$n_markers * x$n_celltypes > x$n_genes) {
      return("The marker gene blocks do not fit into n_genes.")
    }
  ),
  label = "synthetic CellSweep params",
  hint = paste(
    "n_empty must be at least 30, marker_weight above 1, and the",
    "fractions within their respective ranges."
  ),
  fields = list(
    n_real = p_int(
      600L,
      range = "[1,)",
      doc = paste(
        "Number of real barcodes. Cell types are assigned round-robin",
        "over them."
      )
    ),
    n_empty = p_int(
      2000L,
      range = "[30,)",
      doc = paste(
        "Number of empty droplets. The ambient profile is estimated",
        "off these, so at least 30 and preferably a lot more."
      )
    ),
    n_genes = p_int(200L, range = "[1,)", doc = "Number of genes."),
    n_celltypes = p_int(3L, range = "[1,)", doc = "Number of cell types."),
    n_markers = p_int(
      20L,
      range = "[1,)",
      doc = paste(
        "Width of each cell type's marker block. The blocks are",
        "contiguous and disjoint, so `n_markers * n_celltypes` has to",
        "fit into `n_genes`."
      )
    ),
    marker_weight = p_dbl(
      25,
      range = "(1,)",
      doc = paste(
        "Enrichment of a marker gene over background in its own cell",
        "type's profile. Must exceed 1."
      )
    ),
    ambient_dominance = p_dbl(
      0.6,
      range = "[0,1]",
      doc = paste(
        "Fraction of the soup coming from the first cell type. The",
        "remainder is flat background."
      )
    ),
    alpha_mean = p_dbl(
      0.3,
      range = "[0.02,0.85]",
      doc = "Mean planted ambient fraction across real barcodes."
    ),
    alpha_sd = p_dbl(
      0.12,
      range = "[0,)",
      doc = "Spread of the planted ambient fraction."
    ),
    real_lib_size = p_int(
      3000L,
      range = "[1,)",
      doc = "Expected library size of a real barcode."
    ),
    empty_lib_size = p_int(
      120L,
      range = "[1,)",
      doc = "Expected library size of an empty droplet."
    )
  )
)

spec_sc_mtx_io <- param_spec(
  name = "sc_mtx_io",
  title = "Wrapper function to provide data for mtx-based loading",
  checker = "ScMtxIO",
  label = "MTX IO params",
  hint = "cells_as_rows and has_hdr must be booleans.",
  extra_ctor = quote({
    checkmate::assertFileExists(path_mtx)
    checkmate::assertFileExists(path_obs)
    checkmate::assertFileExists(path_var)
    path_mtx <- path.expand(path_mtx)
    path_obs <- path.expand(path_obs)
    path_var <- path.expand(path_var)
  }),
  extra_check = quote({
    files_ok <- purrr::map_lgl(
      c("path_mtx", "path_obs", "path_var"),
      \(n) checkmate::testFileExists(x[[n]])
    )
    if (!all(files_ok)) {
      return(paste(
        "Some of the files specified in the config for mtx ingest do not",
        "exist. Please check the provided params."
      ))
    }
  }),
  fields = list(
    path_mtx = p_free(doc = "Path to the .mtx file"),
    path_obs = p_free(doc = "Path to the file containing cell/barcode info."),
    path_var = p_free(doc = "Path to the file containing gene/variable info."),
    cells_as_rows = p_lgl(doc = "Do cells represent the rows or columns."),
    has_hdr = p_lgl(doc = "Do the plain text files have headers.")
  )
)

spec_sc_min_quality <- param_spec(
  name = "sc_min_quality",
  title = paste(
    "Wrapper function to generate QC metric params for single",
    "cell"
  ),
  checker = "ScMinQC",
  label = "single cell QC params",
  hint = paste(
    "min_unique_genes, min_lib_size and min_cells must be",
    "integers; target_size must be a float."
  ),
  fields = list(
    min_unique_genes = p_int(
      100L,
      doc = "Minimum number of unique genes per cell/spot to be included."
    ),
    min_lib_size = p_int(
      250L,
      doc = "Minimum library size per cell/spot to be included."
    ),
    min_cells = p_int(
      10L,
      doc = paste(
        "Minimum number of cells a gene has to be expressed to be",
        "included."
      )
    ),
    target_size = p_dbl(10000, doc = "The target size for the normalisation.")
  )
)

spec_sc_hvg <- param_spec(
  name = "sc_hvg",
  title = "Wrapper function for HVG detection parameters.",
  checker = "ScHvg",
  label = "HVG params",
  hint = "loess_span must be in [0.1, 1]; num_bin must be an integer.",
  fields = list(
    method = p_choice(
      "vst",
      c("vst", "meanvarbin", "dispersion", "residual"),
      doc = paste(
        "`\"residual\"` ranks genes by the residual variance of a",
        "model fitted with [bixverse::fit_residuals_sc()], and needs",
        "that fit on the object first. It also treats `hvg_no` as a",
        "per-group count and returns the union across groups, so a",
        "grouped fit can select more than `hvg_no` genes."
      )
    ),
    loess_span = p_dbl(
      0.3,
      range = "[0.1, 1]",
      doc = paste(
        "The span parameter for the loess function that is used to",
        "standardise the variance for `method = \"vst\"`."
      )
    ),
    num_bin = p_int(20L, doc = "Not yet implemented."),
    bin_method = p_choice(
      "equal_width",
      c("equal_width", "equal_freq"),
      doc = "The binning method."
    )
  )
)

spec_sc_pca <- param_spec(
  name = "sc_pca",
  title = "Wrapper for PCA specifically designed for single cells",
  checker = "ScPca",
  label = "single cell PCA params",
  hint = paste(
    "mean_center, normalise_variance, randomised and clr must be",
    "single booleans; size_factor must be a single numeric."
  ),
  fields = list(
    mean_center = p_lgl(TRUE, doc = "Shall the data be mean centred"),
    normalise_variance = p_lgl(
      TRUE,
      doc = "Shall the data have normalised variance"
    ),
    randomised = p_lgl(
      TRUE,
      doc = "Shall fast, approximate randomised SVD be used."
    ),
    clr = p_lgl(
      FALSE,
      doc = paste(
        "Shall the CLR-type `PFlogPF` be applied, see Booeshaghi, et",
        "al."
      )
    ),
    size_factor = p_dbl(
      10000,
      doc = paste(
        "The used size factor during I/O. It needs to be the same as",
        "during I/O to have correct results when using the `PFlogPF`",
        "transformation."
      )
    )
  )
)

spec_sc_sctransform <- param_spec(
  name = "sc_sctransform",
  title = "Wrapper function for scTransform (v2) parameters",
  description = paste(
    "Defaults are sctransform's own with `vst.flavor = \"v2\"`",
    "applied. Only the step-1 fit scales with `n_genes` and",
    "`n_cells`: those bound the subsample the negative binomial",
    "models are fitted on, and every later pass streams gene by",
    "gene. Raising them costs fitting time, not memory."
  ),
  references = "Choudhary and Satija, Genome Biology, 2022.",
  checker = "ScSctransform",
  label = "scTransform params",
  hint = paste(
    "n_genes, n_cells and min_cells must be positive integers;",
    "bw_adjust, outlier_th and poisson_diff_theta must be",
    "positive; clip_min and clip_max must be single numerics or",
    "NULL."
  ),
  extra_ctor = quote(assert_clip_range(clip_min, clip_max)),
  extra_check = quote({
    res <- check_clip_pair(x, "scTransform params")
    if (!isTRUE(res)) {
      return(res)
    }
  }),
  fields = list(
    n_genes = p_int(
      2000L,
      range = "[1,)",
      doc = "Genes in the step-1 subsample."
    ),
    n_cells = p_int(
      2000L,
      range = "[1,)",
      doc = "Cells in the step-1 subsample."
    ),
    min_cells = p_int(
      5L,
      range = "[1,)",
      doc = paste(
        "Minimum number of cells a gene must be detected in to be",
        "modelled."
      )
    ),
    bw_adjust = p_dbl(
      3,
      range = "(0,)",
      doc = paste(
        "Bandwidth multiplier for the kernel regression that",
        "regularises the parameters."
      )
    ),
    gmean_eps = p_dbl(1, range = "[0,)", doc = "Offset in the geometric mean."),
    outlier_th = p_dbl(
      10,
      range = "(0,)",
      doc = paste(
        "Threshold, in median absolute deviations, past which a",
        "step-1 fit is treated as an outlier."
      )
    ),
    poisson_diff_theta = p_dbl(
      0.001,
      range = "(0,)",
      doc = paste(
        "Below this, the fitted dispersion is taken as the Poisson",
        "limit."
      )
    ),
    clip_min = p_dbl(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "Lower residual clipping bound. `NULL` uses `-sqrt(n_cells)`.",
        "Must be given together with `clip_max`."
      )
    ),
    clip_max = p_dbl(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "Upper residual clipping bound. `NULL` uses `sqrt(n_cells)`.",
        "Must be given together with `clip_min`."
      )
    )
  )
)

spec_sc_apr <- param_spec(
  name = "sc_apr",
  title = "Wrapper function for analytic Pearson residual parameters",
  description = paste(
    "The closed-form alternative to scTransform: one shared",
    "dispersion instead of a fitted model per gene. Much cheaper,",
    "and on most data sets it ranks genes about as well."
  ),
  references = "Lause, Berens and Kobak, Genome Biology, 2021.",
  checker = "ScApr",
  label = "analytic Pearson params",
  hint = paste(
    "theta must be positive (Inf is allowed); min_cells must be a",
    "non-negative integer; clip_min and clip_max must be single",
    "numerics or NULL."
  ),
  extra_ctor = quote(assert_clip_range(clip_min, clip_max)),
  extra_check = quote({
    res <- check_clip_pair(x, "analytic Pearson params")
    if (!isTRUE(res)) {
      return(res)
    }
  }),
  fields = list(
    theta = p_dbl(
      100,
      range = "(0,]",
      doc = paste(
        "The shared negative binomial dispersion. `Inf` gives the",
        "Poisson limit."
      )
    ),
    min_cells = p_int(
      5L,
      range = "[0,)",
      doc = paste(
        "Minimum number of cells a gene must be detected in to be",
        "retained. `0L` keeps everything."
      )
    ),
    clip_min = p_dbl(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "Lower residual clipping bound. `NULL` uses `-sqrt(n_cells)`.",
        "Must be given together with `clip_max`."
      )
    ),
    clip_max = p_dbl(
      NULL,
      null_ok = TRUE,
      doc = paste(
        "Upper residual clipping bound. `NULL` uses `sqrt(n_cells)`.",
        "Must be given together with `clip_min`."
      )
    )
  )
)

spec_sc_knn <- param_spec(
  name = "sc_knn",
  title = "Parameters for single cell kNN searches",
  checker = "ScKnn",
  label = "kNN params",
  hint = paste(
    "k, n_trees, m, ef_construction and ef_search must be",
    "positive integers; delta must be a positive numeric;",
    "diversify_prob must be a numeric in [0, 1]; extract_knn must",
    "be a single boolean; search_budget, ef_budget, n_list and",
    "n_probe must be NULL or positive integers."
  ),
  fields = list(
    k = p_int(15L, range = "[1,)", doc = "Number of neighbours."),
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
      c("cosine", "euclidean"),
      doc = "Distance metric to use."
    ),
    n_trees = p_int(50L, range = "[1,)", doc = "Annoy param: number of trees."),
    search_budget = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "Annoy param: optional search budget per tree. If `NULL`,",
        "defaults to `n_trees * k * 20L` internally."
      )
    ),
    delta = p_dbl(
      0.001,
      range = "(0,)",
      doc = "NNDescent param: early termination criterion."
    ),
    diversify_prob = p_dbl(
      0,
      range = "[0,1]",
      doc = paste(
        "NNDescent param: diversification probability applied at the",
        "end of index construction."
      )
    ),
    ef_budget = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "NNDescent param: optional query budget. Higher values",
        "improve recall at the cost of speed."
      )
    ),
    extract_knn = p_lgl(
      FALSE,
      doc = paste(
        "NNDescent param: hand back the graph the descent already",
        "built instead of beam searching it. Skips the query pass",
        "entirely, so it is much faster, at the cost of some recall.",
        "Rows the descent never filled come back padded with",
        "duplicate edges. Ignored by every other method."
      )
    ),
    m = p_int(
      16L,
      range = "[1,)",
      doc = "HNSW param: number of connections between layers."
    ),
    ef_construction = p_int(
      200L,
      range = "[1,)",
      doc = paste(
        "HNSW param: size of the dynamic candidate list during",
        "construction."
      )
    ),
    ef_search = p_int(
      100L,
      range = "[1,)",
      doc = paste(
        "HNSW param: size of the candidate list at query time. Higher",
        "values improve recall at the cost of speed."
      )
    ),
    n_list = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "IVF param: number of clusters to generate. If `NULL`,",
        "defaults to `sqrt(n)` internally."
      )
    ),
    n_probe = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      doc = paste(
        "IVF param: number of clusters to query. If `NULL`, defaults",
        "to `sqrt(n_list)` internally."
      )
    )
  )
)

spec_sc_dsb <- param_spec(
  name = "sc_dsb",
  title = "Default parameters for DSB ADT normalisation",
  checker = "ScDsb",
  label = "DSB params",
  hint = paste(
    "denoise_counts and use_isotype_controls must be single",
    "logicals; pseudocount must be a single numeric > 0;",
    "quantile_low must be NULL or in [0, 1); quantile_high must",
    "be NULL or in (0, 1]."
  ),
  extra_ctor = quote({
    if (xor(is.null(quantile_low), is.null(quantile_high))) {
      stop("quantile_low and quantile_high must both be provided or both NULL.")
    }
    if (
      !is.null(quantile_low) &&
        !is.null(quantile_high) &&
        quantile_low >= quantile_high
    ) {
      stop("quantile_low must be strictly less than quantile_high.")
    }
  }),
  extra_check = quote({
    if (xor(is.null(x$quantile_low), is.null(x$quantile_high))) {
      return(
        "quantile_low and quantile_high must both be provided or both NULL."
      )
    }
    if (
      !is.null(x$quantile_low) &&
        !is.null(x$quantile_high) &&
        x$quantile_low >= x$quantile_high
    ) {
      return("quantile_low must be strictly less than quantile_high.")
    }
  }),
  fields = list(
    denoise_counts = p_lgl(
      TRUE,
      doc = "Run Step II (cell-to-cell technical noise removal)."
    ),
    use_isotype_controls = p_lgl(
      TRUE,
      doc = paste(
        "Include isotype controls in the noise matrix in Step II.",
        "Requires `isotype_indices` to be passed at call time."
      )
    ),
    pseudocount = p_dbl(
      10,
      range = "(0,)",
      doc = paste(
        "Pseudocount added before the log transform. The DSB paper",
        "recommends `10` with empty droplets and `1` without."
      )
    ),
    quantile_low = p_dbl(
      NULL,
      range = "[0,1)",
      null_ok = TRUE,
      doc = paste(
        "Optional numeric in `[0, 1)`. Lower quantile for per-protein",
        "output clipping. If `NULL` (and `quantile_high` is also",
        "`NULL`), no clipping is applied."
      )
    ),
    quantile_high = p_dbl(
      NULL,
      range = "(0,1]",
      null_ok = TRUE,
      doc = paste(
        "Optional numeric in `(0, 1]`. Upper quantile for per-protein",
        "output clipping. If `NULL` (and `quantile_low` is also",
        "`NULL`), no clipping is applied."
      )
    )
  )
)

spec_sctype_cells <- param_spec(
  name = "sctype_cells",
  title = "Parameters for the per-cell ScType assignment",
  description = paste(
    "Controls the per-cell path of [assign_sc_type()]: how the",
    "raw ScType scores are rescaled, how hard the scores get",
    "smoothed over the sNN graph, and where the cut-offs for an",
    "Unknown call and for a mixed cluster sit."
  ),
  references = "Zhou et al., NIPS, 2004.",
  checker = "SctypeCell",
  label = "ScType cell params",
  class_tag = "params_sctype_cells",
  fields = list(
    alpha = p_dbl(
      0.5,
      range = "[0,1]",
      doc = paste(
        "Numeric in `[0, 1]`. Self-retention during smoothing. Each",
        "iteration computes `alpha * original + (1 - alpha) *",
        "neighbour_average`."
      )
    ),
    iterations = p_int(
      2L,
      range = "[0,)",
      doc = paste(
        "Integer >= 0. Number of smoothing iterations. `0` disables",
        "smoothing."
      )
    ),
    tolerance = p_dbl(
      1e-04,
      range = "(0,)",
      doc = "Numeric > 0. Convergence tolerance for the smoothing."
    ),
    calibration = p_choice(
      "none",
      c("none", "column_z"),
      doc = paste(
        "`\"column_z\"` standardises each cell type's score column",
        "across cells, which removes the bias towards cell types",
        "whose marker sets happen to produce larger scores."
      )
    ),
    score_floor = p_dbl(
      0.25,
      range = "[0,)",
      doc = paste(
        "Numeric >= 0. Minimum score for a cell to get a call instead",
        "of `NA`."
      )
    ),
    purity_threshold = p_dbl(
      0.9,
      range = "[0,1]",
      doc = paste(
        "Numeric in `[0, 1]`. Cluster purity above which the hybrid",
        "assignment keeps the cluster-level call."
      )
    )
  )
)

spec_symphony_map <- param_spec(
  name = "symphony_map",
  title = "Default parameters for Symphony query mapping",
  checker = "SymphonyMap",
  label = "Symphony map params",
  hint = "sigma and lambda must be non-negative numerics.",
  class_tag = "params_symphony_map",
  fields = list(
    sigma = p_dbl(
      0.1,
      range = "[0,)",
      doc = paste(
        "Soft-clustering fuzziness for query -> reference centroid",
        "assignment. Symphony R default is 0.1."
      )
    ),
    lambda = p_dbl(
      1,
      range = "[0,)",
      doc = paste(
        "Ridge penalty on batch coefficients. Symphony R hardcodes",
        "1.0."
      )
    )
  )
)

spec_ligand_target <- param_spec(
  name = "ligand_target",
  title = "Parameters for ligand to target influence computation",
  checker = "LigandTarget",
  label = "ligand-target params",
  hint = paste(
    "lr_sig_hub, gr_hub, ltf_cutoff and damping_factor must be single",
    "numerics in [0, 1]; tol must be a single positive numeric; max_iter",
    "must be a single positive integer; topology_correction and",
    "secondary_targets must be single booleans."
  ),
  extra_ctor = quote(max_iter <- as.integer(max_iter)),
  fields = list(
    lr_sig_hub = p_dbl(
      0,
      range = "[0,1]",
      doc = paste(
        "Numeric in `[0, 1]`. Hub correction strength for the",
        "ligand-receptor / signalling layer. 0 disables correction."
      )
    ),
    gr_hub = p_dbl(
      0,
      range = "[0,1]",
      doc = paste(
        "Numeric in `[0, 1]`. Hub correction strength for the gene",
        "regulatory layer. 0 disables correction."
      )
    ),
    ltf_cutoff = p_dbl(
      0.99,
      range = "[0,1]",
      doc = paste(
        "Numeric in `[0, 1]`. Quantile cutoff applied to the",
        "intermediate ligand-to-TF matrix."
      )
    ),
    damping_factor = p_dbl(
      0.5,
      range = "[0,1]",
      doc = "Numeric in `[0, 1]`. PageRank-style damping factor."
    ),
    tol = p_dbl(
      1e-06,
      range = "(0,)",
      doc = "Numeric > 0. Convergence tolerance for the propagation step."
    ),
    max_iter = p_int(
      1000L,
      range = "[1,)",
      integerish = TRUE,
      doc = "Integer >= 1. Maximum iterations for the propagation step."
    ),
    topology_correction = p_lgl(FALSE, doc = "Apply topology correction."),
    secondary_targets = p_lgl(
      FALSE,
      doc = "Run a second round through targets."
    )
  )
)
