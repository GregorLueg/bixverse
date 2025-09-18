# correlation, covariance, equality tests --------------------------------------

## simple versions -------------------------------------------------------------

set.seed(123)
mat <- matrix(data = rnorm(12 * 14), nrow = 12, ncol = 14)
rownames(mat) <- sprintf("sample_%i", 1:12)
colnames(mat) <- sprintf("feature_%i", 1:14)

# Pearson
expect_equivalent(
  current = rs_cor(mat, spearman = FALSE),
  target = cor(mat),
  info = "Correlation equivalence test Rust <> R"
)
# Spearman
expect_equivalent(
  current = rs_cor(mat, spearman = TRUE),
  target = cor(mat, method = "spearman"),
  info = "Spearman Correlation equivalence test Rust <> R"
)
# Co-variance
expect_equivalent(
  current = rs_covariance(mat),
  target = cov(mat),
  info = "Covariance equivalence test Rust <> R"
)

if (requireNamespace("coop", quietly = TRUE)) {
  expect_equivalent(
    current = rs_cos(mat),
    target = coop::cosine(mat),
    info = "Cosne equivalence test Rust <> R"
  )
}

## correlation two matrices ----------------------------------------------------

set.seed(246)
mat_2 <- matrix(data = rnorm(12 * 8), nrow = 12, ncol = 8)
rownames(mat_2) <- sprintf("sample_%i", 1:12)
colnames(mat_2) <- sprintf("feature_%i", 1:8)

# Pearson - two matrices
expect_equivalent(
  current = rs_cor2(mat, mat_2, spearman = FALSE),
  target = cor(mat, mat_2),
  info = "Correlation equivalence test Rust <> R (two matrices)"
)
# Spearman
expect_equivalent(
  current = rs_cor2(mat, mat_2, spearman = TRUE),
  target = cor(mat, mat_2, method = "spearman"),
  info = "Spearman Correlation equivalence test Rust <> R (two matrices)"
)

## upper triangle versions -----------------------------------------------------

# Check if the upper triangle class behaves as expected
cor_data <- rs_cor_upper_triangle(mat, spearman = FALSE, shift = 1L)

cor_class <- bixverse:::upper_triangular_sym_mat$new(
  values = cor_data,
  features = colnames(mat),
  shift = 1L
)

expect_equal(
  current = cor_class$get_sym_matrix(.verbose = FALSE),
  target = cor(mat),
  info = "Upper triangle class test Rust <> R"
)

## covariance to cor -----------------------------------------------------------

cov_data <- rs_covariance(mat)

expect_equivalent(
  current = rs_cov2cor(cov_data),
  target = cor(mat),
  info = "cov2cor equivalence test Rust <> R"
)

## mutual information ----------------------------------------------------------

if (requireNamespace("infotheo", quietly = TRUE)) {
  # equal width strategy

  mi_rust_result_1.1 <- rs_mutual_info(
    mat,
    n_bins = NULL,
    strategy = "equal_width",
    normalise = FALSE
  )

  # ensure that the same discretisation is used
  infotheo_res_1.1 <- infotheo::mutinformation(infotheo::discretize(
    mat,
    disc = "equalwidth",
    nbins = sqrt(nrow(mat))
  ))

  mi_rust_result_1.2 <- rs_mutual_info(
    mat_2,
    n_bins = NULL,
    strategy = "equal_width",
    normalise = FALSE
  )

  # ensure that the same discretisation is used
  infotheo_res_1.2 <- infotheo::mutinformation(infotheo::discretize(
    mat_2,
    disc = "equalwidth",
    nbins = sqrt(nrow(mat))
  ))

  expect_equivalent(
    current = mi_rust_result_1.1,
    target = infotheo_res_1.1,
    info = "mutual information rest Rust <> R - equal width - matrix 1"
  )

  expect_equivalent(
    current = mi_rust_result_1.2,
    target = infotheo_res_1.2,
    info = "mutual information rest Rust <> R - equal width - matrix 2"
  )

  # equal frequency strategy

  mi_rust_result_2.1 <- rs_mutual_info(
    mat,
    n_bins = NULL,
    strategy = "equal_freq",
    normalise = FALSE
  )

  infotheo_res_2.1 <- infotheo::mutinformation(infotheo::discretize(
    mat,
    disc = "equalfreq",
    nbins = sqrt(nrow(mat))
  ))

  mi_rust_result_2.2 <- rs_mutual_info(
    mat_2,
    n_bins = NULL,
    strategy = "equal_freq",
    normalise = FALSE
  )

  infotheo_res_2.2 <- infotheo::mutinformation(infotheo::discretize(
    mat_2,
    disc = "equalfreq",
    nbins = sqrt(nrow(mat))
  ))

  expect_equivalent(
    current = mi_rust_result_2.1,
    target = infotheo_res_2.1,
    info = "mutual information rest Rust <> R - equal distance - matrix 1"
  )

  expect_equivalent(
    current = mi_rust_result_2.2,
    target = infotheo_res_2.2,
    info = "mutual information rest Rust <> R - equal distance - matrix 2"
  )
}

## pointwise mutual information ------------------------------------------------

set.seed(123)
mat_lgl <- matrix(
  data = as.logical(rbinom(n = 10 * 12, size = 1, prob = 0.35)),
  nrow = 10,
  ncol = 12
)
rownames(mat_lgl) <- sprintf("sample_%i", 1:10)
colnames(mat_lgl) <- sprintf("feature_%i", 1:12)

expected_res <- qs2::qs_read("./test_data/npmi_res.qs")
expected_res_norm <- qs2::qs_read("./test_data/npmi_res_norm.qs")

res <- rs_pointwise_mutual_info(mat_lgl, normalise = FALSE)
res_normalised <- rs_pointwise_mutual_info(mat_lgl, normalise = TRUE)

expect_equal(
  current = res,
  target = expected_res,
  info = "nPMI calculations - non normalised"
)

expect_equal(
  current = res_normalised,
  target = expected_res_norm,
  info = "nPMI calculations - normalised"
)

## distances -------------------------------------------------------------------

### euclidean ------------------------------------------------------------------

rust_res <- rs_dist(mat, distance_type = "euclidean")
r_res <- as.matrix(dist(t(mat), method = "euclidean"))

expect_equivalent(
  current = rust_res,
  target = r_res,
  info = "Distance - Euclidean"
)

### manhattan ------------------------------------------------------------------

rust_res <- rs_dist(mat, distance_type = "manhattan")
r_res <- as.matrix(dist(t(mat), method = "manhattan"))

expect_equivalent(
  current = rust_res,
  target = r_res,
  info = "Distance - Manhattan"
)

### canberra -------------------------------------------------------------------

rust_res <- rs_dist(mat, distance_type = "canberra")
r_res <- as.matrix(dist(t(mat), method = "canberra"))

expect_equivalent(
  current = rust_res,
  target = r_res,
  info = "Distance - Canberra"
)

### cosine distance ------------------------------------------------------------

if (requireNamespace("coop", quietly = TRUE)) {
  rust_res <- rs_dist(mat, distance_type = "cosine")
  r_res <- 1 - abs(coop::cosine(mat))

  expect_equivalent(
    current = rust_res,
    target = r_res,
    info = "Distance - Cosine"
  )
}

# hypergeom distributions ------------------------------------------------------

m <- 10
n <- 7
k <- 8
x <- 0:(k + 1)

rust_vals <- purrr::map_dbl(
  x,
  ~ {
    rs_phyper(
      q = .x,
      m = m,
      n = n,
      k = k
    )
  }
)

r_vals <- purrr::map_dbl(
  x,
  ~ {
    phyper(
      q = .x,
      m = m,
      n = n,
      k = k,
      lower.tail = F
    )
  }
)

expect_equal(
  current = rust_vals,
  target = r_vals,
  info = "Hypergeometric test values for Rust <> R."
)

# pca --------------------------------------------------------------------------

r_pca_res <- prcomp(mat)

rs_pca_res <- rs_prcomp(mat, scale = FALSE)

# Absolute as signs can flip
expect_equivalent(
  current = abs(rs_pca_res$scores),
  target = abs(r_pca_res$x),
  info = "Factors for PCA (via SVD) for Rust <> R."
)

# The last remaining PCs are approaching numerical null space and are not
# stable in the implementations anymore... They can be ignored. At most the
# first 11 PCs should be the same with this data
expect_equivalent(
  current = abs(rs_pca_res$v[, 1:11]),
  target = abs(r_pca_res$rotation[, 1:11]),
  info = "Loadings/rotation for PCA (via SVD) for Rust <> R."
)

expect_equal(
  current = rs_pca_res$s,
  target = r_pca_res$sdev,
  info = "Standard deviation of the eigenvalues for PCA (via SVD) for Rust <> R."
)

# graph stuff ------------------------------------------------------------------

## data ------------------------------------------------------------------------

edge_dt <- data.table::data.table(
  from = c("a", "b", "c", "d", "d"),
  to = c("b", "c", "d", "a", "e")
)

edge_dt_weighted <- data.table::data.table(
  from = c("a", "b", "c", "d", "d"),
  to = c("b", "c", "d", "a", "e"),
  weight = c(1, 1, 0.5, 0.4, 0.25)
)

unique_nodes <- unique(c(edge_dt$from, edge_dt$to))

personalised_v1 <- c(1, 0, 0, 0, 0)
personalised_v2 <- rep(1, 5) / 5

## tests -----------------------------------------------------------------------

if (requireNamespace("igraph", quietly = TRUE)) {
  # graphs
  g_undir <- igraph::graph_from_data_frame(edge_dt, directed = FALSE)
  g_dir <- igraph::graph_from_data_frame(edge_dt, directed = TRUE)
  g_weighted <- igraph::graph_from_data_frame(
    edge_dt_weighted,
    directed = FALSE
  )

  # version 1 - igraph
  igraph_res_undir_v1 <- igraph::page_rank(
    graph = g_undir,
    personalized = personalised_v1
  )$vector

  igraph_res_dir_v1 <- igraph::page_rank(
    graph = g_dir,
    personalized = personalised_v1
  )$vector

  igraph_res_weighted_v1 <- igraph::page_rank(
    graph = g_weighted,
    personalized = personalised_v1
  )$vector

  # version 2 - igraph
  igraph_res_undir_v2 <- igraph::page_rank(
    graph = g_undir,
    personalized = personalised_v2
  )$vector

  igraph_res_dir_v2 <- igraph::page_rank(
    graph = g_dir,
    personalized = personalised_v2
  )$vector

  igraph_res_weighted_v2 <- igraph::page_rank(
    graph = g_weighted,
    personalized = personalised_v2
  )$vector

  # version 1 - rust
  rs_res_undir_v1 <- rs_page_rank(
    node_names = unique_nodes,
    from = edge_dt$from,
    to = edge_dt$to,
    weights = NULL,
    personalised = personalised_v1,
    undirected = TRUE
  )

  rs_res_dir_v1 <- rs_page_rank(
    node_names = unique_nodes,
    from = edge_dt$from,
    to = edge_dt$to,
    weights = NULL,
    personalised = personalised_v1,
    undirected = FALSE
  )

  rs_res_weighted_v1 <- rs_page_rank(
    node_names = unique_nodes,
    from = edge_dt_weighted$from,
    to = edge_dt_weighted$to,
    weights = edge_dt_weighted$weight,
    personalised = personalised_v1,
    undirected = TRUE
  )

  # version 2 - rust
  rs_res_undir_v2 <- rs_page_rank(
    node_names = unique_nodes,
    from = edge_dt$from,
    to = edge_dt$to,
    weights = NULL,
    personalised = personalised_v2,
    undirected = TRUE
  )

  rs_res_dir_v2 <- rs_page_rank(
    node_names = unique_nodes,
    from = edge_dt$from,
    to = edge_dt$to,
    weights = NULL,
    personalised = personalised_v2,
    undirected = FALSE
  )

  rs_res_weighted_v2 <- rs_page_rank(
    node_names = unique_nodes,
    from = edge_dt_weighted$from,
    to = edge_dt_weighted$to,
    weights = edge_dt_weighted$weight,
    personalised = personalised_v2,
    undirected = TRUE
  )

  # version 1
  cor_undir_v1 <- cor(igraph_res_undir_v1, rs_res_undir_v1)
  cor_dir_v1 <- cor(igraph_res_dir_v1, rs_res_dir_v1)
  cor_weighted_v1 <- cor(igraph_res_weighted_v1, rs_res_weighted_v1)

  expect_true(
    cor_undir_v1 > 0.99,
    info = "Rust personalsied Page Rank implementation undirected network (v1)."
  )

  expect_true(
    cor_dir_v1 > 0.99,
    info = "Rust personalsied Page Rank implementation directed network (v1)."
  )

  expect_true(
    cor_weighted_v1 > 0.99,
    info = "Rust personalsied Page Rank implementation weighted network (v1)."
  )

  # version 2
  cor_undir_v2 <- cor(igraph_res_undir_v2, rs_res_undir_v2)
  cor_dir_v2 <- cor(igraph_res_dir_v2, rs_res_dir_v2)
  cor_weighted_v2 <- cor(igraph_res_weighted_v2, rs_res_weighted_v2)

  expect_true(
    cor_undir_v2 > 0.99,
    info = "Rust personalsied Page Rank implementation undirected network (v2)."
  )

  expect_true(
    cor_dir_v2 > 0.99,
    info = "Rust personalsied Page Rank implementation directed network (v2)."
  )

  expect_true(
    cor_weighted_v2 > 0.99,
    info = "Rust personalsied Page Rank implementation weighted network (v2)."
  )
}

# set similarities -------------------------------------------------------------

## data ------------------------------------------------------------------------

set_a <- letters[1:5]
set_b <- letters[2:7]

jaccard <- length(intersect(set_a, set_b)) / length(union(set_a, set_b))
overlap_coef <- length(intersect(set_a, set_b)) /
  min(c(length(set_a), length(set_b)))

## results ---------------------------------------------------------------------

rs_jaccard <- rs_set_similarity(set_a, set_b, overlap_coefficient = FALSE)
rs_overlap_coef <- rs_set_similarity(set_a, set_b, overlap_coefficient = TRUE)

expect_equal(
  current = jaccard,
  target = rs_jaccard,
  info = "Jaccard similarity Rust <> R"
)

expect_equal(
  current = overlap_coef,
  target = rs_overlap_coef,
  info = "Overlap coefficient Rust <> R"
)

## matrix version --------------------------------------------------------------

### data -----------------------------------------------------------------------

set_c <- letters[5:11]
set_d <- letters[3:7]
set_e <- letters[20:25]

list_a <- list(set_a = set_a, set_b = set_b)
list_b <- list(set_c = set_c, set_d = set_d, set_e = set_e)

expected_matrix_jaccard <- matrix(
  data = c(0.09090909, 0.4285714, 0, 0.3, 0.8333333, 0),
  nrow = 2,
  ncol = 3,
  byrow = TRUE
)

expected_matrix_overlap <- matrix(
  data = c(0.2, 0.6, 0, 0.5, 1, 0),
  nrow = 2,
  ncol = 3,
  byrow = TRUE
)

### tests ----------------------------------------------------------------------

jaccard_mat <- rs_set_similarity_list(
  s_1_list = list_a,
  s_2_list = list_b,
  overlap = FALSE
)

overlap_mat <- rs_set_similarity_list(
  s_1_list = list_a,
  s_2_list = list_b,
  overlap = TRUE
)

expect_equivalent(
  current = jaccard_mat,
  target = expected_matrix_jaccard,
  tolerance = 1e-7
)

expect_equivalent(
  current = overlap_mat,
  target = expected_matrix_overlap,
  tolerance = 1e-7
)

# effect sizes -----------------------------------------------------------------

## data ------------------------------------------------------------------------

set.seed(42L)
mat_a <- matrix(rnorm(12, mean = 5, sd = 1), nrow = 3, ncol = 4)
mat_b <- matrix(rnorm(12, mean = 7, sd = 1), nrow = 3, ncol = 4)

expected_es <- c(-0.7859511, -0.4890762, -0.1603383, -0.2796240)
expected_se <- c(0.8474333, 0.8286131, 0.8178075, 0.8204770)

expected_es_no_cor <- c(-1.2032370, -0.7487419, -0.2454669, -0.4280850)
expected_se_no_cor <- c(0.8873077, 0.8446209, 0.8195656, 0.8257954)

## tests -----------------------------------------------------------------------

rs_results <- rs_hedges_g(
  mat_a = mat_a,
  mat_b = mat_b,
  small_sample_correction = TRUE
)

rs_results_no_cor <- rs_hedges_g(
  mat_a = mat_a,
  mat_b = mat_b,
  small_sample_correction = FALSE
)

expect_equal(
  current = rs_results$effect_sizes,
  target = expected_es,
  tolerance = 1e-7,
  info = paste("Hedge effect size with correction")
)

expect_equal(
  current = rs_results$standard_errors,
  target = expected_se,
  tolerance = 1e-7,
  info = paste("Hedge standard error with correction")
)

expect_equal(
  current = rs_results_no_cor$effect_sizes,
  target = expected_es_no_cor,
  tolerance = 1e-7,
  info = paste("Hedge effect size with correction")
)

expect_equal(
  current = rs_results_no_cor$standard_errors,
  target = expected_se_no_cor,
  tolerance = 1e-7,
  info = paste("Hedge standard error with correction")
)

# loess ------------------------------------------------------------------------

x <- as.numeric(seq_len(100))
y <- 2 * x + rnorm(20, 0, 0.1)

r_fit <- loess(
  y ~ x,
  span = 0.2,
  degree = 1L,
  normalise = FALSE
)
rs_predictions <- rs_2d_loess(x, y, 0.2, 1)

# r predictions
r_predictions <- predict(r_fit, x)

expect_equal(
  current = rs_predictions$predicted,
  target = r_predictions,
  tolerance = 1e-4,
  info = "Rust version of a Loess implementation - predictions"
)

# r residuals
r_residuals <- residuals(r_fit)

expect_true(
  current = cor(rs_predictions$residuals, r_residuals) > 0.99,
  info = "Rust version of a Loess implementation - residuals"
)
