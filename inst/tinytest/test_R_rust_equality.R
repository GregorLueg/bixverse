# correlation, covariance, equality tests --------------------------------------

## simple versions -------------------------------------------------------------

set.seed(123)
mat <- matrix(data = rnorm(100), nrow = 10, ncol = 10)
rownames(mat) <- sprintf("sample_%i", 1:10)
colnames(mat) <- sprintf("feature_%i", 1:10)

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
mat_2 <- matrix(data = rnorm(100), nrow = 10, ncol = 10)
rownames(mat_2) <- sprintf("sample_%i", 1:10)
colnames(mat_2) <- sprintf("feature_%i", 1:10)

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

expect_equivalent(
  current = abs(rs_pca_res$v),
  target = abs(r_pca_res$rotation),
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

unique_nodes <- unique(c(edge_dt$from, edge_dt$to))

personalised_v1 <- c(1, 0, 0, 0, 0)
personalised_v2 <- rep(1, 5) / 5

## tests -----------------------------------------------------------------------

if (requireNamespace("igraph", quietly = TRUE)) {
  g_undir <- igraph::graph_from_data_frame(edge_dt, directed = FALSE)
  g_dir <- igraph::graph_from_data_frame(edge_dt, directed = TRUE)

  igraph_res_undir_v1 <- igraph::page_rank(
    graph = g_undir,
    personalized = personalised_v1
  )$vector

  igraph_res_dir_v1 <- igraph::page_rank(
    graph = g_dir,
    personalized = personalised_v1
  )$vector

  igraph_res_undir_v2 <- igraph::page_rank(
    graph = g_undir,
    personalized = personalised_v2
  )$vector

  igraph_res_dir_v2 <- igraph::page_rank(
    graph = g_dir,
    personalized = personalised_v2
  )$vector

  rs_res_undir_v1 <- rs_page_rank(
    node_names = unique_nodes,
    from = edge_dt$from,
    to = edge_dt$to,
    personalised = personalised_v1,
    undirected = TRUE
  )

  rs_res_dir_v1 <- rs_page_rank(
    node_names = unique_nodes,
    from = edge_dt$from,
    to = edge_dt$to,
    personalised = personalised_v1,
    undirected = FALSE
  )

  rs_res_undir_v2 <- rs_page_rank(
    node_names = unique_nodes,
    from = edge_dt$from,
    to = edge_dt$to,
    personalised = personalised_v2,
    undirected = TRUE
  )

  rs_res_dir_v2 <- rs_page_rank(
    node_names = unique_nodes,
    from = edge_dt$from,
    to = edge_dt$to,
    personalised = personalised_v2,
    undirected = FALSE
  )

  cor_undir_v1 <- cor(igraph_res_undir_v1, rs_res_undir_v1)

  cor_dir_v1 <- cor(igraph_res_dir_v1, rs_res_dir_v1)

  cor_undir_v2 <- cor(igraph_res_undir_v2, rs_res_undir_v2)

  cor_dir_v2 <- cor(igraph_res_dir_v2, rs_res_dir_v2)

  expect_true(
    cor_undir_v1 > 0.99,
    info = "Rust personalsied Page Rank implementation undirected network (v1)."
  )

  expect_true(
    cor_dir_v1 > 0.99,
    info = "Rust personalsied Page Rank implementation directed network (v1)."
  )

  expect_true(
    cor_undir_v2 > 0.99,
    info = "Rust personalsied Page Rank implementation undirected network (v2)."
  )

  expect_true(
    cor_dir_v2 > 0.99,
    info = "Rust personalsied Page Rank implementation directed network (v2)."
  )
} else {
  exit_file("igraph package not available for comparison tests")
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
