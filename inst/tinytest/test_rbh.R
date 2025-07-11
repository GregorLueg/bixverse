# rbh test ---------------------------------------------------------------------

## set similiarity -------------------------------------------------------------

### generate data --------------------------------------------------------------

set.seed(123)

modules_set_a <- purrr::map(
  1:3,
  ~ {
    sample(letters, 10)
  }
)
names(modules_set_a) <- LETTERS[1:3]
data_a <- data.table::setDT(stack(modules_set_a))[, `:=`(
  origin = "set_A",
  ind = as.character(ind)
)]

modules_set_b <- purrr::map(
  1:2,
  ~ {
    sample(letters, 8)
  }
)
names(modules_set_b) <- letters[1:2]
data_b <- data.table::setDT(stack(modules_set_b))[, `:=`(
  origin = "set_B",
  ind = as.character(ind)
)]

full_data <- data.table::rbindlist(list(data_a, data_b))

### expected similarities ------------------------------------------------------

overlap_coef_sims <- c(0.625, 0.375)
jaccard_sim <- c(0.3846154, 0.2)
origin <- c("set_A_C", "set_A_A")
targets <- c("set_B_a", "set_B_b")

### test class -----------------------------------------------------------------

#### overlap coefficients ------------------------------------------------------

object <- rbh_graph(
  full_data,
  rbh_type = "set",
  dataset_col = "origin",
  module_col = "ind",
  value_col = "values"
)

expect_error(
  current = rbh_graph(
    full_data,
    rbh_type = "cor",
    dataset_col = "origin",
    module_col = "ind",
    value_col = "values"
  ),
  info = "RBH class error with wrong input (v1)"
)

object <- generate_rbh_graph(
  object,
  minimum_similarity = 0,
  overlap_coefficient = TRUE,
  .debug = FALSE
)

overlap_res <- get_rbh_res(object) %>% data.table::setorder(-similiarity)

expect_equivalent(
  current = overlap_res$combined_origin,
  target = origin,
  info = "Origin test for RBH graph"
)
expect_equivalent(
  current = overlap_res$combined_target,
  target = targets,
  info = "Targets test for RBH graph"
)
expect_equivalent(
  current = overlap_res$similiarity,
  target = overlap_coef_sims,
  info = "Overlap coefficients test for RBH graph"
)

#### jaccard -------------------------------------------------------------------

object <- generate_rbh_graph(
  object,
  minimum_similarity = 0,
  overlap_coefficient = FALSE,
  .debug = FALSE
)

overlap_res <- get_rbh_res(object) %>% data.table::setorder(-similiarity)

expect_equivalent(
  current = overlap_res$similiarity,
  target = jaccard_sim,
  info = "Overlap coefficients test for RBH graph",
  tolerance = 10e-6
)

## correlations ----------------------------------------------------------------

### generate data --------------------------------------------------------------

set.seed(42)

matrix_a <- matrix(data = rnorm(16), nrow = 4, ncol = 4)
rownames(matrix_a) <- sprintf("gene_%i", 1:4)
colnames(matrix_a) <- sprintf("ic_%i", 1:4)

matrix_b <- matrix(data = rnorm(15), nrow = 5, ncol = 3)
rownames(matrix_b) <- sprintf("gene_%i", 1:5)
colnames(matrix_b) <- sprintf("ic_%i", 1:3)

matrix_c <- matrix(data = rnorm(20), nrow = 4, ncol = 5)
rownames(matrix_c) <- sprintf("gene_%i", 10:13)
colnames(matrix_c) <- sprintf("ic_%i", 1:5)

expected_cor_sim <- data.table(
  origin = "origin_1",
  target = "origin_2",
  origin_modules = c("ic_2", "ic_3", "ic_4"),
  target_modules = c("ic_3", "ic_1", "ic_2"),
  similiarity = c(0.8420369, 0.8487127, 0.9582737)
) %>%
  .[, `:=`(
    combined_origin = paste(origin, origin_modules, sep = "_"),
    combined_target = paste(target, target_modules, sep = "_")
  )]
expected_cor_sim_spearman <- data.table(
  origin = "origin_1",
  target = "origin_2",
  origin_modules = c("ic_2", "ic_3", "ic_4"),
  target_modules = c("ic_3", "ic_1", "ic_2"),
  similiarity = 1
) %>%
  .[, `:=`(
    combined_origin = paste(origin, origin_modules, sep = "_"),
    combined_target = paste(target, target_modules, sep = "_")
  )]


full_data <- list(
  origin_1 = matrix_a,
  origin_2 = matrix_b
)

### test class -----------------------------------------------------------------

#### pearson correlations ------------------------------------------------------

object <- rbh_graph(
  full_data,
  rbh_type = "cor"
)

expect_error(
  current = rbh_graph(
    full_data,
    rbh_type = "set"
  ),
  info = "RBH class error with wrong input (v2)"
)

object <- generate_rbh_graph(
  object,
  minimum_similarity = 0,
  spearman = FALSE
)

cor_res_pearson <- get_rbh_res(object) %>% data.table::setorder(combined_origin)

expect_equivalent(
  current = cor_res_pearson,
  target = expected_cor_sim,
  info = "RBH correlation result (pearson)",
  tolerance = 1e-7
)

#### spearman correlations -----------------------------------------------------

object <- generate_rbh_graph(
  object,
  minimum_similarity = 0,
  spearman = TRUE
)

cor_res_spearman <- get_rbh_res(object) %>%
  data.table::setorder(combined_origin)

expect_equivalent(
  current = cor_res_spearman,
  target = expected_cor_sim_spearman,
  info = "RBH correlation result (Spearman)",
  tolerance = 1e-7
)

### version with no feature overlap --------------------------------------------

# no overlapping features...
full_data_empty <- list(
  origin_1 = matrix_a,
  origin_2 = matrix_c
)

object <- rbh_graph(
  full_data_empty,
  rbh_type = "cor"
)

object <- generate_rbh_graph(
  object,
  minimum_similarity = 0,
  spearman = FALSE
)

expect_true(
  current = nrow(get_rbh_res(object)) == 0,
  info = "RBH correlation result (no overlapping features)"
)
