# rbh test ---------------------------------------------------------------------

## generate data ---------------------------------------------------------------

set.seed(123)

modules_set_a <- purrr::map(
  1:5,
  ~ {
    sample(letters, 10)
  }
)
names(modules_set_a) <- LETTERS[1:5]
data_a <- data.table::setDT(stack(modules_set_a))[, `:=`(
  origin = "set_A",
  ind = as.character(ind)
)]

modules_set_b <- purrr::map(
  1:3,
  ~ {
    sample(letters, 8)
  }
)
names(modules_set_b) <- LETTERS[1:3]
data_b <- data.table::setDT(stack(modules_set_b))[, `:=`(
  origin = "set_B",
  ind = as.character(ind)
)]

full_data <- data.table::rbindlist(list(data_a, data_b))

## expected similarities -------------------------------------------------------

overlap_coef_sims <- c(0.5, 0.375, 0.375, 0.375, 0.375)
jaccard_sim <- c(0.2857143, 0.2, 0.2, 0.2, 0.2)
origin <- c("set_A_A", "set_A_B", "set_A_C", "set_A_C", "set_A_E")
targets <- c("set_B_B", "set_B_A", "set_B_A", "set_B_C", "set_B_C")

## test class ------------------------------------------------------------------

### overlap coefficients -------------------------------------------------------

object <- rbh_graph(
  full_data,
  dataset_col = "origin",
  module_col = "ind",
  value_col = "values"
)

object <- generate_rbh_graph(
  object,
  minimum_similarity = 0,
  overlap_coefficient = TRUE,
  .debug = FALSE
)

overlap_res <- get_rbh_res(object)

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

### jaccard --------------------------------------------------------------------

object <- generate_rbh_graph(
  object,
  minimum_similarity = 0,
  overlap_coefficient = FALSE,
  .debug = FALSE
)

overlap_res <- get_rbh_res(object)

expect_equivalent(
  current = overlap_res$similiarity,
  target = jaccard_sim,
  info = "Overlap coefficients test for RBH graph",
  tolerance = 10e-6
)
