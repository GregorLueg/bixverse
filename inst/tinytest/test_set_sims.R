# rbh test ---------------------------------------------------------------------

## generate data ---------------------------------------------------------------

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

## expected similarities -------------------------------------------------------

overlap_coef_sims <- c(0.625, 0.375)
jaccard_sim <- c(0.3846154, 0.2)
origin <- c("set_A_C", "set_A_A")
targets <- c("set_B_a", "set_B_b")

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

### jaccard --------------------------------------------------------------------

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
