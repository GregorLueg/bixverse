devtools::load_all()
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

object <- rbh_graph(
  full_data,
  dataset_col = "origin",
  module_col = "ind",
  value_col = "values"
)

object <- generate_rbh_graph(
  object,
  minimum_similarity = .1,
  overlap_coefficient = TRUE,
  .debug = FALSE
)

overlap_res <- get_rbh_res(object)

list_of_list <- object@module_data

rs_rbh_sets(
  list_of_list,
  overlap_coefficient = FALSE,
  min_similarity = 0.1,
  TRUE
)

rextendr::document()

length(intersect(list_of_list$set_A$A, list_of_list$set_B$b)) /
  length(union(list_of_list$set_A$A, list_of_list$set_B$b))


length(intersect(list_of_list$set_A$C, list_of_list$set_B$a)) /
  length(union(list_of_list$set_A$C, list_of_list$set_B$a))
