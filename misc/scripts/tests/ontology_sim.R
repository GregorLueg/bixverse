test_onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c", "g"),
  child = c("b", "c", "d", "e", "f", "h")
)

library(ontologyIndex)
library(magrittr)
library(data.table)

rextendr::document()

devtools::load_all()

data("hpo")

calculate_wang_sim(test_onto)

test_class <- ontology(test_onto, .verbose = FALSE)

test_class <- calculate_wang_sim_onto(test_class)

get_sim_matrix(test_class, as_data_table = TRUE)
