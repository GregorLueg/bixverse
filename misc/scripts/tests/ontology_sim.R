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

hpo_data <- as.data.table(stack(hpo$children)) %>%
  setnames(old = c("values", "ind"), new = c("child", "parent"))

test_class <- ontology(hpo_data, .verbose = FALSE)

test_class <- calculate_wang_sim_onto(test_class)

test_class <- filter_similarities(object = test_class, alpha = 0.001)

get_results(test_class)

test_class@final_results

get_params(test_class)

devtools::load_all()

critval <- calculate_critical_value(test_class, alpha = 0.001)

devtools::load_all()

data <- test_class@sim_mat$get_cor_vector()

filtered_data <- rs_filter_onto_sim(
  sim_vals = data$cor_data,
  names = data$features,
  threshold = critval
)

length(filtered_data$t1)

rextendr::document()

features <- test_class@sim_mat
