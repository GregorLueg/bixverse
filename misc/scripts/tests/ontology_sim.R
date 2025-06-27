test_onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c", "g"),
  child = c("b", "c", "d", "e", "f", "h")
)

library(ontologyIndex)
library(magrittr)
library(data.table)

rextendr::document()

data("hpo")

hpo_data <- as.data.table(stack(hpo$children)) %>%
  setnames(old = c("values", "ind"), new = c("child", "parent"))

c(ancestors, descendants) %<-% get_ontology_ancestry(test_onto)

ic_data <- calculate_information_content(descendants)

tictoc::tic()
combined <- calculate_semantic_sim(
  similarity_type = "combined",
  terms = sort(names(ancestors)),
  ancestor_list = ancestors,
  ic_list = ic_data
)
tictoc::toc()


combined[1:5, 1:5]

tictoc::tic()
sim_wang <- rs_onto_sim_wang(test_onto$parent, test_onto$child, w = 0.8)
tictoc::toc()

sim_wang_mat <- sim_wang$sym_mat %>%
  `colnames<-`(sim_wang$names) %>%
  `rownames<-`(sim_wang$names)

sim_wang_mat[rownames(combined)[1:5], colnames(combined)[1:5]]

sim_wang$names

summary(c(sim_wang))

sim_wang[1:5, 1:5]

#           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
# [1,] 1.0000000 0.4807122 0.6212121 0.3825911 0.8145401 0.4807122
# [2,] 0.4807122 1.0000000 0.7641509 0.4767442 0.5901639 0.5901639
# [3,] 0.6212121 0.7641509 1.0000000 0.6428571 0.7641509 0.7641509
# [4,] 0.3825911 0.4767442 0.6428571 1.0000000 0.4767442 0.4767442
# [5,] 0.8145401 0.5901639 0.7641509 0.4767442 1.0000000 0.5901639
# [6,] 0.4807122 0.5901639 0.7641509 0.4767442 0.5901639 1.0000000
