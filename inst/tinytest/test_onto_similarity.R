# library ----------------------------------------------------------------------

library(magrittr)
library(zeallot)

# test ontology semantic similarity --------------------------------------------

## expected values -------------------------------------------------------------

test_onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f")
)

expected_ancestors <- list(
  a = "a",
  b = c("b", "a"),
  c = c("c", "b", "a"),
  d = c("d", "b", "a"),
  e = c("e", "b", "a"),
  f = c("f", "c", "b", "a")
)

expected_descendants <- list(
  "a" = c("b", "c", "d", "e", "f", "a"),
  "b" = c("b", "c", "d", "e", "f"),
  "c" = c("c", "f"),
  "d" = c("d"),
  "e" = c("e"),
  "f" = c("f")
)

expected_ic_data <- list(
  a = 0,
  b = 0.1823216,
  c = 1.098612,
  d = 1.791759,
  e = 1.791759,
  f = 1.791759
)

expected_resnik <- readRDS("./test_data/semantic_sim_resnik.rds")
expected_lin <- readRDS("./test_data/semantic_sim_lin.rds")
expected_combined <- readRDS("./test_data/semantic_sim_combined.rds")

expected_sim_filtered <- data.table::data.table(
  term1 = "c",
  term2 = "f",
  filtered_sim = 1.098612
)

expected_critval <- 1.098612

## separate functions ----------------------------------------------------------

### ancestry ------------------------------------------------------------------

c(ancestors, descendants) %<-% get_ontology_ancestry(test_onto)

expect_equivalent(
  current = ancestors[names(expected_ancestors)],
  target = expected_ancestors,
  info = "Ontology similarity ancestors test"
)

expect_equivalent(
  current = descendants[names(expected_descendants)],
  target = expected_descendants,
  info = "Ontology similarity descendants test"
)

### information content --------------------------------------------------------

ic_data <- calculate_information_content(descendants)

expect_equivalent(
  current = ic_data[names(expected_ic_data)],
  target = expected_ic_data,
  info = "Ontology similarity test for information content.",
  tolerance = 1e-5
)

### semantic similarity --------------------------------------------------------

resnik <- calculate_semantic_sim(
  similarity_type = "resnik",
  terms = sort(names(ancestors)),
  ancestor_list = ancestors,
  ic_list = ic_data
)

lin <- calculate_semantic_sim(
  similarity_type = "lin",
  terms = sort(names(ancestors)),
  ancestor_list = ancestors,
  ic_list = ic_data
)

combined <- calculate_semantic_sim(
  similarity_type = "combined",
  terms = sort(names(ancestors)),
  ancestor_list = ancestors,
  ic_list = ic_data
)

expect_equivalent(
  current = resnik,
  target = expected_resnik,
  info = "Ontology similarity test for semantic semilarity (Resnik).",
  tolerance = 1e-5
)

expect_equivalent(
  current = lin,
  target = expected_lin,
  info = "Ontology similarity test for semantic semilarity (Lin).",
  tolerance = 1e-5
)

expect_equivalent(
  current = combined,
  target = expected_combined,
  info = "Ontology similarity test for semantic semilarity (combined type).",
  tolerance = 1e-5
)

## class -----------------------------------------------------------------------

test_class <- ontology(test_onto, .verbose = FALSE)

test_class <- calculate_semantic_sim_onto(test_class, sim_type = "resnik")

crit_val <- get_params(test_class)$semantic_similarity$critval

semantic_similarity_filtered <- get_semantic_similarities(test_class)

expect_equivalent(
  current = semantic_similarity_filtered,
  target = expected_sim_filtered,
  info = "Ontology similarity test for class with filtering on crit value - similarities",
  tolerance = 1e-5
)

expect_equivalent(
  current = crit_val,
  target = expected_critval,
  info = "Ontology similarity test for class with filtering on crit value - critical value",
  tolerance = 1e-5
)
